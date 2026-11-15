function [cand, N, G] = pupil_hybrid_candidates(C, maskUNet, opts)
% All pupil candidates of one ROI crop with their features (for labelling and for choosing the scoring rule).
% Candidates: (a) every U-Net blob, (b) for every chain of nested dark regions (threshold ladder) the innermost region
% with an edge, grown outward while it stays round and does not jump in area (as pupil_hybrid_mask).
% Filters here are LOOSE on purpose (the final rule is chosen from labelled frames).
%   cand(i): .m (mask) .src ('unet'|'dark') .area .diam .cx .cy .circ .sol .ecc .con (ring median contrast)
%            .sec (8 sector contrasts, sorted ascending) .dark (median N inside) .stab .dglint (glint distance / radius)
%            .dctr (distance to crop centre / crop diagonal) .t (threshold index) .unet_ov (fraction inside U-Net mask)
if nargin < 3, opts = struct(); end
o = defaults(opts);
C = double(C); [H, W] = size(C);
G = C >= prctile(C(:), o.glint_pct);
G = imdilate(G, strel('disk', o.glint_dilate));
Cf = C; if any(G(:)) && nnz(G) < 0.3 * numel(G), Cf = regionfill(C, G); end
Cb = imgaussfilt(Cf, o.smooth_sigma);
lo = prctile(Cb(:), 1); hi = prctile(Cb(:), 99);
N = min(max((Cb - lo) / max(1, hi - lo), 0), 1);
minA = pi * (o.min_diam / 2) ^ 2; maxD = o.max_diam_frac * min(H, W);
GL = regionprops(imdilate(C >= prctile(C(:), 99.5), strel('disk', 1)), 'Centroid', 'Area');
gl = zeros(0, 2); if ~isempty(GL), gl = cat(1, GL([GL.Area] >= 3).Centroid); end
U = false(H, W); if ~isempty(maskUNet), U = logical(maskUNet); end
cand = struct('m', {}, 'src', {}, 'area', {}, 'diam', {}, 'cx', {}, 'cy', {}, 'circ', {}, 'sol', {}, 'ecc', {}, 'con', {}, ...
    'sec', {}, 'dark', {}, 'stab', {}, 'dglint', {}, 'dctr', {}, 't', {}, 'unet_ov', {});
% ---- (a) U-Net blobs
if any(U(:))
    CC = bwconncomp(U);
    for i = 1:CC.NumObjects
        m = false(H, W); m(CC.PixelIdxList{i}) = true;
        if nnz(m) < minA, continue; end
        q = feats(m, N, G, gl, H, W); if q.diam > maxD, continue; end
        q.src = 'unet'; q.stab = NaN; q.t = NaN; q.unet_ov = 1;
        cand(end+1) = orderfields(q, cand); %#ok<AGROW>
    end
end
% ---- (b) dark chains
qs = prctile(N(:), o.ladder); nT = numel(qs);
reg = struct('m', {}, 't', {}, 'area', {}, 'diam', {}, 'circ', {}, 'sol', {}, 'con', {});
for t = 1:nT
    B = imfill(imopen(N <= qs(t), strel('disk', 1)), 'holes');
    CC = bwconncomp(B);
    R = regionprops(CC, 'Area', 'Perimeter', 'Solidity', 'EquivDiameter', 'BoundingBox');
    for i = 1:numel(R)
        if R(i).Area < minA || R(i).EquivDiameter > maxD, continue; end
        bb = R(i).BoundingBox;
        if bb(1) <= 1.5 || bb(2) <= 1.5 || bb(1) + bb(3) >= W - 0.5 || bb(2) + bb(4) >= H - 0.5, continue; end
        m = false(H, W); m(CC.PixelIdxList{i}) = true;
        circ = min(1, 4 * pi * R(i).Area / max(1, R(i).Perimeter) ^ 2);
        reg(end+1) = struct('m', m, 't', t, 'area', R(i).Area, 'diam', R(i).EquivDiameter, 'circ', circ, ...
            'sol', R(i).Solidity, 'con', ring_med(N, m, R(i).EquivDiameter, G)); %#ok<AGROW>
    end
end
nR = numel(reg); parent = zeros(1, nR);
for i = 1:nR
    for j = 1:nR
        if reg(j).t == reg(i).t + 1 && nnz(reg(i).m & reg(j).m) >= 0.9 * reg(i).area, parent(i) = j; break, end
    end
end
isLeaf = true(1, nR); isLeaf(parent(parent > 0)) = false;
for i = find(isLeaf)
    path = i; k = i;
    while parent(k) > 0, k = parent(k); path(end+1) = k; end %#ok<AGROW>
    sel = 0;
    for p = 1:numel(path)
        q = reg(path(p));
        good = q.circ >= o.min_circ && q.sol >= o.min_solidity && q.con >= o.min_contrast;
        if sel == 0
            if good, sel = p; end
        else
            prev = reg(path(sel));
            if good && q.area <= o.max_growth * prev.area && q.con >= 0.6 * prev.con, sel = p; else, break; end
        end
    end
    if sel == 0, continue; end
    m = imclose(reg(path(sel)).m, strel('disk', 1));
    % duplicate of a candidate already in the list (chains that merged)?
    dup = false;
    for c = 1:numel(cand)
        if strcmp(cand(c).src, 'dark') && nnz(cand(c).m & m) / nnz(cand(c).m | m) > 0.8, dup = true; break, end
    end
    if dup, continue; end
    q = feats(m, N, G, gl, H, W);
    s = 1; d0 = reg(path(sel)).diam;
    for p = sel - 1:-1:1, if reg(path(p)).diam >= 0.7 * d0, s = s + 1; else, break; end, end
    for p = sel + 1:numel(path), if reg(path(p)).diam <= 1.3 * d0, s = s + 1; else, break; end, end
    q.src = 'dark'; q.stab = s; q.t = reg(path(sel)).t; q.unet_ov = nnz(m & U) / max(1, nnz(m));
    cand(end+1) = orderfields(q, cand); %#ok<AGROW>
end
end

function q = feats(m, N, G, gl, H, W)
R = regionprops(m, 'Area', 'Centroid', 'Perimeter', 'Solidity', 'EquivDiameter', 'Eccentricity');
[~, i] = max([R.Area]); R = R(i);
q = struct('m', m, 'src', '', 'area', R.Area, 'diam', R.EquivDiameter, 'cx', R.Centroid(1), 'cy', R.Centroid(2), ...
    'circ', min(1, 4 * pi * R.Area / max(1, R.Perimeter) ^ 2), 'sol', R.Solidity, 'ecc', R.Eccentricity, ...
    'con', ring_med(N, m, R.EquivDiameter, G), 'sec', sector_con(N, m, R, G), 'dark', median(N(m & ~G)), ...
    'stab', NaN, 'dglint', NaN, 'dctr', hypot(R.Centroid(1) - W / 2, R.Centroid(2) - H / 2) / hypot(W, H), 't', NaN, 'unet_ov', NaN);
if isempty(q.dark), q.dark = median(N(m)); end
if ~isempty(gl), q.dglint = min(hypot(gl(:,1) - R.Centroid(1), gl(:,2) - R.Centroid(2))) / max(R.EquivDiameter / 2, 4); end
end

function c = ring_med(N, m, d, G)
r = max(2, round(0.3 * d));
ring = imdilate(m, strel('disk', r)) & ~m & ~G; inside = m & ~G;
if nnz(ring) < 5 || nnz(inside) < 5, c = 0; return; end
c = median(N(ring)) - median(N(inside));
end

function s = sector_con(N, m, R, G)
% contrast in 8 angular sectors of the ring around the region (sorted ascending): low values = a side of the region
% that is NOT bordered by brighter tissue (eyelid margin, eye corner, another dark area)
d = R.EquivDiameter; r = max(2, round(0.3 * d));
ring = imdilate(m, strel('disk', r)) & ~m & ~G; inside = m & ~G;
s = zeros(1, 8);
if nnz(ring) < 8 || nnz(inside) < 5, return; end
[yy, xx] = find(ring); ang = atan2(yy - R.Centroid(2), xx - R.Centroid(1));
sec = min(8, floor((ang + pi) / (2 * pi) * 8) + 1);
vi = median(N(inside)); vr = N(ring);
for k = 1:8
    v = vr(sec == k);
    if isempty(v), s(k) = 0; else, s(k) = median(v) - vi; end
end
s = sort(s);
end

function o = defaults(o)
d = struct('glint_pct', 98, 'glint_dilate', 2, 'smooth_sigma', 1.0, 'min_diam', 5, 'max_diam_frac', 0.45, ...
    'ladder', [1 2 3 4 5 7 9 12 15 18 22 26 30 35], 'min_circ', 0.5, 'min_solidity', 0.75, 'min_contrast', 0.08, ...
    'max_growth', 1.6);
f = fieldnames(d);
for i = 1:numel(f), if ~isfield(o, f{i}), o.(f{i}) = d.(f{i}); end, end
end
