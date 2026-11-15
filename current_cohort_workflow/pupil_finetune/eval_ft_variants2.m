function eval_ft_variants2(diagDir, candDir, fitDir, modelFile, pkgPupilDir, ftDir, mode)
% Held-out check of the blob CHOICE when several U-Net blobs pass max p >= 0.9 (fine-tuned fold models).
% Variants: maxp (current), mass (sum of p), round (circularity), glint (closest to the brightest spot),
% tie_round (near-ties of max p broken by circularity), oracle (any confident blob correct = upper bound).
addpath(pkgPupilDir);
M = load(modelFile); S = M.INPUT;
[X, T, Tm, Cr] = build_dataset(diagDir, candDir, fitDir, M, S);
FOLDS = {{'9669m_black', '9670f_green1'}, {'9669m_green1', '9670f_white'}, {'9669m_green2', '9670f_black'}, {'9669m_white', '9670f_green2'}};
if strcmp(mode, 'folds'), names = arrayfun(@(f) sprintf('unet_ft_fold%d.mat', f), 1:4, 'UniformOutput', false);
else, names = {'unet_ft_test_day1.mat', 'unet_ft_test_day2.mat'}; end
V = {'maxp', 'mass', 'round', 'glint', 'tie_round', 'oracle'};
out = table();
for f = 1:numel(names)
    N = load(fullfile(ftDir, names{f}));
    if strcmp(mode, 'folds'), test = ismember(T.mouse, FOLDS{f}); else, test = startsWith(T.session, sprintf('day%d_', f)); end
    idx = find(test);
    P = predict(N.net, X(:,:,:,idx), 'MiniBatchSize', 32);
    for q = 1:numel(idx)
        i = idx(q); h = T.cropH(i); w = T.cropW(i);
        pm = imresize(P(:,:,2,q), [h w], 'bilinear');
        m0 = imresize(P(:,:,2,q) > 0.5, [h w], 'nearest');
        m0 = imopen(m0, strel('disk', 1)); m0 = imclose(m0, strel('disk', 2)); m0 = imfill(m0, 'holes'); m0 = bwareaopen(m0, 50);
        CC = bwconncomp(m0);
        mx = cellfun(@(ix) max(pm(ix)), CC.PixelIdxList); ok = mx >= 0.9;
        row = table(T.session(i), T.frame(i), T.label(i), nnz(ok), 'VariableNames', {'session', 'frame', 'label', 'n_conf'});
        g = glint_xy(Cr{i});
        for v = 1:numel(V)
            if strcmp(V{v}, 'oracle')
                res = "miss"; if ~any(ok) && T.label(i) == -1, res = "none_ok"; end
                cand = find(ok);
                if isempty(cand), res = judge(false(h, w), Tm{i}, T.label(i)); end
                for c = cand
                    mm = false(h, w); mm(CC.PixelIdxList{c}) = true; r = judge(mm, Tm{i}, T.label(i));
                    if r == "correct", res = r; break; else, res = r; end
                end
                row.(V{v}) = res;
            else
                b = choose(CC, pm, ok, mx, V{v}, g, h, w);
                mm = false(h, w); if b > 0, mm(CC.PixelIdxList{b}) = true; end
                row.(V{v}) = judge(mm, Tm{i}, T.label(i));
            end
        end
        out = [out; row]; %#ok<AGROW>
    end
end
writetable(out, fullfile(ftDir, sprintf('variants2_%s.csv', mode)));
fprintf('frames with >1 confident blob: %d of %d\n', nnz(out.n_conf > 1), height(out));
fprintf('%-10s %7s %6s %5s %6s %8s | recall  bad-pick rate\n', 'variant', 'correct', 'wrong', 'miss', 'false', 'none-ok');
for v = 1:numel(V)
    r = out.(V{v});
    nc = nnz(r == "correct"); nw = nnz(r == "wrong"); nm = nnz(r == "miss"); nf = nnz(r == "false"); no = nnz(r == "none_ok");
    fprintf('%-10s %7d %6d %5d %6d %8d | %5.1f%%  %5.1f%%\n', V{v}, nc, nw, nm, nf, no, 100 * nc / max(1, nnz(out.label >= 1)), ...
        100 * (nw + nf) / max(1, nc + nw + nf));
end
end

function b = choose(CC, pm, ok, mx, how, g, h, w)
b = 0; if ~any(ok), return; end
cand = find(ok);
switch how
    case 'maxp', [~, k] = max(mx(cand));
    case 'mass', s = cellfun(@(ix) sum(pm(ix)), CC.PixelIdxList(cand)); [~, k] = max(s);
    case 'round', cr = circ(CC, cand, h, w); [~, k] = max(cr);
    case 'glint'
        if any(isnan(g)), [~, k] = max(mx(cand));
        else
            d = zeros(size(cand));
            for j = 1:numel(cand), [yy, xx] = ind2sub([h w], CC.PixelIdxList{cand(j)}); d(j) = hypot(mean(xx) - g(1), mean(yy) - g(2)); end
            [~, k] = min(d);
        end
    case 'tie_round'
        near = mx(cand) >= max(mx(cand)) - 0.01; cr = circ(CC, cand, h, w); cr(~near) = -Inf; [~, k] = max(cr);
end
b = cand(k);
end

function cr = circ(CC, cand, h, w)
cr = zeros(size(cand));
for j = 1:numel(cand)
    mm = false(h, w); mm(CC.PixelIdxList{cand(j)}) = true;
    R = regionprops(mm, 'Area', 'Perimeter');
    cr(j) = 4 * pi * R(1).Area / max(1, R(1).Perimeter) ^ 2;
end
end

function g = glint_xy(C)
% centre of the largest very bright spot (>= 90 % of the crop maximum and >= the 99.5th percentile)
C = double(C); g = [NaN NaN];
thr = max(0.9 * max(C(:)), prctile(C(:), 99.5));
B = C >= thr; CCg = bwconncomp(B);
if CCg.NumObjects == 0, return; end
np = cellfun(@numel, CCg.PixelIdxList); [~, k] = max(np);
[yy, xx] = ind2sub(size(C), CCg.PixelIdxList{k}); g = [mean(xx) mean(yy)];
end

function res = judge(m, tgt, lab)
if lab >= 1
    if ~any(m(:)), res = "miss"; return; end
    iou = nnz(m & tgt) / max(1, nnz(m | tgt));
    R = regionprops(m, 'Centroid'); c = round(R(1).Centroid);
    inside = tgt(min(end, max(1, c(2))), min(end, max(1, c(1))));
    if iou >= 0.3 || inside, res = "correct"; else, res = "wrong"; end
elseif lab == -1
    if any(m(:)), res = "false"; else, res = "none_ok"; end
else
    if any(m(:)), res = "blob_unknown"; else, res = "none_unknown"; end
end
end

function [X, T, Tm, Cr] = build_dataset(diagDir, candDir, fitDir, M, S)
L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string');
meta = jsondecode(fileread(fullfile(diagDir, 'meta.json')));
X = zeros(S, S, 1, 0, 'uint8'); T = table(); Tm = {}; Cr = {};
for k = 1:numel(meta)
    key = string(meta(k).key); Sd = load(fullfile(diagDir, [char(key) '.mat'])); Xs = Sd.samples;
    r0 = double(Sd.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
    F = arrayfun(@(j) Xs(:, :, j), 1:size(Xs, 3), 'UniformOutput', false)';
    [~, ~, crops, preps] = pupil_mask_frame(F, rr, M, 'IrisCorrection', true);
    C = load(fullfile(candDir, [char(key) '_cands.mat'])); allc = C.allc;
    Lk = L(L.session == key, :); mouse = regexprep(key, '^day\d+_', '');
    for j = 1:numel(crops)
        row = Lk(Lk.frame == j - 1, :);
        if isempty(row) || isnan(row.label(1)), continue; end
        lab = row.label(1); [h, w] = size(crops{j});
        if lab >= 1, m = logical(allc{j}(lab).m); elseif lab == -1, m = false(h, w); else, m = []; end
        X(:, :, 1, end + 1) = imresize(preps{j}, [S S], 'bilinear'); %#ok<AGROW>
        T = [T; table(key, mouse, j - 1, lab, h, w, 'VariableNames', {'session', 'mouse', 'frame', 'label', 'cropH', 'cropW'})]; %#ok<AGROW>
        Tm{end + 1, 1} = m; %#ok<AGROW>
        Cr{end + 1, 1} = crops{j}; %#ok<AGROW>
    end
end
end
