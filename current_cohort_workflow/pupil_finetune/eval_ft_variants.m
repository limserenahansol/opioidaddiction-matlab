function eval_ft_variants(diagDir, candDir, fitDir, modelFile, pkgPupilDir, ftDir, mode)
% Held-out evaluation of saved fine-tuned U-Nets with several mask post-processing variants.
% mode 'folds': unet_ft_fold1..4.mat (leave-2-mice-out); 'split': unet_ft_test_day1/day2.mat (leave-day-out).
addpath(pkgPupilDir);
M = load(modelFile); S = M.INPUT;
[X, ~, T, Tm] = build_dataset(diagDir, candDir, fitDir, M, S);
FOLDS = {{'9669m_black', '9670f_green1'}, {'9669m_green1', '9670f_white'}, {'9669m_green2', '9670f_black'}, {'9669m_white', '9670f_green2'}};
if strcmp(mode, 'folds'), names = arrayfun(@(f) sprintf('unet_ft_fold%d.mat', f), 1:4, 'UniformOutput', false);
else, names = {'unet_ft_test_day1.mat', 'unet_ft_test_day2.mat'}; end
V = {'pkg_fallback', 'no_fallback', 'maxp_0.9', 'meanp_0.75', 'meanp_0.85'};
out = table();
for f = 1:numel(names)
    if ~isfile(fullfile(ftDir, names{f})), fprintf('missing %s\n', names{f}); continue; end
    N = load(fullfile(ftDir, names{f}));
    if strcmp(mode, 'folds'), test = ismember(T.mouse, FOLDS{f}); else, test = startsWith(T.session, sprintf('day%d_', f)); end
    idx = find(test);
    P = predict(N.net, X(:,:,:,idx), 'MiniBatchSize', 32);
    for q = 1:numel(idx)
        i = idx(q); pm = imresize(P(:,:,2,q), [T.cropH(i) T.cropW(i)], 'bilinear');
        row = table(T.session(i), T.frame(i), T.label(i), 'VariableNames', {'session', 'frame', 'label'});
        for v = 1:numel(V)
            row.(matlab.lang.makeValidName(V{v})) = judge(post(P(:,:,2,q), pm, T.cropH(i), T.cropW(i), V{v}), Tm{i}, T.label(i));
        end
        out = [out; row]; %#ok<AGROW>
    end
end
writetable(out, fullfile(ftDir, sprintf('variants_%s.csv', mode)));
fprintf('%-14s %7s %6s %5s %6s %8s | recall  bad-pick rate\n', 'variant', 'correct', 'wrong', 'miss', 'false', 'none-ok');
for v = 1:numel(V)
    r = out.(matlab.lang.makeValidName(V{v}));
    nc = nnz(r == "correct"); nw = nnz(r == "wrong"); nm = nnz(r == "miss"); nf = nnz(r == "false"); no = nnz(r == "none_ok");
    fprintf('%-14s %7d %6d %5d %6d %8d | %5.1f%%  %5.1f%%\n', V{v}, nc, nw, nm, nf, no, 100 * nc / max(1, nnz(out.label >= 1)), ...
        100 * (nw + nf) / max(1, nc + nw + nf));
end
end

function m = post(p128, pm, h, w, variant)
switch variant
    case 'pkg_fallback'
        ms = p128 > 0.5;
        if nnz(ms) < 50, thr = max(0.40, graythresh(p128)); ms = p128 >= thr; end
        m = imresize(ms, [h w], 'nearest');
    otherwise
        m = imresize(p128 > 0.5, [h w], 'nearest');
end
m = imopen(m, strel('disk', 1)); m = imclose(m, strel('disk', 2)); m = imfill(m, 'holes'); m = bwareaopen(m, 50);
if ~any(m(:)), return; end
CC = bwconncomp(m); np = cellfun(@numel, CC.PixelIdxList);
mx = cellfun(@(ix) max(pm(ix)), CC.PixelIdxList); mn = cellfun(@(ix) mean(pm(ix)), CC.PixelIdxList);
switch variant
    case 'maxp_0.9', ok = mx >= 0.9; [~, b] = max(mx .* ok);
    case 'meanp_0.75', ok = mn >= 0.75; [~, b] = max(np .* ok);
    case 'meanp_0.85', ok = mn >= 0.85; [~, b] = max(np .* ok);
    otherwise, ok = true(size(np)); [~, b] = max(np);
end
m = false(h, w);
if any(ok), m(CC.PixelIdxList{b}) = true; end
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

function [X, Y, T, Tm] = build_dataset(diagDir, candDir, fitDir, M, S)
L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string');
meta = jsondecode(fileread(fullfile(diagDir, 'meta.json')));
X = zeros(S, S, 1, 0, 'uint8'); Y = zeros(S, S, 1, 0, 'uint8'); T = table(); Tm = {};
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
        lab = row.label(1); fit = string(row.fit(1)); [h, w] = size(crops{j});
        if lab >= 1, m = logical(allc{j}(lab).m); kind = fit; if kind ~= "good", kind = "contour_" + fit; end
        elseif lab == -1, m = false(h, w); kind = "none";
        else, m = []; kind = "nocontour"; end
        X(:, :, 1, end + 1) = imresize(preps{j}, [S S], 'bilinear'); %#ok<AGROW>
        if isempty(m), Y(:, :, 1, end + 1) = 0; else, Y(:, :, 1, end + 1) = uint8(imresize(m, [S S], 'nearest')); end %#ok<AGROW>
        T = [T; table(key, mouse, j - 1, lab, kind, h, w, 'VariableNames', {'session', 'mouse', 'frame', 'label', 'kind', 'cropH', 'cropW'})]; %#ok<AGROW>
        Tm{end + 1, 1} = m; %#ok<AGROW>
    end
end
end
