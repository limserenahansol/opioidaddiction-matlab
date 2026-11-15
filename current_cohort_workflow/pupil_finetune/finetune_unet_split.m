function finetune_unet_split(diagDir, candDir, fitDir, modelFile, pkgPupilDir, outDir, varargin)
% Same as finetune_unet_cv, but the held-out set is chosen by 'Split': 'day1' = train on day 2, test on day 1;
% 'day2' = train on day 1, test on day 2 (same mice, other day = how the model is used on days 3-18).
% Fine-tune the package's general pupil U-Net on labelled current-rig frames (same preprocessing, input size, loss and
% augmentation as pupil_train_general.m). Leave-mice-out: fold f holds out 2 mice (both days); 'Fold', 0 = train on all
% mice (final model). Evaluation per held-out frame with the package's mask post-processing (pupil_mask_frame steps
% 4-5) and the largest blob:
%   pupil labelled on a contour -> correct if IoU >= 0.3 with that contour or the blob centroid lies inside it
%   pupil not visible (-1)      -> correct if no blob
%   pupil visible, no contour (0) -> reported separately (location unknown)
p = inputParser;
p.addParameter('Fold', 1); p.addParameter('Split', 'day2'); p.addParameter('Epochs', 25); p.addParameter('LR', 1e-4); p.addParameter('MiniBatch', 16);
p.addParameter('Threads', 0); p.parse(varargin{:}); o = p.Results;
if o.Threads > 0, maxNumCompThreads(o.Threads); end
addpath(pkgPupilDir); if ~exist(outDir, 'dir'), mkdir(outDir); end
rng(21);
M = load(modelFile); S = M.INPUT;
[X, Y, T, Tm] = build_dataset(diagDir, candDir, fitDir, M, S);
FOLDS = {{'9669m_black', '9670f_green1'}, {'9669m_green1', '9670f_white'}, {'9669m_green2', '9670f_black'}, {'9669m_white', '9670f_green2'}};
test = startsWith(T.session, string(o.Split) + "_"); o.Fold = 9; %#ok<NASGU>
useTr = ~test & (T.kind == "good" | T.kind == "none");
fprintf('[fold %d] train %d frames (%d with pupil), test %d frames\n', o.Fold, nnz(useTr), nnz(useTr & T.kind == "good"), nnz(test));
classes = ["background", "pupil"];
dsT = transform(combine(arrayDatastore(X(:,:,:,useTr), 'IterationDimension', 4), ...
    arrayDatastore(Y(:,:,:,useTr), 'IterationDimension', 4)), @(d) aug(d, classes, true));
lg = layerGraph(M.net);
opts = trainingOptions('adam', 'InitialLearnRate', o.LR, 'MaxEpochs', o.Epochs, 'MiniBatchSize', o.MiniBatch, ...
    'Shuffle', 'every-epoch', 'L2Regularization', 1e-4, 'ExecutionEnvironment', 'cpu', 'Verbose', true, ...
    'VerboseFrequency', 20, 'Plots', 'none');
t0 = tic; net = trainNetwork(dsT, lg, opts); fprintf('[fold %d] trained in %.1f min\n', o.Fold, toc(t0) / 60);
INPUT = S; PREPROC = M.PREPROC; %#ok<NASGU>
TRAIN = struct('base', modelFile, 'fold', o.Fold, 'heldOut', {{}}, 'epochs', o.Epochs, 'lr', o.LR, ...
    'nTrain', nnz(useTr), 'when', char(datetime('now'))); %#ok<NASGU>
TRAIN.heldOut = {o.Split};
save(fullfile(outDir, sprintf('unet_ft_test_%s.mat', o.Split)), 'net', 'INPUT', 'PREPROC', 'TRAIN', '-v7.3');

% ---- evaluation on the held-out mice: fine-tuned vs original U-Net
idx = find(test); R = table();
Pn = predict(net, X(:,:,:,idx), 'MiniBatchSize', 32); Po = predict(M.net, X(:,:,:,idx), 'MiniBatchSize', 32);
for q = 1:numel(idx)
    i = idx(q); tgt = Tm{i};
    [okN, dN] = judge(post(Pn(:,:,2,q), T.cropH(i), T.cropW(i)), tgt, T.label(i));
    [okO, dO] = judge(post(Po(:,:,2,q), T.cropH(i), T.cropW(i)), tgt, T.label(i));
    R = [R; table(T.session(i), T.frame(i), T.label(i), T.kind(i), okN, okO, dN, dO, ...
        'VariableNames', {'session', 'frame', 'label', 'kind', 'ft', 'orig', 'ft_diam_ratio', 'orig_diam_ratio'})]; %#ok<AGROW>
end
writetable(R, fullfile(outDir, sprintf('eval_test_%s.csv', o.Split)));
for nm = ["ft", "orig"]
    v = R.(nm);
    fprintf('[fold %d] %-4s  pupil on contour: correct %d / wrong %d / miss %d | not visible: none-ok %d / false %d | visible-no-contour: blob %d / none %d\n', ...
        o.Fold, nm, nnz(v == "correct"), nnz(v == "wrong"), nnz(v == "miss"), nnz(v == "none_ok"), nnz(v == "false"), ...
        nnz(v == "blob_unknown"), nnz(v == "none_unknown"));
end
end

% ======================================================================================================================
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
    Lk = L(L.session == key, :);
    mouse = regexprep(key, '^day\d+_', '');
    for j = 1:numel(crops)
        row = Lk(Lk.frame == j - 1, :);
        if isempty(row) || isnan(row.label(1)), continue; end
        lab = row.label(1); fit = string(row.fit(1)); [h, w] = size(crops{j});
        if lab >= 1
            m = logical(allc{j}(lab).m); kind = fit; if kind ~= "good", kind = "contour_" + fit; end
        elseif lab == -1
            m = false(h, w); kind = "none";
        else
            m = []; kind = "nocontour";
        end
        X(:, :, 1, end + 1) = imresize(preps{j}, [S S], 'bilinear'); %#ok<AGROW>
        if isempty(m), Y(:, :, 1, end + 1) = 0; else, Y(:, :, 1, end + 1) = uint8(imresize(m, [S S], 'nearest')); end %#ok<AGROW>
        T = [T; table(key, mouse, j - 1, lab, kind, h, w, 'VariableNames', {'session', 'mouse', 'frame', 'label', 'kind', 'cropH', 'cropW'})]; %#ok<AGROW>
        Tm{end + 1, 1} = m; %#ok<AGROW>
    end
end
end

function m = post(pm, h, w)
% pupil_mask_frame steps 4-5, then the largest blob
ms = pm > 0.5;
if nnz(ms) < 50, thr = max(0.40, graythresh(pm)); ms = pm >= thr; end
m = imresize(ms, [h w], 'nearest');
m = imopen(m, strel('disk', 1)); m = imclose(m, strel('disk', 2)); m = imfill(m, 'holes'); m = bwareaopen(m, 50);
if any(m(:))
    CC = bwconncomp(m); [~, b] = max(cellfun(@numel, CC.PixelIdxList)); m = false(h, w); m(CC.PixelIdxList{b}) = true;
end
end

function [res, dr] = judge(m, tgt, lab)
dr = NaN;
if lab >= 1
    if ~any(m(:)), res = "miss"; return; end
    iou = nnz(m & tgt) / max(1, nnz(m | tgt));
    R = regionprops(m, 'Centroid'); c = round(R(1).Centroid);
    inside = tgt(min(end, max(1, c(2))), min(end, max(1, c(1))));
    if iou >= 0.3 || inside, res = "correct"; else, res = "wrong"; end
    dr = sqrt(nnz(m) / max(1, nnz(tgt)));
elseif lab == -1
    if any(m(:)), res = "false"; else, res = "none_ok"; end
else
    if any(m(:)), res = "blob_unknown"; else, res = "none_unknown"; end
end
end

function out = aug(d, classes, doAug)
% = pupil_train_general.m aug (copied unchanged)
I = d{1}; L = d{2};
if doAug
    tf = randomAffine2d('Rotation', [-15 15], 'Scale', [0.85 1.2], 'XTranslation', [-10 10], 'YTranslation', [-10 10]);
    rv = affineOutputView(size(I), tf, 'BoundsStyle', 'centerOutput');
    I = imwarp(I, tf, 'OutputView', rv, 'FillValues', median(I(:)));
    L = imwarp(L, tf, 'nearest', 'OutputView', rv, 'FillValues', 0);
    J = im2double(I);
    J = J*(0.75 + 0.5*rand) + (rand-0.5)*0.2;          % contrast, brightness
    J = max(0, min(1, J)) .^ (0.7 + 0.6*rand);          % gamma
    if rand < 0.3, J = imgaussfilt(J, 0.5 + rand); end  % blur
    if rand < 0.3, J = J + 0.03*randn(size(J)); end     % noise
    I = im2uint8(max(0, min(1, J)));
end
out = {I, categorical(L, [0 1], classes)};
end
