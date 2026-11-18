function teach_finetune(setDir, baseModel, pkgPupilDir, outDir, varargin)
% teach_finetune  The NEW pupil model of this cohort (2026-11-17): fine-tune the package's U-Net on the teach set made
% by tools/pupil_teach_dataset.py = the user's own drawings (hard frames) + the package's reliable normal frames.
% Same preprocessing (pupil_mask_frame: glint, CLAHE, 128x128), loss, augmentation and options as finetune_unet_cv.m.
%   'Fold', 0  -> train on every session (the final model)
%   'Fold', f  -> hold out the mice of FOLDS{f} (both days) and evaluate on their samples with the package's mask
%                 post-processing and its 'maxp' choice (the blob with the highest pupil probability, >= 0.9):
%                 user drawings: correct if IoU >= 0.3 or the blob centre lies inside the drawn pupil
%                 package frames: agreement (IoU) with the package's own pupil
% 'UserRepeat' copies of every user sample go into training (each copy augmented differently).
p = inputParser;
p.addParameter('Fold', 0); p.addParameter('Epochs', 30); p.addParameter('LR', 1e-4); p.addParameter('MiniBatch', 16);
p.addParameter('UserRepeat', 6); p.addParameter('Tag', 'teach'); p.addParameter('Threads', 0);
p.parse(varargin{:}); o = p.Results;
if o.Threads > 0, maxNumCompThreads(o.Threads); end
addpath(pkgPupilDir); if ~exist(outDir, 'dir'), mkdir(outDir); end
rng(23 + o.Fold);
M = load(baseModel); S = M.INPUT;
[X, Y, T, Tm] = build(setDir, M, S);
FOLDS = {{'9669m_green1'}, {'9669m_green2'}, {'9669m_white', '9670f_black'}, {'9670f_green1', '9669m_black'}, ...
    {'9670f_green2', '9670f_white'}};   % fold 5 added 2026-11-17 when the 9670f drawings came in
if o.Fold > 0, test = ismember(T.mouse, FOLDS{o.Fold}); else, test = false(height(T), 1); end
tr = find(~test);
isUser = T.kind == "user" | T.kind == "user_pos" | T.kind == "none";
tr = [tr; repmat(tr(isUser(tr)), o.UserRepeat - 1, 1)];
fprintf('[%s fold %d] train %d samples (%d user x%d, %d package), test %d (%d user)\n', o.Tag, o.Fold, numel(tr), ...
    nnz(isUser(~test)), o.UserRepeat, nnz(~isUser & ~test), nnz(test), nnz(test & isUser));
classes = ["background", "pupil"];
dsT = transform(combine(arrayDatastore(X(:,:,:,tr), 'IterationDimension', 4), ...
    arrayDatastore(Y(:,:,:,tr), 'IterationDimension', 4)), @(d) aug(d, classes, true));
lg = layerGraph(M.net);
opts = trainingOptions('adam', 'InitialLearnRate', o.LR, 'MaxEpochs', o.Epochs, 'MiniBatchSize', o.MiniBatch, ...
    'Shuffle', 'every-epoch', 'L2Regularization', 1e-4, 'ExecutionEnvironment', 'cpu', 'Verbose', true, ...
    'VerboseFrequency', 50, 'Plots', 'none');
t0 = tic; net = trainNetwork(dsT, lg, opts); fprintf('[%s fold %d] trained in %.1f min\n', o.Tag, o.Fold, toc(t0) / 60);
INPUT = S; PREPROC = M.PREPROC; %#ok<NASGU>
TRAIN = struct('base', baseModel, 'set', setDir, 'fold', o.Fold, 'heldOut', {{}}, 'epochs', o.Epochs, 'lr', o.LR, ...
    'userRepeat', o.UserRepeat, 'nTrain', numel(tr), 'nUser', nnz(isUser(~test)), 'when', char(datetime('now')), ...
    'note', 'teach set 2026-11-17: user drawings (hard frames) + package reliable normal frames'); %#ok<NASGU>
if o.Fold > 0, TRAIN.heldOut = FOLDS{o.Fold}; end
save(fullfile(outDir, sprintf('%s_fold%d.mat', o.Tag, o.Fold)), 'net', 'INPUT', 'PREPROC', 'TRAIN', '-v7.3');
if o.Fold == 0, return; end
% ---- held-out evaluation: new vs base model, package maxp choice
idx = find(test); R = table();
Pn = predict(net, X(:,:,:,idx), 'MiniBatchSize', 32); Po = predict(M.net, X(:,:,:,idx), 'MiniBatchSize', 32);
for q = 1:numel(idx)
    i = idx(q); tgt = Tm{i};
    [okN, dN, iN] = judge(post_maxp(Pn(:,:,2,q), T.cropH(i), T.cropW(i)), tgt, T.kind(i));
    [okO, dO, iO] = judge(post_maxp(Po(:,:,2,q), T.cropH(i), T.cropW(i)), tgt, T.kind(i));
    R = [R; table(T.session(i), T.frame(i), T.kind(i), okN, okO, dN, dO, iN, iO, 'VariableNames', ...
        {'session', 'frame', 'kind', 'new', 'base', 'new_diam_ratio', 'base_diam_ratio', 'new_iou', 'base_iou'})]; %#ok<AGROW>
end
writetable(R, fullfile(outDir, sprintf('%s_eval_fold%d.csv', o.Tag, o.Fold)));
for nm = ["new", "base"]
    u = R.(nm)(R.kind ~= "package"); pk = R.(nm)(R.kind == "package");
    fprintf('[%s fold %d] %-4s  user drawings: correct %d / wrong %d / miss %d | package frames: agree %d / not %d\n', ...
        o.Tag, o.Fold, nm, nnz(u == "correct"), nnz(u == "wrong"), nnz(u == "miss"), nnz(pk == "correct"), ...
        nnz(pk ~= "correct"));
end
end

% ======================================================================================================================
function [X, Y, T, Tm] = build(setDir, M, S)
meta = jsondecode(fileread(fullfile(setDir, 'meta.json')));
X = zeros(S, S, 1, 0, 'uint8'); Y = zeros(S, S, 1, 0, 'uint8'); T = table(); Tm = {};
for k = 1:numel(meta)
    key = string(meta(k).key); Sd = load(fullfile(setDir, [char(key) '.mat']));
    r0 = double(Sd.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
    F = arrayfun(@(j) Sd.samples(:, :, j), 1:size(Sd.samples, 3), 'UniformOutput', false)';
    [~, ~, crops, preps] = pupil_mask_frame(F, rr, M, 'IrisCorrection', true);
    kinds = string(Sd.kind); frames = double(Sd.frames(:));
    mouse = regexprep(key, '^day\d+_', '');
    for j = 1:numel(crops)
        [h, w] = size(crops{j}); m = logical(Sd.masks(:, :, j));
        if ~isequal(size(m), [h w])
            error('%s sample %d: mask %dx%d but the crop is %dx%d', key, j, size(m, 1), size(m, 2), h, w);
        end
        X(:, :, 1, end + 1) = imresize(preps{j}, [S S], 'bilinear'); %#ok<AGROW>
        Y(:, :, 1, end + 1) = uint8(imresize(m, [S S], 'nearest')); %#ok<AGROW>
        T = [T; table(key, mouse, frames(j), strtrim(kinds(j)), h, w, ...
            'VariableNames', {'session', 'mouse', 'frame', 'kind', 'cropH', 'cropW'})]; %#ok<AGROW>
        Tm{end + 1, 1} = m; %#ok<AGROW>
    end
end
end

function m = post_maxp(pm, h, w)
% pupil_mask_frame steps 4-5, then the package's maxp choice: the blob with the highest pupil probability (>= 0.9)
ms = pm > 0.5;
if nnz(ms) < 50, thr = max(0.40, graythresh(pm)); ms = pm >= thr; end
m = imresize(ms, [h w], 'nearest'); P = imresize(pm, [h w], 'bilinear');
m = imopen(m, strel('disk', 1)); m = imclose(m, strel('disk', 2)); m = imfill(m, 'holes'); m = bwareaopen(m, 50);
out = false(h, w);
if any(m(:))
    CC = bwconncomp(m); mp = cellfun(@(ix) max(P(ix)), CC.PixelIdxList);
    [best, b] = max(mp);
    if best >= 0.9, out(CC.PixelIdxList{b}) = true; end
end
m = out;
end

function [res, dr, iou] = judge(m, tgt, kind)
dr = NaN; iou = NaN;
if kind == "none"
    if any(m(:)), res = "false"; else, res = "none_ok"; end
    return
end
if ~any(m(:)), res = "miss"; return; end
iou = nnz(m & tgt) / max(1, nnz(m | tgt));
R = regionprops(m, 'Centroid'); c = round(R(1).Centroid);
inside = tgt(min(end, max(1, c(2))), min(end, max(1, c(1))));
thr = 0.3; if kind == "package", thr = 0.5; end
if iou >= thr || (inside && kind ~= "package"), res = "correct"; else, res = "wrong"; end
dr = sqrt(nnz(m) / max(1, nnz(tgt)));
end

function out = aug(d, classes, doAug)
% = pupil_train_general.m / finetune_unet_cv.m aug (copied unchanged)
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
