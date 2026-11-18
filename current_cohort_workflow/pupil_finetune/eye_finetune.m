function eye_finetune(setDir, baseModel, outModel, varargin)
% eye_finetune  ONE-TIME whole-eye model for a mouse whose pupil cannot be seen (9670f black: eye damage since day
% 1; user 2026-11-17: "we can just use whole eye as alternate"). Fine-tunes the package's U-Net to segment the
% visible eye between the eyelids, from the user's outlines (tools/wholeeye_build_set.py -> setDir).
% Same network, loss, options and augmentation as finetune_unet_cv.m / teach_finetune.m; input = eye_prep (no glint
% erasing). Class 2 keeps the base network's name "pupil" but means "eye" here.
p = inputParser;
p.addParameter('Epochs', 30); p.addParameter('LR', 2e-4); p.addParameter('MiniBatch', 16);
p.addParameter('Repeat', 40); p.addParameter('Threads', 6);
p.addParameter('Note', '9670f black whole eye, user outlines of 2026-11-17');   % 2026-11-18: saved in TRAIN.note
p.parse(varargin{:}); o = p.Results;
if o.Threads > 0, maxNumCompThreads(o.Threads); end
rng(17);
M = load(baseModel); S = M.INPUT;
meta = jsondecode(fileread(fullfile(setDir, 'meta.json')));
n = numel(meta);
X = zeros(S, S, 1, n, 'uint8'); Y = zeros(S, S, 1, n, 'uint8');
for k = 1:n
    I = imread(fullfile(setDir, meta(k).img));
    Mk = imread(fullfile(setDir, meta(k).mask)) > 0;
    X(:, :, 1, k) = imresize(eye_prep(I), [S S], 'bilinear');
    Y(:, :, 1, k) = uint8(imresize(Mk, [S S], 'nearest'));
end
idx = repmat((1:n)', o.Repeat, 1);
classes = ["background", "pupil"];
ds = transform(combine(arrayDatastore(X(:, :, :, idx), 'IterationDimension', 4), ...
    arrayDatastore(Y(:, :, :, idx), 'IterationDimension', 4)), @(d) aug(d, classes));
opts = trainingOptions('adam', 'InitialLearnRate', o.LR, 'MaxEpochs', o.Epochs, 'MiniBatchSize', o.MiniBatch, ...
    'Shuffle', 'every-epoch', 'L2Regularization', 1e-4, 'ExecutionEnvironment', 'cpu', 'Verbose', true, ...
    'VerboseFrequency', 25, 'Plots', 'none');
fprintf('[eye] %d outlines x%d, %d epochs\n', n, o.Repeat, o.Epochs);
t0 = tic; net = trainNetwork(ds, layerGraph(M.net), opts);
fprintf('[eye] trained in %.1f min\n', toc(t0) / 60);
INPUT = S; PREP = 'eye_prep.m: adapthisteq of the grey eye crop (no glint erasing), then imresize to INPUT x INPUT'; %#ok<NASGU>
TRAIN = struct('base', baseModel, 'set', setDir, 'n', n, 'repeat', o.Repeat, 'epochs', o.Epochs, 'lr', o.LR, ...
    'when', char(datetime('now')), 'note', o.Note); %#ok<NASGU>
if ~exist(fileparts(outModel), 'dir'), mkdir(fileparts(outModel)); end
save(outModel, 'net', 'INPUT', 'PREP', 'TRAIN', '-v7.3');
Pm = predict(net, X, 'MiniBatchSize', 16);
for k = 1:n   % how well it fits the outlines it learned from (not a test)
    m = Pm(:, :, 2, k) > 0.5; t = Y(:, :, 1, k) > 0;
    fprintf('[eye] outline %d (day %d frame %d): IoU %.2f\n', k, meta(k).day, meta(k).frame, nnz(m & t) / max(1, nnz(m | t)));
end
end

function out = aug(d, classes)
% = pupil_train_general.m / finetune_unet_cv.m / teach_finetune.m augmentation (copied unchanged)
I = d{1}; L = d{2};
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
out = {I, categorical(L, [0 1], classes)};
end
