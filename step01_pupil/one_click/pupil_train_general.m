function pupil_train_general(pairCsv, modelOut, varargin)
% pupil_train_general  Train ONE U-Net for all mice / days from the per-video pseudo-labels (2026-11-12).
%
%   pupil_train_general(pairCsv, modelOut)                 e.g. train_pairs.csv -> general_pupil_unet_holdout3.mat
%   pupil_train_general(..., 'Epochs', 12, 'TestCsv', test_pairs.csv)
%
% Images in the per-video runs are ROI crops already preprocessed as in step1 (glint suppression + CLAHE), masks 0/255.
% Every pair is resized to INPUT x INPUT (default 128) - ROIs are 95-150 px - so one network fits every ROI;
% at tracking time the mask is resized back to the ROI, so diameters stay in ROI pixels as before.
% Augmentation: random rotation / scale / shift, brightness-contrast-gamma, blur, noise (labels warped with nearest).
% Saves  net, INPUT, preprocessing settings and the training log in modelOut (.mat); if TestCsv is given, also the Dice per
% held-out pair and per run (compare with that run's own per-video validation Dice).

p = inputParser;
p.addParameter('Epochs', 12);
p.addParameter('Input', 128);
p.addParameter('MiniBatch', 32);
p.addParameter('LR', 1e-3);
p.addParameter('TestCsv', '');
p.addParameter('ValFrac', 0.05);
p.parse(varargin{:});
o = p.Results; S = o.Input;
rng(7);

T = readtable(pairCsv, 'TextType', 'string', 'Delimiter', ',');
fprintf('[train] %d pairs from %s\n', height(T), pairCsv);
[X, Y] = load_pairs(T, S);

% validation = 5 % of the training RUNS (whole runs, so neighbouring frames do not leak)
runs = unique(T.run); nv = max(2, round(o.ValFrac*numel(runs)));
vr = runs(randperm(numel(runs), nv)); isV = ismember(T.run, vr);
fprintf('[train] train %d pairs / val %d pairs (%d runs)\n', nnz(~isV), nnz(isV), nv);

classes = ["background","pupil"];
dsT = transform(combine(arrayDatastore(X(:,:,:,~isV),'IterationDimension',4), arrayDatastore(Y(:,:,:,~isV),'IterationDimension',4)), @(d) aug(d, classes, true));
dsV = transform(combine(arrayDatastore(X(:,:,:,isV),'IterationDimension',4),  arrayDatastore(Y(:,:,:,isV),'IterationDimension',4)),  @(d) aug(d, classes, false));

lg = unetLayers([S S 1], 2, 'EncoderDepth', 4, 'NumFirstEncoderFilters', 32);
frac = mean(Y(:)>0); w = 1./sqrt([1-frac frac]); w = w/mean(w);
pxName = lg.Layers(end).Name;
lg = replaceLayer(lg, pxName, dicePixelClassificationLayer('Name','dice','Classes',classes));
fprintf('[train] pupil pixel fraction %.3f (class weights would be %.2f / %.2f); loss = generalized Dice\n', frac, w);

itPerEpoch = floor(nnz(~isV)/o.MiniBatch);
opts = trainingOptions('adam', 'InitialLearnRate', o.LR, 'MaxEpochs', o.Epochs, 'MiniBatchSize', o.MiniBatch, ...
    'Shuffle','every-epoch', 'ValidationData', dsV, 'ValidationFrequency', max(50, round(itPerEpoch/2)), ...
    'LearnRateSchedule','piecewise', 'LearnRateDropPeriod', max(1, round(o.Epochs/3)), 'LearnRateDropFactor', 0.3, ...
    'L2Regularization', 1e-4, 'ExecutionEnvironment','auto', 'Verbose', true, 'VerboseFrequency', 100, 'Plots','none', ...
    'OutputNetwork','best-validation-loss');
t0 = tic;
[net, info] = trainNetwork(dsT, lg, opts);
fprintf('[train] done in %.1f min\n', toc(t0)/60);

INPUT = S; PREPROC = struct('glint_pct', 98, 'glint_dilate', 2, 'clahe', true); %#ok<NASGU>
TRAIN = struct('pairCsv', pairCsv, 'nPairs', height(T), 'epochs', o.Epochs, 'when', char(datetime('now')), ...
               'finalValLoss', info.FinalValidationLoss); %#ok<NASGU>
save(modelOut, 'net', 'INPUT', 'PREPROC', 'TRAIN', '-v7.3');
fprintf('[train] saved %s\n', modelOut);

if ~isempty(o.TestCsv)
    E = readtable(o.TestCsv, 'TextType', 'string', 'Delimiter', ',');
    [Xe, Ye] = load_pairs(E, S);
    Ypred = semanticseg(Xe, net, 'OutputType','uint8', 'MiniBatchSize', 128) == 2;   % H x W x N, class 2 = pupil
    d = zeros(height(E),1);
    for i = 1:height(E)
        a = Ypred(:,:,i); b = squeeze(Ye(:,:,1,i)) > 0;   % compare at INPUT size
        d(i) = 2*nnz(a&b) / max(1, nnz(a)+nnz(b));
    end
    E.dice_general = d;
    [g, runName] = findgroups(E.run);
    perRun = table(runName, splitapply(@median, d, g), 'VariableNames', {'run','median_dice_general'});
    [~, base] = fileparts(modelOut);
    writetable(E, fullfile(fileparts(modelOut), [base '_test_pairs_dice.csv']));
    writetable(perRun, fullfile(fileparts(modelOut), [base '_test_runs_dice.csv']));
    fprintf('[test] held-out pairs: median Dice %.3f (IQR %.3f-%.3f), runs: median of run medians %.3f\n', ...
        median(d), prctile(d,25), prctile(d,75), median(perRun.median_dice_general));
end
end

function [X, Y] = load_pairs(T, S)
n = height(T); X = zeros(S, S, 1, n, 'uint8'); Y = zeros(S, S, 1, n, 'uint8');
for i = 1:n
    I = imread(T.image(i)); if size(I,3) > 1, I = rgb2gray(I); end
    L = imread(T.label(i)); if size(L,3) > 1, L = rgb2gray(L); end
    X(:,:,1,i) = imresize(I, [S S], 'bilinear');
    Y(:,:,1,i) = uint8(imresize(L > 0, [S S], 'nearest'));
    if mod(i, 5000) == 0, fprintf('  loaded %d / %d\n', i, n); end
end
end

function out = aug(d, classes, doAug)
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
