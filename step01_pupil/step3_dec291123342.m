% ============================================================
% STEP-3 (Final): Apply trained U-Net using Step-1/2 settings
% - Uses explicit run folder: run_20250807_131311
% - Matches Step-1 preprocessing & ROI
% - Robust model pick if Step-2 config missing
% - Writes visualization video + CSV (+ clean CSV)
% ============================================================
clear; clc; close all; rng default;
% BASE      = 'I:\addictionAug21\nov13day10\2705\';
% VIDEO     = fullfile(BASE,'blue','run_20251113_141204.avi'); 
% RUN_ROOT_PARENT = fullfile(BASE,'training_data1');                 % where run_* live
% RUN_FOLDER_NAME = 'run_20251113_132549';     



% ------------ USER ROOT ------------
% VIDEO = fullfile(BASE,'side','Fr1_6872_suc_d2_side.avi');
% VIDEO = fullfile(BASE,'side','Fr1_6872_suc_d2_side.avi');
BASE      = 'H:\addiction_dec29\day18\6099\';
VIDEO     = fullfile(BASE,'black','6099BBBB0039 26-01-22 11-01-38.avi');
   % <-- change per video
   % <-- change per video
   % <-- change per video
   % <-- change per video
RUN_ROOT_PARENT = fullfile(BASE,'training_data1');                 % where run_* live
RUN_FOLDER_NAME = 'run_20260122_112425';                               % <- your timestamp run
RUN_ROOT        = fullfile(RUN_ROOT_PARENT, RUN_FOLDER_NAME);
MODEL_DIR       = fullfile(BASE,'trainpupil');                     % where .mat models live
% ---- Per-video trackinH:\addiction_dec29\day3\6099\training_data1\run_20260101_121355H:\addiction_dec29\day3\6099\training_data1\run_20260101_121355H:\addiction_dec29\day3\6099\training_data1\run_20260101_121355g cleanup (no logic changes) ----
if ~exist('PROFILE','var'), PROFILE = 'default'; end

if strcmpi(PROFILE,'dilated_pupil')
    minBlobArea         = 600;   %run_20260106_101046run_20260106_101046 ignore tiny/glint blobs
    searchRadiusPx      = 45;    % larger search window
    maxDiamChangeFrac   = 0.35;  % allow bigger day-to-day changes
    fallback.minPix     = 800;   % don’t accept microscopic fallbacks
end

% Outputs
paths.out.dir   = fullfile(BASE,'pupil_outputs_6099black');
if ~exist(paths.out.dir,'dir'), mkdir(paths.out.dir); end
paths.out.video = fullfile(paths.out.dir,'pupil_tracking_FINAL_U-Net_visualization2.mp4');
paths.out.csv   = fullfile(paths.out.dir,'pupil_diameter_FINAL_U-Net_results2.csv');

fprintf('RUN_ROOT: %s\n', RUN_ROOT);
fprintf('MODEL_DIR: %s\n', MODEL_DIR);
fprintf('OUTPUT DIR: %s\n', paths.out.dir);

% ------------ Load Step-1 config (required) ------------
step1_cfg = fullfile(RUN_ROOT,'step1_config.mat');
assert(exist(step1_cfg,'file')==2, 'Missing %s', step1_cfg);
S1all = load(step1_cfg);             % contains S1, roi, RUN_ROOT, DIRS, counts
S1    = S1all.S1;

% ROI file from Step-1 (or fallback path)
% (Compat change) Prefer ROI saved inside step1_config.mat; else fall back to roi.mat path
if isfield(S1all,'ROI_FILE')
    ROI_FILE = S1all.ROI_FILE;
else
    ROI_FILE = fullfile(BASE,'newoutput','roi.mat');
end
fprintf('ROI_FILE: %s\n', ROI_FILE);

% ------------ Try to load Step-2 config (optional) ------------
step2_cfg = fullfile(RUN_ROOT,'step2_config.mat');
modelFile = '';
if exist(step2_cfg,'file')==2
    S2all = load(step2_cfg);                 % has stableModel / modelFileTS (if saved)
    if isfield(S2all,'stableModel') && exist(S2all.stableModel,'file')
        modelFile = S2all.stableModel;
    elseif isfield(S2all,'modelFileTS') && exist(S2all.modelFileTS,'file')
        modelFile = S2all.modelFileTS;
    end
else
    fprintf('NOTE: step2_config.mat not found in %s. Will auto-pick model from %s\n', RUN_ROOT, MODEL_DIR);
end

% (Compat change) If still empty, try the copy saved by Step-2 in RUN_ROOT
if isempty(modelFile)
    runCopy = fullfile(RUN_ROOT,'trainedPupilNet.mat');
    if exist(runCopy,'file')==2
        modelFile = runCopy;
    end
end

% If still empty, pick newest trainedPupilNet_*.mat in MODEL_DIR
if isempty(modelFile)
    m = dir(fullfile(MODEL_DIR,'trainedPupilNet_*.mat'));
    if ~isempty(m)
        [~,im] = max([m.datenum]);
        modelFile = fullfile(m(im).folder, m(im).name);
    else
        % final fallback: stable name
        modelFile = fullfile(MODEL_DIR,'trainedPupilNet_PseudoLabeled.mat');
    end
end
assert(exist(modelFile,'file')==2, 'Could not find any model .mat in %s', MODEL_DIR);
fprintf('Model: %s\n', modelFile);

% ------------ Match Step-1 preprocessing exactly ------------
PREPROC_MODE.clahe      = S1.USE_CLAHE;
PREPROC_MODE.stretchlim = S1.USE_STRETCHLIM;
PREPROC_MODE.none       = ~(PREPROC_MODE.clahe || PREPROC_MODE.stretchlim);
if PREPROC_MODE.stretchlim
    stretchTail = S1.ADAPTIVE_CLIP;  % same tail as Step-1
else
    stretchTail = 0.01;
end

% Tracking/cleanup params
fallback.minPix  = 50;
fallback.baseThr = 0.40;
fallback.useOtsu = true;

minBlobArea = 50;
SE_open  = strel('disk',1);
SE_close = strel('disk',2);

warmupFrames       = 30;
searchRadiusPx     = 30;
alphaCenter        = 0.18;
alphaDiam          = 0.12;
maxCenterJumpPx    = 18;
maxDiamChangeFrac  = 0.18;
holdLastIfRejected = true;

displayZoomInterp  = 'bilinear';
drawMaskOutline    = true;
reportEvery        = 300;

% ------------ Load model ------------
fprintf('Loading trained model...\n');
Snet = load(modelFile); net = [];
if isfield(Snet,'net'), net = Snet.net;
else
    fn = fieldnames(Snet);
    for k=1:numel(fn)
        obj = Snet.(fn{k});
        if isa(obj,'DAGNetwork') || isa(obj,'SeriesNetwork') || isa(obj,'dlnetwork')
            net = obj; break;
        end
    end
end
assert(~isempty(net),'No network object found in %s', modelFile);

% Input size and channels
inSize     = localGetInputSize(net);      % [H W C]
targetSize = inSize(1:2);
needRGB    = (numel(inSize)>=3 && inSize(3)==3);
fprintf('Model input size: [%d %d %d]\n', inSize(1), inSize(2), max(1,inSize(3)));

% ------------ Video & ROI ------------
vr = VideoReader(VIDEO);
% (Compat change) Prefer ROI loaded directly from step1_config.mat if present
if isfield(S1all,'roi') && numel(S1all.roi)==4
    roi = double(S1all.roi);
    disp('Loaded ROI from step1_config.mat.');
elseif exist(ROI_FILE,'file')
    R = load(ROI_FILE); roi = double(R.roi); disp('Loaded ROI from file.');
else
    f0 = readFrame(vr); figure; imshow(f0); title('Select ROI used in Step-1/2');
    roi = round(getrect()); close;
    assert(~isempty(roi) && all(roi(3:4)>0),'Invalid ROI.'); save(ROI_FILE,'roi');
end
vr = VideoReader(VIDEO); % rewind

% Writer (pad odd dims if needed)
first = read(vr,1);
needH = mod(size(first,1),2); needW = mod(size(first,2),2);
padVec = [needH needW 0];
if needH || needW, padFcn = @(im) padarray(im, padVec, 'replicate', 'post');
else,               padFcn = @(im) im;
end
try
    vw = VideoWriter(paths.out.video,'MPEG-4');
catch
    warning('MPEG-4 not available; falling back to Motion JPEG AVI.');
    [p,f,~] = fileparts(paths.out.video);
    paths.out.video = fullfile(p,[f '.avi']);
    vw = VideoWriter(paths.out.video,'Motion JPEG AVI');
end
vw.FrameRate = vr.FrameRate; vw.Quality = 95; open(vw);

% ------------ Main loop ------------
pupilDiameters  = [];
frameTimestamps = [];
centerX = []; centerY = [];
usedFallbackProb = [];

lastGoodDiameter = NaN;
lastGoodCenter   = [NaN NaN];
anchorCenter     = [NaN NaN];
warmCenters      = [];

frameIdx = 0; t0 = tic;
while hasFrame(vr)
    frameIdx = frameIdx + 1;
    frame = readFrame(vr);

    % ROI + grayscale
    eyeROI = imcrop(frame, roi);
    if isempty(eyeROI), continue; end
    if size(eyeROI,3)==3, grayROI = rgb2gray(eyeROI); else, grayROI = eyeROI; end
    if islogical(grayROI), grayROI = uint8(grayROI)*255; end
    if ~isa(grayROI,'uint8'), grayROI = im2uint8(mat2gray(grayROI)); end

    % Preprocess to MATCH Step-1
    imForNet = preproc(grayROI, PREPROC_MODE, stretchTail);

    % Net input
    netInput = imresize(imForNet, targetSize, 'nearest');
    if needRGB, netInput = repmat(netInput,[1 1 3]); end

    % Segmentation (with optional score-map fallback)
    usedFallback = false;
    [maskSmall, scoreOK] = localSemanticMask(netInput, net);
    if (~scoreOK || nnz(maskSmall) < fallback.minPix)
        sc = localScoreMap(netInput, net);
        if ~isempty(sc)
            thr = fallback.baseThr;
            if fallback.useOtsu, thr = max(thr, graythresh(sc)); end
            maskSmall = (sc >= thr);
            usedFallback = true;
        end
    end

    % Resize to ROI, morph cleanup
    mask = imresize(maskSmall, size(grayROI), 'nearest');
    mask = imopen(mask, SE_open);
    mask = imclose(mask, SE_close);
    mask = imfill(mask, 'holes');
    if minBlobArea>0, mask = bwareaopen(mask, minBlobArea); end

    % Anchor & search window
    statsAll = regionprops(mask, 'Area','Centroid','EquivDiameter');
    if ~isempty(statsAll)
        [~,kmax] = max([statsAll.Area]); c0 = statsAll(kmax).Centroid;
        if numel(warmCenters) < warmupFrames
            warmCenters(end+1,:) = c0; %#ok<AGROW>
        end
        if numel(warmCenters) >= max(5, round(0.5*warmupFrames))
            anchorCenter = median(warmCenters,1,'omitnan');
        end
    end
    refCenter = lastGoodCenter;
    if any(isnan(refCenter))
        refCenter = anchorCenter;
        if any(isnan(refCenter)) && ~isempty(statsAll)
            refCenter = statsAll(kmax).Centroid;
        end
    end
    if all(isfinite(refCenter))
        W = diskWindow(size(mask), refCenter, searchRadiusPx);
        mask = mask & W;
        statsAll = regionprops(mask, 'Area','Centroid','EquivDiameter');
    end

    % Choose blob nearest to prev/anchor; plausibility checks
    diameter = NaN; centerROI = [NaN NaN];
    if ~isempty(statsAll)
        centers = cat(1, statsAll.Centroid);
        if all(isfinite(lastGoodCenter))
            d = hypot(centers(:,1)-lastGoodCenter(1), centers(:,2)-lastGoodCenter(2));
            [~,k] = min(d);
        elseif all(isfinite(anchorCenter))
            d = hypot(centers(:,1)-anchorCenter(1), centers(:,2)-anchorCenter(2));
            [~,k] = min(d);
        else
            [~,k] = max([statsAll.Area]);
        end
        centerROI = statsAll(k).Centroid;
        diameter  = statsAll(k).EquivDiameter;

        % Jump/size guards
        accept = true;
        if all(isfinite(lastGoodCenter))
            if hypot(centerROI(1)-lastGoodCenter(1), centerROI(2)-lastGoodCenter(2)) > maxCenterJumpPx
                accept = false;
            end
        end
        if accept && ~isnan(lastGoodDiameter) && ~isempty(maxDiamChangeFrac)
            if abs(diameter-lastGoodDiameter) > maxDiamChangeFrac*max(1,lastGoodDiameter)
                accept = false;
            end
        end
        if ~accept
            if holdLastIfRejected && ~isnan(lastGoodDiameter)
                centerROI = lastGoodCenter; diameter = lastGoodDiameter;
            else
                centerROI = [NaN NaN]; diameter = NaN;
            end
        end
    else
        if holdLastIfRejected && ~isnan(lastGoodDiameter)
            centerROI = lastGoodCenter; diameter = lastGoodDiameter;
        end
    end

    % EMA smoothing
    if ~any(isnan(centerROI))
        if all(isfinite(lastGoodCenter))
            centerROI = (1-alphaCenter)*lastGoodCenter + alphaCenter*centerROI;
        end
        lastGoodCenter = centerROI;
    end
    if ~isnan(diameter)
        if ~isnan(lastGoodDiameter)
            diameter = (1-alphaDiam)*lastGoodDiameter + alphaDiam*diameter;
        end
        lastGoodDiameter = diameter;
    end

    % Save per-frame results
    pupilDiameters  = [pupilDiameters; lastGoodDiameter];
    frameTimestamps = [frameTimestamps; vr.CurrentTime];
    if all(isfinite(lastGoodCenter))
        centerX = [centerX; lastGoodCenter(1)];
        centerY = [centerY; lastGoodCenter(2)];
    else
        centerX = [centerX; NaN]; centerY = [centerY; NaN];
    end
    usedFallbackProb = [usedFallbackProb; usedFallback];

    % --- Visualization (robust sizes) ---
    imDisp = imForNet; 
    if ~isa(imDisp,'uint8'), imDisp = im2uint8(mat2gray(imDisp)); end
    if size(imDisp,3)==1, imDisp = repmat(imDisp,[1 1 3]); end

    zoomPanel = imDisp;
    if drawMaskOutline
        per = bwperim(mask);
        if size(zoomPanel,3)==1, zoomPanel = repmat(zoomPanel,[1 1 3]); end
        zoomPanel(:,:,3) = max(zoomPanel(:,:,3), uint8(255*per)); % blue outline
    end
    if ~isnan(lastGoodDiameter) && all(isfinite(lastGoodCenter))
        zoomPanel = insertShape(zoomPanel,'Circle',...
            [lastGoodCenter lastGoodDiameter/2], 'LineWidth',2,'Color','cyan');
    end

    % Annotated frame
    annotated = frame;
    if ~isnan(lastGoodDiameter) && all(isfinite(lastGoodCenter))
        centerFrame = [roi(1)+lastGoodCenter(1), roi(2)+lastGoodCenter(2)];
        annotated = insertShape(annotated,'Circle',...
            [centerFrame lastGoodDiameter/2],'LineWidth',2,'Color','yellow');
        annotated = insertText(annotated,[10 10], ...
            sprintf('Diameter: %.2f px', lastGoodDiameter), ...
            'FontSize',12,'BoxColor','black','TextColor','white');
    end

    % === FORCE MATCHED HEIGHTS ===
    H = size(annotated,1);
    zoomPanel = imresize(zoomPanel, [H NaN], 'Method', displayZoomInterp);  % exact height
    if size(zoomPanel,3)==1, zoomPanel = repmat(zoomPanel,[1 1 3]); end
    if ~isa(zoomPanel,'uint8'), zoomPanel = im2uint8(mat2gray(zoomPanel)); end

    % Separator with same height & 3ch
    sep = zeros(H, 10, 3, 'like', annotated);

    % Now safe to concatenate horizontally
    outFrame = [annotated, sep, zoomPanel];

    % Pad to even dims for MPEG-4 if needed
    outFrame = padFcn(outFrame);
    writeVideo(vw, outFrame);

    if reportEvery>0 && mod(frameIdx, reportEvery)==0
        fprintf('Frame %d | t=%.1fs | diam=%.2f px | fallback=%d\n', ...
            frameIdx, vr.CurrentTime, lastGoodDiameter, usedFallback);
    end
end

% Save CSV + quick plots
close(vw); disp('--- Final processing done. ---');

Tout = table(frameTimestamps, pupilDiameters, centerX, centerY, usedFallbackProb, ...
    'VariableNames',{'Timestamp_seconds','Diameter_pixels','CenterX_px','CenterY_px','UsedScoreFallback'});
writetable(Tout, paths.out.csv);
disp(['CSV -> ', paths.out.csv]);

figure('Name','Pupil Diameter'); 
raw = pupilDiameters; if all(isnan(raw)), warning('All diameters are NaN'); raw(:)=0; end
try
    cleaned  = fillmissing(raw,'linear'); smoothed = sgolayfilt(cleaned,3,11);
catch
    cleaned  = fillmissing(raw,'linear'); smoothed = movmean(cleaned,11);
end
plot(frameTimestamps, raw,'Color',[0.8 0.8 0.8],'DisplayName','Raw'); hold on;
plot(frameTimestamps, smoothed,'b','LineWidth',1.5,'DisplayName','Smoothed'); grid on;
xlabel('Time (s)'); ylabel('Diameter (px)'); legend('show');
saveas(gcf, fullfile(paths.out.dir,'Final_U-Net_TimeSeries.png'));

% write clean CSV (adds smoothed column)
Tclean = Tout; Tclean.SmoothedDiameter_px = smoothed;
writetable(Tclean, fullfile(paths.out.dir,'pupil_diameter_FINAL_U-Net_results_CLEAN.csv'));

fprintf('\nOutputs written to:\n  %s\n  %s\n  %s\n', ...
    paths.out.video, ...
    paths.out.csv, ...
    fullfile(paths.out.dir,'pupil_diameter_FINAL_U-Net_results_CLEAN.csv'));

% ---------------------- Local functions ----------------------
function Iout = preproc(grayROI,mode,stretchTail)
    if islogical(grayROI), grayROI = uint8(grayROI)*255; end
    if ~isa(grayROI,'uint8'), grayROI = im2uint8(mat2gray(grayROI)); end
    if mode.clahe
        Iout = adapthisteq(grayROI);
    elseif mode.stretchlim
        Iout = imadjust(grayROI, stretchlim(grayROI, stretchTail), []);
    else
        Iout = grayROI;
    end
end

function W = diskWindow(sz, center, radius)
    [X,Y] = meshgrid(1:sz(2), 1:sz(1));
    W = ((X-center(1)).^2 + (Y-center(2)).^2) <= radius^2;
end

function inSize = localGetInputSize(net)
    inSize = [];
    try
        inSize = net.Layers(1).InputSize;     % DAG/Series
    catch
        try
            if isprop(net,'InputSizes') && ~isempty(net.InputSizes)
                s = net.InputSizes{1};
                if isstruct(s) && isfield(s,'Size'), inSize = s.Size;
                elseif isnumeric(s), inSize = s;
                end
            end
        catch
        end
    end
    assert(~isempty(inSize),'Could not determine network input size.');
    if numel(inSize)==2, inSize = [inSize 1]; end
end

function [maskSmall, scoreOK] = localSemanticMask(Iin, net)
    scoreOK = true;
    try
        [C, ~] = semanticseg(Iin, net);
    catch
        C = semanticseg(Iin, net);
        scoreOK = false;
    end
    maskSmall = (C=='pupil');
end

function sc = localScoreMap(Iin, net)
    sc = [];
    try
        [C, scores] = semanticseg(Iin, net);
        if isempty(scores), return; end
        K = size(scores,3);
        if K<=1
            sc = im2double(scores(:,:,1));
        else
            labs = categories(C);
            if numel(labs)==K
                idx = find(strcmpi(labs,'pupil'),1);
                if isempty(idx), idx = min(2,K); end
            else
                idx = min(2,K);
            end
            sc = im2double(scores(:,:,idx));
        end
    catch
        sc = [];
    end
end
