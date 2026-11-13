function Q = pupil_track_session(videoFile, roi, model, outDir, varargin)
% pupil_track_session  Pupil diameter for one video with the GENERAL U-Net - no clicks, no training (2026-11-12).
%
%   Q = pupil_track_session(videoFile, roi, model, outDir)
%   Q = pupil_track_session(..., 'MaxFrames', 3000, 'Batch', 128, 'QCVideo', false)
%
% videoFile  the pupil camera .avi (the one whose frame count matches the camera TTLs)
% roi        [x y w h] eye ROI in video pixels (roi.mat of that mouse-day)
% model      path to the general model .mat (pupil_train_general) or the loaded struct
% outDir     output folder; writes the SAME files and columns as step3_dec291123342.m:
%              pupil_diameter_FINAL_U-Net_results2.csv   Timestamp_seconds, Diameter_pixels, CenterX_px, CenterY_px, UsedScoreFallback
%              pupil_diameter_FINAL_U-Net_results_CLEAN.csv   (+ SmoothedDiameter_px = fillmissing + sgolay(3,11))
%              Final_U-Net_TimeSeries.png, qc_frames.png (12 frames with the mask outline), qc.json
% Per frame (as step 1 made the training images): crop ROI -> gray -> glint suppression (98th pct, dilate 2) -> CLAHE
%   -> resize to the net input -> pupil probability; mask = p > 0.5 (fallback: p >= max(0.40, Otsu) when < 50 px).
% Then step 3's tracking UNCHANGED: open/close/fill/area filter, warm-up anchor (30 frames), search radius 30 px,
% reject centre jumps > 18 px or diameter changes > 18 % (hold last), EMA alpha centre 0.18 / diameter 0.12.
% Difference from step 3: every frame is processed (step 3 skipped frame 1), so rows = video frames = camera TTLs.

p = inputParser;
p.addParameter('MaxFrames', Inf);
p.addParameter('Batch', 128);
p.addParameter('QCVideo', false);
p.parse(varargin{:});
o = p.Results;

if ischar(model) || isstring(model), M = load(model); else, M = model; end
net = M.net; S = M.INPUT; PP = M.PREPROC;
if ~exist(outDir, 'dir'), mkdir(outDir); end
roi = double(roi(:)');

% ---- step 3 tracking parameters (default profile) ----
fallback.minPix = 50; fallback.baseThr = 0.40;
minBlobArea = 50; SE_open = strel('disk',1); SE_close = strel('disk',2);
warmupFrames = 30; searchRadiusPx = 30; alphaCenter = 0.18; alphaDiam = 0.12;
maxCenterJumpPx = 18; maxDiamChangeFrac = 0.18; holdLastIfRejected = true;

vr = VideoReader(videoFile);
nEst = floor(vr.Duration * vr.FrameRate);
nMax = min(nEst + 10, o.MaxFrames);
ts = nan(nMax,1); dia = nan(nMax,1); cx = nan(nMax,1); cy = nan(nMax,1); fb = false(nMax,1);
rawFound = false(nMax,1); rejected = false(nMax,1);
qcIdx = round(linspace(1, max(1, min(nEst, nMax)), 12)); qcImg = {}; qcMask = {};

lastGoodCenter = [NaN NaN]; lastGoodDiameter = NaN; anchorCenter = [NaN NaN]; warmCenters = zeros(0,2);
roiSz = [];   % = size of the imcrop result (imcrop gives h+1 x w+1 for [x y w h]); step 3 also used size(grayROI)
k = 0; t0 = tic;
X = zeros(S, S, 1, o.Batch, 'uint8'); G = cell(o.Batch,1); tb = zeros(o.Batch,1); nb = 0;
while hasFrame(vr) && k + nb < nMax
    tb(nb+1) = vr.CurrentTime;
    F = readFrame(vr);
    C = imcrop(F, roi); if size(C,3) == 3, C = rgb2gray(C); end
    if ~isa(C, 'uint8'), C = im2uint8(mat2gray(C)); end
    N = preprocess_like_step1(C, PP);
    if isempty(roiSz), roiSz = size(N); [XX, YY] = meshgrid(1:roiSz(2), 1:roiSz(1)); end
    nb = nb + 1; G{nb} = N; X(:,:,1,nb) = imresize(N, [S S], 'bilinear');
    if nb == o.Batch || ~hasFrame(vr) || k + nb >= nMax
        P = predict(net, X(:,:,:,1:nb), 'MiniBatchSize', nb);   % S x S x 2 x nb softmax
        for j = 1:nb
            k = k + 1; ts(k) = tb(j);
            pm = P(:,:,2,j);
            ms = pm > 0.5; usedFallback = false;
            if nnz(ms) < fallback.minPix
                thr = max(fallback.baseThr, graythresh(pm)); ms = pm >= thr; usedFallback = true;
            end
            mask = imresize(ms, roiSz, 'nearest');
            mask = imopen(mask, SE_open); mask = imclose(mask, SE_close); mask = imfill(mask, 'holes');
            mask = bwareaopen(mask, minBlobArea);
            % ---------------- step 3 tracking (unchanged) ----------------
            statsAll = regionprops(mask, 'Area', 'Centroid', 'EquivDiameter');
            if ~isempty(statsAll)
                [~,kmax] = max([statsAll.Area]); c0 = statsAll(kmax).Centroid;
                if size(warmCenters,1) < warmupFrames, warmCenters(end+1,:) = c0; end %#ok<AGROW>
                if size(warmCenters,1) >= max(5, round(0.5*warmupFrames)), anchorCenter = median(warmCenters,1,'omitnan'); end
            end
            refCenter = lastGoodCenter;
            if any(isnan(refCenter))
                refCenter = anchorCenter;
                if any(isnan(refCenter)) && ~isempty(statsAll), refCenter = statsAll(kmax).Centroid; end
            end
            if all(isfinite(refCenter))
                mask = mask & (((XX-refCenter(1)).^2 + (YY-refCenter(2)).^2) <= searchRadiusPx^2);
                statsAll = regionprops(mask, 'Area', 'Centroid', 'EquivDiameter');
            end
            rawFound(k) = ~isempty(statsAll);
            if ~isempty(statsAll)
                centers = cat(1, statsAll.Centroid);
                if all(isfinite(lastGoodCenter)), [~,ii] = min(hypot(centers(:,1)-lastGoodCenter(1), centers(:,2)-lastGoodCenter(2)));
                elseif all(isfinite(anchorCenter)), [~,ii] = min(hypot(centers(:,1)-anchorCenter(1), centers(:,2)-anchorCenter(2)));
                else, [~,ii] = max([statsAll.Area]);
                end
                centerROI = statsAll(ii).Centroid; diameter = statsAll(ii).EquivDiameter; accept = true;
                if all(isfinite(lastGoodCenter)) && hypot(centerROI(1)-lastGoodCenter(1), centerROI(2)-lastGoodCenter(2)) > maxCenterJumpPx
                    accept = false;
                end
                if accept && ~isnan(lastGoodDiameter) && abs(diameter-lastGoodDiameter) > maxDiamChangeFrac*max(1,lastGoodDiameter)
                    accept = false;
                end
                if ~accept
                    rejected(k) = true;
                    if holdLastIfRejected && ~isnan(lastGoodDiameter), centerROI = lastGoodCenter; diameter = lastGoodDiameter;
                    else, centerROI = [NaN NaN]; diameter = NaN;
                    end
                end
            else
                centerROI = [NaN NaN]; diameter = NaN;
                if holdLastIfRejected && ~isnan(lastGoodDiameter), centerROI = lastGoodCenter; diameter = lastGoodDiameter; end
            end
            if ~any(isnan(centerROI))
                if all(isfinite(lastGoodCenter)), centerROI = (1-alphaCenter)*lastGoodCenter + alphaCenter*centerROI; end
                lastGoodCenter = centerROI;
            end
            if ~isnan(diameter)
                if ~isnan(lastGoodDiameter), diameter = (1-alphaDiam)*lastGoodDiameter + alphaDiam*diameter; end
                lastGoodDiameter = diameter;
            end
            dia(k) = lastGoodDiameter; cx(k) = lastGoodCenter(1); cy(k) = lastGoodCenter(2); fb(k) = usedFallback;
            % --------------------------------------------------------------
            if any(qcIdx == k), qcImg{end+1} = G{j}; qcMask{end+1} = mask; end %#ok<AGROW>
        end
        nb = 0;
        if mod(k, 3000) < o.Batch
            fprintf('  frame %d / ~%d  (%.0f fps)  diam %.2f px\n', k, nEst, k/toc(t0), lastGoodDiameter);
        end
    end
end
ts = ts(1:k); dia = dia(1:k); cx = cx(1:k); cy = cy(1:k); fb = fb(1:k); rawFound = rawFound(1:k); rejected = rejected(1:k);

Tout = table(ts, dia, cx, cy, double(fb), 'VariableNames', ...
    {'Timestamp_seconds','Diameter_pixels','CenterX_px','CenterY_px','UsedScoreFallback'});
writetable(Tout, fullfile(outDir, 'pupil_diameter_FINAL_U-Net_results2.csv'));
raw = dia; if all(isnan(raw)), raw(:) = 0; end
cleaned = fillmissing(raw, 'linear');
try, smoothed = sgolayfilt(cleaned, 3, 11); catch, smoothed = movmean(cleaned, 11); end
Tc = Tout; Tc.SmoothedDiameter_px = smoothed;
cleanCsv = fullfile(outDir, 'pupil_diameter_FINAL_U-Net_results_CLEAN.csv');
writetable(Tc, cleanCsv);

fig = figure('Visible','off','Color','w','Position',[100 100 1100 380]);
plot(ts, raw, 'Color', [0.8 0.8 0.8], 'DisplayName', 'Diameter (EMA)'); hold on;
plot(ts, smoothed, 'b', 'LineWidth', 1.2, 'DisplayName', 'Smoothed'); grid on;
xlabel('Video time (s)'); ylabel('Pupil diameter (ROI px)'); legend('Location','best');
[~, vn] = fileparts(videoFile); title(sprintf('%s  (n = %d frames)', vn, k), 'Interpreter', 'none');
exportgraphics(fig, fullfile(outDir, 'Final_U-Net_TimeSeries.png'), 'Resolution', 150); close(fig);

if ~isempty(qcImg)
    th = 160; tw = round(th * roi(3) / roi(4));
    tiles = repmat({zeros(th, tw, 3, 'uint8')}, 1, 12);           % 2 rows x 6, black where a frame is missing
    for i = 1:min(12, numel(qcImg))
        I3 = repmat(qcImg{i}, [1 1 3]); per = bwperim(qcMask{i});
        I3(:,:,1) = max(I3(:,:,1), uint8(255*per)); I3(:,:,2) = min(I3(:,:,2), uint8(255*~per));
        tiles{i} = imresize(I3, [th tw]);
    end
    imwrite([cat(2, tiles{1:6}); cat(2, tiles{7:12})], fullfile(outDir, 'qc_frames.png'));
end

Q = struct('video', char(videoFile), 'roi', roi, 'nFrames', k, 'fps', k/toc(t0), 'seconds', toc(t0), ...
    'fracNoBlob', mean(~rawFound), 'fracRejected', mean(rejected), 'fracFallback', mean(fb), ...
    'medianDiameter', median(dia, 'omitnan'), 'iqrDiameter', prctile(dia, 75) - prctile(dia, 25), ...
    'fracNaN', mean(isnan(dia)), 'cleanCsv', char(cleanCsv));
fid = fopen(fullfile(outDir, 'qc.json'), 'w'); fprintf(fid, '%s', jsonencode(Q, 'PrettyPrint', true)); fclose(fid);
fprintf('  done: %d frames in %.0f s (%.0f fps) | no blob %.1f %% | rejected %.1f %% | median diam %.2f px\n', ...
    k, Q.seconds, Q.fps, 100*Q.fracNoBlob, 100*Q.fracRejected, Q.medianDiameter);
end

function N = preprocess_like_step1(G, PP)
% = preprocess_for_step1 in step1_2dec29_nove222222222.m (SUPPRESS_GLINT true, USE_CLAHE true)
Gd = G;
Gl = Gd >= prctile(double(Gd(:)), PP.glint_pct);
Gl = imdilate(Gl, strel('disk', PP.glint_dilate));
if any(Gl(:)), Gd(Gl) = uint8(median(double(Gd(~Gl)), 'omitnan')); end
if PP.clahe, N = adapthisteq(Gd); else, N = Gd; end
end
