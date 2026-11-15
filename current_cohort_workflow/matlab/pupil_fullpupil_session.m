function outDir = pupil_fullpupil_session(sessionDir, varargin)
% pupil_fullpupil_session  Full-pupil trace of one session whose pupil fills the eye opening (fully dilated pupil).
%
%   outDir = pupil_fullpupil_session('<pupil_out_root>\dayN\<cage>\<mouse>')
%   ... = pupil_fullpupil_session(..., 'MaxFrames', 600)       % quick test
%
% Rule (from the user, 2026-11-15): when the whole eye is dark, the pupil is the WHOLE dark area inside the eyelids.
% The package's tracker cannot measure that: its U-Net takes a small dark part, and its step 3 keeps only the part of
% the mask within 30 px of the last centre, so no pupil can be wider than 60 px (a full pupil here is ~120 px).
% This function reads the same video and ROI as the package (pupil_outputs_general\qc.json) and, per frame:
%   1) the ROI crop exactly as pupil_mask_frame makes it (imcrop, rgb2gray, uint8)
%   2) fullpupil_mask (pupil_finetune\fullpupil_mask.m): the whole dark area inside the lids
%   3) a leak into dark fur or shadow: the region touches the ROI edge on more than BorderMax (40) pixels AND is more
%      than 5 % wider than the last accepted value -> rejected (hold the last value). An eye that is only cut by the
%      ROI edge (same or smaller size) is measured as it is.
%      Then the package's step-3 rules on that one region, unchanged except the search-radius cut: reject a centre
%      jump > 18 px or a diameter change > 18 % (hold the last value), EMA (centre 0.18, diameter 0.12), forget the
%      last value after ReacquireFrames frames without an accepted one (qc.json reacquireFrames, else 15)
%   4) the package's cleaning (fillmissing + sgolay(3, 11) -> SmoothedDiameter_px), but no values are made up:
%      gaps longer than MaxGap (30 frames, ~1 s) and the ends stay empty (NaN) instead of being interpolated or
%      extrapolated. The diameter is the equivalent diameter of the whole mask (all its pieces).
% Output in <sessionDir>\pupil_outputs_fullpupil\ (the package's folders are not touched; a test run with MaxFrames
% writes to pupil_outputs_fullpupil_test\ unless OutName is given):
%   pupil_fullpupil_results2.csv / pupil_fullpupil_CLEAN.csv  same columns as the package's results2 / CLEAN files
%   fullpupil_raw.csv   per frame: raw diameter, centre, Otsu level, accepted / rejected / no region, re-locks, and
%                       darkest_px = the black region around the package tracker's pupil (its centre from the package's
%                       results2.csv): ~20-30 px for a normal pupil, > 60 px when the pupil fills the eye. The 1-minute
%                       switch of fullpupil_sessions.py (mouse option) is decided on it.
%   fullpupil_qc.json   numbers;  fullpupil_sheet.png  24 frames over the session with the outline
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'pupil_finetune'));          % fullpupil_mask.m
p = inputParser; p.addParameter('MaxFrames', Inf); p.addParameter('OutName', '');
p.addParameter('ReacquireFrames', NaN); p.addParameter('BorderMax', 40); p.addParameter('LeakGain', 0.05);
p.addParameter('MaxGap', 30); p.parse(varargin{:}); o = p.Results;
Q = jsondecode(fileread(fullfile(sessionDir, 'pupil_outputs_general', 'qc.json')));
roi = double(Q.roi(:)');
if isempty(o.OutName)
    if isfinite(o.MaxFrames), o.OutName = 'pupil_outputs_fullpupil_test'; else, o.OutName = 'pupil_outputs_fullpupil'; end
end
if isnan(o.ReacquireFrames)   % the same value the package used for this session
    if isfield(Q, 'reacquireFrames') && isfinite(double(Q.reacquireFrames)), o.ReacquireFrames = double(Q.reacquireFrames);
    else, o.ReacquireFrames = 15; end
end
if isfield(Q, 'frameResized') && logical(Q.frameResized)
    error('%s: the package resized the frames of this session; not supported here', sessionDir);
end
outDir = fullfile(sessionDir, o.OutName); if ~exist(outDir, 'dir'), mkdir(outDir); end
pkR = readtable(fullfile(sessionDir, 'pupil_outputs_general', 'pupil_diameter_FINAL_U-Net_results2.csv'), ...
    'VariableNamingRule', 'preserve');           % the package tracker's pupil centre per frame (seed of darkest_px)
seedX = pkR.CenterX_px; seedY = pkR.CenterY_px;
vr = VideoReader(Q.video);
nEst = floor(vr.Duration * vr.FrameRate); nMax = min(nEst + 10, o.MaxFrames);
darkD = nan(nMax, 1);
alphaCenter = 0.18; alphaDiam = 0.12; maxCenterJumpPx = 18; maxDiamChangeFrac = 0.18;
ts = nan(nMax, 1); dia = nan(nMax, 1); cx = nan(nMax, 1); cy = nan(nMax, 1);
rawD = nan(nMax, 1); rawX = nan(nMax, 1); rawY = nan(nMax, 1); thr = nan(nMax, 1);
state = zeros(nMax, 1, 'int8');                           % 1 accepted, 2 rejected (held), 3 leak at the ROI edge, 0 none
border = nan(nMax, 1);
reacq = zeros(0, 1);
qcIdx = round(linspace(1, max(1, min(nEst, nMax)), 24)); qcImg = {}; qcMask = {}; qcK = [];
lastGoodCenter = [NaN NaN]; lastGoodDiameter = NaN; missRun = 0; k = 0; t0 = tic;
while hasFrame(vr) && k < nMax
    t = vr.CurrentTime; fr = readFrame(vr); k = k + 1; ts(k) = t;
    C = imcrop(fr, roi); if size(C, 3) == 3, C = rgb2gray(C); end
    if ~isa(C, 'uint8'), C = im2uint8(mat2gray(C)); end
    [m, fi] = fullpupil_mask(C); thr(k) = fi.thr;
    if any(m(:)) && k <= numel(seedX) && isfinite(seedX(k)) && isfinite(seedY(k))
        darkD(k) = darkest_part(fi.S, m, fi.G, seedX(k), seedY(k));
    end
    if any(m(:))
        border(k) = nnz(m(1, :)) + nnz(m(end, :)) + nnz(m(:, 1)) + nnz(m(:, end));
        R = regionprops(m, 'Centroid', 'Area');                   % centre of the largest piece,
        [~, big] = max([R.Area]);                                 % diameter of the whole mask (all its pieces)
        centerROI = R(big).Centroid; diameter = 2 * sqrt(nnz(m) / pi);
        rawD(k) = diameter; rawX(k) = centerROI(1); rawY(k) = centerROI(2);
    end
    if any(m(:)) && border(k) > o.BorderMax && ~isnan(lastGoodDiameter) && diameter > (1 + o.LeakGain) * lastGoodDiameter
        state(k) = 3; centerROI = lastGoodCenter; diameter = lastGoodDiameter;   % leak at the ROI edge
    elseif any(m(:))
        accept = true;
        if all(isfinite(lastGoodCenter)) && hypot(centerROI(1) - lastGoodCenter(1), centerROI(2) - lastGoodCenter(2)) > maxCenterJumpPx
            accept = false;
        end
        if accept && ~isnan(lastGoodDiameter) && abs(diameter - lastGoodDiameter) > maxDiamChangeFrac * max(1, lastGoodDiameter)
            accept = false;
        end
        if accept, state(k) = 1;
        else, state(k) = 2; centerROI = lastGoodCenter; diameter = lastGoodDiameter;
        end
    else
        centerROI = lastGoodCenter; diameter = lastGoodDiameter;
    end
    if ~any(isnan(centerROI))
        if all(isfinite(lastGoodCenter)), centerROI = (1 - alphaCenter) * lastGoodCenter + alphaCenter * centerROI; end
        lastGoodCenter = centerROI;
    end
    if ~isnan(diameter)
        if ~isnan(lastGoodDiameter), diameter = (1 - alphaDiam) * lastGoodDiameter + alphaDiam * diameter; end
        lastGoodDiameter = diameter;
    end
    dia(k) = lastGoodDiameter; cx(k) = lastGoodCenter(1); cy(k) = lastGoodCenter(2);
    if o.ReacquireFrames > 0
        if state(k) == 1, missRun = 0; else, missRun = missRun + 1; end
        if missRun >= o.ReacquireFrames
            lastGoodCenter = [NaN NaN]; lastGoodDiameter = NaN; missRun = 0; reacq(end + 1, 1) = k; %#ok<AGROW>
        end
    end
    if any(qcIdx == k), qcImg{end + 1} = C; qcMask{end + 1} = m; qcK(end + 1) = k; end %#ok<AGROW>
    if mod(k, 5000) == 0, fprintf('  frame %d / ~%d (%.0f fps) full pupil %.1f px\n', k, nEst, k / toc(t0), lastGoodDiameter); end
end
n = k; ts = ts(1:n); dia = dia(1:n); cx = cx(1:n); cy = cy(1:n);
rawD = rawD(1:n); rawX = rawX(1:n); rawY = rawY(1:n); thr = thr(1:n); state = state(1:n); border = border(1:n);
darkD = darkD(1:n);
Tout = table(ts, dia, cx, cy, zeros(n, 1), 'VariableNames', ...
    {'Timestamp_seconds', 'Diameter_pixels', 'CenterX_px', 'CenterY_px', 'UsedScoreFallback'});
writetable(Tout, fullfile(outDir, 'pupil_fullpupil_results2.csv'));
% the package's cleaning without made-up values: short gaps (<= MaxGap frames) are filled linearly as in the package;
% longer gaps and the ends before the first / after the last value stay NaN
raw = dia;
if all(isnan(raw))
    smoothed = nan(n, 1);
else
    filled = fillmissing(raw, 'linear', 'EndValues', 'nearest');          % only for the filter (no NaN inside)
    try, smoothed = sgolayfilt(filled, 3, 11); catch, smoothed = movmean(filled, 11); end
    gap = isnan(fillmissing(raw, 'linear', 'EndValues', 'none', 'MaxGap', o.MaxGap));
    smoothed(gap) = NaN;
end
nanRuns = diff([0; isnan(smoothed); 0]); runLen = find(nanRuns == -1) - find(nanRuns == 1);
Tc = Tout; Tc.SmoothedDiameter_px = smoothed;
writetable(Tc, fullfile(outDir, 'pupil_fullpupil_CLEAN.csv'));
Traw = table((1:n)', ts, rawD, rawX, rawY, thr, border, double(state), ismember((1:n)', reacq), darkD, 'VariableNames', ...
    {'frame', 't_s', 'raw_diam_px', 'raw_cx_px', 'raw_cy_px', 'otsu_level', 'roi_edge_px', 'state_1acc_2rej_3edge_0none', ...
     'reacquired', 'darkest_px'});
writetable(Traw, fullfile(outDir, 'fullpupil_raw.csv'));
pk = readtable(fullfile(sessionDir, 'pupil_outputs_general', 'pupil_diameter_FINAL_U-Net_results_CLEAN.csv'), ...
    'VariableNamingRule', 'preserve');
info = struct('session', sessionDir, 'video', Q.video, 'roi', roi, 'nFrames', n, 'packageFrames', height(pk), ...
    'fracAccepted', mean(state == 1), 'fracRejected', mean(state == 2), 'fracEdgeLeak', mean(state == 3), ...
    'fracNoRegion', mean(state == 0), 'borderMax', o.BorderMax, 'leakGain', o.LeakGain, ...
    'reacquireFrames', o.ReacquireFrames, 'maxGap', o.MaxGap, 'fracNaN', mean(isnan(smoothed)), ...
    'longestNaNRun', max([0; runLen(:)]), ...
    'nReacquired', numel(reacq), 'medianDiameter', median(dia, 'omitnan'), 'iqrDiameter', iqr(dia), ...
    'medianDiameterPackage', median(pk.SmoothedDiameter_px, 'omitnan'), 'when', char(datetime('now')), ...
    'rule', 'whole dark area inside the eyelids (fullpupil_mask.m); package step-3 rules without the 30-px search radius');
fid = fopen(fullfile(outDir, 'fullpupil_qc.json'), 'w'); fwrite(fid, jsonencode(info, 'PrettyPrint', true)); fclose(fid);
% check sheet: 24 frames over the session
fig = figure('Color', 'w', 'Position', [50 50 1600 1100], 'Visible', 'off');
tl = tiledlayout(fig, 4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
for q = 1:numel(qcImg)
    nexttile(tl); I = double(qcImg{q}); lo = prctile(I(:), 0.5); hi = prctile(I(:), 99.5);
    imshow(min(1, max(0, (I - lo) / max(1, hi - lo))) .^ 0.7); hold on;
    if any(qcMask{q}(:)), contour(qcMask{q}, [0.5 0.5], 'r', 'LineWidth', 1.1); end
    title(sprintf('frame %d: %.0f px (U-Net %.0f)', qcK(q), dia(qcK(q)), pk.SmoothedDiameter_px(min(height(pk), qcK(q)))), 'FontSize', 8);
end
title(tl, sprintf('%s - full pupil (whole dark area inside the lids), red outline', strrep(sessionDir, '\', '\\')), 'Interpreter', 'none');
exportgraphics(fig, fullfile(outDir, 'fullpupil_sheet.png'), 'Resolution', 100); close(fig);
fprintf('full pupil: %d frames (package %d) | accepted %.1f %% rejected %.1f %% edge leak %.1f %% none %.1f %% | median %.1f px (U-Net %.1f px) | %.0f s\n', ...
    n, height(pk), 100 * info.fracAccepted, 100 * info.fracRejected, 100 * info.fracEdgeLeak, 100 * info.fracNoRegion, info.medianDiameter, ...
    info.medianDiameterPackage, toc(t0));
end

% ======================================================================================================================
function d = darkest_part(S, A, G, sx, sy)
% Diameter (px) of the black region around the tracker's pupil (seed sx, sy = the package's centre, ROI pixels): inside
% the dark area A, pixels darker than halfway between its darkest 5 % and its 75th percentile (glints count as dark),
% opened with a disk of 3 px, the piece at the seed (or the nearest piece), holes filled. Same rule as darkest_part in
% fullpupil_sessions.py. A normal pupil gives ~20-30 px; a pupil that fills the eye gives the large black region.
d = NaN;
v = S(A & ~G);
if numel(v) < 50, return; end
dk = prctile(v, 5); md = prctile(v, 75);
B = ((S < dk + 0.5 * (md - dk)) | G) & A;
B = imopen(B, strel('disk', 3));
CCb = bwconncomp(B);
if CCb.NumObjects == 0, return; end
L = labelmatrix(CCb);
r = min(max(round(sy), 1), size(B, 1)); c = min(max(round(sx), 1), size(B, 2));
l = L(r, c);
if l == 0
    [~, idx] = bwdist(B); l = L(idx(r, c));
end
d = 2 * sqrt(nnz(imfill(L == l, 'holes')) / pi);
end
