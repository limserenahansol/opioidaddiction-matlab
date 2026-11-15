function make_fullpupil_labels(diagDir, keys, outDir, pkgPupilDir, modelFile, candDir, fitDir)
% make_fullpupil_labels  Full-pupil masks for sessions whose pupil fills the eye opening (fully dilated pupil).
%
%   make_fullpupil_labels(diagDir, {'day2_9669m_green1', ...}, outDir, pkgPupilDir, modelFile, candDir, fitDir)
%
% Rule (from the user, 2026-11-15): when the whole eye is dark, the pupil is the WHOLE dark area inside the eyelids,
% not only its darkest part. For every sampled frame of the given sessions (the diag samples that the fine-tuning
% uses, pupil_diag_extract3.py) this writes the mask of that dark area, on exactly the ROI crop the U-Net sees:
%   1) ROI crop as pupil_mask_frame makes it (imcrop convention); smoothed (Gaussian sigma 1.5)
%   2) glints = the brightest 1 % of the crop, dilated 3 px
%   3) Otsu threshold between the dark eye and the skin (glint pixels left out), dark = below it
%   4) open with a disk of 4 px (cuts thin dark lines: lid margin hairs, the spout's shadow), keep the dark
%      component with the largest area in the central half of the crop, fill holes (glints inside the eye)
%   5) glints on the edge of the eye: their pixels are added where a closing (disk 6) of the mask covers them;
%      a glint outside the eye (below the lid) is not added
% Output per session: <key>_fullpupil.mat (masks{j} for sample j-1, thr, diam_px) and <key>_fullpupil_sheet.png
% (outline red = full-pupil mask, yellow = the old consensus label, if any) for checking by eye before training.
if ischar(keys) || isstring(keys), keys = cellstr(keys); end
addpath(pkgPupilDir); if ~exist(outDir, 'dir'), mkdir(outDir); end
M = load(modelFile);
L = table();
if nargin >= 7 && ~isempty(fitDir), L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string'); end
for k = 1:numel(keys)
    key = keys{k};
    Sd = load(fullfile(diagDir, [key '.mat'])); Xs = Sd.samples;
    r0 = double(Sd.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
    F = arrayfun(@(j) Xs(:, :, j), 1:size(Xs, 3), 'UniformOutput', false)';
    [~, ~, crops] = pupil_mask_frame(F, rr, M, 'IrisCorrection', true);
    C = []; if nargin >= 6 && ~isempty(candDir), C = load(fullfile(candDir, [key '_cands.mat'])); end
    n = numel(crops); masks = cell(n, 1); thr = nan(n, 1); diam = nan(n, 1);
    fig = figure('Color', 'w', 'Position', [50 50 1600 1100], 'Visible', 'off');
    tl = tiledlayout(fig, 4, ceil(n / 4), 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:n
        I = double(crops{j});
        [m, fi] = fullpupil_mask(crops{j});          % steps 1-5 above (one shared definition)
        thr(j) = fi.thr;
        masks{j} = m; diam(j) = fi.diam;
        nexttile(tl);
        lo = prctile(I(:), 0.5); hi = prctile(I(:), 99.5);
        imshow(min(1, max(0, (I - lo) / max(1, hi - lo))) .^ 0.7); hold on;
        if any(m(:)), contour(m, [0.5 0.5], 'r', 'LineWidth', 1.2); end
        lab = NaN;
        if ~isempty(L)
            row = L(L.session == string(key) & L.frame == j - 1, :);
            if ~isempty(row), lab = row.label(1); end
        end
        if ~isempty(C) && isfinite(lab) && lab >= 1
            contour(logical(C.allc{j}(lab).m), [0.5 0.5], 'y', 'LineWidth', 0.8);
        end
        title(sprintf('#%d  full %.0f px  (old label %g)', j - 1, diam(j), lab), 'FontSize', 8);
    end
    title(tl, sprintf('%s: red = whole dark area inside the lids (full pupil), yellow = old label', strrep(key, '_', ' ')));
    exportgraphics(fig, fullfile(outDir, [key '_fullpupil_sheet.png']), 'Resolution', 110); close(fig);
    save(fullfile(outDir, [key '_fullpupil.mat']), 'masks', 'thr', 'diam', 'key');
    fprintf('%s: %d frames, full-pupil diameter median %.0f px (range %.0f-%.0f)\n', key, n, median(diam, 'omitnan'), ...
        min(diam), max(diam));
end
end
