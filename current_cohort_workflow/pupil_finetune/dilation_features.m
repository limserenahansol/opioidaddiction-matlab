function dilation_features(diagDir, fitDir, pkgPupilDir, modelFile, outCsv)
% dilation_features  Per labelled frame: is the eye fully dark (pupil fills the eye opening) or is there a normal
% pupil inside a lighter iris? Numbers for choosing a rule, on every sampled frame of every session:
%   unet_px   diameter of the tracker's blob (pupil_mask_frame + pupil_maxp_mask, the production model)
%   full_px   diameter of the whole dark area inside the lids (fullpupil_mask)
%   pupil_lv  median level of the blob (or the 5th percentile of the dark area if there is no blob)
%   rest_lv   median level of the dark area outside the blob (dilated 3 px): iris (normal eye) or more pupil
%   skin_lv   median level outside the dark area
%   contrast  (rest_lv - pupil_lv) / (skin_lv - pupil_lv): small = the rest of the eye is as dark as the pupil
%   dark_frac share of the dark area darker than pupil_lv + 0.25 x (skin_lv - pupil_lv)
% Label columns come from fitDir/labels_consensus.csv (how = fullpupil_rule marks the fully dilated frames).
addpath(pkgPupilDir);
M = load(modelFile);
L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string');
meta = jsondecode(fileread(fullfile(diagDir, 'meta.json')));
out = table();
for k = 1:numel(meta)
    key = string(meta(k).key);
    Sd = load(fullfile(diagDir, [char(key) '.mat'])); Xs = Sd.samples;
    r0 = double(Sd.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
    F = arrayfun(@(j) Xs(:, :, j), 1:size(Xs, 3), 'UniformOutput', false)';
    [masks, ~, crops, ~, ~, ~, probs] = pupil_mask_frame(F, rr, M, 'IrisCorrection', true);
    for j = 1:numel(crops)
        [mm, info] = pupil_maxp_mask(masks{j}, probs{j}, 0.9);
        [A, fi] = fullpupil_mask(crops{j});
        S = fi.S; G = fi.G;
        if any(mm(:)), pl = median(S(mm & ~G)); else, pl = prctile(S(A & ~G), 5); end
        rest = A & ~imdilate(mm, strel('disk', 3)) & ~G;
        rl = median(S(rest)); sk = fi.skin;
        con = (rl - pl) / max(1, sk - pl);
        df = mean(S(A & ~G) < pl + 0.25 * (sk - pl));
        % thin ring just outside the tracker's blob (2-6 px), inside the dark area: iris in a normal eye, more pupil
        % when the pupil fills the eye
        ring = imdilate(mm, strel('disk', 6)) & ~imdilate(mm, strel('disk', 2)) & A & ~G;
        if any(mm(:)) && nnz(ring) >= 20
            ringc = (median(S(ring)) - pl) / max(1, sk - pl);
            ringdark = mean(S(ring) < pl + 0.15 * (sk - pl));
        else
            ringc = NaN; ringdark = NaN;
        end
        row = L(L.session == key & L.frame == j - 1, :);
        if isempty(row), lab = NaN; how = ""; else, lab = row.label(1); how = row.how(1); end
        out = [out; table(key, j - 1, lab, how, info.diam, fi.diam, pl, rl, sk, con, df, ringc, ringdark, ...
            'VariableNames', {'session', 'frame', 'label', 'how', 'unet_px', 'full_px', 'pupil_lv', 'rest_lv', ...
            'skin_lv', 'contrast', 'dark_frac', 'ring_contrast', 'ring_dark'})]; %#ok<AGROW>
    end
    fprintf('%s done\n', key);
end
writetable(out, outCsv);
end
