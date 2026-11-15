function eval_fullpupil(diagDir, candDir, fitDir, pkgPupilDir, modelDirs, modelNames, outCsv)
% eval_fullpupil  Held-out check of fold models with the tracker's own detector: pupil_mask_frame (with and without
% IrisCorrection) + pupil_maxp_mask (most probable blob, kept if max p >= 0.9). Labels: labels_consensus.csv of
% fitDir with the candidate masks of candDir (made by add_fullpupil_labels.m: in fully dilated sessions the label
% is the whole dark area inside the lids).
%   modelDirs  = {dir1, dir2, ...} each with unet_ft_fold1..4.mat (leave-mice-out folds of finetune_unet_cv)
%   modelNames = {'new', 'old', ...}
% Per frame: label kind (pupil / none / nocontour), full = label from the full-pupil rule, result (correct, wrong,
% miss, false, none_ok, ...), IoU and diameter ratio (equivalent diameter of the chosen blob / of the label).
addpath(pkgPupilDir);
FOLDS = {{'9669m_black', '9670f_green1'}, {'9669m_green1', '9670f_white'}, {'9669m_green2', '9670f_black'}, {'9669m_white', '9670f_green2'}};
L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string');
meta = jsondecode(fileread(fullfile(diagDir, 'meta.json')));
out = table();
for mi = 1:numel(modelDirs)
    for f = 1:4
        M = load(fullfile(modelDirs{mi}, sprintf('unet_ft_fold%d.mat', f)));
        for k = 1:numel(meta)
            key = string(meta(k).key); mouse = regexprep(key, '^day\d+_', '');
            if ~ismember(mouse, FOLDS{f}), continue; end
            Sd = load(fullfile(diagDir, [char(key) '.mat'])); Xs = Sd.samples;
            r0 = double(Sd.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
            F = arrayfun(@(j) Xs(:, :, j), 1:size(Xs, 3), 'UniformOutput', false)';
            C = load(fullfile(candDir, [char(key) '_cands.mat'])); allc = C.allc;
            for iris = [true false]
                [masks, ~, crops, ~, irisApplied, ~, probs] = pupil_mask_frame(F, rr, M, 'IrisCorrection', iris);
                for j = 1:numel(crops)
                    row = L(L.session == key & L.frame == j - 1, :);
                    if isempty(row) || isnan(row.label(1)), continue; end
                    lab = row.label(1); [h, w] = size(crops{j});
                    if lab >= 1, tgt = logical(allc{j}(lab).m); elseif lab == -1, tgt = false(h, w); else, tgt = []; end
                    [mm, info] = pupil_maxp_mask(masks{j}, probs{j}, 0.9);
                    [res, iou, dr] = judge(mm, tgt, lab);
                    out = [out; table(string(modelNames{mi}), f, key, mouse, j - 1, lab, row.how(1) == "fullpupil_rule", ...
                        row.how(1) == "agree" | row.how(1) == "fullpupil_rule", iris, irisApplied(j), res, iou, dr, ...
                        info.score, info.diam, 'VariableNames', {'model', 'fold', 'session', 'mouse', 'frame', 'label', ...
                        'full', 'agreed_or_full', 'iris', 'iris_applied', 'result', 'iou', 'diam_ratio', 'score', 'diam_px'})]; %#ok<AGROW>
                end
            end
        end
        fprintf('%s fold %d done\n', modelNames{mi}, f);
    end
end
writetable(out, outCsv);
for mi = 1:numel(modelNames)
    for iris = [true false]
        R = out(out.model == modelNames{mi} & out.iris == iris, :);
        Fu = R(R.full, :); Nn = R(~R.full & R.label >= 1, :); Ag = R(~R.full & R.agreed_or_full & R.label >= 1, :);
        Ne = R(~R.full & R.label == -1, :);
        fprintf('\n[%s, IrisCorrection %d]\n', modelNames{mi}, iris);
        fprintf('  full-pupil frames (%d): chosen blob / label diameter median %.2f (IQR %.2f-%.2f), IoU>=0.6 in %.0f %%, miss %d, iris cut %d\n', ...
            height(Fu), median(Fu.diam_ratio, 'omitnan'), prctile(Fu.diam_ratio, 25), prctile(Fu.diam_ratio, 75), ...
            100 * mean(Fu.iou >= 0.6), nnz(Fu.result == "miss"), nnz(Fu.iris_applied));
        fprintf('  other pupil frames (%d): correct %d (%.0f %%), wrong %d, miss %d; agreed (%d): correct %.0f %%, wrong %.0f %% of picks\n', ...
            height(Nn), nnz(Nn.result == "correct"), 100 * mean(Nn.result == "correct"), nnz(Nn.result == "wrong"), ...
            nnz(Nn.result == "miss"), height(Ag), 100 * mean(Ag.result == "correct"), ...
            100 * nnz(Ag.result == "wrong") / max(1, nnz(Ag.result ~= "miss")));
        fprintf('  not-visible frames (%d): false blob %d\n', height(Ne), nnz(Ne.result == "false"));
    end
end
end

function [res, iou, dr] = judge(m, tgt, lab)
iou = NaN; dr = NaN;
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
