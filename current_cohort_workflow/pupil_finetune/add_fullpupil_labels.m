function add_fullpupil_labels(candDir, fitDir, fullDir, outCandDir, outFitDir, exclude)
% add_fullpupil_labels  New label set in which the frames of fully dilated sessions are labelled with the full-pupil
% mask (make_fullpupil_labels.m). The original candidate and label folders are not changed:
%   outCandDir/<key>_cands.mat   = the original candidates, plus the full-pupil mask appended as the LAST candidate of
%                                  every frame of the sessions found in fullDir (other sessions copied unchanged)
%   outFitDir/labels_consensus.csv = the original consensus labels; for those frames label = that last candidate,
%                                  fit = good, how = fullpupil_rule (votes = old label)
%   exclude: {key, frame; ...} frames whose full-pupil mask is wrong (e.g. it leaks into dark fur) -> label NaN
%            (finetune_unet_cv skips them)
if nargin < 6, exclude = cell(0, 2); end
if ~exist(outCandDir, 'dir'), mkdir(outCandDir); end
if ~exist(outFitDir, 'dir'), mkdir(outFitDir); end
L = readtable(fullfile(fitDir, 'labels_consensus.csv'), 'TextType', 'string');
L.votes = string(L.votes); L.how = string(L.how); L.fit = string(L.fit); L.conf = string(L.conf);
files = dir(fullfile(candDir, '*_cands.mat'));
for f = 1:numel(files)
    key = erase(files(f).name, '_cands.mat');
    src = fullfile(candDir, files(f).name); dst = fullfile(outCandDir, files(f).name);
    fp = fullfile(fullDir, [key '_fullpupil.mat']);
    if ~isfile(fp), copyfile(src, dst); continue; end
    C = load(src); Fm = load(fp);
    allc = C.allc;
    tmpl = struct('m', []);                       % field names of the candidates (a frame can have none)
    for j = 1:numel(allc)
        if ~isempty(allc{j}), tmpl = allc{j}(1); break; end
    end
    fn = fieldnames(tmpl);
    for j = 1:numel(allc)
        s = allc{j};
        new = tmpl;
        for q = 1:numel(fn), new.(fn{q}) = []; end
        new.m = Fm.masks{j};
        if isempty(s), allc{j} = new; else, allc{j} = [s(:); new]'; end
        r = find(L.session == string(key) & L.frame == j - 1, 1);
        newLab = numel(allc{j});
        isEx = any(strcmp(exclude(:, 1), key) & cell2mat(exclude(:, 2)) == j - 1);
        if isempty(r)
            continue;
        end
        L.votes(r) = "old " + string(L.label(r));
        if isEx
            L.label(r) = NaN; L.how(r) = "fullpupil_excluded";
        else
            L.label(r) = newLab; L.fit(r) = "good"; L.how(r) = "fullpupil_rule";
        end
    end
    C.allc = allc; save(dst, '-struct', 'C');
    fprintf('%s: full-pupil mask added as the last candidate of %d frames\n', key, numel(allc));
end
writetable(L, fullfile(outFitDir, 'labels_consensus.csv'));
fprintf('labels: %d frames, %d full-pupil, %d excluded\n', height(L), nnz(L.how == "fullpupil_rule"), nnz(L.how == "fullpupil_excluded"));
end
