function run_step9_one_click
% run_step9_one_click — one-click Step 9: EFA (Step7) vs PR-trajectory (Step8) cluster overlap
%
% Requires:
%   1) Latest longitudinal run (same as Step 8): pr_step08_locate_data
%   2) Step 7: one of
%        run_*/figs/addiction_score_efa_mousefit/addiction_score_mousefit_by_mouse.csv
%        run_*/figs/**/AI_mouse_summary.csv
%      (must include mouse_key + AI_mouse or AI_mouse_z for k-means)
%   3) Step 8: run_pr_licking_full_analysis (USE_DEMO = false) — must write
%        step08_pr_trajectory_downstream/output_full/MouseMean_LickRate_by_PRLevel.csv
%      (same mouse_key as EFA / lick MEGA). Or set OVERRIDE_PR_MOUSEMEAN to that CSV.
%      PR clusters match run_pr_trajectory_phenotype_cluster (within-trial, K=2).
%
% Output: out_step9_auto/efa_kmeans_for_step9.csv, pr_trajectory_for_step9.csv,
%         and compare/ (ARI, NMI, confusion) from compare_trajectory_clusterings.m
%
% Edit K_EFA / K_PR / OVERRIDE_RUN below if needed.
%
% Where the inputs live (same run_* for Step 7 EFA, Step 4/8 lick MEGA, and Step 8 tables):
%   <run_*>/figs/lick_MEGA/per_trial_features.csv     ← Step 4/8 lick MEGA
%   <run_*>/figs/.../AI_mouse_summary.csv             ← Step 7 (first match below)
%   <run_*>/ALL_mice_longitudinal.csv
%   step08/.../output_full/MouseMean_LickRate_by_PRLevel.csv  ← written by run_pr_licking_full_analysis

K_EFA = 2;       % EFA: k-means on AI mouse score
K_PR  = 2;       % PR-trajectory k-means (match run_pr_trajectory_phenotype_cluster)
RATE_PR = 'within';
% If nonempty: same run_* you use in Step 8/9 (e.g. run_009 on K:). pr_step08_locate_data() otherwise picks LATEST.
OVERRIDE_RUN = 'K:\addiction_concate_Dec_2025\longitudinal_outputs\run_009';   % e.g. 'K:\addiction_concate_Dec_2025\longitudinal_outputs\run_009'
% If nonempty: full path to MouseMean_LickRate_by_PRLevel.csv (bypasses step08 output_full)
OVERRIDE_PR_MOUSEMEAN = '';
rng(1);

here = fileparts(mfilename('fullpath'));
rootM = fileparts(here);
step08 = fullfile(rootM, 'step08_pr_trajectory_downstream');
addpath(step08);

% --- run folder = same as Step 8 lick MEGA (must match EFA) ---
if isstring(OVERRIDE_RUN) && strlength(OVERRIDE_RUN) > 0
    OVERRIDE_RUN = char(OVERRIDE_RUN);
end
p = pr_step08_locate_data(OVERRIDE_RUN);
runDir = p.runDir;
fprintf('Step9 one-click: runDir = %s\n', runDir);
fprintf('  lick MEGA: %s\n  longitudinal: %s\n', p.trial_csv, p.long_csv);

% --- EFA / addiction index table ---
efaCandidates = {
    fullfile(runDir, 'figs', 'addiction_index_EFA', 'AI_mouse_summary.csv')
    fullfile(runDir, 'figs', 'addiction_score_efa_mousefit', 'addiction_score_mousefit_by_mouse.csv')
    fullfile(runDir, 'figs', 'addiction_MEGA', 'AI_mouse_summary.csv')
    };
efaPath = '';
for i = 1:numel(efaCandidates)
    if exist(efaCandidates{i}, 'file') == 2
        efaPath = efaCandidates{i};
        break
    end
end
if isempty(efaPath) || exist(efaPath, 'file') ~= 2
    error([ 'No EFA / addiction mouse summary found under:\n  %s\\figs\\…\n', ...
        'Run analyze_addiction_score_efa_mousefit (or compute_addiction_index_EFA_mousefit_*) first.' ], runDir);
end
fprintf('EFA source: %s\n', efaPath);
Tefa = readtable(efaPath);
Tefa.mouse_key = cellstr(step09_normalize_mouse_key(Tefa.mouse_key));
if ismember('AI_mouse', Tefa.Properties.VariableNames)
    Xefa = double(Tefa.AI_mouse);
elseif ismember('AI_mouse_z', Tefa.Properties.VariableNames)
    Xefa = double(Tefa.AI_mouse_z);
else
    error('Table needs AI_mouse or AI_mouse_z. Columns: %s', strjoin(Tefa.Properties.VariableNames, ', '));
end
ok = isfinite(Xefa);
Tefa = Tefa(ok, :);
Xefa = Xefa(ok);
if numel(Xefa) < 4
    error('Too few mice with finite AI for k-means.');
end
try
    cEfa = kmeans(Xefa(:), K_EFA, 'Replicates', 20, 'MaxIter', 1000, 'Display', 'off', 'Start', 'plus');
catch
    cEfa = kmeans(Xefa(:), K_EFA, 'Replicates', 20, 'MaxIter', 1000, 'Display', 'off');
end
Eout = table(Tefa.mouse_key, cEfa, 'VariableNames', {'mouse_key', 'cluster'});

% --- PR trajectory: MouseMean from Step 8 (same mouse_key as lick MEGA / EFA cohort) ---
if isstring(OVERRIDE_PR_MOUSEMEAN) && strlength(OVERRIDE_PR_MOUSEMEAN) > 0
    OVERRIDE_PR_MOUSEMEAN = char(OVERRIDE_PR_MOUSEMEAN);
end
if ~isempty(OVERRIDE_PR_MOUSEMEAN)
    mm = char(OVERRIDE_PR_MOUSEMEAN);
    if exist(mm, 'file') ~= 2
        error('OVERRIDE_PR_MOUSEMEAN not found: %s', mm);
    end
else
    mm = fullfile(step08, 'output_full', 'MouseMean_LickRate_by_PRLevel.csv');
    if exist(mm, 'file') ~= 2
        enr = fullfile(step08, 'output_full', 'per_trial_enriched_step08.csv');
        if exist(enr, 'file') == 2
            try
                fprintf(1, 'Step9: building MouseMean from per_trial_enriched (pr_step08_rebuild_mousemean_from_enriched)...\n');
                pr_step08_rebuild_mousemean_from_enriched(enr, fileparts(mm));
            catch ME
                warning('Step9: auto-build of MouseMean failed: %s', ME.message);
            end
        end
    end
    if exist(mm, 'file') ~= 2
        enr = fullfile(step08, 'output_full', 'per_trial_enriched_step08.csv');
        hint = '';
        if exist(enr, 'file') == 2
            hint = sprintf([ '\n  Found %s — you can also run: pr_step08_rebuild_mousemean_from_enriched\n' ...
                '  or run_pr_licking_full_analysis (USE_DEMO = false, OVERRIDE_RUN if needed, SKIP_PLOTS = true).\n' ], enr);
        else
            hint = sprintf([ '\n  (No per_trial_enriched: run run_pr_licking_full_analysis, USE_DEMO = false, ', ...
                'OVERRIDE_RUN matching runDir.)\n' ]);
        end
        error([ 'Missing:\n  %s\n', ...
            '  Build with: pr_step08_rebuild_mousemean_from_enriched or run_pr_licking_full_analysis.', ...
            '  Or set OVERRIDE_PR_MOUSEMEAN to the CSV. Do not use output_demo (demo m001, … vs EFA).\n', ...
            '%s' ], mm, hint);
    end
end
fprintf('PR lick table: %s\n', mm);
S = readtable(mm);
S.mouse_key = cellstr(step09_normalize_mouse_key(S.mouse_key));
[mkPr, idPr] = pr_step09_build_traj_cluster(S, K_PR, RATE_PR);
PrTab = table(mkPr, idPr, 'VariableNames', {'mouse_key', 'trajectory_cluster'});

% --- work dir ---
outA = fullfile(here, 'out_step9_auto');
if exist(outA, 'dir') ~= 7, mkdir(outA); end
efaOut = fullfile(outA, 'efa_kmeans_for_step9.csv');
prOut  = fullfile(outA, 'pr_trajectory_for_step9.csv');
writetable(Eout, efaOut);
writetable(PrTab, prOut);
fprintf('Wrote: %s\nWrote: %s\n', efaOut, prOut);

% --- compare ---
outCmp = fullfile(outA, 'compare');
compare_trajectory_clusterings(efaOut, prOut, 'outDir', outCmp, 'efaLabelVar', 'cluster', 'prLabelVar', 'trajectory_cluster');
fprintf('\nStep 9 one-click done. See: %s\n', outCmp);
end

% ----- local: same k-means layout as run_pr_trajectory_phenotype_cluster.m -----
function [mice, idx] = pr_step09_build_traj_cluster(S, K, RATE_KIND)
mice = cellstr(unique(S.mouse_key));
tr = unique(S.Trial);
nM = numel(mice);
nT = numel(tr);
Xw = nan(nM, nT);
Xb = nan(nM, nT);
for i = 1:nM
    for j = 1:nT
        r = strcmp(S.mouse_key, mice{i}) & S.Trial == tr(j);
        if any(r) && ismember('MeanWithinTrialLickRate', S.Properties.VariableNames)
            Xw(i, j) = S.MeanWithinTrialLickRate(find(r, 1, 'first'));
        end
        if any(r) && ismember('MeanBetweenTrialLickRate', S.Properties.VariableNames)
            Xb(i, j) = S.MeanBetweenTrialLickRate(find(r, 1, 'first'));
        end
    end
end
if strcmpi(RATE_KIND, 'within')
    X = Xw;
elseif strcmpi(RATE_KIND, 'between')
    X = Xb;
    if size(X, 2) > 1
        X = X(:, 2:end);
    end
else
    X = [Xw, Xb(:, 2:end)];
end
mu = mean(X, 1, 'omitnan');
for c = 1:size(X, 2)
    ix = isnan(X(:, c));
    if any(ix)
        X(ix, c) = mu(c);
    end
end
Z = (X - mean(X, 2, 'omitnan')) ./ (std(X, 0, 2, 'omitnan') + eps);
try
    idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 1000, 'Display', 'off', 'Start', 'plus');
catch
    idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 1000, 'Display', 'off');
end
mice = mice(:);
idx = idx(:);
end
