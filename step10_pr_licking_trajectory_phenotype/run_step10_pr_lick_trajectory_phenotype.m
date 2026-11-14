function run_step10_pr_lick_trajectory_phenotype
% run_step10_pr_lick_trajectory_phenotype — PR-based licking trajectory phenotyping (no EFA)
%
% This step benchmarks a paper-style *trajectory* approach using only your PR/licking
% design (no lever-press SA). "Trajectory" = same behavioral features measured per mouse
% across multiple days/sessions. Core idea (adapted to your metrics):
%   - Effort / motivation: PR breakpoint (max PR trial; proxy for how far animals push the ratio)
%   - Consummatory / reward engagement: within-trial lick rate
%   - Persistent / inter-trial seeking: between-trial lick rate (optional ITI metrics later)
%
% A suggested Methods sentence (edit as you like):
%  Because we do not have lever-press self-administration, we adapted group-based
%  trajectory analysis to PR-derived licking measures (breakpoint, within- and
%  between-trial lick rate), identifying mice with distinct longitudinal patterns of
%  motivation- and reward-related licking.
%
% Inputs: Step 8  output  per_trial_enriched_step08.csv  (or OVERRIDE_ENRICHED)
%         Falls back to  MouseMean_LickRate_by_PRLevel.csv  if there is not enough
%         longitudinal day structure (one row per mouse with 3 features).
%
% Output:  out_step10/  (B: full day columns; Passive 6-10 lick rules in X) and, by default,
%          out_step10_sens_excl6_10/  (A: k-means input drops days 6-10; same dsum) for sensitivity

K_TRAJ = 2;              % 2-3 for small N
MIN_DAYS = 2;            % if fewer distinct days, use mouse-level 3-D features instead of day-trajectory
R_CORR_DROP = 0.8;       % if |r| between two features is >= this on mouse×day table, print warning
PR_CAP_W = 50;           % Hz, within-trial
PR_CAP_B = 22;           % between-trial (Hz)
PR_MIN_TTR = 0.02;
% PR sessions start at day 3 in your design; trajectory analysis starts at day 4.
DAY_INDEX_MIN = 4;
% Passive: no licking in Period "During" for these calendar days — trial lick fields -> NaN, and
% mouse×day lick means (within/between) for Passive on these days are forced to NaN for analysis.
% Set to [] to disable. (Lick-based metrics nulled for Passive During on these days.)
% PR breakpoint: also nulled in mouse×day dsum and omitted from k-means on these days for *all* groups,
% so Active PR on 6-10 does not drive clustering while licks are invalid/absent in this window.
% Set to [] to keep PR on all calendar days in X.
PR_EXCLUDE_DAYS = 6:10;
PASSIVE_DURING_LICK_DAYS   = 6:10;    % set to []: no lick nans; use full lick to see if clusters still separate
PASSIVE_DURING_PERIOD_NAME = 'During';
% (A) vs (B) sensitivity: run k-means twice — same dsum/Passive rules, but (A) drops
% all trajectory columns for these calendar days from the clustering design matrix
% (no 5-feature×day columns for 6–10, so no NaN-structure from Passive During 6–10 in X).
% (B) uses the full set of day columns with Passive 6–10 lick features NaN + omit rules.
KMEANS_EXCLUDE_DAYS_SENS = 6:10;   % [] to disable the second (A) output folder
% Clustering features (day-traj: each × calendar day; order = block order in X)
% 1=PR max trial, 2=total lick count, 3=within Hz, 4=between Hz, 5=mean 1/IEI (Hz-like, within-trial IEI)
STEP10_TRAJ_FNS = {'pr_breakpoint', 'total_lick_count', 'mean_within_hz', 'mean_between_hz', 'mean_iti_lick_freq_hz'};

%% paths
here = fileparts(mfilename('fullpath'));
rootM = fileparts(here);
step08 = fullfile(rootM, 'step08_pr_trajectory_downstream');
addpath(step08);
outDirB = fullfile(here, 'out_step10');
if exist(outDirB, 'dir') ~= 7, mkdir(outDirB); end
outDirA = fullfile(here, 'out_step10_sens_excl6_10');

OVERRIDE_RUN = [];
OVERRIDE_ENRICHED = fullfile(step08, 'output_full', 'per_trial_enriched_step08.csv');
OVERRIDE_MOUSEMEAN = fullfile(step08, 'output_full', 'MouseMean_LickRate_by_PRLevel.csv');

% --- load enriched
enrP = pr_step10_resolve_enrich_path(OVERRIDE_ENRICHED, step08);
E = readtable(enrP);
E = pr_step10_coerce_day_index(E);
E = pr_step10_sanitize_enriched(E, PR_MIN_TTR, PR_CAP_W, PR_CAP_B);
E = pr_step10_naan_passive_during_lick_trials(E, PASSIVE_DURING_LICK_DAYS, PASSIVE_DURING_PERIOD_NAME);
if ismember('day_index', E.Properties.VariableNames)
    n0 = height(E);
    E = E(double(E.day_index) >= DAY_INDEX_MIN, :);
    fprintf('Step10: day_index >= %d only (PR from day 3; analysis from day 4): %d -> %d rows.\n', ...
        DAY_INDEX_MIN, n0, height(E));
else
    warning('step10:dayIndex:missing', ['Column day_index not found: cannot filter to day>=%d. ', ...
        'Using all rows (add day_index in Step8 export if you need day-4+ only).'], DAY_INDEX_MIN);
end
if height(E) < 1
    error('No rows after day filter (DAY_INDEX_MIN=%d). Check per_trial_enriched_step08.csv.', DAY_INDEX_MIN);
end
Tlong = table; %#ok<NASGU>
dsum0 = table; %#ok<NASGU>
if ismember('day_index', E.Properties.VariableNames) && ismember('Trial', E.Properties.VariableNames) && height(E) > 0
    dsum0 = pr_step10_aggregate_mouse_day(E);
    n1 = height(dsum0);
    dsum0 = dsum0(isfinite(dsum0.day_index) & (double(dsum0.day_index) >= DAY_INDEX_MIN), :);
    if height(dsum0) < n1
        fprintf('Step10: dsum (mouse x day) dropped %d row(s) with day_index < %d or invalid.\n', n1 - height(dsum0), DAY_INDEX_MIN);
    end
    dsum0 = pr_step10_naan_dsum_passive_lick_by_day(dsum0, E, PASSIVE_DURING_LICK_DAYS);
    dsum0 = pr_step10_naan_dsum_pr_breakpoint_by_day(dsum0, PR_EXCLUDE_DAYS);
    Tlong = dsum0;
end

runList = {outDirB, '[]', ' (B) full 5x days: Passive 6-10 lick NaN/omit; PR 6-10 NaN+omit (all groups)'};
if ~isempty(KMEANS_EXCLUDE_DAYS_SENS)
    if exist(outDirA, 'dir') ~= 7, mkdir(outDirA); end
    runList(end+1, 1:3) = {outDirA, mat2str(KMEANS_EXCLUDE_DAYS_SENS(:)'), ' (A) k-means design matrix excludes days 6-10 (same dsum)'}; %#ok<AGROW>
end
for ir = 1:size(runList, 1)
    outDir = char(runList{ir, 1});
    exclS = strtrim(char(runList{ir, 2}));
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    if strcmp(exclS, '[]') || strlength(exclS) < 1
        exclD = [];
    else
        exclD = str2num(exclS); %#ok<ST2NM> string from mat2str
    end
    dsum = dsum0;
    if height(dsum) < 1
        mks = cellstr(unique(E.mouse_key, 'stable'));
        X = pr_step10_mouse_level_from_enriched(E, mks);
        fLabs = {'maxPR_dge4', 'totLick_dge4', 'meanWithin_dge4', 'meanBetween_dge4', 'itiLickFrq_dge4'};
        mode = 'mouse3';
        fprintf('Run %d %s: dsum empty — mouse-level 5 features from enriched (day>=%d).\n', ir, strtrim(char(runList{ir,3})), DAY_INDEX_MIN);
    else
        writetable(dsum, fullfile(outDir, 'pr_licking_mouse_by_day.csv'));
        pr_step10_feature_correlation_heatmap(dsum, R_CORR_DROP, outDir, DAY_INDEX_MIN);
        dT = sort(unique(double(dsum.day_index(:))));
        dT = dT(isfinite(dT) & dT >= DAY_INDEX_MIN);
        if ~isempty(exclD)
            dT = dT(~ismember(dT, exclD(:).'));
        end
        mode = 'mouse3';
        if numel(dT) >= MIN_DAYS
            [X, mks, fLabs] = pr_step10_build_traj_matrix_from_day(dsum, STEP10_TRAJ_FNS, DAY_INDEX_MIN, exclD);
            if ~isempty(X) && size(X, 1) >= 3
                mode = 'daytraj';
            end
        else
            fprintf('Run %d: only %d day(s) after excl %s; using mouse-level 5 features.\n', ir, numel(dT), mat2str(exclD));
        end
        if strcmp(mode, 'mouse3')
            mks = cellstr(unique(E.mouse_key, 'stable'));
            X = pr_step10_mouse_level_from_enriched(E, mks);
            fLabs = {'maxPR_dge4', 'totLick_dge4', 'meanWithin_dge4', 'meanBetween_dge4', 'itiLickFrq_dge4'};
            fprintf('Run %d %s: mouse-level 5 features from enriched (day>=%d).\n', ir, strtrim(char(runList{ir,3})), DAY_INDEX_MIN);
        else
            if isempty(exclD)
                exmsg = 'none';
            else
                exmsg = mat2str(exclD(:).');
            end
            fprintf('Run %d %s: day-trajectory stack (K-means excl days = %s).\n', ir, strtrim(char(runList{ir,3})), exmsg);
        end
    end
    gMouse = pr_step10_group_per_mouse(E, mks);
    omitX = false(size(X, 1), size(X, 2));
    if exist('fLabs', 'var') && iscell(fLabs) && numel(fLabs) == size(X, 2) && strcmpi(mode, 'daytraj')
        if ~isempty(PASSIVE_DURING_LICK_DAYS)
            omitX = pr_step10_omit_passive_lick_columns(fLabs, mks, gMouse, PASSIVE_DURING_LICK_DAYS);
        end
        if ~isempty(PR_EXCLUDE_DAYS)
            omitX = pr_step10_omit_pr_expday_columns(fLabs, mks, PR_EXCLUDE_DAYS, omitX);
        end
    end
    X = pr_step10_impute_col_means(X, omitX);
    Z = pr_step10_zrows(X);
    if size(Z, 1) < 3
        error('Not enough subjects for k-means (n=%d).', size(Z,1));
    end
    K = pr_step10_choose_k(K_TRAJ, size(Z, 1));
    try
        idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 2000, 'Display', 'off', 'Start', 'plus');
    catch
        idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 2000, 'Display', 'off');
    end
    Out = table(mks(:), idx, gMouse(:), 'VariableNames', {'mouse_key', 'trajectory_phenotype_cluster', 'Group'});
    writetable(Out, fullfile(outDir, 'pr_licking_trajectory_clusters.csv'));
    fprintf('Wrote: %s  (K=%d, n=%d)%s\n', fullfile(outDir, 'pr_licking_trajectory_clusters.csv'), K, numel(mks), strtrim(char(runList{ir,3})));
    pr_step10_crosstab_group_cluster(Out, outDir);
    pr_step10_figure_group_cluster_heatmap(Out, outDir);
    pr_step10_figure_radar_or_bar(X, fLabs, idx, mks, outDir, DAY_INDEX_MIN);
    if height(dsum) >= 1 && exist('Tlong', 'var') && istable(Tlong) && height(Tlong) > 0 && strcmpi(mode, 'daytraj')
        pr_step10_figure_cluster_means_by_day(dsum, Out, outDir, DAY_INDEX_MIN);
        pr_step10_spaghetti_breakpoint(Tlong, Out, outDir, DAY_INDEX_MIN);
    end
    mDoc = fullfile(here, 'docs', 'METHODS_step10.md');
    if exist(mDoc, 'file') == 2
        try
            copyfile(mDoc, fullfile(outDir, 'METHODS_step10.md'), 'f');
            fprintf('Methods: copied to %s\n', fullfile(outDir, 'METHODS_step10.md'));
        catch ME
            warning('step10:methodsCopy', 'Could not copy METHODS for run %d: %s', ir, ME.message);
        end
    end
    fprintf('Step 10 run %d done%s. No EFA used. -> %s\n', ir, strtrim(char(runList{ir,3})), outDir);
end
if isempty(KMEANS_EXCLUDE_DAYS_SENS)
    fprintf('Step 10: primary output %s (sensitivity run disabled: KMEANS_EXCLUDE_DAYS_SENS=[]).\n', outDirB);
else
    fprintf('Step 10: primary (B) %s; sensitivity (A) %s\n', outDirB, outDirA);
end
end

%% -------------------------------------------------------------------------
function g = pr_step10_group_per_mouse(E, mks)
% One label per mouse from Step 8 enriched `Group` (Active / Passive / …).
g = repmat({'Unknown'}, numel(mks), 1);
if ~ismember('Group', E.Properties.VariableNames)
    return
end
for i = 1:numel(mks)
    r = strcmp(E.mouse_key, mks{i});
    if ~any(r)
        g{i} = 'Unknown';
        continue
    end
    gv = E.Group(r);
    % isfinite() is numeric-only; Group is text -> use strlength (0 for "" and for <missing> string in recent MATLAB)
    s = strtrim(string(gv(:)));
    s = s(strlength(s) > 0);
    if isempty(s)
        g{i} = 'Unknown';
        continue
    end
    uq = unique(s, 'stable');
    if numel(uq) > 1
        warning('step10:group:mixed', 'Mouse %s has multiple Group labels; using first: %s', mks{i}, char(uq(1)));
    end
    g{i} = char(strtrim(uq(1)));
end
end

function pr_step10_crosstab_group_cluster(Out, outDir)
if ~ismember('Group', Out.Properties.VariableNames)
    return
end
C = double(Out.trajectory_phenotype_cluster);
uC = unique(C, 'sorted');
% Stable group order: Active, Passive, then the rest
allG = cellstr(string(Out.Group));
uG0 = unique(allG, 'stable');
pref = {'Active', 'Passive'};
uG = {};
for p = 1:numel(pref)
    k = strcmpi(uG0, pref{p});
    if any(k)
        uG{end+1} = uG0{find(k, 1)}; %#ok<AGROW>
    end
end
for j = 1:numel(uG0)
    if ~any(strcmpi(uG, uG0{j}))
        uG{end+1} = uG0{j}; %#ok<AGROW>
    end
end
M = zeros(numel(uC), numel(uG));
for i = 1:numel(uC)
    for j = 1:numel(uG)
        M(i, j) = sum(C == uC(i) & strcmp(allG, uG{j}));
    end
end
rown = arrayfun(@(c) sprintf('cluster_%d', c), uC, 'uni', 0);
Twide = array2table(M, 'VariableNames', matlab.lang.makeValidName(strrep(uG, ' ', '_')), 'RowNames', rown);
writetable(Twide, fullfile(outDir, 'pr_cluster_x_group_matrix.csv'), 'WriteRowNames', true);
fprintf('Cluster x Group (counts, rows=cluster, columns=Group as in CSV variable names):\n');
disp(Twide);
nL = numel(uC) * numel(uG);
cluCol = zeros(nL, 1);
gCol = cell(nL, 1);
cntCol = zeros(nL, 1);
r0 = 0;
for i = 1:numel(uC)
    for j = 1:numel(uG)
        r0 = r0 + 1;
        cluCol(r0) = uC(i);
        gCol{r0} = uG{j};
        cntCol(r0) = M(i, j);
    end
end
Tlong = table(cluCol, gCol, cntCol, 'VariableNames', {'trajectory_phenotype_cluster', 'Group', 'Count'});
writetable(Tlong, fullfile(outDir, 'pr_cluster_x_group_long.csv'));
fprintf('Wrote: %s, %s\n', fullfile(outDir, 'pr_cluster_x_group_matrix.csv'), fullfile(outDir, 'pr_cluster_x_group_long.csv'));
for j = 1:numel(uG)
    gnm = strtrim(uG{j});
    if any(strcmpi(gnm, {'active', 'passive'}))
        idxm = strcmpi(allG, uG{j});
        cc = C(idxm);
        fprintf('%s (n=%d):  ', gnm, numel(cc));
        for k = 1:numel(uC)
            fprintf('cluster %d = %d  ', uC(k), sum(cc == uC(k)));
        end
        fprintf('\n');
    end
end
end

function p = pr_step10_resolve_enrich_path(ovEnr, step08)
p = char(ovEnr);
if exist(p, 'file') ~= 2
    p = fullfile(step08, 'output_full', 'per_trial_enriched_step08.csv');
end
if exist(p, 'file') ~= 2
    error('Missing per_trial_enriched_step08.csv under step08 output_full. Run run_pr_licking_full_analysis.');
end
end

function p = pr_step10_resolve_mm(ov, step08)
p = char(ov);
if exist(p, 'file') ~= 2
    p = fullfile(step08, 'output_full', 'MouseMean_LickRate_by_PRLevel.csv');
end
if exist(p, 'file') ~= 2
    error('Missing MouseMean_LickRate_by_PRLevel.csv. Run run_pr_licking_full_analysis or pr_step08_rebuild_mousemean_from_enriched.');
end
end

function k = pr_step10_choose_k(kWant, n)
if n < 3
    error('Step10: need at least 3 subjects (n=%d).', n);
end
k = min(kWant, n - 1);
k = max(2, k);
k = min(k, n - 1);
end

function X = pr_step10_impute_col_means(X, omit)
% Column means; skip imputation for cells with omit==true (leave NaN for z-scoring to zero later).
if nargin < 2 || isempty(omit)
    omit = false(size(X, 1), size(X, 2));
end
X = double(X);
if ~isequal(size(omit), size(X))
    omit = false(size(X, 1), size(X, 2));
end
for c = 1:size(X, 2)
    col = X(:, c);
    om = omit(:, c);
    take = isfinite(col) & ~om;
    mu = mean(col(take), 'omitnan');
    if ~isfinite(mu)
        mu = 0;
    end
    for r = 1:size(X, 1)
        if ~isfinite(col(r)) && ~om(r)
            col(r) = mu;
        end
    end
    X(:, c) = col;
end
end

function E = pr_step10_sanitize_enriched(E, minTtr, capW, capB)
if ismember('time_to_reward', E.Properties.VariableNames) && ismember('lick_count', E.Properties.VariableNames)
    ttr = max(double(E.time_to_reward), minTtr);
    w = double(E.lick_count) ./ ttr;
elseif ismember('lick_rate', E.Properties.VariableNames)
    w = double(E.lick_rate);
else
    error('enriched: need lick_count+time_to_reward or lick_rate');
end
w = w(:);
w(isfinite(w) & (w < 0 | w > capW)) = NaN;
E.WithinTrialLickRate = w;
if ismember('BetweenTrialLickRate', E.Properties.VariableNames)
    b = double(E.BetweenTrialLickRate(:));
    b(isfinite(b) & (b < 0 | b > capB)) = NaN;
    E.BetweenTrialLickRate = b;
end
E.mouse_key = cellstr(string(E.mouse_key));
if ismember('day_index', E.Properties.VariableNames)
    E.day_index = double(E.day_index);
end
E.Trial = double(E.Trial);
end

function E = pr_step10_coerce_day_index(E)
% Ensure day_index is double so the day>=4 filter cannot silently fail (text / categorical from CSV).
if ~ismember('day_index', E.Properties.VariableNames)
    return
end
v = E.day_index;
if isnumeric(v)
    di = double(v(:));
elseif iscategorical(v)
    di = double(v(:));
elseif iscell(v) || isstring(v)
    di = str2double(string(v(:)));
else
    di = str2double(string(v(:)));
end
E.day_index = di;
end

function E = pr_step10_naan_passive_during_lick_trials(E, dayRng, periodName)
if isempty(dayRng) || ~ismember('day_index', E.Properties.VariableNames) || ~ismember('Group', E.Properties.VariableNames)
    return
end
if ~ismember('Period', E.Properties.VariableNames) && ~ismember('PeriodClean', E.Properties.VariableNames)
    return
end
n0 = 0;
perCol = 'Period';
if ~ismember('Period', E.Properties.VariableNames)
    perCol = 'PeriodClean';
end
for r = 1:height(E)
    if ~strcmpi(strtrim(char(string(E.Group(r)))), 'passive')
        continue
    end
    d = double(E.day_index(r));
    if ~isfinite(d) || ~ismember(d, dayRng(:))
        continue
    end
    p = strtrim(char(string(E.(perCol)(r))));
    if ~strcmpi(p, periodName)
        continue
    end
    n0 = n0 + 1;
    if ismember('WithinTrialLickRate', E.Properties.VariableNames)
        E.WithinTrialLickRate(r) = NaN;
    end
    if ismember('BetweenTrialLickRate', E.Properties.VariableNames)
        E.BetweenTrialLickRate(r) = NaN;
    end
    if ismember('lick_rate', E.Properties.VariableNames)
        E.lick_rate(r) = NaN;
    end
    if ismember('lick_count', E.Properties.VariableNames)
        E.lick_count(r) = NaN;
    end
    if ismember('iei_median', E.Properties.VariableNames)
        E.iei_median(r) = NaN;
    end
end
if n0 > 0
    fprintf('Step10: nanned lick fields on %d trial row(s) (Passive, %s, day in %s).\n', n0, periodName, mat2str(dayRng));
end
end

function dsum = pr_step10_naan_dsum_passive_lick_by_day(dsum, E, dayRng)
% Force mouse×day within/between to NaN for Passive on these calendar days (analysis rule).
if isempty(dayRng) || isempty(dsum) || height(dsum) < 1 || ~ismember('Group', E.Properties.VariableNames)
    return
end
mk2g = pr_step10_mouse_group_map_from_E(E);
n0 = 0;
for r = 1:height(dsum)
    mk = char(strtrim(string(dsum.mouse_key(r))));
    d = double(dsum.day_index(r));
    if ~isKey(mk2g, mk) || ~strcmpi(mk2g(mk), 'passive')
        continue
    end
    if ~isfinite(d) || ~ismember(d, dayRng(:))
        continue
    end
    n0 = n0 + 1;
    dsum.mean_within_hz(r) = NaN;
    dsum.mean_between_hz(r) = NaN;
    if ismember('total_lick_count', dsum.Properties.VariableNames)
        dsum.total_lick_count(r) = NaN;
    end
    if ismember('mean_iti_lick_freq_hz', dsum.Properties.VariableNames)
        dsum.mean_iti_lick_freq_hz(r) = NaN;
    end
end
if n0 > 0
    fprintf('Step10: nanned lick metrics (within/between/total/ITIfreq) on %d mouse×day row(s) (Passive, day in %s).\n', n0, mat2str(dayRng));
end
end

function dsum = pr_step10_naan_dsum_pr_breakpoint_by_day(dsum, dayRng)
% Set pr_breakpoint to NaN for every mouse×day row in dayRng (all groups) for trajectory k-means.
% Pair with pr_step10_omit_pr_expday_columns for the design matrix.
if isempty(dayRng) || isempty(dsum) || height(dsum) < 1
    return
end
if ~ismember('pr_breakpoint', dsum.Properties.VariableNames)
    return
end
n0 = 0;
for r = 1:height(dsum)
    d = double(dsum.day_index(r));
    if ~isfinite(d) || ~ismember(d, dayRng(:))
        continue
    end
    dsum.pr_breakpoint(r) = NaN;
    n0 = n0 + 1;
end
if n0 > 0
    fprintf('Step10: nanned pr_breakpoint on %d mouse×day row(s) (all groups, day in %s).\n', n0, mat2str(dayRng));
end
end

function M = pr_step10_mouse_group_map_from_E(E)
M = containers.Map('KeyType', 'char', 'ValueType', 'char');
if ~ismember('Group', E.Properties.VariableNames)
    return
end
uK = unique(cellstr(string(E.mouse_key)), 'stable');
for i = 1:numel(uK)
    f = find(strcmp(cellstr(string(E.mouse_key)), uK{i}), 1, 'first');
    M(uK{i}) = char(strtrim(string(E.Group(f))));
end
end

function omit = pr_step10_omit_passive_lick_columns(fLabs, mks, gPer, dayRng)
omit = false(numel(mks), numel(fLabs));
if isempty(dayRng) || isempty(fLabs)
    return
end
for j = 1:numel(fLabs)
    f = fLabs{j};
    isW = ~isempty(strfind(f, 'WithinHz')); %#ok<STREMP> 
    isB = ~isempty(strfind(f, 'BetweenHz')); %#ok<STREMP>
    isT = ~isempty(strfind(f, 'TotalL')); %#ok<STREMP>
    isI = ~isempty(strfind(f, 'ITIfr')); %#ok<STREMP>
    if ~(isW || isB || isT || isI) || isempty(strfind(f, 'expDay')) %#ok<STREMP>
        continue
    end
    tok = regexp(f, 'expDay([0-9.]+)', 'tokens', 'once');
    if isempty(tok)
        continue
    end
    c1 = tok{1};
    if iscell(c1) && numel(c1) >= 1
        sday = c1{1};
    else
        sday = c1;
    end
    d = str2double(sday);
    if ~isfinite(d) || ~ismember(d, dayRng(:))
        continue
    end
    for i = 1:min(numel(mks), numel(gPer))
        if strcmpi(strtrim(char(string(gPer{i}))), 'passive')
            omit(i, j) = true;
        end
    end
end
end

function omit = pr_step10_omit_pr_expday_columns(fLabs, mks, dayRng, omit0)
% Mark PR expDay* columns in fLabs for omit imputation: all rows, for calendar days in dayRng.
omit = false(numel(mks), numel(fLabs));
if nargin >= 4 && ~isempty(omit0) && isequal(size(omit0), [numel(mks), numel(fLabs)])
    omit = omit0;
end
if isempty(dayRng) || isempty(fLabs)
    return
end
for j = 1:numel(fLabs)
    f = fLabs{j};
    if isempty(regexp(f, '^PR\s+expDay', 'once'))
        continue
    end
    tok = regexp(f, 'expDay([0-9.]+)', 'tokens', 'once');
    if isempty(tok)
        continue
    end
    c1 = tok{1};
    if iscell(c1) && numel(c1) >= 1
        sday = c1{1};
    else
        sday = c1;
    end
    d = str2double(sday);
    if ~isfinite(d) || ~ismember(d, dayRng(:))
        continue
    end
    for i = 1:numel(mks)
        omit(i, j) = true;
    end
end
end

function dsum = pr_step10_aggregate_mouse_day(E)
G = findgroups(E.mouse_key, E.day_index);
dsum = table;
dsum.mouse_key = cellstr(splitapply(@(x) x(1), E.mouse_key, G));
dsum.day_index = splitapply(@(x) x(1), E.day_index, G);
dsum.pr_breakpoint = splitapply(@(t) max(double(t), [], 'omitnan'), E.Trial, G);
dsum.mean_within_hz = splitapply(@(w) mean(double(w), 'omitnan'), E.WithinTrialLickRate, G);
if ismember('BetweenTrialLickRate', E.Properties.VariableNames)
    dsum.mean_between_hz = splitapply(@(b) mean(double(b), 'omitnan'), E.BetweenTrialLickRate, G);
else
    dsum.mean_between_hz = nan(height(dsum), 1);
end
if ismember('lick_count', E.Properties.VariableNames)
    dsum.total_lick_count = splitapply(@(x) pr_step10_sum_lick_count(x), E.lick_count, G);
else
    dsum.total_lick_count = nan(height(dsum), 1);
end
if ismember('iei_median', E.Properties.VariableNames)
    dsum.mean_iti_lick_freq_hz = splitapply(@(ie) pr_step10_mean_iti_lick_freq_hz(ie), E.iei_median, G);
else
    dsum.mean_iti_lick_freq_hz = nan(height(dsum), 1);
end
end

function s = pr_step10_sum_lick_count(x)
% Sum licks for (mouse,day); if every trial is NaN, return NaN (not 0 from sum('omitnan')).
x = double(x(:));
if isempty(x) || all(isnan(x)) || all(~isfinite(x))
    s = NaN;
else
    s = sum(x, 'omitnan');
end
end

function v = pr_step10_mean_iti_lick_freq_hz(ie)
% Within-trial: mean(1/IEI) in Hz-like units (IEI in s). NaN if no valid IEI.
ie = double(ie(:));
ie = ie(isfinite(ie) & ie > 1e-6);
if isempty(ie)
    v = NaN;
else
    v = mean(1 ./ ie, 'omitnan');
end
end

function [X, mks, fLabs] = pr_step10_build_traj_matrix_from_day(dsum, fns, dayMin, excludeDays)
% Each row = one mouse, columns = stacked [feat1 all days, feat2 all days, ...] in calendar day order.
% dList = only days >= dayMin (again), so the first column is never "experiment day 1" — label uses expDay<N>.
% excludeDays: optional. If non-empty, these calendar day indices are omitted from dList (sensitivity: no
% 5x-feature columns for those days, so Passive-missing-lick structure cannot drive k-means there).
if nargin < 3, dayMin = 4; end
if nargin < 4, excludeDays = []; end
mks = cellstr(unique(dsum.mouse_key, 'stable'));
dList = sort(unique(double(dsum.day_index(:))));
dList = dList(isfinite(dList) & dList >= dayMin);
if ~isempty(excludeDays)
    exd = double(excludeDays(:)).';
    dList = dList(~ismember(dList, exd));
end
if isempty(dList)
    error('step10:buildTraj: no day_index values >= %d in mouse-day table (check CSV day_index).', dayMin);
end
T = numel(dList);
F = numel(fns);
X = nan(numel(mks), T * F);
for im = 1:numel(mks)
    dsub = dsum(strcmp(dsum.mouse_key, mks{im}), :);
    for it = 1:T
        r = abs(double(dsub.day_index) - dList(it)) < 1e-9;
        for ifn = 1:F
            c = fns{ifn};
            v = dsub.(c);
            col = (ifn-1) * T + it;
            if any(r) && ismember(c, dsub.Properties.VariableNames)
                X(im, col) = mean(v(r), 'omitnan');
            end
        end
    end
end
fLabs = cell(1, T*F);
ftag = cell(size(fns));
for ifn = 1:F
    switch fns{ifn}
        case 'pr_breakpoint', ftag{ifn} = 'PR';
        case 'total_lick_count', ftag{ifn} = 'TotalL';
        case 'mean_within_hz', ftag{ifn} = 'WithinHz';
        case 'mean_between_hz', ftag{ifn} = 'BetweenHz';
        case 'mean_iti_lick_freq_hz', ftag{ifn} = 'ITIfr';
        otherwise, ftag{ifn} = fns{ifn};
    end
    for it = 1:T
        % expDay7 = 7th calendar day in the file, NOT "column 7" or "day 1" on the x-axis
        fLabs{(ifn-1)*T+it} = sprintf('%s expDay%g', ftag{ifn}, dList(it));
    end
end
end

function X3 = pr_step10_mouse_level_from_enriched(E, mks)
% One row per mouse: [max PR, total licks, mean within Hz, mean between Hz, mean 1/IEI] over all E rows.
n = numel(mks);
X3 = nan(n, 5);
for i = 1:n
    r = strcmp(E.mouse_key, mks{i});
    X3(i, 1) = max(E.Trial(r), [], 'omitnan');
    if ismember('lick_count', E.Properties.VariableNames)
        lc = double(E.lick_count(r));
        if isempty(lc) || all(isnan(lc)) || all(~isfinite(lc))
            X3(i, 2) = NaN;
        else
            X3(i, 2) = sum(lc, 'omitnan');
        end
    end
    X3(i, 3) = mean(E.WithinTrialLickRate(r), 'omitnan');
    if ismember('BetweenTrialLickRate', E.Properties.VariableNames)
        X3(i, 4) = mean(E.BetweenTrialLickRate(r), 'omitnan');
    end
    if ismember('iei_median', E.Properties.VariableNames)
        X3(i, 5) = pr_step10_mean_iti_lick_freq_hz(E.iei_median(r));
    end
end
end

function X3 = pr_step10_mousemean_three(S, mks)
% Optional: from MouseMean CSV (all days); prefer mouse_level_from_enriched after day filter.
n = numel(mks);
X3 = nan(n, 3);
for i = 1:n
    r = strcmp(S.mouse_key, mks{i});
    X3(i, 1) = max(double(S.Trial(r)), [], 'omitnan');
    X3(i, 2) = mean(S.MeanWithinTrialLickRate(r), 'omitnan');
    if ismember('MeanBetweenTrialLickRate', S.Properties.VariableNames)
        X3(i, 3) = mean(S.MeanBetweenTrialLickRate(r), 'omitnan');
    else
        X3(i, 3) = NaN;
    end
end
end

function Z = pr_step10_zrows(X)
% Column z-score; drop near-constant columns only (do not drop rows / mice)
X = double(X);
v = var(X, 0, 1, 'omitnan');
X = X(:, v >= 1e-12);
if isempty(X)
    error('All features are (near) constant; cannot z-score for clustering.');
end
Z = (X - mean(X, 1, 'omitnan')) ./ (std(X, 0, 1, 'omitnan') + eps);
Z(~isfinite(Z)) = 0;
end

function pr_step10_feature_correlation_heatmap(dsum, rDrop, outDir, dayMin)
if nargin < 4, dayMin = 4; end
fn = {'pr_breakpoint', 'total_lick_count', 'mean_within_hz', 'mean_between_hz', 'mean_iti_lick_freq_hz'};
L  = {'PR max trial', 'Total lick count', 'Within Hz', 'Between Hz', 'Mean 1/IEI (Hz)'};
M = [];
ok = false(1, numel(fn));
for k = 1:numel(fn)
    if ismember(fn{k}, dsum.Properties.VariableNames)
        M = [M, double(dsum.(fn{k}))]; %#ok<AGROW>
        ok(k) = true;
    end
end
L = L(ok);
if isempty(M)
    return
end
R = corr(M, 'rows', 'pairwise');
writetable(array2table(R, 'VariableNames', L, 'RowNames', L), fullfile(outDir, 'pr_feature_corr_mouse_day.csv'), 'WriteRowNames', true);
nF = size(R, 1);
f = figure('Color', 'w', 'Position', [100, 100, 90 + 70 * nF, 80 + 70 * nF]);
imagesc(R);
set(gca, 'XTick', 1:nF, 'XTickLabel', L, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nF, 'YTickLabel', L);
caxis([-1, 1]);
colormap(flipud(parula));
colorbar
title(sprintf('Feature correlation (mouse-day, day>=%d)', dayMin));
for i = 1:nF
    for j = 1:nF
        if i < j
            a = R(i, j);
            if isfinite(a) && abs(a) >= rDrop
                fprintf(2, 'High |r| (>=%.2f): %s vs %s: r=%.2f  (consider keeping one for clustering)\n', rDrop, L{i}, L{j}, a);
            end
        end
    end
end
try
    exportgraphics(f, fullfile(outDir, 'pr_feature_correlation_heatmap.png'), 'Resolution', 200);
    close(f);
catch
    try
        saveas(f, fullfile(outDir, 'pr_feature_correlation_heatmap.png'));
        close(f);
    end
end
end

function pr_step10_figure_radar_or_bar(X, fLabs, idx, mks, outDir, dayMin)
if nargin < 6, dayMin = 4; end
% Mean raw-feature profile per cluster. Column order matches pr_step10_build_traj_matrix_from_day:
% blocks: PR, total licks, within, between, ITIfr — each x calendar day (all days >= dayMin).
% Mouse-3 mode uses five pooled columns (*dge4) instead.
K = max(idx);
nc = min(30, size(X, 2));
fig = figure('Color', 'w', 'Position', [80, 60, 1000, 440]);
tiledlayout(1, K, 'TileSpacing', 'compact');
isPooled5 = (nc >= 5) && numel(fLabs) >= 5;
if isPooled5
    isPooled5 = all(cellfun(@(s) ~isempty(regexpi(char(s), 'dge4', 'once')), fLabs(1:5)));
end
for k = 1:K
    nexttile;
    m = mean(X(idx == k, 1:nc), 1, 'omitnan');
    se = std(X(idx == k, 1:nc), 0, 1, 'omitnan') / max(1, sqrt(sum(idx==k)));
    errorbar(1:nc, m, se, 'o-', 'LineWidth', 1.2, 'CapSize', 0);
    if isPooled5
        set(gca, 'XTick', 1:5, 'XTickLabel', cellstr(fLabs(1:5)), 'XTickLabelRotation', 15, 'FontSize', 8);
        xlabel(sprintf('5 pooled features (day>=%d); Passive lick-missing days are NaN / impute-omitted', dayMin), 'FontSize', 9);
    elseif numel(fLabs) >= nc
        set(gca, 'XTick', 1:nc, 'XTickLabel', fLabs(1:nc), 'XTickLabelRotation', 50, 'FontSize', 6, 'TickLabelInterpreter', 'none');
        xlabel({sprintf('Block order: PR, TotalLicks, Within, Between, ITIfr — each x calendar day'); ...
            sprintf('(day_index>=%d; Passive 6–10 lick features often NaN)', dayMin)}, 'FontSize', 7);
    else
        xlabel(sprintf('Column 1..%d (trajectory feature stack; days>=%d only in Step10)', nc, dayMin));
    end
    ylabel('Raw (pre-z-score) feature');
    title(sprintf('Cluster %d mean profile  n=%d  (day>=%d)', k, sum(idx==k), dayMin));
    grid on
end
try
    exportgraphics(fig, fullfile(outDir, 'pr_cluster_mean_trajectory_shape.png'), 'Resolution', 200);
    close(fig);
catch
    try
        saveas(fig, fullfile(outDir, 'pr_cluster_mean_trajectory_shape.png'));
        close(fig);
    end
end
end

function pr_step10_figure_group_cluster_heatmap(Out, outDir)
% Visual comparison: manifest Active/Passive vs k-means trajectory cluster (mouse counts).
if ~ismember('Group', Out.Properties.VariableNames)
    return
end
C = double(Out.trajectory_phenotype_cluster);
G = cellstr(strtrim(string(Out.Group(:))));
uC = unique(C, 'sorted');
Kc = numel(uC);
M = zeros(2, Kc);
for j = 1:Kc
    M(1, j) = sum(strcmpi(G, 'active') & C == uC(j));
    M(2, j) = sum(strcmpi(G, 'passive') & C == uC(j));
end
fig = figure('Color', 'w', 'Position', [100, 100, 520, 380]);
imagesc(M);
axis tight
colormap(flipud(gray(256)))
cb = colorbar;
cb.Label.String = 'Mouse count';
set(gca, 'YTick', 1:2, 'YTickLabel', {'Active (manifest)', 'Passive (manifest)'});
set(gca, 'XTick', 1:Kc, 'XTickLabel', arrayfun(@(c) sprintf('Trajectory C%d', c), uC, 'uni', 0));
xlabel('k-means cluster (same z-scored PR+within+between stack as Step10)');
mxv = max(M(:));
if ~isfinite(mxv) || mxv < 1
    mxv = 1;
end
for ii = 1:2
    for jj = 1:Kc
        v = M(ii, jj);
        tcol = [0, 0, 0];
        if v / mxv > 0.45
            tcol = [1, 1, 1];
        end
        text(jj, ii, sprintf('%d', v), 'HorizontalAlignment', 'center', 'Color', tcol, 'FontWeight', 'bold', 'FontSize', 12);
    end
end
title({'Manifest group vs trajectory-derived cluster', '(counts; compare to pr_cluster_x_group_matrix.csv)'}, 'FontSize', 10);
try
    exportgraphics(fig, fullfile(outDir, 'pr_group_vs_trajectory_cluster_heatmap.png'), 'Resolution', 200);
    close(fig);
catch
    try
        saveas(fig, fullfile(outDir, 'pr_group_vs_trajectory_cluster_heatmap.png'));
        close(fig);
    end
end
end

function pr_step10_figure_cluster_means_by_day(dsum, Out, outDir, dayMin)
if nargin < 4, dayMin = 4; end
% Same three features that enter k-means: mean mouse-day value per cluster per calendar day (mean ± SEM).
% This shows how the *combined* phenotypes separate over time, not only PR.
mksO = cellstr(Out.mouse_key);
cluId = double(Out.trajectory_phenotype_cluster);
mk2k = containers.Map('KeyType', 'char', 'ValueType', 'double');
for j = 1:numel(mksO)
    mk2k(mksO{j}) = cluId(j);
end
clu = nan(height(dsum), 1);
for r = 1:height(dsum)
    mk = char(string(dsum.mouse_key(r)));
    if isKey(mk2k, mk)
        clu(r) = mk2k(mk);
    end
end
ok = isfinite(clu);
dsum2 = dsum(ok, :);
clu2 = clu(ok);
uDays = sort(unique(double(dsum2.day_index)));
Kc = max(clu2);
cmap = lines(Kc);
feats = {'pr_breakpoint', 'total_lick_count', 'mean_within_hz', 'mean_between_hz', 'mean_iti_lick_freq_hz'};
ftit = {'PR breakpoint (max trial)', 'Total lick count (sum/trial)', 'Within-trial lick rate (Hz)', 'Between-trial lick rate (Hz)', 'Mean 1/IEI (Hz, within-trial)'};
nFt = numel(feats);
fig = figure('Color', 'w', 'Position', [60, 20, 900, 200 + 130 * nFt]);
tl = tiledlayout(nFt, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for fi = 1:nFt
    nexttile(tl);
    hold on
    fn = feats{fi};
    if ~ismember(fn, dsum2.Properties.VariableNames)
        title([ftit{fi} ' (missing in dsum)'], 'Interpreter', 'none');
        continue
    end
    for k = 1:Kc
        mu = nan(size(uDays));
        se = nan(size(uDays));
        for t = 1:numel(uDays)
            d = uDays(t);
            row = (clu2 == k) & (abs(double(dsum2.day_index) - d) < 1e-9);
            v = dsum2.(fn);
            v = v(row);
            v = v(isfinite(v));
            if isempty(v)
                continue
            end
            mu(t) = mean(v, 'omitnan');
            se(t) = std(v, 0, 1, 'omitnan') / max(1, sqrt(numel(v)));
        end
        errorbar(uDays, mu, se, 'o-', 'Color', cmap(k, :), 'LineWidth', 1.4, 'CapSize', 0, 'DisplayName', sprintf('Cluster %d', k));
    end
    ylabel('Mean ± SEM (mouse-day level)', 'Interpreter', 'none');
    title(ftit{fi}, 'Interpreter', 'none');
    grid on
    if fi == 1
        legend('Location', 'northeast', 'Box', 'off');
    end
    if fi == nFt
        xlabel(sprintf('Calendar day index (all >= %d in Step10)', dayMin), 'Interpreter', 'none');
    end
end
sgtitle(fig, '5 features (same as k-means): cluster mean ± SEM vs day (missing lick days stay NaN→0 in z)', 'FontSize', 10);
try
    exportgraphics(fig, fullfile(outDir, 'pr_multifeature_cluster_means_by_day.png'), 'Resolution', 200);
    close(fig);
catch
    try
        saveas(fig, fullfile(outDir, 'pr_multifeature_cluster_means_by_day.png'));
        close(fig);
    end
end
end

function pr_step10_spaghetti_breakpoint(Tlong, Out, outDir, dayMin)
if nargin < 4, dayMin = 4; end
% Per mouse per day: line color = k-means cluster; line style = Active (solid) vs Passive (dashed) from Step 8 Group.
Kc = max(double(Out.trajectory_phenotype_cluster));
cmap = lines(max(Kc, 1));
mksu = cellstr(Out.mouse_key);
clu = double(Out.trajectory_phenotype_cluster);
hasG = ismember('Group', Out.Properties.VariableNames);
fig = figure('Color', 'w', 'Position', [80, 100, 1040, 520]);
hold on
for i = 1:height(Out)
    r = strcmp(Tlong.mouse_key, mksu{i});
    if ~any(r)
        continue
    end
    t = Tlong(r, :);
    [d, o] = sort(t.day_index);
    y = t.pr_breakpoint; y = y(o);
    k = clu(i);
    if hasG
        gstr = strtrim(char(string(Out.Group{i})));
        if any(strcmpi(gstr, {'active', 'passive'}))
            if strcmpi(gstr, 'active')
                ls = '-';
            else
                ls = '--';
            end
        else
            ls = ':';
        end
    else
        ls = '-';
    end
    plot(d, y, 'LineStyle', ls, 'Color', cmap(k, :), 'LineWidth', 0.95, 'HandleVisibility', 'off');
end
xlabel(sprintf('day index (>= %d in Step10)', dayMin), 'Interpreter', 'none');
ylabel('PR breakpoint (max PR trial) per day');
if hasG
    title(sprintf('PR breakpoint vs day (color=trajectory cluster; solid=Active, dashed=Passive; day>=%d)', dayMin), 'Interpreter', 'none');
else
    title(sprintf('PR breakpoint vs day (cluster=color; day>=%d)', dayMin));
end
grid on
box off
ax = gca;
hC = gobjects(Kc, 1);
for k = 1:Kc
    hC(k) = plot(ax, nan, nan, '-', 'Color', cmap(k, :), 'LineWidth', 2);
end
legend(hC, arrayfun(@(k) sprintf('Cluster %d (trajectory)', k), 1:Kc, 'uni', 0), 'Location', 'northeast', 'Box', 'off');
if hasG
    annotation('textbox', [0.18, 0.01, 0.64, 0.06], 'String', 'Group line style: Active = solid, Passive = dashed, other = dotted', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
end
try
    exportgraphics(fig, fullfile(outDir, 'pr_breakpoint_spaghetti_by_cluster.png'), 'Resolution', 200);
    close(fig);
catch
    try
        saveas(fig, fullfile(outDir, 'pr_breakpoint_spaghetti_by_cluster.png'));
        close(fig);
    end
end
% Split: Active | Passive, color = same cluster
if hasG
    fig2 = figure('Color', 'w', 'Position', [90, 80, 1100, 480]);
    tl = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:2
        if s == 1
            want = 'active';
            subt = 'Active (assigned in manifest)';
        else
            want = 'passive';
            subt = 'Passive (assigned in manifest)';
        end
        nexttile(tl, s);
        hold on
        for i = 1:height(Out)
            gstr = strtrim(char(string(Out.Group{i})));
            if ~strcmpi(gstr, want)
                continue
            end
            r2 = strcmp(Tlong.mouse_key, mksu{i});
            if ~any(r2)
                continue
            end
            t2 = Tlong(r2, :);
            [d2, o2] = sort(t2.day_index);
            y2 = t2.pr_breakpoint; y2 = y2(o2);
            k2 = clu(i);
            plot(d2, y2, '-', 'Color', cmap(k2, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        end
        xlabel(sprintf('day index (>= %d)', dayMin), 'Interpreter', 'none');
        ylabel('PR breakpoint');
        title(subt, 'Interpreter', 'none');
        grid on
    end
    try
        exportgraphics(fig2, fullfile(outDir, 'pr_breakpoint_spaghetti_by_cluster_split_Active_Passive.png'), 'Resolution', 200);
        close(fig2);
    catch
        try
            saveas(fig2, fullfile(outDir, 'pr_breakpoint_spaghetti_by_cluster_split_Active_Passive.png'));
            close(fig2);
        end
    end
end
end
