function make_longitudinal_plotsall_statistic_NEWCOHORT_2026()
% Longitudinal plots + stats, updated for NEW COHORT + NEW TIMELINE (days 1–18)
%
% Cohort rules (user-provided):
% - Folder = cage (4-digit), subfolder = mouse color (mouse ID)
% - Sex: f/m
% - Group: a=active, p=passive  (stored as isPassive: 0=active, 1=passive)
% - Pairing (pair_id):
%   1) 6100: black(active) vs orange,red(passive)
%   2) 0911: red(active) vs orange(passive)
%   3) 0911: white(active) vs black(passive)
%   4) 0910: black(active) vs orange,red(passive)
%   5) 6099: orange(active) vs red(passive)
%   6) 6099: black(active) vs white(passive)  [6099 white died day1–13 only]
%
% Timeline (analysis phases; day1-2 habituation ignored):
%   PRE        : day 3–5  (water PR)
%   DURING     : day 6–10 (passive training morphine; active has PR, passive typically not)
%   POST       : day 11–13 (all morphine PR)
%   WITHDRAWAL : day 14–16 (water PR)
%   REEXPOSURE : day 17–18 (morphine PR)
%
% “Less reliable” transition days: configurable.
% Default requested set: [4 6 11 14 17] (edit below if you want).
%
% This script runs TWO variants automatically:
%   (1) INCLUDE transition days
%   (2) EXCLUDE transition days
% and saves outputs into run_*/figs/phases_variant_* and run_*/figs/stats_variant_*

%% ===================== CONFIG =====================
cfg = struct();

% Where the longitudinal_outputs live
cfg.tryPath = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Habituation days to ignore entirely (drop from analysis)
cfg.ignoreDays = [1 2];

% Less-reliable days (you can edit)
cfg.transitionDays = [4 6 11 14 17];

% If you still want the old pupil-only exclusion, keep it here (else set empty)
cfg.pupilOnlyExclude = ["7597_black"];  % set [] to disable

% Phase definition
cfg.phaseNames = ["PRE","DURING","POST","WITHDRAWAL","REEXPOSURE"];

% Save per-variant day-table CSV
cfg.saveDayTableCSV = true;

%% ===================== COHORT METADATA =====================
% Canonical mouse_key format in THIS script: "cage_color" (lowercase), e.g., "6100_red"
C = cohortTable_NEWCOHORT_2026();  % returns table with mouse_key_norm, sex, isPassive, pair_id

%% ===================== locate latest run and read CSV =====================
tryPath = cfg.tryPath;
if ~exist(tryPath,'dir')
    here = pwd; cand = here;
    for up=1:6
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), tryPath = p; break; end
        cand = fileparts(cand);
    end
end
if ~exist(tryPath,'dir'), error('Cannot find longitudinal_outputs folder.'); end

d = dir(fullfile(tryPath,'run_*')); assert(~isempty(d),'No run_* folders found.');
[~,idx] = max([d.datenum]); runDir = fullfile(d(idx).folder, d(idx).name);

csvPath = fullfile(runDir, 'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0, 'Missing %s', csvPath);

fprintf('Reading: %s\n', csvPath);
T0 = readtable(csvPath, 'VariableNamingRule','preserve');

% Ensure types
T0 = ensureString(T0,'day_name');
T0 = ensureString(T0,'mouse_key');
if hasVar(T0,'Session_Paradigm'), T0 = ensureString(T0,'Session_Paradigm'); end
if hasVar(T0,'sex'), T0 = ensureString(T0,'sex'); end
if hasVar(T0,'pair_id'), T0 = ensureString(T0,'pair_id'); end

if hasVar(T0,'isPassive') && ~isnumeric(T0.isPassive), T0.isPassive = double(T0.isPassive); end
if ~isnumeric(T0.day_index),   T0.day_index   = double(T0.day_index);   end
if hasVar(T0,'session_idx') && ~isnumeric(T0.session_idx)
    T0.session_idx = double(T0.session_idx);
end

% Normalize mouse_key and merge cohort metadata
T0.mouse_key_norm = normalizeMouseKey(T0.mouse_key);

% Join cohort: overwrite isPassive if cohort provides it (and warn on mismatch)
T0 = joinCohortIntoFrameTable(T0, C);

% Pupil-only exclusion list
if hasVar(T0,'Diameter_px') && ~isempty(cfg.pupilOnlyExclude)
    mk = string(T0.mouse_key_norm);
    for k=1:numel(cfg.pupilOnlyExclude)
        mask = (mk == lower(cfg.pupilOnlyExclude(k)));
        if any(mask)
            fprintf('Pupil-only exclude: setting Diameter_px=NaN for %d rows (%s)\n', nnz(mask), cfg.pupilOnlyExclude(k));
            T0.Diameter_px(mask) = NaN;
        end
    end
end

% Drop frames from TTL-only "extra" trials
if hasVar(T0,'TrialRequirement')
    trr = T0.TrialRequirement; if ~isstring(trr), trr = string(trr); end
    extraMask = false(height(T0),1);
    if hasVar(T0,'Trial'), extraMask = (T0.Trial>0) & strcmpi(strtrim(trr),"extra");
    else, extraMask = strcmpi(strtrim(trr),"extra");
    end
    if any(extraMask)
        fprintf('Excluding %d frames with TrialRequirement="extra" (%.2f%%).\n', ...
            nnz(extraMask), 100*nnz(extraMask)/height(T0));
        T0(extraMask,:) = [];
    end
end

% Create base fig folders
figDir = fullfile(runDir, 'figs'); if ~exist(figDir,'dir'), mkdir(figDir); end
fprintf('Base figs folder: %s\n', figDir);

%% ===================== RUN TWO VARIANTS =====================
variants = { ...
    struct('name','INCLUDE_TRANSITION_DAYS', 'excludeTransitions', false), ...
    struct('name','EXCLUDE_TRANSITION_DAYS', 'excludeTransitions', true)};

for v=1:numel(variants)
    V = variants{v};
    fprintf('\n================ VARIANT: %s ================\n', V.name);

    % Filter by ignoreDays and (optionally) transitionDays
    T = T0;
    T = T(~ismember(T.day_index, cfg.ignoreDays), :);
    if V.excludeTransitions
        T = T(~ismember(T.day_index, cfg.transitionDays), :);
    end

    % Add phase label at the frame level (useful for per-trial/peri-event splits if needed)
    T.phase = phaseFromDay(T.day_index);

    % Keep only rows that fall into defined phases (3–18)
    T = T(~isundefined(T.phase), :);

    % Build output folders
    figDirVar   = fullfile(figDir, ['phases_variant_' V.name]);
    statsOutDir = fullfile(figDir, ['stats_variant_'  V.name]);
    if ~exist(figDirVar,'dir'), mkdir(figDirVar); end
    if ~exist(statsOutDir,'dir'), mkdir(statsOutDir); end

    % ---- Compute timebase and normalize TTLs
    tb = pickTimebase(T);
    hasLick   = hasVar(T,'Lick_TTL');
    hasReward = hasVar(T,'Injector_TTL');
    hasPupil  = hasVar(T,'Diameter_px');

    if hasLick
        T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
        T.Lick_TTL = T.Lick_TTL > 0.5;
    end
    if hasReward
        T.Injector_TTL(isnan(T.Injector_TTL)) = 0;
        T.Injector_TTL = T.Injector_TTL > 0.5;
    end

    % ---- Per-day table D (computed from this variant’s T)
    [D, TT_all] = computeDayAndTrialTables(T, tb, hasLick, hasReward, hasPupil);

    % Save D if requested
    if cfg.saveDayTableCSV
        outCSV = fullfile(statsOutDir, 'D_day_table.csv');
        try
            writetable(D, outCSV);
        catch
            warning('Could not write %s', outCSV);
        end
    end

    % ---- PLOTS (per-day spaghetti + group means) ----
    colsAll = lines(numel(unique(D.mouse_key_norm,'stable')));

    % 0) Requirement
    if hasVar(D,'RequirementLast')
        f = figure('Color','w'); hold on
        plotAllMiceLines(D, 'RequirementLast', 'Requirement (last in JSONL)', colsAll, '-');
        title(['Requirement by day — ' V.name]);
        savepng(f, fullfile(figDirVar,'req_by_day_all_mice.png')); close(f);
    end

    % 1) TTL metrics overview
    if hasLick || hasReward
        f = figure('Color','w'); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

        ax=nexttile; hold(ax,'on');
        if hasVar(D,'lick_freq_per_min'), plotAllMiceLines(D,'lick_freq_per_min','Lick frequency (/min)',colsAll,'-'); end
        if hasVar(D,'rew_freq_per_min'),  plotAllMiceLines(D,'rew_freq_per_min','Reward frequency (/min)',colsAll,'--'); end
        title(ax,'Frequencies');

        ax=nexttile; hold(ax,'on');
        if hasVar(D,'lick_meanDur_s'), plotAllMiceLines(D,'lick_meanDur_s','Mean duration (s)',colsAll,'-'); end
        if hasVar(D,'rew_meanDur_s'),  plotAllMiceLines(D,'rew_meanDur_s','Mean duration (s)',colsAll,'--'); end
        title(ax,'Mean durations');

        ax=nexttile; hold(ax,'on');
        if hasVar(D,'lick_totalDur_s'), plotAllMiceLines(D,'lick_totalDur_s','Total duration (s/session)',colsAll,'-'); end
        if hasVar(D,'rew_totalDur_s'),  plotAllMiceLines(D,'rew_totalDur_s','Total duration (s/session)',colsAll,'--'); end
        title(ax,'Totals');

        ax=nexttile; hold(ax,'on');
        if hasVar(D,'lick_medianIEI_s'), plotAllMiceLines(D,'lick_medianIEI_s','Median IEI / IRI (s)',colsAll,'-'); end
        if hasVar(D,'rew_medianIRI_s'),  plotAllMiceLines(D,'rew_medianIRI_s','Median IEI / IRI (s)',colsAll,'--'); end
        title(ax,'Intervals');

        title(tl, ['TTL-derived metrics across days — ' V.name]);
        savepng(f, fullfile(figDirVar,'lick_reward_metrics_by_day.png')); close(f);
    end

    % 2) Group means ± SEM vs day (Active vs Passive)
    if hasVar(D,'isPassive') && any(ismember(D.isPassive,[0 1]))
        plotGroupVsDay_split(D, figDirVar, ['_' V.name]);
    end

    % 3) Phase-wise summary (boxplots per phase; Active vs Passive)
    plotPhaseSummaries(D, figDirVar, V.name);

    % 4) Pupil correlations (if present)
    if hasPupil && (hasVar(D,'lick_freq_per_min') || hasVar(D,'rew_freq_per_min'))
        plotPupilCorrelations(D, fullfile(figDirVar,'pupil_correlations_overview.png'), V.name);
    end

    % 5) Per-trial trajectories
    if ~isempty(TT_all) && hasVar(TT_all,'Trial')
        plotPerTrialAll(TT_all, fullfile(figDirVar,'per_trial_all_mice.png'));
        if hasVar(TT_all,'isPassive')
            plotPerTrialByGroup(TT_all, figDirVar);
        end
    end

    % 6) Peri-event pupil Δ (if present)
    if hasPupil && (hasLick || hasReward)
        win   = [-2 6]; baseW = [-2 0]; dt = medianDT(tb);
        if hasLick,   periEventPlot(T, tb, 'Lick_TTL',     win, baseW, dt, figDirVar, 'pupil_peri_lick_active_vs_passive.png'); end
        if hasReward, periEventPlot(T, tb, 'Injector_TTL', win, baseW, dt, figDirVar, 'pupil_peri_reward_active_vs_passive.png'); end
    end

    % ---- STATS SUITE updated to PHASE bins ----
    metrics_for_stats = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
                         'lick_meanDur_s','rew_meanDur_s','lick_totalDur_s','rew_totalDur_s', ...
                         'lick_medianIEI_s','rew_medianIRI_s','pupil_mean'};
    labels_for_stats  = {'Requirement','Lick freq (/min)','Reward freq (/min)', ...
                         'Lick mean dur (s)','Reward mean dur (s)','Lick total dur (s)','Reward total dur (s)', ...
                         'Lick median IEI (s)','Reward median IRI (s)','Pupil mean (px)'};

    runStatsComparisons_and_Correlations_PHASES(D, metrics_for_stats, labels_for_stats, statsOutDir);
    withinGroup_PhaseANOVA_and_Pairs(D, 0, metrics_for_stats, labels_for_stats, statsOutDir); % Active-only mice
    withinGroup_PhaseANOVA_and_Pairs(D, 1, metrics_for_stats, labels_for_stats, statsOutDir); % Had-passive mice
    computeAndPlot_AP_Index_PHASES(D, metrics_for_stats, labels_for_stats, statsOutDir);
    activePassive_byPhase(D, metrics_for_stats, labels_for_stats, statsOutDir);

    fprintf('Variant done: %s\n', V.name);
end

fprintf('\nALL DONE. Outputs are under:\n  %s\n', figDir);
end

%% ===================== COHORT TABLE =====================
function C = cohortTable_NEWCOHORT_2026()
% mouse_key_norm must be lowercase "cage_color"
mouse_key_norm = [ ...
    "6100_red"
    "6100_orange"
    "6100_black"
    "0911_red"
    "0911_orange"
    "0911_black"
    "0911_white"
    "0910_red"
    "0910_orange"
    "0910_black"
    "6099_red"
    "6099_orange"
    "6099_black"
    "6099_white"   % included for completeness (died day13; data ends early)
];

sex = [ ...
    "f"
    "f"
    "f"
    "f"
    "f"
    "f"
    "f"
    "m"
    "m"
    "m"
    "m"
    "m"
    "m"
    "m"
];

% isPassive: 0=active, 1=passive
isPassive = [ ...
    1  % 6100_red passive
    1  % 6100_orange passive
    0  % 6100_black active
    0  % 0911_red active
    1  % 0911_orange passive
    1  % 0911_black passive
    0  % 0911_white active
    1  % 0910_red passive
    1  % 0910_orange passive
    0  % 0910_black active
    1  % 6099_red passive
    0  % 6099_orange active
    0  % 6099_black active
    1  % 6099_white passive (died day13)
];

% Pair IDs (same for one active vs one/two passives)
pair_id = [ ...
    "PAIR_6100"   % red
    "PAIR_6100"   % orange
    "PAIR_6100"   % black
    "PAIR_0911_A" % 0911 red–orange
    "PAIR_0911_A"
    "PAIR_0911_B" % 0911 black–white
    "PAIR_0911_B"
    "PAIR_0910"   % 0910 black vs orange/red
    "PAIR_0910"
    "PAIR_0910"
    "PAIR_6099_A" % 6099 orange–red
    "PAIR_6099_A"
    "PAIR_6099_B" % 6099 black–white
    "PAIR_6099_B"
];

C = table(mouse_key_norm, sex, isPassive, pair_id);
end

%% ===================== PHASE LABELS =====================
function ph = phaseFromDay(day_index)
di = double(day_index(:));
lab = strings(size(di));

lab(di>=3  & di<=5 )  = "PRE";
lab(di>=6  & di<=10)  = "DURING";
lab(di>=11 & di<=13)  = "POST";
lab(di>=14 & di<=16)  = "WITHDRAWAL";
lab(di>=17 & di<=18)  = "REEXPOSURE";

ph = categorical(lab, ["PRE","DURING","POST","WITHDRAWAL","REEXPOSURE"]);
end

%% ===================== CORE COMPUTATION =====================
function [D, TT_all] = computeDayAndTrialTables(T, tb, hasLick, hasReward, hasPupil)
% Per-session table S -> per-day table D (median across sessions), plus per-trial table.

% ----- per-session keys -----
needCols = {'mouse_key_norm','day_index','day_name'};
if hasVar(T,'session_idx'), needCols{end+1}='session_idx'; else, T.session_idx = zeros(height(T),1); needCols{end+1}='session_idx'; end
sessKeys = unique(T(:,needCols),'rows','stable');

S = table();
S.mouse_key_norm = sessKeys.mouse_key_norm;
S.day_index      = sessKeys.day_index;
S.day_name       = sessKeys.day_name;
S.session_idx    = sessKeys.session_idx;

S.phase          = phaseFromDay(S.day_index);

S.RequirementLast = nan(height(S),1);
S.isPassive       = nan(height(S),1);
S.SessionMinutes  = nan(height(S),1);
[S.lick_n,S.lick_freq_per_min,S.lick_meanDur_s,S.lick_totalDur_s,S.lick_medianIEI_s] = deal(nan(height(S),1));
[S.rew_n,S.rew_freq_per_min,S.rew_meanDur_s,S.rew_totalDur_s,S.rew_medianIRI_s]       = deal(nan(height(S),1));
S.pupil_mean = nan(height(S),1);

% Keep cohort fields if present
S.sex    = strings(height(S),1);
S.pair_id= strings(height(S),1);

for i=1:height(sessKeys)
    mk  = sessKeys.mouse_key_norm(i);
    di  = sessKeys.day_index(i);
    ss  = sessKeys.session_idx(i);

    rows = (T.mouse_key_norm==mk) & (T.day_index==di) & (T.session_idx==ss);
    if ~any(rows), continue; end

    % cohort fields (string-safe mode)
    if hasVar(T,'sex'),     S.sex(i)     = modeString(T.sex(rows)); end
    if hasVar(T,'pair_id'), S.pair_id(i) = modeString(T.pair_id(rows)); end

    tb_s = tb(rows);
    sessionDur_s = finiteRange(tb_s);
    if ~(isfinite(sessionDur_s) && sessionDur_s>0) && hasVar(T,'Frame')
        fr = double(T.Frame(rows)); if any(isfinite(fr)), sessionDur_s = (max(fr)-min(fr))/30; end
    end
    S.SessionMinutes(i) = sessionDur_s/60;

    if hasVar(T,'RequirementLast')
        S.RequirementLast(i) = mean(double(T.RequirementLast(rows)),'omitnan');
    end
    if hasVar(T,'isPassive')
        ip = double(T.isPassive(rows)); ip = ip(isfinite(ip));
        if ~isempty(ip), S.isPassive(i) = mode(round(ip)); end
    end
    if hasPupil
        S.pupil_mean(i) = mean(double(T.Diameter_px(rows)), 'omitnan');
    end

    if hasLick
        [n,meanDur,totalDur,medianIEI] = eventMetrics(tb_s, T.Lick_TTL(rows));
        S.lick_n(i)=n; S.lick_meanDur_s(i)=meanDur; S.lick_totalDur_s(i)=totalDur; S.lick_medianIEI_s(i)=medianIEI;
        if S.SessionMinutes(i)>0, S.lick_freq_per_min(i) = n/S.SessionMinutes(i); end
    end
    if hasReward
        [n,meanDur,totalDur,medianIRI] = eventMetrics(tb_s, T.Injector_TTL(rows));
        S.rew_n(i)=n; S.rew_meanDur_s(i)=meanDur; S.rew_totalDur_s(i)=totalDur; S.rew_medianIRI_s(i)=medianIRI;
        if S.SessionMinutes(i)>0, S.rew_freq_per_min(i) = n/S.SessionMinutes(i); end
    end
end

% ----- per-day collapse (median across sessions) -----
dayKeys = unique(S(:,{'mouse_key_norm','day_index','day_name'}),'rows','stable');
D = dayKeys;
D.phase = phaseFromDay(D.day_index);

agg = @(x) median(x,'omitnan');

coreCols = {'RequirementLast','isPassive','SessionMinutes', ...
           'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
           'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s','pupil_mean'};

for c = 1:numel(coreCols), D.(coreCols{c}) = nan(height(D),1); end

% propagate cohort labels per day (mode across sessions) using string-safe mode
D.sex     = strings(height(D),1);
D.pair_id = strings(height(D),1);

for i=1:height(D)
    mk = D.mouse_key_norm(i); di = D.day_index(i);
    r  = (S.mouse_key_norm==mk) & (S.day_index==di);
    if ~any(r), continue; end
    for c = 1:numel(coreCols)
        nm = coreCols{c};
        D.(nm)(i) = agg(S.(nm)(r));
    end
    if any(strlength(S.sex(r))>0),     D.sex(i)     = modeString(S.sex(r)); end
    if any(strlength(S.pair_id(r))>0), D.pair_id(i) = modeString(S.pair_id(r)); end
end

% ----- per-trial table -----
TT_all = table();
if hasVar(T,'Trial')
    TT_all = computePerTrial(T, tb, hasLick, hasReward, hasPupil);
    if ~isempty(TT_all)
        TT_all.phase = phaseFromDay(TT_all.day_index);
    end
end
end

%% ===================== PHASE SUMMARY PLOTS =====================
function plotPhaseSummaries(D, figDir, tag)
% For each metric, show Active vs Passive within each phase (per-mouse phase medians)
if ~hasVar(D,'phase') || ~hasVar(D,'isPassive'), return; end

metrics = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
           'lick_totalDur_s','rew_totalDur_s','pupil_mean'};
labels  = {'Requirement','Lick freq (/min)','Reward freq (/min)', ...
           'Lick total dur (s)','Reward total dur (s)','Pupil mean (px)'};

phases = categories(D.phase);

for i=1:numel(metrics)
    y = metrics{i};
    if ~hasVar(D,y), continue; end
    sub = D(isfinite(D.(y)) & ismember(D.isPassive,[0 1]) & ~isundefined(D.phase), ...
            {'mouse_key_norm','phase','isPassive',y});
    if isempty(sub), continue; end

    % per-mouse, per-phase medians
    M = groupsummary(sub, {'mouse_key_norm','phase','isPassive'}, 'median', y);
    M.Properties.VariableNames(end) = {'y_med'};

    fig = figure('Color','w','Position',[120 120 1160 420]);
    tl  = tiledlayout(1,numel(phases),'TileSpacing','compact','Padding','compact');

    csvRows = {};
    for p=1:numel(phases)
        nexttile; hold on
        Mp = M(string(M.phase)==string(phases{p}), :);
        if isempty(Mp), axis off; title(phases{p}); continue; end

        Xa = Mp.y_med(Mp.isPassive==0); Xa = Xa(isfinite(Xa));
        Xp = Mp.y_med(Mp.isPassive==1); Xp = Xp(isfinite(Xp));

        plotTwoGroupsDots(Xa, Xp, labels{i}, phases{p});
        [pval,~] = safeRanksum(Xa, Xp);
        yl = ylim; ybar = yl(2) - 0.05*range(yl);
        drawSigBar(gca, 1, 2, ybar, pval);

        csvRows(end+1,:) = {y, string(phases{p}), numel(Xa), numel(Xp), nanmedian(Xa), nanmedian(Xp), pval}; %#ok<AGROW>
    end
    title(tl, sprintf('%s — Active vs Passive by PHASE (per-mouse phase medians) — %s', labels{i}, tag));
    savepng(fig, fullfile(figDir, sprintf('phase_ACTIVEvsPASSIVE_%s.png', y)));
    close(fig);

    if ~isempty(csvRows)
        Tcsv = cell2table(csvRows, 'VariableNames', ...
            {'metric','phase','n_active','n_passive','median_active','median_passive','p_ranksum'});
        try, writetable(Tcsv, fullfile(figDir, sprintf('phase_ACTIVEvsPASSIVE_%s.csv', y))); catch, end
    end
end
end

%% ===================== UPDATED STATS (PHASES) =====================
function runStatsComparisons_and_Correlations_PHASES(D, metrics, labels, outDir)
% (1) Any-Passive (mouse-level) vs Active-only (mouse-level) using per-mouse medians over ALL included days
% (2) Same comparison within each PHASE using per-mouse phase medians
% (3) Pearson correlation heatmaps (overall + per phase)

if ~hasVar(D,'mouse_key_norm') || ~hasVar(D,'isPassive')
    warning('mouse_key_norm/isPassive missing; skipping stats.'); return;
end

% mouse-level group assignment: had any passive day?
Gmouse = groupsummary(D, 'mouse_key_norm', @(x) any(x==1), 'isPassive');
Gmouse.Properties.VariableNames(end) = {'hasPassive'}; % logical

keyD = buildKey(D, {'mouse_key_norm'});
keyG = buildKey(Gmouse, {'mouse_key_norm'});
[tf,loc] = ismember(keyD, keyG);
D.mouse_hasPassive = nan(height(D),1);
D.mouse_hasPassive(tf) = double(Gmouse.hasPassive(loc(tf)));

% keep only metrics that exist
keep = metrics(ismember(metrics, D.Properties.VariableNames));
labels = labels(ismember(metrics, keep));
metrics = keep;

csvRows = {};

% (1) ALL-days per-mouse medians
for i = 1:numel(metrics)
    y = metrics{i};

    M = groupsummary(D(:,{'mouse_key_norm','mouse_hasPassive',y}), ...
                     {'mouse_key_norm','mouse_hasPassive'}, 'median', y);
    if isempty(M), continue; end
    M.Properties.VariableNames(end) = {'y_med'};

    Xa = M.y_med(M.mouse_hasPassive==0);
    Xp = M.y_med(M.mouse_hasPassive==1);

    [p,statsOut] = safeRanksum(Xa, Xp);

    fig = figure('Color','w','Position',[100 100 560 420]); hold on
    plotTwoGroupsDots(Xa, Xp, labels{i}, 'ALL included days (per-mouse median)');
    ylims = ylim; ybar = ylims(2) - 0.05*range(ylims);
    drawSigBar(gca, 1, 2, ybar, p);
    title(sprintf('%s — Any Passive vs None (ALL included days)', labels{i}));
    savepng(fig, fullfile(outDir, sprintf('stat_ANYvsNONE_ALL_%s.png', y))); close(fig);

    csvRows(end+1,:) = {y,'ALL',numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),p,statsOut.delta}; %#ok<AGROW>
end

% (2) Per-phase per-mouse medians
if hasVar(D,'phase') && ~all(isundefined(D.phase))
    phases = categories(D.phase);
    for i = 1:numel(metrics)
        y = metrics{i};
        fig = figure('Color','w','Position',[80 80 1160 420]);
        tl  = tiledlayout(1,numel(phases),'TileSpacing','compact','Padding','compact');

        for pIdx = 1:numel(phases)
            ph = phases{pIdx};
            nexttile; hold on

            sub = D(string(D.phase)==string(ph), {'mouse_key_norm','mouse_hasPassive',y});
            if isempty(sub), axis off; title(ph); continue; end

            Mp = groupsummary(sub, {'mouse_key_norm','mouse_hasPassive'}, 'median', y);
            if isempty(Mp), axis off; title(ph); continue; end
            Mp.Properties.VariableNames(end) = {'y_med'};

            Xa = Mp.y_med(Mp.mouse_hasPassive==0);
            Xp = Mp.y_med(Mp.mouse_hasPassive==1);

            [pval,~] = safeRanksum(Xa, Xp);
            plotTwoGroupsDots(Xa, Xp, labels{i}, ph);
            ylims = ylim; ybar = ylims(2) - 0.05*range(ylims);
            drawSigBar(gca, 1, 2, ybar, pval);

            csvRows(end+1,:) = {y,ph,numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),pval,NaN}; %#ok<AGROW>
        end

        title(tl, sprintf('%s — Any Passive vs None (per PHASE; per-mouse medians)', labels{i}));
        savepng(fig, fullfile(outDir, sprintf('stat_PHASES_%s.png', y))); close(fig);
    end
end

% Write CSV
if ~isempty(csvRows)
    Tcsv = cell2table(csvRows, 'VariableNames', ...
        {'metric','level','n_activeOnly','n_hadPassive','median_activeOnly','median_hadPassive','p_ranksum','cliffs_delta'});
    outCSV = fullfile(outDir, 'stats_pvalues_ANYvsNONE_and_PHASES.csv');
    try, writetable(Tcsv, outCSV); catch, warning('Could not write %s', outCSV); end
end

% (3) Correlation matrices (overall + per phase)
corrVarsAll = {'RequirementLast','lick_freq_per_min','rew_freq_per_min','lick_totalDur_s','rew_totalDur_s','pupil_mean'};
corrVarsAll = corrVarsAll(ismember(corrVarsAll, D.Properties.VariableNames));
if numel(corrVarsAll) >= 2
    plotCorrMatrix(D, corrVarsAll, fullfile(outDir,'corr_matrix_ALL.png'), 'ALL rows');
    plotCorrMatrix(D(D.mouse_hasPassive==0,:), corrVarsAll, fullfile(outDir,'corr_matrix_ACTIVEONLY.png'), 'Active-only');
    plotCorrMatrix(D(D.mouse_hasPassive==1,:), corrVarsAll, fullfile(outDir,'corr_matrix_HADPASSIVE.png'), 'Had-passive');

    if hasVar(D,'phase') && ~all(isundefined(D.phase))
        phases = categories(D.phase);
        for pIdx=1:numel(phases)
            ph = phases{pIdx};
            sub = D(string(D.phase)==string(ph), :);
            if ~isempty(sub)
                plotCorrMatrix(sub, corrVarsAll, fullfile(outDir, sprintf('corr_matrix_%s.png', ph)), ['Phase ' ph]);
            end
        end
    end
end
end

function withinGroup_PhaseANOVA_and_Pairs(D, groupFlag, metrics, labels, outDir)
% Within a group of mice (Active-only=0, Had-passive=1), test across PHASES.
% Uses per-mouse medians within each phase. KW omnibus + FDR pairwise heatmap.

if ~hasVar(D,'isPassive') || ~hasVar(D,'phase'), return; end

gname = ternary(groupFlag==0, 'ACTIVEONLY', 'HADPASSIVE');
phases = categories(D.phase);

% mouse-level group membership
Gmouse = groupsummary(D, 'mouse_key_norm', @(x) any(x==1), 'isPassive');
Gmouse.Properties.VariableNames(end) = {'hasPassive'};

if groupFlag==0
    miceKeep = string(Gmouse.mouse_key_norm(Gmouse.hasPassive==0));
else
    miceKeep = string(Gmouse.mouse_key_norm(Gmouse.hasPassive==1));
end

for i=1:numel(metrics)
    y = metrics{i};
    if ~hasVar(D,y), continue; end
    lbl = labels{i};

    sub = D(ismember(D.mouse_key_norm, miceKeep) & ~isundefined(D.phase), {'mouse_key_norm','phase',y});
    M = groupsummary(sub, {'mouse_key_norm','phase'}, 'median', y);
    if isempty(M), continue; end
    M.Properties.VariableNames(end) = {'y_med'};

    % Assemble per-phase vectors
    groups = cell(1,numel(phases));
    for p=1:numel(phases)
        xp = M.y_med(string(M.phase)==string(phases{p}));
        groups{p} = xp(isfinite(xp));
    end
    p_kw = kruskal_p_fallback(groups);

    % Pairwise tests (paired on mouse when possible)
    pairs   = nchoosek(1:numel(phases),2);
    p_pair  = nan(size(pairs,1),1);
    n_used  = zeros(size(pairs,1),1);

    for k=1:size(pairs,1)
        p1 = pairs(k,1); p2 = pairs(k,2);

        M1 = M(string(M.phase)==string(phases{p1}), {'mouse_key_norm','y_med'});
        M2 = M(string(M.phase)==string(phases{p2}), {'mouse_key_norm','y_med'});
        J  = innerjoin(M1,M2,'Keys','mouse_key_norm', ...
                       'LeftVariables',{'mouse_key_norm','y_med'}, ...
                       'RightVariables','y_med');

        if ~isempty(J)
            x1 = J.y_med_M1; x2 = J.y_med_M2;
            good = isfinite(x1) & isfinite(x2);
            x1 = x1(good); x2 = x2(good);
            n_used(k) = numel(x1);

            if n_used(k) >= 3
                if all(x1 == x2), p_pair(k)=1;
                else
                    try, p_pair(k) = signrank(x1,x2);
                    catch, p_pair(k) = safeRanksum(x1,x2);
                    end
                end
                continue;
            end
        end

        % fallback unpaired
        p_pair(k) = safeRanksum(groups{p1}, groups{p2});
        n_used(k) = min(numel(groups{p1}), numel(groups{p2}));
    end

    [p_adj, sig_mask] = fdr_bh_local(p_pair, 0.05);

    % Heatmap of -log10(p_fdr)
    H = nan(numel(phases));
    for k=1:size(pairs,1)
        i1 = pairs(k,1); i2 = pairs(k,2);
        val = -log10(max(p_adj(k), eps));
        H(i1,i2)=val; H(i2,i1)=val;
    end

    fig = figure('Color','w','Position',[120 120 520 460]);
    imagesc(H); axis square; colorbar; colormap(parula);
    title(sprintf('%s — %s (KW p=%.3g)', gname, lbl, p_kw));
    set(gca,'XTick',1:numel(phases),'XTickLabel',phases,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(phases),'YTickLabel',phases);

    for r=1:numel(phases)
        for c=1:numel(phases)
            if r==c, text(c,r,'—','HorizontalAlignment','center'); continue; end
            kk = find((pairs(:,1)==min(r,c)) & (pairs(:,2)==max(r,c)), 1);
            if isempty(kk), continue; end
            star = sigStarFromP(p_adj(kk));
            if ~isempty(star)
                text(c,r,star,'HorizontalAlignment','center','FontWeight','bold','Color','k');
            end
        end
    end
    savepng(fig, fullfile(outDir, sprintf('withinPHASE_%s_%s.png', gname, y))); close(fig);

    % CSV dump
    pair_lbl = strings(size(pairs,1),1);
    for k=1:size(pairs,1)
        pair_lbl(k) = sprintf('%s vs %s', phases{pairs(k,1)}, phases{pairs(k,2)});
    end
    Ctab = table(pair_lbl, p_pair, p_adj, sig_mask, n_used, ...
        'VariableNames', {'pair','p_raw','p_fdr','is_sig','n'});
    try, writetable(Ctab, fullfile(outDir, sprintf('withinPHASE_%s_%s_pairs.csv', gname, y))); catch, end
end
end

function computeAndPlot_AP_Index_PHASES(D, metrics, labels, outDir)
% Phase-level AP index:
%   For each day d that has both Active and Passive, compute I(d)=(A-P)/(A+P)
%   Then compare I(d) distributions across phases.

if ~hasVar(D,'isPassive') || ~hasVar(D,'phase'), return; end
phases = categories(D.phase);

for i=1:numel(metrics)
    y = metrics{i};
    if ~hasVar(D,y), continue; end
    lbl = labels{i};

    days = sort(unique(D.day_index));
    I = nan(size(days));
    Ph = strings(size(days));

    for k=1:numel(days)
        di = days(k);
        A = D{D.day_index==di & D.isPassive==0, y}; A = A(isfinite(A));
        P = D{D.day_index==di & D.isPassive==1, y}; P = P(isfinite(P));
        if ~isempty(A) && ~isempty(P)
            a = median(A,'omitnan'); p = median(P,'omitnan');
            if isfinite(a) && isfinite(p) && (a+p)~=0
                I(k) = (a - p) / (a + p);
                Ph(k) = modeString(string(D.phase(D.day_index==di)));
            end
        end
    end

    good = isfinite(I) & strlength(Ph)>0;
    if nnz(good) < 3, continue; end

    % Plot I by day
    fig = figure('Color','w','Position',[120 120 700 420]); hold on
    plot(days(good), I(good), '-o','LineWidth',1.6,'MarkerSize',5,'Color',[0.2 0.2 0.2]);
    yline(0,'-','Color',[0.6 0.6 0.6]); grid on; box on;
    xlabel('Day'); ylabel(sprintf('I=(A-P)/(A+P) — %s', lbl));
    title(sprintf('Active–Passive Index by day — %s', lbl));
    savepng(fig, fullfile(outDir, sprintf('AP_index_by_day_%s.png', y))); close(fig);

    % Compare I across phases
    groups = cell(1,numel(phases));
    for pIdx=1:numel(phases)
        groups{pIdx} = I(good & Ph==string(phases{pIdx}));
    end
    p_kw = kruskal_p_fallback(groups);

    pairs = nchoosek(1:numel(phases),2);
    p_pair = nan(size(pairs,1),1);
    for k=1:size(pairs,1)
        p1=pairs(k,1); p2=pairs(k,2);
        p_pair(k) = safeRanksum(groups{p1}, groups{p2});
    end
    [p_adj, ~] = fdr_bh_local(p_pair, 0.05);

    H = nan(numel(phases));
    for k=1:size(pairs,1)
        i1=pairs(k,1); i2=pairs(k,2);
        val = -log10(max(p_adj(k), eps));
        H(i1,i2)=val; H(i2,i1)=val;
    end

    fig = figure('Color','w','Position',[120 120 520 460]);
    imagesc(H); axis square; colorbar; colormap(parula);
    title(sprintf('AP index — phase tests (KW p=%.3g): %s', p_kw, lbl));
    set(gca,'XTick',1:numel(phases),'XTickLabel',phases,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(phases),'YTickLabel',phases);

    for r=1:numel(phases)
        for c=1:numel(phases)
            if r==c, text(c,r,'—','HorizontalAlignment','center'); continue; end
            i1=min(r,c); i2=max(r,c);
            kk = find(pairs(:,1)==i1 & pairs(:,2)==i2, 1);
            if isempty(kk), continue; end
            st = sigStarFromP(p_adj(kk));
            if ~isempty(st), text(c,r,st,'HorizontalAlignment','center','FontWeight','bold','Color','k'); end
        end
    end
    savepng(fig, fullfile(outDir, sprintf('AP_index_phases_heatmap_%s.png', y))); close(fig);

    pair_lbl = strings(size(pairs,1),1);
    for k=1:size(pairs,1)
        pair_lbl(k) = sprintf('%s vs %s', phases{pairs(k,1)}, phases{pairs(k,2)});
    end
    Tcsv = table(repmat(string(lbl),numel(p_adj),1), pair_lbl, p_pair, p_adj, ...
        'VariableNames', {'metric','pair','p_raw','p_fdr'});
    try, writetable(Tcsv, fullfile(outDir, sprintf('AP_index_phase_pairs_%s.csv', y))); catch, end
end
end

function activePassive_byPhase(D, metrics, labels, outDir)
% Active vs Passive within each phase using per-mouse phase medians.
if ~hasVar(D,'isPassive') || ~hasVar(D,'phase'), return; end
phases = categories(D.phase);

for i=1:numel(metrics)
    y = metrics{i};
    if ~hasVar(D,y), continue; end
    lbl = labels{i};

    sub = D(~isundefined(D.phase) & isfinite(D.(y)) & ismember(D.isPassive,[0 1]), ...
            {'mouse_key_norm','phase','isPassive',y});
    if isempty(sub), continue; end

    M = groupsummary(sub, {'mouse_key_norm','phase','isPassive'}, 'median', y);
    M.Properties.VariableNames(end) = {'y_med'};

    fig = figure('Color','w','Position',[80 80 1160 420]);
    tl  = tiledlayout(1,numel(phases),'TileSpacing','compact','Padding','compact');

    csvRows = {};
    for pIdx=1:numel(phases)
        ph = phases{pIdx};
        nexttile; hold on
        Mp = M(string(M.phase)==string(ph), :);
        if isempty(Mp), axis off; title(ph); continue; end

        Xa = Mp.y_med(Mp.isPassive==0); Xa = Xa(isfinite(Xa));
        Xp = Mp.y_med(Mp.isPassive==1); Xp = Xp(isfinite(Xp));

        plotTwoGroupsDots(Xa, Xp, lbl, ph);
        [pval,~] = safeRanksum(Xa, Xp);
        yl = ylim; ybar = yl(2) - 0.05*range(yl);
        drawSigBar(gca, 1, 2, ybar, pval);

        csvRows(end+1,:) = {y,ph,numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),pval}; %#ok<AGROW>
    end

    title(tl, sprintf('%s — Active vs Passive (per PHASE; per-mouse medians)', lbl));
    savepng(fig, fullfile(outDir, sprintf('stat_ACTIVEvsPASSIVE_PHASES_%s.png', y)));
    close(fig);

    if ~isempty(csvRows)
        Tcsv = cell2table(csvRows, 'VariableNames', ...
            {'metric','phase','n_active','n_passive','median_active','median_passive','p_ranksum'});
        try, writetable(Tcsv, fullfile(outDir, sprintf('stat_ACTIVEvsPASSIVE_PHASES_%s.csv', y))); catch, end
    end
end
end

%% ===================== COHORT JOIN HELPERS =====================
function T = joinCohortIntoFrameTable(T, C)
% Adds: sex, pair_id, and a cohort-derived isPassive_cohort.
% If T already has isPassive, compare and overwrite with cohort where known.

% Map from mouse_key_norm -> cohort fields
keyT = string(T.mouse_key_norm);
keyC = string(C.mouse_key_norm);

[tf, loc] = ismember(keyT, keyC);

% Add fields
T.sex    = strings(height(T),1);
T.pair_id= strings(height(T),1);
T.isPassive_cohort = nan(height(T),1);

T.sex(tf)     = C.sex(loc(tf));
T.pair_id(tf) = C.pair_id(loc(tf));
T.isPassive_cohort(tf) = double(C.isPassive(loc(tf)));

% Overwrite isPassive if cohort has it
if hasVar(T,'isPassive')
    % warn if mismatch (only where both finite)
    both = tf & isfinite(T.isPassive) & isfinite(T.isPassive_cohort);
    if any(both)
        mismatch = both & (round(T.isPassive) ~= round(T.isPassive_cohort));
        if any(mismatch)
            fprintf('WARNING: isPassive mismatch for %d rows. Overwriting with cohort isPassive.\n', nnz(mismatch));
        end
    end
    % overwrite where cohort known
    T.isPassive(tf) = T.isPassive_cohort(tf);
else
    T.isPassive = T.isPassive_cohort;
end
end

function mk = normalizeMouseKey(mouse_key_col)
% Converts typical forms like "6100 red", "6100_red", "6100-red", "6100/red" into "6100_red" (lower)
mk = string(mouse_key_col);
mk = lower(strtrim(mk));
mk = regexprep(mk,'[\/\\]+','_');
mk = regexprep(mk,'[\s\-]+','_');
mk = regexprep(mk,'__+','_');

% If someone has "6100red" without separator, try to split into digits + letters
% (only if it starts with 4 digits)
for i=1:numel(mk)
    s = mk(i);
    if ~contains(s,'_')
        m = regexp(s,'^(\d{4})([a-z]+)$','tokens','once');
        if ~isempty(m)
            mk(i) = string(m{1}) + "_" + string(m{2});
        end
    end
end
end

%% ===================== PLOTTING / METRICS HELPERS =====================
function out = firstMatch(T, pat)
v = T.Properties.VariableNames;
hit = find(contains(v,pat),1,'first');
if isempty(hit), out = ''; else, out = v{hit}; end
end

function tf = hasVar(T,var)
tf = ismember(var, T.Properties.VariableNames);
end

function r = finiteRange(x)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x)-min(x); end
end

function dt = medianDT(tb)
x = double(tb(:)); d = diff(x); d(~isfinite(d)) = [];
dt = median(d,'omitnan'); if ~isfinite(dt) || dt<=0, dt = 1/30; end
end

function tb = pickTimebase(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
tb = nan(height(T),1);
for i=1:numel(cands)
    nm = cands{i};
    if hasVar(T,nm)
        v = double(T.(nm));
        if any(isfinite(v)), tb = v; return; end
    end
end
if hasVar(T,'Frame'), tb = double(T.Frame)/30; end
end

function [n, meanDur, totalDur, medianIEI] = eventMetrics(tb_s, ttl)
ttl = logical(ttl(:)); tb_s = double(tb_s(:));
if numel(tb_s) < 2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
d = diff([false; ttl; false]); onIdx = find(d==1); offIdx = find(d==-1)-1;
n = numel(onIdx); if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
dt = [diff(tb_s); median(diff(tb_s),'omitnan')]; segDur = zeros(n,1);
for k=1:n, ii = onIdx(k):offIdx(k); segDur(k) = sum(dt(ii)); end
meanDur  = mean(segDur,'omitnan'); totalDur = sum(segDur,'omitnan');
onTimes = tb_s(onIdx); if numel(onTimes)>=2, medianIEI = median(diff(onTimes),'omitnan'); else, medianIEI = NaN; end
end

function plotAllMiceLines(D, ycol, ylab, cols, ls)
if nargin<5, ls='-'; end
hold on
mice = unique(D.mouse_key_norm,'stable');
for i=1:numel(mice)
    sub = D(D.mouse_key_norm==mice(i),:); [~,ord] = sort(double(sub.day_index)); sub = sub(ord,:);
    if ~hasVar(sub,ycol), continue; end
    xx = double(sub.day_index); yy = double(sub.(ycol));
    good = isfinite(xx) & isfinite(yy); if ~any(good), continue; end
    plot(xx(good), yy(good), 'LineStyle',ls, 'Marker','o', 'LineWidth',0.7, 'MarkerSize',3.0, ...
         'Color',[0.7 0.7 0.7]);
end
% mean ± SEM across mice by day
ds = groupsummary(D, 'day_index', {'mean'}, ycol);
es = groupsummary(D, 'day_index', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ycol);
x = double(ds.day_index); y = double(ds.("mean_"+ycol)); e = double(es.("fun1_"+ycol));
if ~isempty(x)
    xv = [x; flipud(x)]; yv = [y-e; flipud(y+e)];
    patch('XData',xv,'YData',yv,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.6);
    plot(x, y, 'k-', 'LineWidth',1.8);
end
xlabel('Day'); ylabel(ylab); grid on; box on; applyRobustYLim(gca);
end

function savepng(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',150); catch, print(fh, fn, '-dpng','-r150'); end
end

function plotGroupVsDay_split(D, figDir, tag)
metrics = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
           'lick_meanDur_s','rew_meanDur_s','lick_totalDur_s','rew_totalDur_s', ...
           'lick_medianIEI_s','rew_medianIRI_s','pupil_mean'};
labels  = {'Requirement','Lick freq (/min)','Reward freq (/min)', ...
           'Lick mean dur (s)','Reward mean dur (s)','Lick total dur (s)','Reward total dur (s)', ...
           'Lick median IEI (s)','Reward median IRI (s)','Pupil mean (px)'};

keep = metrics(ismember(metrics, D.Properties.VariableNames));

% ACTIVE
sub = D(D.isPassive==0,:);
if ~isempty(sub)
    f = figure('Color','w');
    tl = tiledlayout(ceil(numel(keep)/2),2,'TileSpacing','compact','Padding','compact');
    for i=1:numel(keep)
        nexttile; hold on
        ycol = keep{i};
        ds = groupsummary(sub, 'day_index', {'mean'}, ycol);
        es = groupsummary(sub, 'day_index', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ycol);
        x = ds.day_index; y = ds.("mean_"+ycol); e = es.("fun1_"+ycol);
        shadeRibbon(x,y,e); xlabel('Day'); ylabel(labels{i}); grid on; box on; applyRobustYLim(gca);
    end
    title(tl,['ACTIVE: means ± SEM across days ' tag]);
    savepng(f, fullfile(figDir,['group_active_by_day' tag '.png']));
    close(f);
end

% PASSIVE
sub = D(D.isPassive==1,:);
if ~isempty(sub)
    f = figure('Color','w');
    tl = tiledlayout(ceil(numel(keep)/2),2,'TileSpacing','compact','Padding','compact');
    for i=1:numel(keep)
        nexttile; hold on
        ycol = keep{i};
        ds = groupsummary(sub, 'day_index', {'mean'}, ycol);
        es = groupsummary(sub, 'day_index', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ycol);
        x = ds.day_index; y = ds.("mean_"+ycol); e = es.("fun1_"+ycol);
        shadeRibbon(x,y,e); xlabel('Day'); ylabel(labels{i}); grid on; box on; applyRobustYLim(gca);
    end
    title(tl,['PASSIVE: means ± SEM across days ' tag]);
    savepng(f, fullfile(figDir,['group_passive_by_day' tag '.png']));
    close(f);
end
end

function shadeRibbon(x,y,e)
x = double(x); y = double(y); e = double(e); e(~isfinite(e))=0;
fill([x; flipud(x)], [y-e; flipud(y+e)], [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.5);
plot(x,y,'-o','LineWidth',1.6,'MarkerSize',5,'Color',[0.2 0.2 0.2]);
applyRobustYLim(gca);
end

function applyRobustYLim(ax)
if nargin<1, ax=gca; end
ln = findobj(ax,'Type','line');
vals = [];
for i=1:numel(ln)
    yi = get(ln(i),'YData'); yi = yi(:); yi = yi(isfinite(yi));
    vals = [vals; yi]; %#ok<AGROW>
end
if numel(vals) < 2, return; end
lo = prctile(vals,2); hi = prctile(vals,98);
if ~isfinite(lo) || ~isfinite(hi) || hi<=lo, return; end
pad = 0.08*(hi-lo);
try, ylim(ax, [lo - pad, hi + pad]); catch, end
end

function plotPupilCorrelations(D, outPath, tag)
f = figure('Color','w'); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Pupil vs Lick
nexttile; hold on
if hasVar(D,'pupil_mean') && hasVar(D,'lick_freq_per_min')
    good = isfinite(D.pupil_mean) & isfinite(D.lick_freq_per_min);
    scatter(D.pupil_mean(good), D.lick_freq_per_min(good), 12, 'filled','MarkerFaceAlpha',0.6);
    if any(good)
        R = corrcoef(D.pupil_mean(good), D.lick_freq_per_min(good)); r1 = R(1,2);
        ttl1 = sprintf('%s: Pupil vs Lick (r=%.2f)', tag, r1);
    else
        ttl1 = sprintf('%s: Pupil vs Lick', tag);
    end
    xlabel('Mean pupil (px)'); ylabel('Lick freq (/min)'); title(ttl1); grid on; box on; applyRobustYLim(gca);
else
    axis off
end

% Pupil vs Reward
nexttile; hold on
if hasVar(D,'pupil_mean') && hasVar(D,'rew_freq_per_min')
    good = isfinite(D.pupil_mean) & isfinite(D.rew_freq_per_min);
    scatter(D.pupil_mean(good), D.rew_freq_per_min(good), 12, 'filled','MarkerFaceAlpha',0.6);
    if any(good)
        R = corrcoef(D.pupil_mean(good), D.rew_freq_per_min(good)); r2 = R(1,2);
        ttl2 = sprintf('%s: Pupil vs Reward (r=%.2f)', tag, r2);
    else
        ttl2 = sprintf('%s: Pupil vs Reward', tag);
    end
    xlabel('Mean pupil (px)'); ylabel('Reward freq (/min)'); title(ttl2); grid on; box on; applyRobustYLim(gca);
else
    axis off
end

savepng(f, outPath); close(f);
end

%% ---------------- Per-trial + peri-event ----------------
function TT = computePerTrial(T, tb, hasLick, hasReward, hasPupil)
g = unique(T(:,{'mouse_key_norm','day_index','day_name','session_idx'}),'rows','stable');
out = {};
for i=1:height(g)
    mk = g.mouse_key_norm(i); di = g.day_index(i); ss = g.session_idx(i);
    rs = (T.mouse_key_norm==mk & T.day_index==di & T.session_idx==ss);
    if ~any(rs) || ~hasVar(T,'Trial'), continue; end
    Trial = double(T.Trial(rs)); good  = isfinite(Trial) & Trial>0; if ~any(good), continue; end

    tb_s  = tb(rs); tb_s = tb_s(good); tr = Trial(good);

    passFlag = NaN;
    if hasVar(T,'isPassive')
        ip = double(T.isPassive(rs)); ip = ip(good);
        if any(isfinite(ip)), passFlag = mode(round(ip(isfinite(ip)))); end
    end

    if hasLick,   lick = T.Lick_TTL(rs);     lick = lick(good); else, lick=false(size(tr)); end
    if hasReward, rew  = T.Injector_TTL(rs); rew  = rew(good);  else, rew=false(size(tr)); end
    if hasPupil,  pd   = double(T.Diameter_px(rs)); pd = pd(good); else, pd = nan(size(tr)); end

    trials = sort(unique(tr));
    for k = 1:numel(trials)
        mask = (tr==trials(k)); tbk = tb_s(mask); lk=lick(mask); rw=rew(mask); pdK = pd(mask);
        [ln,~,~,lIEI] = eventMetrics(tbk, lk);
        [rn,~,rDur,~] = eventMetrics(tbk, rw);
        trialDur = finiteRange(tbk);
        out(end+1,:) = {mk, di, g.day_name(i), ss, trials(k), passFlag, trialDur, ln, lIEI, rn, rDur, mean(pdK,'omitnan')}; %#ok<AGROW>
    end
end

if isempty(out)
    TT = table('Size',[0 12], ...
        'VariableTypes',{'string','double','string','double','double','double','double','double','double','double','double','double'}, ...
        'VariableNames',{'mouse_key_norm','day_index','day_name','session_idx','Trial','isPassive','trial_duration_s','lick_count','lick_medianIEI_s','rew_count','rew_dur_s','pupil_mean'});
    return;
end

TT = cell2table(out, 'VariableNames',{'mouse_key_norm','day_index','day_name','session_idx','Trial','isPassive','trial_duration_s','lick_count','lick_medianIEI_s','rew_count','rew_dur_s','pupil_mean'});
end

function plotPerTrialAll(TT, outPath)
if isempty(TT), return; end
M = groupsummary(TT, 'Trial', 'median', {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
E = groupsummary(TT, 'Trial', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
f = figure('Color','w'); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
plotErr(nexttile, M.Trial, M.median_trial_duration_s,  E.fun1_trial_duration_s,  'Trial duration (s)');
plotErr(nexttile, M.Trial, M.median_lick_count,        E.fun1_lick_count,        'Licks per trial');
plotErr(nexttile, M.Trial, M.median_rew_count,         E.fun1_rew_count,         'Rewards per trial');
plotErr(nexttile, M.Trial, M.median_rew_dur_s,         E.fun1_rew_dur_s,         'Reward duration per trial (s)');
plotErr(nexttile, M.Trial, M.median_pupil_mean,        E.fun1_pupil_mean,        'Pupil mean per trial (px)');
title(tl,'Per-trial trajectories (median ± SEM across mice)'); savepng(f, outPath); close(f);
end

function plotPerTrialByGroup(TT, figDir)
if ~ismember('isPassive', TT.Properties.VariableNames) || isempty(TT), return; end
labels = {'ACTIVE','PASSIVE'};
for gp = [0 1]
    sub = TT(TT.isPassive==gp, :);
    if isempty(sub), continue; end
    M = groupsummary(sub, 'Trial', 'median', {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
    E = groupsummary(sub, 'Trial', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
    f = figure('Color','w'); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    plotErr(nexttile, M.Trial, M.median_trial_duration_s,  E.fun1_trial_duration_s,  'Trial duration (s)');
    plotErr(nexttile, M.Trial, M.median_lick_count,        E.fun1_lick_count,        'Licks per trial');
    plotErr(nexttile, M.Trial, M.median_rew_count,         E.fun1_rew_count,         'Rewards per trial');
    plotErr(nexttile, M.Trial, M.median_rew_dur_s,         E.fun1_rew_dur_s,         'Reward duration per trial (s)');
    plotErr(nexttile, M.Trial, M.median_pupil_mean,        E.fun1_pupil_mean,        'Pupil mean per trial (px)');
    title(tl, ['Per-trial trajectories (median ± SEM) — ' labels{gp+1}]);
    savepng(f, fullfile(figDir, ['per_trial_' lower(labels{gp+1}) '.png'])); close(f);
end
end

function plotErr(ax, x, y, e, yl)
axes(ax); hold on
if isempty(x) || isempty(y), return; end
x = double(x); y = double(y); e = double(e); e(~isfinite(e)) = 0;
fill([x; flipud(x)], [y-e; flipud(y+e)], [0.9 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.45);
plot(x, y, '-o', 'LineWidth',1.2, 'MarkerSize',4, 'Color',[0.2 0.2 0.2]);
grid on; box on; xlabel('Trial'); ylabel(yl); applyRobustYLim(gca);
end

function periEventPlot(T, tb, ttlName, win, baseW, dt, figDir, outName)
if ~hasVar(T,ttlName) || ~hasVar(T,'Diameter_px'), return; end
f = figure('Color','w'); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gp = 0:1
    if hasVar(T,'isPassive'), rowsG = (double(T.isPassive)==gp); else, rowsG = true(height(T),1); end
    sk = unique(T(rowsG,{'mouse_key_norm','day_index','session_idx'}),'rows','stable');
    allM = []; tAxis = [];
    for i=1:height(sk)
        mk=sk.mouse_key_norm(i); di=sk.day_index(i); ss=sk.session_idx(i);
        rr = rowsG & (T.mouse_key_norm==mk) & (T.day_index==di) & (T.session_idx==ss);
        if ~any(rr), continue; end
        t = tb(rr); y = double(T.Diameter_px(rr)); ttl = logical(T.(ttlName)(rr));
        if ~any(isfinite(y)), continue; end
        onTimes = ttlOnsets(t, ttl, 0.5);
        if isempty(onTimes), continue; end
        [tAx, ~, ~, ~, matSess] = periEventAverage(t, y, onTimes, win, baseW, dt);
        if isempty(tAx), continue; end
        tAxis = tAx; allM = [allM; matSess.']; %#ok<AGROW>
    end
    ax=nexttile; hold(ax,'on')
    if ~isempty(allM)
        mu  = mean(allM,1,'omitnan');
        se  = std(allM,0,1,'omitnan')./sqrt(sum(isfinite(allM),1));
        fill([tAxis; flipud(tAxis)], [mu'-se'; flipud(mu'+se')], [0.9 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.45);
        plot(tAxis, mu, 'LineWidth',1.8, 'Color',[0.2 0.2 0.2]);
        yline(0,'-','Color',[0.6 0.6 0.6]); xlabel('Time from event (s)'); ylabel('\Delta pupil (px)');
        title(sprintf('%s — %s', ttlPretty(ttlName), ternary(gp,'Passive','Active'))); grid on; box on; applyRobustYLim(ax);
    else
        axis(ax,'off'); text(0.5,0.5,'No events','Units','normalized','HorizontalAlignment','center');
    end
end
title(tl, sprintf('Peri-event pupil \\Delta  [%s]', ttlPretty(ttlName)));
savepng(f, fullfile(figDir,outName)); close(f);
end

function lab = ttlPretty(nm)
if strcmpi(nm,'Lick_TTL'), lab='Lick onsets';
elseif strcmpi(nm,'Injector_TTL'), lab='Reward onsets';
else, lab=nm; end
end

function on = ttlOnsets(t, ttl, refractory_s)
d = diff([false; ttl(:); false]); onIdx = find(d==1); on = double(t(onIdx));
if isempty(on), return; end
on = sort(on);
if refractory_s>0
    keep = true(size(on));
    for k=2:numel(on)
        if on(k)-on(k-1)<refractory_s, keep(k)=false; end
    end
    on = on(keep);
end
end

function [tAxis, m, e, n, M] = periEventAverage(t, pupil, eventTimes, win, baseW, dt)
pre = abs(win(1)); post = win(2); tAxis = (-pre:dt:post).';
if isempty(tAxis), m=[]; e=[]; n=0; M=[]; return; end
good = isfinite(t) & isfinite(pupil); if ~any(good), m=[]; e=[]; n=0; M=[]; return; end
ti = t(good); yi = pupil(good); F = griddedInterpolant(ti, yi, 'linear', 'nearest');
tmin = min(ti); tmax = max(ti);
keepEv = eventTimes + win(1) >= tmin & eventTimes + win(2) <= tmax; ev = eventTimes(keepEv);
if isempty(ev), m=[]; e=[]; n=0; M=[]; return; end
nEv = numel(ev); M = nan(numel(tAxis), nEv); bMask = (tAxis>=baseW(1) & tAxis<baseW(2));
for k=1:nEv
    tr = ev(k) + tAxis; yk = F(tr); b = mean(yk(bMask),'omitnan'); M(:,k) = yk - b;
end
m = mean(M,2,'omitnan'); e = std(M,0,2,'omitnan') ./ sqrt(sum(isfinite(M),2)); n = nEv;
end

function out = ternary(cond, a, b)
if cond, out=a; else, out=b; end
end

%% ===================== STATS HELPERS =====================
function drawSigBar(ax, x1, x2, yTop, p)
if ~isfinite(p) || p>=0.05, return; end
hold(ax,'on');
yl = ylim(ax); dy = 0.04*range(yl);
plot(ax, [x1 x1 x2 x2],[yTop yTop+dy yTop+dy yTop], 'k-','LineWidth',1.2);
text(ax, mean([x1 x2]), yTop+dy*1.05, sigStarFromP(p), 'HorizontalAlignment','center','FontWeight','bold');
end

function s = sigStarFromP(p)
if ~isfinite(p), s=''; return; end
if p<0.001, s='***';
elseif p<0.01, s='**';
elseif p<0.05, s='*';
else, s='';
end
end

function [p, out] = safeRanksum(xa, xp)
xa = xa(isfinite(xa)); xp = xp(isfinite(xp));
p = NaN; delta = NaN;
if isempty(xa) || isempty(xp)
    out = struct('delta',delta,'es_str','n.s.'); return;
end
try
    p = ranksum(xa, xp);
catch
    try
        [~,p] = ttest2(xa, xp, 'Vartype','unequal');
    catch
        p = NaN;
    end
end
delta = cliffsDelta(xa, xp);
es_str = sprintf('\\Delta=%.2f', delta);
if isfinite(p)
    if p < 0.001, es_str = [es_str, ', ***'];
    elseif p < 0.01, es_str = [es_str, ', **'];
    elseif p < 0.05, es_str = [es_str, ', *'];
    end
end
out = struct('delta',delta,'es_str',es_str);
end

function d = cliffsDelta(a, b)
a = a(:); b = b(:);
a = a(isfinite(a)); b = b(isfinite(b));
if isempty(a) || isempty(b), d = NaN; return; end
cnt = 0; gt = 0; lt = 0;
for i=1:numel(a)
    gt = gt + sum(a(i) > b);
    lt = lt + sum(a(i) < b);
    cnt = cnt + numel(b);
end
d = (gt - lt) / cnt;
end

function [p_fdr, sig_mask] = fdr_bh_local(pvals, alpha)
if nargin<2, alpha=0.05; end
p = pvals(:); n = numel(p); p(isnan(p)) = 1;
[ps, idx] = sort(p);
thresh = (1:n)'/n * alpha;
sig_sorted = ps <= thresh;
if any(sig_sorted)
    k = find(sig_sorted,1,'last');
    cutoff = thresh(k);
else
    cutoff = 0;
end
sig_mask = p <= cutoff;
adj = nan(size(p));
adj(idx) = min( cummin( (n./(1:n)') .* ps ), 1);
p_fdr = adj(:);
end

function p = kruskal_p_fallback(groups)
try
    allx = []; grp = [];
    for i=1:numel(groups)
        xi = groups{i}; xi = xi(isfinite(xi));
        allx = [allx; xi(:)];
        grp  = [grp; repmat(i,numel(xi),1)];
    end
    if isempty(allx) || numel(unique(grp))<2, p=NaN; return; end
    p = kruskalwallis(allx, grp, 'off');
catch
    k = numel(groups);
    allx = []; grp = [];
    for i=1:k
        xi = groups{i}; xi = xi(isfinite(xi));
        allx = [allx; xi(:)];
        grp  = [grp; repmat(i,numel(xi),1)];
    end
    if isempty(allx) || numel(unique(grp))<2, p=NaN; return; end
    [r, ~] = tiedrank(allx);
    n = numel(allx);
    N_i = accumarray(grp,1,[k,1]);
    R_i = accumarray(grp,r,[k,1]);
    H = (12/(n*(n+1))) * sum( (R_i.^2) ./ N_i ) - 3*(n+1);
    df = k-1;
    p = 1 - chi2cdf(H, df);
end
end

function plotTwoGroupsDots(Xa, Xp, ylab, xlabRight)
if nargin<4, xlabRight=''; end
hold on
try
    swarmchart(ones(size(Xa))*1, Xa, 16, 'filled','MarkerFaceAlpha',0.6);
    swarmchart(ones(size(Xp))*2, Xp, 16, 'filled','MarkerFaceAlpha',0.6);
catch
    jitter = @(n) (rand(n,1)-0.5)*0.25;
    scatter(1 + jitter(numel(Xa)), Xa, 16, 'filled','MarkerFaceAlpha',0.6);
    scatter(2 + jitter(numel(Xp)), Xp, 16, 'filled','MarkerFaceAlpha',0.6);
end
m1 = mean(Xa,'omitnan'); s1 = std(Xa,'omitnan')/sqrt(max(1,sum(isfinite(Xa))));
m2 = mean(Xp,'omitnan'); s2 = std(Xp,'omitnan')/sqrt(max(1,sum(isfinite(Xp))));
errorbar([1 2],[m1 m2],[s1 s2],'k_','LineWidth',1.4);
xlim([0.5 2.5]); grid on; box on
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel(ylab); if ~isempty(xlabRight), xlabel(xlabRight); end
applyRobustYLim(gca);
end

function plotCorrMatrix(Tin, vars, outPath, tag)
X = []; names = {};
for i=1:numel(vars)
    if ismember(vars{i}, Tin.Properties.VariableNames)
        v = double(Tin.(vars{i}));
        if any(isfinite(v))
            X = [X, v(:)]; %#ok<AGROW>
            names{end+1} = vars{i}; %#ok<AGROW>
        end
    end
end
if size(X,2) < 2, return; end
good = all(isfinite(X),2); X = X(good,:); if size(X,1) < 3, return; end

R = corrcoef(X);
n = size(X,1);
P = ones(size(R));
for i=1:size(R,1)
    for j=1:size(R,2)
        if i==j, P(i,j)=0; continue; end
        r = max(min(R(i,j),0.999999),-0.999999);
        t = r * sqrt( (n-2) / max(1e-12, 1 - r^2) );
        cdf = tcdf_local(abs(t), n-2);
        P(i,j) = 2*(1 - cdf);
    end
end

fig = figure('Color','w','Position',[100 100 560 520]);
imagesc(R, [-1 1]); axis square; colorbar
title(sprintf('Pearson r — %s', tag));
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',45);
set(gca,'YTick',1:numel(names),'YTickLabel',names);
for i=1:numel(names)
    for j=1:numel(names)
        if i==j, continue; end
        p = P(i,j);
        text(j,i, sprintf('%.2f%s', R(i,j), sigStarFromP(p)), ...
            'HorizontalAlignment','center','Color','k','FontSize',9,'FontWeight','bold');
    end
end
savepng(fig, outPath); close(fig);
end

function F = tcdf_local(t, nu)
t = abs(t);
x = nu ./ (nu + t.^2);
F = 1 - 0.5*betainc(x, nu/2, 0.5);
end

function T = ensureString(T, vname)
if ~hasVar(T,vname), return; end
if ~isstring(T.(vname)), T.(vname) = string(T.(vname)); end
end

function k = buildKey(T, keys)
k = repmat("", height(T), 1);
for i = 1:numel(keys)
    v = T.(keys{i});
    if ~isstring(v), v = string(v); end
    if i == 1, k = v;
    else, k = k + "|" + v;
    end
end
end

%% ===================== STRING-SAFE MODE (FIX) =====================
function s = modeString(x)
% Return most frequent (non-missing) string. Tie -> first (stable).
x = string(x);
x = x(strlength(x)>0 & ~ismissing(x));
if isempty(x)
    s = "";
    return;
end
[u, ~, ic] = unique(x, 'stable');
cnt = accumarray(ic, 1);
[~, k] = max(cnt);
s = u(k);
end
