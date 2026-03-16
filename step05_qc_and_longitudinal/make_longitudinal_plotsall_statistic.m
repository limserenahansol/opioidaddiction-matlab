function make_longitudinal_plotsall_statistic() 
% All-mice longitudinal plots with Week/Group splits + STATS, BIN TESTS & INDEX.
% Reads latest longitudinal_outputs/run_*/ALL_mice_longitudinal.csv
% Outputs to run_*/figs, run_*/figs/weeks, and run_*/figs/stats

%% -------- locate latest run and read CSV --------
tryPath = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(tryPath,'dir')
    here = pwd; cand = here;
    for up=1:5
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
T = readtable(csvPath, 'VariableNamingRule','preserve');
% --- PUPIL-ONLY exclusion: 7597 black ---
if hasVar(T,'mouse_key') && hasVar(T,'Diameter_px')
    mk      = string(T.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));  % normalize "7597_black" etc.
    mask7597 = contains(mk_norm,"7597") & contains(mk_norm,"black");
    if any(mask7597)
        fprintf('Pupil-only exclude: setting Diameter_px=NaN for %d rows (7597 black).\n', nnz(mask7597));
        T.Diameter_px(mask7597) = NaN;   % only pupil; lick/reward untouched
    end
end

% Ensure types
T = ensureString(T,'day_name');
T = ensureString(T,'mouse_key');
if hasVar(T,'Session_Paradigm'), T = ensureString(T,'Session_Paradigm'); end
if hasVar(T,'isPassive') && ~isnumeric(T.isPassive), T.isPassive = double(T.isPassive); end
if ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if hasVar(T,'session_idx') && ~isnumeric(T.session_idx)
    T.session_idx = double(T.session_idx);
end

% Create figs folders
figDir      = fullfile(runDir, 'figs'); if ~exist(figDir,'dir'), mkdir(figDir); end
figDirWeeks = fullfile(figDir, 'weeks'); if ~exist(figDirWeeks,'dir'), mkdir(figDirWeeks); end
fprintf('Saving figs to:\n  %s\n  %s\n', figDir, figDirWeeks);

%% -------- drop frames from TTL-only "extra" trials --------
if hasVar(T,'TrialRequirement')
    trr = T.TrialRequirement; if ~isstring(trr), trr = string(trr); end
    extraMask = false(height(T),1);
    if hasVar(T,'Trial'), extraMask = (T.Trial>0) & strcmpi(strtrim(trr),"extra");
    else, extraMask = strcmpi(strtrim(trr),"extra");
    end
    if any(extraMask)
        fprintf('Excluding %d frames with TrialRequirement="extra" (%.2f%% of rows).\n', ...
            nnz(extraMask), 100*nnz(extraMask)/height(T));
        T(extraMask,:) = [];
    end
end

%% -------- helpers & TTL normalization --------
tb = pickTimebase(T);
hasLick   = hasVar(T,'Lick_TTL');
hasReward = hasVar(T,'Injector_TTL');
hasPupil  = hasVar(T,'Diameter_px');

if hasLick
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0;  T.Lick_TTL = T.Lick_TTL > 0.5;
end
if hasReward
    T.Injector_TTL(isnan(T.Injector_TTL)) = 0;  T.Injector_TTL = T.Injector_TTL > 0.5;
end

%% -------- per-session metrics --------
sessKeys = unique(T(:,{'mouse_key','day_index','day_name','session_idx'}),'rows','stable');

S = table();
S.mouse_key   = sessKeys.mouse_key;
S.day_index   = sessKeys.day_index;
S.day_name    = sessKeys.day_name;
S.session_idx = sessKeys.session_idx;

S.RequirementLast = nan(height(S),1);
S.isPassive       = nan(height(S),1);
S.SessionMinutes  = nan(height(S),1);
[S.lick_n,S.lick_freq_per_min,S.lick_meanDur_s,S.lick_totalDur_s,S.lick_medianIEI_s] = deal(nan(height(S),1));
[S.rew_n,S.rew_freq_per_min,S.rew_meanDur_s,S.rew_totalDur_s,S.rew_medianIRI_s]       = deal(nan(height(S),1));
S.pupil_mean = nan(height(S),1);

for i=1:height(sessKeys)
    mk  = sessKeys.mouse_key(i); di  = sessKeys.day_index(i); ss  = sessKeys.session_idx(i);
    rows = T.mouse_key==mk & T.day_index==di & T.session_idx==ss; if ~any(rows), continue; end
    tb_s = tb(rows);
    sessionDur_s = finiteRange(tb_s);
    if ~(isfinite(sessionDur_s) && sessionDur_s>0) && hasVar(T,'Frame')
        fr = double(T.Frame(rows)); if any(isfinite(fr)), sessionDur_s = (max(fr)-min(fr))/30; end
    end
    S.SessionMinutes(i) = sessionDur_s/60;

    if hasVar(T,'RequirementLast'), S.RequirementLast(i) = mean(double(T.RequirementLast(rows)),'omitnan'); end
    if hasVar(T,'isPassive')
        ip = double(T.isPassive(rows)); ip = ip(isfinite(ip));
        if ~isempty(ip), S.isPassive(i) = mode(round(ip)); end
    end
    if hasPupil, S.pupil_mean(i) = mean(double(T.Diameter_px(rows)), 'omitnan'); end

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

%% -------- collapse per-session → per-day (median across sessions) --------
dayKeys = unique(T(:,{'mouse_key','day_index','day_name'}),'rows','stable');
D = dayKeys; agg = @(x) median(x,'omitnan');

coreCols = {'RequirementLast','isPassive','SessionMinutes', ...
           'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
           'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s','pupil_mean'};

for c = 1:numel(coreCols), D.(coreCols{c}) = nan(height(D),1); end
for i=1:height(D)
    mk = D.mouse_key(i); di = D.day_index(i);
    r  = S.mouse_key==mk & S.day_index==di; if ~any(r), continue; end
    for c = 1:numel(coreCols), nm = coreCols{c}; D.(nm)(i) = agg(S.(nm)(r)); end
end

% ------- bring Immersion/TST/HOT (median per day across sessions) --------
D = addDayScalarFromT(D, T, 'Immersion_Latency_s');        % adds if present
D = addDayPatternFromT(D, T, 'TST_Pct_');                  % any/all TST % columns
D = addDayPatternFromT(D, T, 'HOT_Pct_');                  % any/all HOT % columns

% ------- Week & Group helpers -------
D.week = zeros(height(D),1);
D.week(D.day_index>=1 & D.day_index<=5)  = 1;
D.week(D.day_index>=6 & D.day_index<=11) = 2;

%% ===================== PLOTS (All days, then Week×Group) =====================

% ---------- 0) Requirement & TTL metrics ----------
cols = lines(numel(unique(D.mouse_key))); % per-mouse lines (muted)
if hasVar(D,'RequirementLast')
    f = figure('Color','w'); hold on
    plotAllMiceLines(D, 'RequirementLast', 'Requirement (last in JSONL)', cols, '-');
    title('Requirement by day (per mouse)'); savepng(f, fullfile(figDir,'req_by_day_all_mice.png')); close(f);
end

if hasLick || hasReward
    f = figure('Color','w'); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile; hold(ax,'on');
    plotAllMiceLines(D, 'lick_freq_per_min', 'Lick frequency (/min)', cols, '-');
    plotAllMiceLines(D, 'rew_freq_per_min',  'Reward frequency (/min)', cols, '--'); title(ax,'Frequencies');
    legend(ax,{'Lick (per-mouse)','Reward (per-mouse)'},'Location','southoutside','NumColumns',2,'Box','off');

    ax=nexttile; hold(ax,'on');
    plotAllMiceLines(D, 'lick_meanDur_s', 'Mean duration (s)', cols, '-');
    plotAllMiceLines(D, 'rew_meanDur_s',  'Mean duration (s)', cols, '--'); title(ax,'Mean durations');

    ax=nexttile; hold(ax,'on');
    plotAllMiceLines(D, 'lick_totalDur_s', 'Total duration (s/session)', cols, '-');
    plotAllMiceLines(D, 'rew_totalDur_s',  'Total duration (s/session)', cols, '--'); title(ax,'Totals');

    ax=nexttile; hold(ax,'on');
    plotAllMiceLines(D, 'lick_medianIEI_s', 'Median IEI / IRI (s)', cols, '-');
    plotAllMiceLines(D, 'rew_medianIRI_s',  'Median IEI / IRI (s)', cols, '--'); title(ax,'Intervals');

    title(tl, 'TTL-derived metrics across days (solid=lick, dashed=reward)');
    savepng(f, fullfile(figDir,'lick_reward_metrics_by_day.png')); close(f);
end
% LICK-only and REWARD-only panels
plotMetricsByDaySeparate(D, figDir);

% ---------- 1) Immersion/TST/HOT (if present) ----------
immName = firstMatch(D,'Immersion_Latency_s');
tstName = firstMatch(D,'TST_Pct_');
hotName = firstMatch(D,'HOT_Pct_');
if ~isempty(immName) || ~isempty(tstName) || ~isempty(hotName)
    f = figure('Color','w'); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    if ~isempty(immName), ax=nexttile; hold(ax,'on'); plotAllMiceLines(D, immName, 'Immersion latency (s)', cols, '-'); end
    if ~isempty(tstName), ax=nexttile; hold(ax,'on'); plotAllMiceLines(D, tstName, ['TST % ' stripPrefix(tstName)], cols, '-'); end
    if ~isempty(hotName), ax=nexttile; hold(ax,'on'); plotAllMiceLines(D, hotName, ['HOT % ' stripPrefix(hotName)], cols, '-'); end
    savepng(f, fullfile(figDir,'imm_tst_hot_by_day.png')); close(f);
end

% ---------- 2) Per-trial trajectories ----------
TT_all = table();
if hasVar(T,'Trial')
    TT_all = computePerTrial(T, tb, hasLick, hasReward, hasPupil);
    plotPerTrialAll(TT_all, fullfile(figDir,'per_trial_all_mice.png'));
end
% split by group (ACTIVE vs PASSIVE)
if ~isempty(TT_all) && ismember('isPassive', TT_all.Properties.VariableNames)
    plotPerTrialByGroup(TT_all, figDir);
end

% ---------- 3) Group means ± SEM vs day (active vs passive, all days) ----------
if any(D.isPassive==0 | D.isPassive==1)
    plotGroupVsDay_split(D, figDir, '');  % all days
end

% ---------- 4) Pupil ↔ lick/reward correlations (all days) ----------
if hasPupil && (hasLick || hasReward)
    plotPupilCorrelations(D, fullfile(figDir,'pupil_correlations_overview.png'), 'ALL DAYS');
end

% ---------- 5) Peri-event pupil Δ (all days) ----------
if hasPupil && (hasLick || hasReward)
    win   = [-2 6]; baseW = [-2 0]; dt = medianDT(tb);
    if hasLick,   periEventPlot(T, tb, 'Lick_TTL',     win, baseW, dt, figDir, 'pupil_peri_lick_active_vs_passive.png'); end
    if hasReward, periEventPlot(T, tb, 'Injector_TTL', win, baseW, dt, figDir, 'pupil_peri_reward_active_vs_passive.png'); end
    if hasLick,   periEventPlot_All(T, tb, 'Lick_TTL',     win, baseW, dt, figDir, 'pupil_peri_lick_ALL.png'); end
    if hasReward, periEventPlot_All(T, tb, 'Injector_TTL', win, baseW, dt, figDir, 'pupil_peri_reward_ALL.png'); end
end

% ----- OPTIONAL: per-parameter separate daily plots: ALL, ACTIVE, PASSIVE -----
metrics = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
           'lick_meanDur_s','rew_meanDur_s','lick_totalDur_s','rew_totalDur_s', ...
           'lick_medianIEI_s','rew_medianIRI_s','pupil_mean'};
labels  = {'Requirement','Lick freq (/min)','Reward freq (/min)', ...
           'Lick mean dur (s)','Reward mean dur (s)','Lick total dur (s)','Reward total dur (s)', ...
           'Lick median IEI (s)','Reward median IRI (s)','Pupil mean (px)'};
colsAll = lines(numel(unique(D.mouse_key)));

for i = 1:numel(metrics)
    y = metrics{i}; if ~ismember(y, D.Properties.VariableNames), continue; end
    lbl = labels{i};
    f = figure('Color','w'); hold on
    plotAllMiceLines(D, y, lbl, colsAll, '-');
    title([lbl ' — ALL mice (daily)']); savepng(f, fullfile(figDir, ['daily_' y '_ALL.png'])); close(f);
    if any(D.isPassive==0)
        sub = D(D.isPassive==0,:);
        f = figure('Color','w'); hold on
        plotAllMiceLines(sub, y, lbl, colsAll, '-');
        title([lbl ' — ACTIVE only (daily)']); savepng(f, fullfile(figDir, ['daily_' y '_ACTIVE.png'])); close(f);
    end
    if any(D.isPassive==1)
        sub = D(D.isPassive==1,:);
        f = figure('Color','w'); hold on
        plotAllMiceLines(sub, y, lbl, colsAll, '-');
        title([lbl ' — PASSIVE only (daily)']); savepng(f, fullfile(figDir, ['daily_' y '_PASSIVE.png'])); close(f);
    end
end

%% ===================== WEEK × GROUP (4 subsets) =====================
scopes = { ...
    struct('name','week1_active',  'week',1, 'pass',0), ...
    struct('name','week1_passive', 'week',1, 'pass',1), ...
    struct('name','week2_active',  'week',2, 'pass',0), ...
    struct('name','week2_passive', 'week',2, 'pass',1)};

% (A) TTL single-metric plots → 8 metrics × 4 scopes
metricList = { ...
    'lick_freq_per_min','Lick freq (/min)'; ...
    'rew_freq_per_min', 'Reward freq (/min)'; ...
    'lick_meanDur_s',   'Lick mean dur (s)'; ...
    'rew_meanDur_s',    'Reward mean dur (s)'; ...
    'lick_totalDur_s',  'Lick total dur (s/session)'; ...
    'rew_totalDur_s',   'Reward total dur (s/session)'; ...
    'lick_medianIEI_s', 'Lick median IEI (s)'; ...
    'rew_medianIRI_s',  'Reward median IRI (s)'};

for sc = 1:numel(scopes)
    subD = D( D.week==scopes{sc}.week & D.isPassive==scopes{sc}.pass, : );
    if isempty(subD), continue; end
    for m=1:size(metricList,1)
        v = metricList{m,1}; lbl = metricList{m,2};
        if ~hasVar(subD,v), continue; end
        f = figure('Color','w'); hold on
        plotAllMiceLines(subD, v, lbl, cols, '-');   % spaghetti
        ttl = sprintf('%s — %s', upper(strrep(scopes{sc}.name,'_',' ')), lbl);
        title(ttl);
        savepng(f, fullfile(figDirWeeks, sprintf('%s_%s.png', scopes{sc}.name, v))); close(f);
    end
end

% (B) Pupil mean by Week×Group
if hasVar(D,'pupil_mean')
    for sc = 1:numel(scopes)
        subD = D( D.week==scopes{sc}.week & D.isPassive==scopes{sc}.pass, : );
        if isempty(subD), continue; end
        f = figure('Color','w'); hold on
        plotAllMiceLines(subD, 'pupil_mean', 'Mean pupil (px)', cols, '-');
        title(sprintf('%s — Mean pupil (px)', upper(strrep(scopes{sc}.name,'_',' '))));
        savepng(f, fullfile(figDirWeeks, sprintf('%s_pupil_mean.png', scopes{sc}.name))); close(f);
    end
end

% (C) Immersion/TST/HOT by Week×Group (if present)
plotImmTstHot_ByWeekGroup(D, cols, figDirWeeks, scopes, immName, tstName, hotName);

% (D) Per-trial trajectories by Week×Group
if ~isempty(TT_all) && hasVar(TT_all,'Trial')
    for sc=1:numel(scopes)
        wk = scopes{sc}.week; ps = scopes{sc}.pass;
        TTg = TT_all;
        TTg.week = zeros(height(TTg),1);
        for i=1:height(TTg)
            di = TTg.day_index(i);
            if di>=1 && di<=5, TTg.week(i)=1; elseif di>=6 && di<=11, TTg.week(i)=2; else, TTg.week(i)=0; end
        end
        TTg = TTg(TTg.week==wk & TTg.isPassive==ps,:);
        if isempty(TTg), continue; end
        nameTag = sprintf('%s',strrep(scopes{sc}.name,'_',' '));
        plotPerTrialSimple(TTg, fullfile(figDirWeeks,sprintf('per_trial_%s.png', scopes{sc}.name)), upper(nameTag));
    end
end

% (E) Correlations: pupil vs selected metrics
if hasVar(D,'pupil_mean')
    corrTargets = {'lick_freq_per_min','lick_n','rew_n','rew_freq_per_min','RequirementLast'};
    corrLabels  = {'Lick freq (/min)','Lick count','Reward count','Reward freq (/min)','Requirement'};
    for sc=1:numel(scopes)
        subD = D( D.week==scopes{sc}.week & D.isPassive==scopes{sc}.pass, : );
        if isempty(subD), continue; end
        f = figure('Color','w'); tl = tiledlayout(1,numel(corrTargets),'TileSpacing','compact','Padding','compact');
        for j=1:numel(corrTargets)
            yv = corrTargets{j};
            if ~hasVar(subD,yv), nexttile; axis off; continue; end
            nexttile; hold on
            good = isfinite(subD.pupil_mean) & isfinite(subD.(yv));
            scatter(subD.pupil_mean(good), subD.(yv)(good), 10,'filled','MarkerFaceAlpha',0.6);
            if any(good)
                try
                    pfit = polyfit(subD.pupil_mean(good), subD.(yv)(good), 1);
                    xx = linspace(min(subD.pupil_mean(good)), max(subD.pupil_mean(good)), 100);
                    yy = polyval(pfit, xx); plot(xx,yy,'k-','LineWidth',1.2);
                catch
                end
                R = corrcoef(subD.pupil_mean(good), subD.(yv)(good));
                r = R(1,2);
                ttl = sprintf('%s (r=%.2f)', corrLabels{j}, r);
            else
                ttl = corrLabels{j};
            end
            xlabel('Mean pupil (px)'); ylabel(corrLabels{j}); title(ttl);
            grid on; box on; applyRobustYLim(gca);
        end
        title(tl, sprintf('Pupil correlations — %s', upper(strrep(scopes{sc}.name,'_',' '))));
        savepng(f, fullfile(figDirWeeks, sprintf('pupil_corr_%s.png', scopes{sc}.name))); close(f);
    end
end

%% ===================== NEW: STATS & CORRELATIONS SUITE =====================
statsOutDir = fullfile(figDir, 'stats'); 
if ~exist(statsOutDir,'dir'), mkdir(statsOutDir); end

metrics_for_stats = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
                     'lick_meanDur_s','rew_meanDur_s','lick_totalDur_s','rew_totalDur_s', ...
                     'lick_medianIEI_s','rew_medianIRI_s','pupil_mean'};

labels_for_stats  = {'Requirement','Lick freq (/min)','Reward freq (/min)', ...
                     'Lick mean dur (s)','Reward mean dur (s)','Lick total dur (s)','Reward total dur (s)', ...
                     'Lick median IEI (s)','Reward median IRI (s)','Pupil mean (px)'};

% (A) Any-Passive vs None (+ stars)
runStatsComparisons_and_Correlations(D, metrics_for_stats, labels_for_stats, statsOutDir);

% (B) Within-group BIN tests (Active-only & Had-passive) with FDR pairwise & heatmaps
withinGroup_BinANOVA_and_Pairs(D, 0, metrics_for_stats, labels_for_stats, statsOutDir); % Active-only mice
withinGroup_BinANOVA_and_Pairs(D, 1, metrics_for_stats, labels_for_stats, statsOutDir); % Had-passive mice

% (C) Active–Passive Index I=(A-P)/(A+P): per-day trend + bin stats
computeAndPlot_AP_Index(D, metrics_for_stats, labels_for_stats, statsOutDir);
% (A2) NEW: Active vs Passive (day-level) within each time bin
activePassive_byBin(D, metrics_for_stats, labels_for_stats, statsOutDir);

% (C2) NEW: Explicit AP-index bin-to-bin Δ plots + heatmap
plot_AP_index_bin_diffs(D, metrics_for_stats, labels_for_stats, statsOutDir);

fprintf('Done.\n');
end  % ===== end main =====


%% ================= helpers =================
function out = firstMatch(T, pat)
v = T.Properties.VariableNames;
hit = find(contains(v,pat),1,'first');
if isempty(hit), out = ''; else, out = v{hit}; end
end

function s = stripPrefix(nm)
s = regexprep(nm,'^(TST|HOT)_Pct_','');
s = strrep(s,'_',' ');
end

function tf = hasVar(T,var)
tf = ismember(var, T.Properties.VariableNames);
end

function r = finiteRange(x)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x)-min(x); end
end
function activePassive_byBin(D, metrics, labels, outDir)
% Day-level state comparison: Active (isPassive==0) vs Passive (isPassive==1)
% within each time bin, using per-mouse median within the bin.
if ~ismember('isPassive', D.Properties.VariableNames), return; end
D.bin5 = dayBin5(D.day_index);
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};

for i=1:numel(metrics)
    y = metrics{i}; if ~ismember(y, D.Properties.VariableNames), continue; end
    lbl = labels{i};

    % Build per-mouse, per-bin, per-state medians
    sub = D(~isundefined(D.bin5) & isfinite(D.(y)) & ismember(D.isPassive,[0 1]), ...
            {'mouse_key','bin5','isPassive',y});
    if isempty(sub), continue; end
    M = groupsummary(sub, {'mouse_key','bin5','isPassive'}, 'median', y);
    M.Properties.VariableNames(end) = {'y_med'};
    
    % Figure: panels = bins
    fig = figure('Color','w','Position',[80 80 1160 420]);
    tl  = tiledlayout(1,numel(binNames),'TileSpacing','compact','Padding','compact');
    csvRows = {};
    for b=1:numel(binNames)
        nexttile; hold on
        Mb = M(string(M.bin5)==binNames{b},:);
        if isempty(Mb), axis off; title(binNames{b}); continue; end
        Xa = Mb.y_med(Mb.isPassive==0); Xa = Xa(isfinite(Xa));
        Xp = Mb.y_med(Mb.isPassive==1); Xp = Xp(isfinite(Xp));
        plotTwoGroupsDots(Xa, Xp, lbl, binNames{b});
        [p,~] = safeRanksum(Xa, Xp);
        yl = ylim; ybar = yl(2) - 0.05*range(yl);
        drawSigBar(gca, 1, 2, ybar, p);
        csvRows(end+1,:) = {y,binNames{b},numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),p}; %#ok<AGROW>
    end
    title(tl, sprintf('%s — Active vs Passive (per bin; per-mouse median)', lbl));
    savepng(fig, fullfile(outDir, sprintf('stat_ACTIVEvsPASSIVE_BINS_%s.png', y)));
    close(fig);

    % CSV dump
    if ~isempty(csvRows)
        Tcsv = cell2table(csvRows, 'VariableNames', ...
            {'metric','bin','n_active','n_passive','median_active','median_passive','p_ranksum'});
        try, writetable(Tcsv, fullfile(outDir, sprintf('stat_ACTIVEvsPASSIVE_BINS_%s.csv', y))); catch, end
    end
end
end
function plot_AP_index_bin_diffs(D, metrics, labels, outDir)
% Build AP index per day, aggregate by bins,
% show ΔI between bins and an FDR-corrected heatmap with stars.

if ~ismember('isPassive', D.Properties.VariableNames), return; end
D.bin5 = dayBin5(D.day_index);
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};

for i=1:numel(metrics)
    y = metrics{i}; if ~ismember(y, D.Properties.VariableNames), continue; end
    lbl = labels{i};

    % Compute daily AP index I(d)
    days = sort(unique(D.day_index));
    I  = nan(size(days));
    for k=1:numel(days)
        di = days(k);
        A = D{D.day_index==di & D.isPassive==0, y}; A = A(isfinite(A));
        P = D{D.day_index==di & D.isPassive==1, y}; P = P(isfinite(P));
        if ~isempty(A) && ~isempty(P)
            a = median(A,'omitnan'); p = median(P,'omitnan');
            if isfinite(a) && isfinite(p) && (a+p)~=0
                I(k) = (a - p) / (a + p);
            end
        end
    end
    good = isfinite(I);
    if nnz(good) < 3, continue; end

    di_used = days(good); Ivec = I(good);
    binLab  = arrayfun(@(d) string(dayToBin(d)), di_used);
    % Distributions per bin + medians
    groups = cell(1,numel(binNames));
    medBin = nan(1,numel(binNames));
    for b=1:numel(binNames)
        xi = Ivec(binLab==binNames{b});
        xi = xi(isfinite(xi));
        groups{b} = xi;
        medBin(b) = median(xi,'omitnan');
    end

    % Pairwise tests on daily I + FDR; store Δ medians
    pairs = nchoosek(1:numel(binNames),2);
    p_pair = nan(size(pairs,1),1);
    dlt    = nan(size(pairs,1),1);
    for k=1:size(pairs,1)
        b1 = pairs(k,1); b2 = pairs(k,2);
        p_pair(k) = safeRanksum(groups{b1}, groups{b2});
        dlt(k)    = medBin(b2) - medBin(b1);
    end
    [p_adj, ~] = fdr_bh_local(p_pair, 0.05);

    % Heatmap of ΔI with stars
    H = nan(numel(binNames));
    for k=1:size(pairs,1)
        i1 = pairs(k,1); i2 = pairs(k,2);
        H(i1,i2) = dlt(k); H(i2,i1) = -dlt(k);
    end
    fig = figure('Color','w','Position',[140 140 560 480]);
    imagesc(H); axis square; colorbar; colormap(parula);
    caxis(max(abs(H(:))).*[-1 1]); % symmetric color scaling
    title(sprintf('\\Delta I between bins — %s', lbl));
    set(gca,'XTick',1:numel(binNames),'XTickLabel',binNames,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(binNames),'YTickLabel',binNames);
    for r=1:numel(binNames)
        for c=1:numel(binNames)
            if r==c, text(c,r,'—','HorizontalAlignment','center'); continue; end
            i1 = min(r,c); i2 = max(r,c);
            k  = find(pairs(:,1)==i1 & pairs(:,2)==i2,1);
            if isempty(k), continue; end
            st = sigStarFromP(p_adj(k));
            if ~isempty(st), text(c,r,st,'HorizontalAlignment','center','FontWeight','bold'); end
        end
    end
    savepng(fig, fullfile(outDir, sprintf('AP_index_binDelta_heatmap_%s.png', y)));
    close(fig);

    % CSV of ΔI and p-values
    pair_lbl = strings(size(pairs,1),1);
    for k=1:size(pairs,1)
        pair_lbl(k) = sprintf('%s vs %s', binNames{pairs(k,1)}, binNames{pairs(k,2)});
    end
    Tcsv = table(repmat(string(lbl),numel(p_adj),1), pair_lbl, dlt, p_pair, p_adj, ...
        'VariableNames', {'metric','pair','delta_median_I','p_raw','p_fdr'});
    try, writetable(Tcsv, fullfile(outDir, sprintf('AP_index_bin_pairs_%s.csv', y))); catch, end
end
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
% Thin gray spaghetti per mouse; bold black mean ± SEM
if nargin<5, ls='-'; end
hold on
mice = unique(D.mouse_key,'stable');
for i=1:numel(mice)
    sub = D(D.mouse_key==mice(i),:); [~,ord] = sort(double(sub.day_index)); sub = sub(ord,:);
    if ~hasVar(sub,ycol), continue; end
    xx = double(sub.day_index); yy = double(sub.(ycol));
    good = isfinite(xx) & isfinite(yy); if ~any(good), continue; end
    plot(xx(good), yy(good), 'LineStyle',ls, 'Marker','o', 'LineWidth',0.7, 'MarkerSize',3.0, ...
         'Color',[0.7 0.7 0.7]);
end
% group mean ± SEM across mice (by day)
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

function TT = computePerTrial(T, tb, hasLick, hasReward, hasPupil)
g = unique(T(:,{'mouse_key','day_index','day_name','session_idx'}),'rows','stable');
rowsPer = @(mk,di,ss) (T.mouse_key==mk & T.day_index==di & T.session_idx==ss);
out = {};
for i=1:height(g)
    mk = g.mouse_key(i); di = g.day_index(i); ss = g.session_idx(i);
    rs = rowsPer(mk,di,ss); if ~any(rs) || ~hasVar(T,'Trial'), continue; end
    Trial = double(T.Trial(rs)); good  = isfinite(Trial) & Trial>0; if ~any(good), continue; end
    tb_s  = tb(rs); tb_s = tb_s(good); tr = Trial(good);
    passFlag = NaN;
    if hasVar(T,'isPassive'), ip = double(T.isPassive(rs)); ip = ip(good);
        if any(isfinite(ip)), passFlag = mode(round(ip(isfinite(ip)))); end
    end
    if hasLick,   lick = T.Lick_TTL(rs);     lick = lick(good); else, lick=false(size(tr)); end
    if hasReward, rew  = T.Injector_TTL(rs); rew  = rew(good);  else, rew=false(size(tr)); end
    if hasPupil,  pd   = double(T.Diameter_px(rs)); pd = pd(good); else, pd = nan(size(tr)); end
    trials = sort(unique(tr));
    for k = 1:numel(trials)
        mask = (tr==trials(k)); tbk = tb_s(mask); lk=lick(mask); rw=rew(mask); pdK = pd(mask);
        [ln,~,~,lIEI] = eventMetrics(tbk, lk); [rn,~,rDur,~] = eventMetrics(tbk, rw);
        trialDur = finiteRange(tbk);
        out(end+1,:) = {mk, di, g.day_name(i), ss, trials(k), passFlag, trialDur, ln, lIEI, rn, rDur, mean(pdK,'omitnan')}; %#ok<AGROW>
    end
end
if isempty(out)
    TT = table('Size',[0 12], ...
        'VariableTypes',{'string','double','string','double','double','double','double','double','double','double','double','double'}, ...
        'VariableNames',{'mouse_key','day_index','day_name','session_idx','Trial','isPassive','trial_duration_s','lick_count','lick_medianIEI_s','rew_count','rew_dur_s','pupil_mean'});
    return;
end
TT = cell2table(out, 'VariableNames',{'mouse_key','day_index','day_name','session_idx','Trial','isPassive','trial_duration_s','lick_count','lick_medianIEI_s','rew_count','rew_dur_s','pupil_mean'});
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

function plotPerTrialSimple(TT, outPath, nameTag)
M = groupsummary(TT, 'Trial', 'median', {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
E = groupsummary(TT, 'Trial', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
f = figure('Color','w'); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
plotErr(nexttile, M.Trial, M.median_trial_duration_s,  E.fun1_trial_duration_s,  [nameTag ': Trial duration (s)']);
plotErr(nexttile, M.Trial, M.median_lick_count,        E.fun1_lick_count,        'Licks per trial');
plotErr(nexttile, M.Trial, M.median_rew_count,         E.fun1_rew_count,         'Rewards per trial');
plotErr(nexttile, M.Trial, M.median_rew_dur_s,         E.fun1_rew_dur_s,         'Reward duration per trial (s)');
plotErr(nexttile, M.Trial, M.median_pupil_mean,        E.fun1_pupil_mean,        'Pupil mean per trial (px)');
title(tl,['Per-trial (median ± SEM) — ' nameTag]); savepng(f, outPath); close(f);
end

function plotErr(ax, x, y, e, yl)
axes(ax); hold on
if isempty(x) || isempty(y), return; end
x = double(x); y = double(y); e = double(e); e(~isfinite(e)) = 0;
fill([x; flipud(x)], [y-e; flipud(y+e)], [0.9 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.45);
plot(x, y, '-o', 'LineWidth',1.2, 'MarkerSize',4, 'Color',[0.2 0.2 0.2]);
grid on; box on; xlabel('Trial'); ylabel(yl); applyRobustYLim(gca);
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
sub = D(D.isPassive==0,:); if ~isempty(sub)
    f = figure('Color','w'); tl = tiledlayout(ceil(numel(keep)/2),2,'TileSpacing','compact','Padding','compact');
    for i=1:numel(keep)
        nexttile; hold on
        ycol = keep{i};
        ds = groupsummary(sub, 'day_index', {'mean'}, ycol);
        es = groupsummary(sub, 'day_index', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ycol);
        x = ds.day_index; y = ds.("mean_"+ycol); e = es.("fun1_"+ycol);
        shadeRibbon(x,y,e); xlabel('Day'); ylabel(labels{i}); grid on; box on; applyRobustYLim(gca);
    end
    title(tl,['ACTIVE: group means ± SEM across days ' tag]); savepng(f, fullfile(figDir,['group_active_by_day' tag '.png'])); close(f);
end

% PASSIVE
sub = D(D.isPassive==1,:); if ~isempty(sub)
    f = figure('Color','w'); tl = tiledlayout(ceil(numel(keep)/2),2,'TileSpacing','compact','Padding','compact');
    for i=1:numel(keep)
        nexttile; hold on
        ycol = keep{i};
        ds = groupsummary(sub, 'day_index', {'mean'}, ycol);
        es = groupsummary(sub, 'day_index', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ycol);
        x = ds.day_index; y = ds.("mean_"+ycol); e = es.("fun1_"+ycol);
        shadeRibbon(x,y,e); xlabel('Day'); ylabel(labels{i}); grid on; box on; applyRobustYLim(gca);
    end
    title(tl,['PASSIVE: group means ± SEM across days ' tag]); savepng(f, fullfile(figDir,['group_passive_by_day' tag '.png'])); close(f);
end
end

function shadeRibbon(x,y,e)
x = double(x); y = double(y); e = double(e); e(~isfinite(e))=0;
fill([x; flipud(x)], [y-e; flipud(y+e)], [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.5);
plot(x,y,'-o','LineWidth',1.6,'MarkerSize',5,'Color',[0.2 0.2 0.2]);
applyRobustYLim(gca);
end

function applyRobustYLim(ax)
% Robust 2–98% clamp with padding; safety for NaN/singletons
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
try
    ylim(ax, [lo - pad, hi + pad]);
catch
end
end

function plotPupilCorrelations(D, outPath, tag)
f = figure('Color','w'); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% Pupil vs Lick
nexttile; hold on
if hasVar(D,'lick_freq_per_min')
    good = isfinite(D.pupil_mean) & isfinite(D.lick_freq_per_min);
    scatter(D.pupil_mean(good), D.lick_freq_per_min(good), 12, 'filled','MarkerFaceAlpha',0.6);
    if any(good)
        R = corrcoef(D.pupil_mean(good), D.lick_freq_per_min(good)); r1 = R(1,2);
        ttl1 = sprintf('%s: Pupil vs Lick (r=%.2f)', tag, r1);
    else, ttl1 = sprintf('%s: Pupil vs Lick', tag); end
    xlabel('Mean pupil (px)'); ylabel('Lick freq (/min)'); title(ttl1); grid on; box on; applyRobustYLim(gca);
else
    axis off
end
% Pupil vs Reward
nexttile; hold on
if hasVar(D,'rew_freq_per_min')
    good = isfinite(D.pupil_mean) & isfinite(D.rew_freq_per_min);
    scatter(D.pupil_mean(good), D.rew_freq_per_min(good), 12, 'filled','MarkerFaceAlpha',0.6);
    if any(good)
        R = corrcoef(D.pupil_mean(good), D.rew_freq_per_min(good)); r2 = R(1,2);
        ttl2 = sprintf('%s: Pupil vs Reward (r=%.2f)', tag, r2);
    else, ttl2 = sprintf('%s: Pupil vs Reward', tag); end
    xlabel('Mean pupil (px)'); ylabel('Reward freq (/min)'); title(ttl2); grid on; box on; applyRobustYLim(gca);
else
    axis off
end
savepng(f, outPath); close(f);
end

function periEventPlot(T, tb, ttlName, win, baseW, dt, figDir, outName)
if ~hasVar(T,ttlName) || ~hasVar(T,'Diameter_px'), return; end
f = figure('Color','w'); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gp = 0:1
    if hasVar(T,'isPassive'), rowsG = (double(T.isPassive)==gp); else, rowsG = true(height(T),1); end
    sk = unique(T(rowsG,{'mouse_key','day_index','session_idx'}),'rows','stable');
    allM = []; tAxis = [];
    for i=1:height(sk)
        mk=sk.mouse_key(i); di=sk.day_index(i); ss=sk.session_idx(i);
        rr = rowsG & T.mouse_key==mk & T.day_index==di & T.session_idx==ss; if ~any(rr), continue; end
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
    for k=2:numel(on), if on(k)-on(k-1)<refractory_s, keep(k)=false; end, end
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

%% -------- bring Immersion/TST/HOT into day table --------
function D = addDayScalarFromT(D, T, col)
if ~hasVar(T,col), return; end
G = groupsummary(T, {'mouse_key','day_index'}, 'median', col);
D = outerJoinInto(D, G, {'mouse_key','day_index'}, ['median_' col], col);
end

function D = addDayPatternFromT(D, T, prefix)
namesT = T.Properties.VariableNames;
hit = namesT(contains(namesT,prefix));
for i=1:numel(hit)
    G = groupsummary(T, {'mouse_key','day_index'}, 'median', hit{i});
    D = outerJoinInto(D, G, {'mouse_key','day_index'}, ['median_' hit{i}], hit{i});
end
end

function D = outerJoinInto(D, G, keys, srcName, dstName)
    if isempty(G), return; end
    if ~ismember(dstName, D.Properties.VariableNames)
        D.(dstName) = nan(height(D),1);
    end
    keyD = buildKey(D, keys);
    keyG = buildKey(G, keys);
    [tf, loc] = ismember(keyD, keyG);
    vals = G.(srcName);
    D.(dstName)(tf) = vals(loc(tf));
end

function plotImmTstHot_ByWeekGroup(D, cols, figDirWeeks, scopes, immName, tstName, hotName)
for sc = 1:numel(scopes)
    wk  = scopes{sc}.week;
    ps  = scopes{sc}.pass;
    tag = upper(strrep(scopes{sc}.name,'_',' '));

    subD = D(D.week==wk & D.isPassive==ps, :);
    if isempty(subD), continue; end

    f = figure('Color','w'); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    ax = nexttile; hold(ax,'on');
    if ~isempty(immName) && ismember(immName, subD.Properties.VariableNames)
        plotAllMiceLines(subD, immName, 'Immersion latency (s)', cols, '-');
        title(ax, [tag ' — Immersion']);
    else
        axis(ax,'off'); text(0.5,0.5,'No Immersion data','Units','normalized','HorizontalAlignment','center');
    end

    ax = nexttile; hold(ax,'on');
    if ~isempty(tstName) && ismember(tstName, subD.Properties.VariableNames)
        plotAllMiceLines(subD, tstName, ['TST % ' stripPrefix(tstName)], cols, '-');
        title(ax, [tag ' — TST']);
    else
        axis(ax,'off'); text(0.5,0.5,'No TST % data','Units','normalized','HorizontalAlignment','center');
    end

    ax = nexttile; hold(ax,'on');
    if ~isempty(hotName) && ismember(hotName, subD.Properties.VariableNames)
        plotAllMiceLines(subD, hotName, ['HOT % ' stripPrefix(hotName)], cols, '-');
        title(ax, [tag ' — HOT']);
    else
        axis(ax,'off'); text(0.5,0.5,'No HOT % data','Units','normalized','HorizontalAlignment','center');
    end

    savepng(f, fullfile(figDirWeeks, sprintf('imm_tst_hot_%s.png', lower(strrep(scopes{sc}.name,'_','_')))));
    close(f);
end
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
        if i == 1
            k = v;
        else
            k = k + "|" + v;
        end
    end
end

function plotMetricsByDaySeparate(D, figDir)
cols = lines(numel(unique(D.mouse_key,'stable')));
lickVars = {'lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s'};
lickLbls = {'Lick frequency (/min)','Lick mean duration (s)', ...
            'Lick total duration (s/session)','Lick median IEI (s)'};
hasLick = any(ismember(lickVars, D.Properties.VariableNames));
if hasLick
    f = figure('Color','w'); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for k=1:numel(lickVars)
        if ~ismember(lickVars{k}, D.Properties.VariableNames), nexttile; axis off; continue; end
        ax = nexttile; hold(ax,'on');
        plotAllMiceLines(D, lickVars{k}, lickLbls{k}, cols, '-');
        title(ax, lickLbls{k}); applyRobustYLim(ax);
    end
    title(tl,'LICK metrics across days'); 
    savepng(f, fullfile(figDir,'lick_metrics_by_day.png')); close(f);
end

rewVars = {'rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s'};
rewLbls = {'Reward frequency (/min)','Reward mean duration (s)', ...
           'Reward total duration (s/session)','Reward median IRI (s)'};
hasRew = any(ismember(rewVars, D.Properties.VariableNames));
if hasRew
    f = figure('Color','w'); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for k=1:numel(rewVars)
        if ~ismember(rewVars{k}, D.Properties.VariableNames), nexttile; axis off; continue; end
        ax = nexttile; hold(ax,'on');
        plotAllMiceLines(D, rewVars{k}, rewLbls{k}, cols, '-');
        title(ax, rewLbls{k}); applyRobustYLim(ax);
    end
    title(tl,'REWARD metrics across days');
    savepng(f, fullfile(figDir,'reward_metrics_by_day.png')); close(f);
end
end

function plotPerTrialByGroup(TT, figDir)
if ~ismember('isPassive', TT.Properties.VariableNames) || isempty(TT), return; end
labels = {'ACTIVE','PASSIVE'};
for gp = [0 1]
    sub = TT(TT.isPassive==gp, :);
    if isempty(sub), continue; end
    M = groupsummary(sub, 'Trial', 'median', ...
        {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
    E = groupsummary(sub, 'Trial', @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), ...
        {'trial_duration_s','lick_count','rew_count','rew_dur_s','pupil_mean'});
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

function periEventPlot_All(T, tb, ttlName, win, baseW, dt, figDir, outName)
if ~ismember('Diameter_px', T.Properties.VariableNames) || ~ismember(ttlName, T.Properties.VariableNames), return; end
f = figure('Color','w'); ax = axes(f); hold(ax,'on');

sk = unique(T(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
allM = []; tAxis = [];
for i=1:height(sk)
    mk=sk.mouse_key(i); di=sk.day_index(i); ss=sk.session_idx(i);
    rr = T.mouse_key==mk & T.day_index==di & T.session_idx==ss;
    if ~any(rr), continue; end
    t  = tb(rr);
    y  = double(T.Diameter_px(rr));
    ttl= logical(T.(ttlName)(rr));
    if ~any(isfinite(y)), continue; end
    onTimes = ttlOnsets(t, ttl, 0.5);
    if isempty(onTimes), continue; end
    [tAx, ~, ~, ~, matSess] = periEventAverage(t, y, onTimes, win, baseW, dt);
    if isempty(tAx), continue; end
    tAxis = tAx; allM = [allM; matSess.']; %#ok<AGROW>
end

if isempty(allM)
    axis(ax,'off'); text(0.5,0.5,'No events','Units','normalized','HorizontalAlignment','center');
else
    mu  = mean(allM,1,'omitnan');
    se  = std(allM,0,1,'omitnan') ./ sqrt(sum(isfinite(allM),1));
    fill([tAxis; flipud(tAxis)], [mu'-se'; flipud(mu'+se')], [0.9 0.9 1.0], 'EdgeColor','none','FaceAlpha',0.45);
    plot(tAxis, mu, 'k-', 'LineWidth',1.8);
    yline(0,'-','Color',[0.6 0.6 0.6]);
    xlabel('Time from event (s)'); ylabel('\Delta pupil (px)');
    title(['Peri-event pupil \Delta  [' ttlPretty(ttlName) ' — ALL mice]']);
    grid on; box on; applyRobustYLim(ax);
end
savepng(f, fullfile(figDir, outName)); close(f);
end

%% ===================== STATS / BIN TESTS / INDEX HELPERS =====================
function runStatsComparisons_and_Correlations(D, metrics, labels, outDir)
% (1) Any-Passive (mouse-level) vs Active-only stats & plots (+stars)
% (2) Day-bin stats (per-mouse medians)
% (3) Pearson correlation heatmaps

if ~ismember('mouse_key', D.Properties.VariableNames) || ~ismember('isPassive', D.Properties.VariableNames)
    warning('mouse_key/isPassive missing; skipping stats.'); return;
end

% mouse-level group assignment: had any passive day?
Gmouse = groupsummary(D, 'mouse_key', @(x) any(x==1), 'isPassive');
Gmouse.Properties.VariableNames(end) = {'hasPassive'}; % logical
keyD = buildKey(D, {'mouse_key'}); keyG = buildKey(Gmouse, {'mouse_key'});
[tf,loc] = ismember(keyD, keyG);
D.mouse_hasPassive = nan(height(D),1);
D.mouse_hasPassive(tf) = double(Gmouse.hasPassive(loc(tf)));

D.bin5 = dayBin5(D.day_index);
keep = metrics(ismember(metrics, D.Properties.VariableNames));
labels = labels(ismember(metrics, keep));
metrics = keep;

% CSV collector
csvRows = {};

% (1) ALL-days per-mouse medians
for i = 1:numel(metrics)
    y = metrics{i};
    M = groupsummary(D(:,{'mouse_key','mouse_hasPassive',y}), {'mouse_key','mouse_hasPassive'}, 'median', y);
    if isempty(M), continue; end
    M.Properties.VariableNames(end) = {'y_med'};
    Xa = M.y_med(M.mouse_hasPassive==0);
    Xp = M.y_med(M.mouse_hasPassive==1);
    [p,statsOut] = safeRanksum(Xa, Xp);

    fig = figure('Color','w','Position',[100 100 560 420]); hold on
    plotTwoGroupsDots(Xa, Xp, labels{i}, 'ALL days (per-mouse median)');
    ylims = ylim; ybar = ylims(2) - 0.05*range(ylims);
    drawSigBar(gca, 1, 2, ybar, p);
    title(sprintf('%s — Any Passive vs None (ALL days)', labels{i}));
    savepng(fig, fullfile(outDir, sprintf('stat_ANYvsNONE_ALL_%s.png', y))); close(fig);

    csvRows(end+1,:) = {y,'ALL',numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),p,statsOut.delta}; %#ok<AGROW>
end

% (2) per-bin per-mouse medians (two-group, with stars)
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};
for i = 1:numel(metrics)
    y = metrics{i};
    fig = figure('Color','w','Position',[80 80 1160 420]); tl = tiledlayout(1,numel(binNames),'TileSpacing','compact','Padding','compact');
    for b = 1:numel(binNames)
        nexttile; hold on
        thisBin = binNames{b};
        sub = D(strcmp(string(D.bin5), thisBin), {'mouse_key','mouse_hasPassive',y});
        if isempty(sub), axis off; title(thisBin); continue; end
        Mb = groupsummary(sub, {'mouse_key','mouse_hasPassive'}, 'median', y);
        if isempty(Mb), axis off; title(thisBin); continue; end
        Mb.Properties.VariableNames(end) = {'y_med'};
        Xa = Mb.y_med(Mb.mouse_hasPassive==0);
        Xp = Mb.y_med(Mb.mouse_hasPassive==1);

        [p,~] = safeRanksum(Xa, Xp);
        plotTwoGroupsDots(Xa, Xp, labels{i}, thisBin);
        ylims = ylim; ybar = ylims(2) - 0.05*range(ylims);
        drawSigBar(gca, 1, 2, ybar, p);

        csvRows(end+1,:) = {y,thisBin,numel(Xa),numel(Xp),nanmedian(Xa),nanmedian(Xp),p,NaN}; %#ok<AGROW>
    end
    title(tl, sprintf('%s — Any Passive vs None (per bin; per-mouse medians)', labels{i}));
    savepng(fig, fullfile(outDir, sprintf('stat_BINS_%s.png', y))); close(fig);
end

% write CSV
if ~isempty(csvRows)
    Tcsv = cell2table(csvRows, 'VariableNames', ...
        {'metric','level','n_activeOnly','n_hadPassive','median_activeOnly','median_hadPassive','p_ranksum','cliffs_delta'});
    outCSV = fullfile(outDir, 'stats_pvalues_ANYvsNONE_and_BINS.csv');
    try, writetable(Tcsv, outCSV); catch, warning('Could not write %s', outCSV); end
end

% (3) Pearson correlation matrices
corrVarsAll = {'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
               'lick_totalDur_s','rew_totalDur_s','pupil_mean'};
corrVarsAll = corrVarsAll(ismember(corrVarsAll, D.Properties.VariableNames));
if numel(corrVarsAll) >= 2
    plotCorrMatrix(D, corrVarsAll, fullfile(outDir,'corr_matrix_ALL.png'), 'ALL rows');
    plotCorrMatrix(D(D.mouse_hasPassive==0,:), corrVarsAll, fullfile(outDir,'corr_matrix_ACTIVEONLY.png'), 'Active-only');
    plotCorrMatrix(D(D.mouse_hasPassive==1,:), corrVarsAll, fullfile(outDir,'corr_matrix_HADPASSIVE.png'), 'Had-passive');
    binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};
    for b = 1:numel(binNames)
        sub = D(strcmp(string(D.bin5), binNames{b}), :);
        if ~isempty(sub)
            plotCorrMatrix(sub, corrVarsAll, fullfile(outDir, sprintf('corr_matrix_%s.png', binNames{b})), ['Bin ' binNames{b}]);
        end
    end
end
end

function withinGroup_BinANOVA_and_Pairs(D, groupFlag, metrics, labels, outDir)
% Within a group of mice (Active-only=0, Had-passive=1), test across bins.
% Uses per-mouse medians in each bin. KW omnibus + FDR pairwise heatmap.
if ~ismember('isPassive', D.Properties.VariableNames), return; end

% label strings
gname = ternary(groupFlag==0, 'ACTIVEONLY', 'HADPASSIVE');
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};

% Determine group membership at the mouse level
Gmouse = groupsummary(D, 'mouse_key', @(x) any(x==1), 'isPassive');
Gmouse.Properties.VariableNames(end) = {'hasPassive'};
if groupFlag==0
    miceKeep = string(Gmouse.mouse_key(Gmouse.hasPassive==0));
else
    miceKeep = string(Gmouse.mouse_key(Gmouse.hasPassive==1));
end

D.bin5 = dayBin5(D.day_index);

for i=1:numel(metrics)
    y = metrics{i}; if ~ismember(y, D.Properties.VariableNames), continue; end
    lbl = labels{i};

    % per-mouse median per bin
    sub = D(ismember(D.mouse_key, miceKeep) & ~isundefined(D.bin5), {'mouse_key','bin5',y});
    M = groupsummary(sub, {'mouse_key','bin5'}, 'median', y);
    if isempty(M), continue; end
    M.Properties.VariableNames(end) = {'y_med'};

    % Assemble data by bins and clean NaNs
    groups = cell(1,numel(binNames));
    for b=1:numel(binNames)
        xb = M.y_med(string(M.bin5)==binNames{b});
        groups{b} = xb(isfinite(xb));
    end
    % Omnibus KW (fallback safe)
    p_kw = kruskal_p_fallback(groups);

    % Pairwise tests (paired if same mouse in both bins) — NaN-robust
    pairs   = nchoosek(1:numel(binNames),2);
    p_pair  = nan(size(pairs,1),1);
    pair_lbl= strings(size(pairs,1),1);
    n_used  = zeros(size(pairs,1),1);

    for k=1:size(pairs,1)
        b1 = pairs(k,1); b2 = pairs(k,2);
        pair_lbl(k) = sprintf('%s vs %s', binNames{b1}, binNames{b2});

        Mb1 = M(string(M.bin5)==binNames{b1}, {'mouse_key','y_med'});
        Mb2 = M(string(M.bin5)==binNames{b2}, {'mouse_key','y_med'});
        J = innerjoin(Mb1,Mb2,'Keys','mouse_key', ...
                      'LeftVariables',{'mouse_key','y_med'}, ...
                      'RightVariables','y_med');

        if ~isempty(J)
            x1 = J.y_med_Mb1;  x2 = J.y_med_Mb2;
            good = isfinite(x1) & isfinite(x2);
            x1 = x1(good);  x2 = x2(good);
            n_used(k) = numel(x1);

            if n_used(k) >= 3
                if all(x1 == x2)
                    p_pair(k) = 1;  % no change
                else
                    try
                        p_pair(k) = signrank(x1, x2);
                    catch
                        p_pair(k) = safeRanksum(x1, x2); % last resort
                    end
                end
                continue
            end
        end

        % Not enough paired data → unpaired on per-bin vectors
        xa = groups{b1}; xb = groups{b2};
        p_pair(k) = safeRanksum(xa, xb);
        n_used(k) = min(numel(xa), numel(xb));
    end

    % FDR correct
    [p_adj, sig_mask] = fdr_bh_local(p_pair, 0.05);

    % Visual: compact heatmap of -log10(p_adj) with stars
    H = nan(numel(binNames)); 
    for k=1:size(pairs,1)
        i1 = pairs(k,1); i2 = pairs(k,2);
        val = -log10(max(p_adj(k), eps));
        H(i1,i2) = val; H(i2,i1) = val;
    end
    fig = figure('Color','w','Position',[120 120 520 460]);
    imagesc(H); axis square; colorbar; colormap(parula);
    title(sprintf('%s — %s (KW p=%.3g)', gname, lbl, p_kw));
    set(gca,'XTick',1:numel(binNames),'XTickLabel',binNames,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(binNames),'YTickLabel',binNames);
    for r=1:numel(binNames)
        for c=1:numel(binNames)
            if r==c, text(c,r,'—','HorizontalAlignment','center'); continue; end
            k = find( (pairs(:,1)==min(r,c)) & (pairs(:,2)==max(r,c)), 1 );
            if isempty(k), continue; end
            star = sigStarFromP(p_adj(k));
            if ~isempty(star), text(c,r,star,'HorizontalAlignment','center','FontWeight','bold','Color','k'); end
        end
    end
    savepng(fig, fullfile(outDir, sprintf('withinBin_%s_%s.png', gname, y))); close(fig);

    % CSV dump
    C = table(pair_lbl, p_pair, p_adj, sig_mask, n_used, 'VariableNames', ...
        {'pair','p_raw','p_fdr','is_sig','n'});
    try, writetable(C, fullfile(outDir, sprintf('withinBin_%s_%s_pairs.csv', gname, y))); catch, end
end
end





function computeAndPlot_AP_Index(D, metrics, labels, outDir)
% For each day d with both Active (isPassive=0) and Passive (isPassive=1) mice,
% compute I(d) = (A - P) / (A + P), where A/P are medians across mice for that day.
% Plot over days with Spearman trend; also run bin-wise KW + FDR pairwise tests.

if ~ismember('isPassive', D.Properties.VariableNames), return; end
D.bin5 = dayBin5(D.day_index);
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};

for i=1:numel(metrics)
    y = metrics{i}; if ~ismember(y, D.Properties.VariableNames), continue; end
    lbl = labels{i};

    % ------ build daily AP index ------
    days = sort(unique(D.day_index));
    I  = nan(size(days));
    for k=1:numel(days)
        di = days(k);
        % NOTE: use brace indexing {} to extract numeric vectors
        A = D{D.day_index==di & D.isPassive==0, y};
        P = D{D.day_index==di & D.isPassive==1, y};
        A = A(isfinite(A)); P = P(isfinite(P));
        if ~isempty(A) && ~isempty(P)
            a = median(A,'omitnan'); p = median(P,'omitnan');
            if isfinite(a) && isfinite(p) && (a+p)~=0
                I(k) = (a - p) / (a + p);
            end
        end
    end

    good = isfinite(I);
    if nnz(good) < 3, continue; end

    % ------ trend across days ------
    dvec = double(days(good));
    Ivec = I(good);
    try
        [rho,p_sp] = corr(dvec(:), Ivec(:), 'Type','Spearman','Rows','complete');
    catch
        rho = NaN; p_sp = NaN;
    end

    fig = figure('Color','w','Position',[120 120 700 420]); hold on
    plot(days(good), I(good), '-o','LineWidth',1.6,'MarkerSize',5,'Color',[0.2 0.2 0.2]);
    yline(0,'-','Color',[0.6 0.6 0.6]);
    grid on; box on; xlabel('Day');
    ylabel(sprintf('Index I = (A-P)/(A+P) — %s', lbl));
    title(sprintf('Active–Passive Index across days — %s (\\rho=%.2f, p=%.3g)', lbl, rho, p_sp));
    yl = ylim; ybar = yl(2) - 0.07*range(yl);
    if isfinite(p_sp) && p_sp<0.05
        text(mean(xlim), ybar, sigStarFromP(p_sp), 'HorizontalAlignment','center','FontWeight','bold');
    end
    savepng(fig, fullfile(outDir, sprintf('AP_index_by_day_%s.png', y)));
    close(fig);

    % ------ bin-wise KW + pairwise FDR on daily index ------
    di_used = days(good);
    binLab  = arrayfun(@(d) string(dayToBin(d)), di_used);
groups  = cell(1,numel(binNames));
for b = 1:numel(binNames)
    groups{b} = Ivec(binLab==binNames{b});  % store numeric vector in cell
end
    p_kw = kruskal_p_fallback(groups);

    pairs   = nchoosek(1:numel(binNames),2);
    p_pair  = nan(size(pairs,1),1);
    for k=1:size(pairs,1)
        b1 = pairs(k,1); b2 = pairs(k,2);
        p_pair(k) = safeRanksum(groups{b1}, groups{b2}); % unpaired on daily I
    end
    [p_adj, ~] = fdr_bh_local(p_pair, 0.05);

    % heatmap of -log10(FDR p)
    H = nan(numel(binNames));
    for k=1:size(pairs,1)
        i1 = pairs(k,1); i2 = pairs(k,2);
        val = -log10(max(p_adj(k), eps));
        H(i1,i2) = val; H(i2,i1) = val;
    end
    fig = figure('Color','w','Position',[120 120 520 460]);
    imagesc(H); axis square; colorbar; colormap(parula);
    title(sprintf('AP index — bin-wise tests (KW p=%.3g): %s', p_kw, lbl));
    set(gca,'XTick',1:numel(binNames),'XTickLabel',binNames,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(binNames),'YTickLabel',binNames);
    for r=1:numel(binNames)
        for c=1:numel(binNames)
            if r==c, text(c,r,'—','HorizontalAlignment','center'); continue; end
            i1 = min(r,c); i2 = max(r,c);
            k  = find(pairs(:,1)==i1 & pairs(:,2)==i2,1);
            if isempty(k), continue; end
            st = sigStarFromP(p_adj(k));
            if ~isempty(st), text(c,r,st,'HorizontalAlignment','center','FontWeight','bold','Color','k'); end
        end
    end
    savepng(fig, fullfile(outDir, sprintf('AP_index_bins_heatmap_%s.png', y)));
    close(fig);

    % CSV of pairwise AP-index p-values
    pair_lbl = strings(size(pairs,1),1);
    for k=1:size(pairs,1), pair_lbl(k) = sprintf('%s vs %s', binNames{pairs(k,1)}, binNames{pairs(k,2)}); end
    Tcsv = table(repmat(string(lbl),numel(p_adj),1), pair_lbl, p_pair, p_adj, ...
                 'VariableNames', {'metric','pair','p_raw','p_fdr'});
    try, writetable(Tcsv, fullfile(outDir, sprintf('AP_index_pairs_%s.csv', y))); catch, end
end
end






function drawSigBar(ax, x1, x2, yTop, p)
% Draws a significance bar with stars between x1 and x2 at yTop (in data units)
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
% Two-sample test + Cliff's delta (unpaired). Returns p; out.delta if requested.
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
% Benjamini–Hochberg FDR (two-sided). Returns adjusted p and significance mask.
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
% adjusted p-values (Benjamini-Hochberg)
adj = nan(size(p));
adj(idx) = min( cummin( (n./(1:n)') .* ps ), 1);
p_fdr = adj(:);
end

function p = kruskal_p_fallback(groups)
% groups: cell array of vectors
try
    % use stats toolbox if available
    allx = []; grp = [];
    for i=1:numel(groups)
        xi = groups{i}; xi = xi(isfinite(xi));
        allx = [allx; xi(:)];
        grp  = [grp; repmat(i,numel(xi),1)];
    end
    if isempty(allx) || numel(unique(grp))<2, p=NaN; return; end
    p = kruskalwallis(allx, grp, 'off');
catch
    % manual KW on ranks (approximate)
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

function c = dayBin5(day_index)
di = double(day_index);
lab = strings(size(di));
lab(di>=3  & di<=5 )  = "D3-5";
lab(di>=6  & di<=8 )  = "D6-8";
lab(di>=9  & di<=11)  = "D9-11";
lab(di>=12 & di<=14)  = "D12-14";
lab(di>=15 & di<=16)  = "D15-16";
c = categorical(lab, ["D3-5","D6-8","D9-11","D12-14","D15-16"]);
end

function b = dayToBin(di)
% returns categorical label for a single day index
if di>=3 && di<=5, b="D3-5";
elseif di>=6 && di<=8, b="D6-8";
elseif di>=9 && di<=11, b="D9-11";
elseif di>=12 && di<=14, b="D12-14";
elseif di>=15 && di<=16, b="D15-16";
else, b="<undef>";
end
end

function plotTwoGroupsDots(Xa, Xp, ylab, xlabRight)
% Dot swarm + mean±SEM for two groups: 1=Active-only, 2=Had-passive
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
set(gca,'XTick',[1 2],'XTickLabel',{'Active-only','Had-passive'});
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
% compute two-tailed p from r using Student's t with df=n-2
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
        text(j,i, sprintf('%.2f%s', R(i,j), sigStarFromP(p)), 'HorizontalAlignment','center','Color','k','FontSize',9,'FontWeight','bold');
    end
end
savepng(fig, outPath); close(fig);
end

function F = tcdf_local(t, nu)
t = abs(t);
x = nu ./ (nu + t.^2);
F = 1 - 0.5*betainc(x, nu/2, 0.5);
end
