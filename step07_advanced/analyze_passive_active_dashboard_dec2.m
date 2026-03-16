function analyze_passive_active_dashboard_dec2()  
% analyze_passive_active_dashboard_dec2
% One-stop cohort visualization + stats for addiction_concate longitudinal outputs
%
% REVISED (FULL, FUNCTIONAL):
%   - Transition days fixed: [3 6 11 14]
%   - Generates BOTH RAW and Δ-from-Pre versions for RM-ANOVA / posthoc plots
%   - Removes ALL grid lines (grid OFF everywhere)
%   - Forces SAME Y-axis for paired panels (Passive vs Active) where applicable
%   - Removes duplicated posthoc comparison bars (keeps unique A–B once)
%   - Adds NEW metrics (per session -> per day -> all plots/stats):
%       * lick_freq_per_10s
%       * lick duration (mean/total already) + bout metrics:
%           bout_n, bout_meanDur_s, bout_totalDur_s
%           bout_freq_per_10s, bout_freq_per_min
%       * PR/Requirement:
%           RequirementLast (existing)
%           Requirement_cum (cumulative across days per mouse)
%           Requirement_speed_per_day (ΔRequirementLast / Δday)
%           Requirement_speed_per_min (RequirementLast / SessionMinutes)
%
% FIXED ERROR:
%   - table vertcat mismatch in allStats: now we harmonize columns before concatenation.

%% ---------- locate latest run + load ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*')); assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix]   = max([d.datenum]);
runDir   = fullfile(d(ix).folder,d(ix).name);
csvPath  = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% Normalize types
T = ensureString(T,'mouse_key');
T = ensureString(T,'day_name');
T = ensureString(T,'Session_Paradigm');

if ~ismember('isPassive',T.Properties.VariableNames)
    T.isPassive = nan(height(T),1);
elseif ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% TTL numerics (keep raw; threshold later)
if ismember('Lick_TTL',T.Properties.VariableNames),     T.Lick_TTL(isnan(T.Lick_TTL)) = 0; end
if ismember('Injector_TTL',T.Properties.VariableNames), T.Injector_TTL(isnan(T.Injector_TTL)) = 0; end

%% ---------- optional legacy exclusions (keep as-is) ----------
% exclude only 6872 black + 8606 forange
if ismember('mouse_key',T.Properties.VariableNames)
    mk      = string(T.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    drop    = (contains(mk_norm,"6872") & contains(mk_norm,"black")) | ...
              (contains(mk_norm,"8606") & contains(mk_norm,"forange"));
    if any(drop)
        fprintf('Excluding %d rows from: %s\n', nnz(drop), strjoin(unique(cellstr(mk(drop)))', ', '));
        T(drop,:) = [];
    end
end

%% ---------- per-session metrics, then per-day medians ----------
[S, D] = fast_session_day_metrics(T, runDir);

% PUPIL-ONLY exclusion: 7597 black
if ismember('mouse_key', D.Properties.VariableNames) && ismember('pupil_mean', D.Properties.VariableNames)
    mk      = string(D.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    mask7597 = contains(mk_norm,"7597") & contains(mk_norm,"black");
    if any(mask7597)
        fprintf('Pupil-only exclude: setting pupil_mean=NaN for %d rows from 7597 black.\n', nnz(mask7597));
        D.pupil_mean(mask7597) = NaN;
    end
end

%% ---------- APPLY NEW COHORT MAP (overrides group assignment) ----------
cohort = buildNewCohortMap();                      % explicit mapping
[S, D] = applyCohortMapToSD(S, D, cohort);         % adds Cage, Color, Sex, GroupType, PairID

% If any mice are NOT in the map, fall back to previous passive detection
if any(ismissing(string(D.GroupType)))
    fprintf('Warning: %d D-rows have missing GroupType (mouse not in cohort map). Falling back to auto detection.\n', ...
        nnz(ismissing(string(D.GroupType))));
    D.HadPassive = classifyHadPassive(S, D);
    D.GroupType(ismissing(string(D.GroupType))) = string(categorical(D.HadPassive(ismissing(string(D.GroupType))), [false true], {'Active','Passive'}));
end

% Ensure D.Group is correct categorical
D.Group = categorical(string(D.GroupType), {'Active','Passive'});

% Exclude habituation days (1-2) from analysis
if ismember('day_index', D.Properties.VariableNames)
    D = D(double(D.day_index) >= 3, :);
end

%% ---------- PR/Requirement handling for PASSIVE during day 6-10 ----------
% Passive mice get replay/forced delivery during 6-10; PR "RequirementLast" is not comparable.
% Set RequirementLast = NaN for (Passive AND day 6-10)
if ismember('RequirementLast', D.Properties.VariableNames)
    mask = (string(D.Group)=="Passive") & (double(D.day_index) >= 6) & (double(D.day_index) <= 10);
    if any(mask)
        fprintf('Setting RequirementLast=NaN for Passive mice on days 6-10 (%d rows).\n', nnz(mask));
        D.RequirementLast(mask) = NaN;
    end
end

%% ---------- period definitions + transition-day option ----------
% FIXED as you requested:
defaultTransitionDays = [3 6 11 14];

runOneDashboard(D, runDir, cohort, true,  defaultTransitionDays);  % with transition days
runOneDashboard(D, runDir, cohort, false, defaultTransitionDays);  % without transition days

fprintf('Done.\n');
end

%% ========================================================================
% RUN ONE VERSION (with or without transition days)
function runOneDashboard(Din, runDir, cohort, includeTransitionDays, transitionDays) %#ok<INUSD>

D = Din;

tag = 'with_transition_days';
if ~includeTransitionDays
    tag = 'without_transition_days';
    if ~isempty(transitionDays)
        D = D(~ismember(double(D.day_index), double(transitionDays)), :);
    end
end

% Assign periods
D.Period = periodOfDay(double(D.day_index));

% Drop undefined
D = D(~isundefined(D.Period), :);

%% ---------- metric list (dynamic) ----------
baseMetrics = { ...
    'RequirementLast','Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min', ...
    'lick_freq_per_min','lick_freq_per_10s', ...
    'lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
    'bout_n','bout_freq_per_10s','bout_freq_per_min','bout_meanDur_s','bout_totalDur_s', ...
    'rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s','pupil_mean'};

baseLabels  = { ...
    'Requirement (PR)','Requirement cumulative','Requirement speed (/day)','Requirement speed (/min)', ...
    'Lick freq (/min)','Lick freq (/10s)', ...
    'Lick mean dur (s)','Lick total dur (s)','Lick median IEI (s)', ...
    'Bout count','Bout freq (/10s)','Bout freq (/min)','Bout mean dur (s)','Bout total dur (s)', ...
    'Reward freq (/min)','Reward mean dur (s)', ...
    'Reward total dur (s)','Reward median IRI (s)','Pupil mean (px)'};

% Optional extras: Immersion; prefer *Nonmoving* for TST/HOT if present.
extraM = {}; extraL = {};
immCol = colLike(D,'Immersion_Latency_s');
if ~isempty(immCol), extraM{end+1}=immCol; extraL{end+1}='Immersion latency (s)'; end %#ok<AGROW>
tstCols = D.Properties.VariableNames(contains(D.Properties.VariableNames,'TST_Pct'));
hotCols = D.Properties.VariableNames(contains(D.Properties.VariableNames,'HOT_Pct'));
tstPick = pickPctColumn(tstCols, 'Nonmoving'); if ~isempty(tstPick), extraM{end+1}=tstPick; extraL{end+1}='TST nonmoving (%)'; end %#ok<AGROW>
hotPick = pickPctColumn(hotCols, 'Nonmoving'); if ~isempty(hotPick), extraM{end+1}=hotPick; extraL{end+1}='HOT nonmoving (%)'; end %#ok<AGROW>

metrics = [baseMetrics, extraM];
labels  = [baseLabels,  extraL];
keepIdx = ismember(metrics, D.Properties.VariableNames);
metrics = metrics(keepIdx); labels = labels(keepIdx);

% Add alternate column names if present (avoids missing desired metrics)
[metrics, labels] = addMetricIfExists(metrics, labels, D, ...
    {'lick_dur_s','lick_duration_s'}, 'Lick duration (s)');
[metrics, labels] = addMetricIfExists(metrics, labels, D, ...
    {'bout_duration_s','lick_bout_duration_s'}, 'Lick bout duration (s)');
[metrics, labels] = addMetricIfExists(metrics, labels, D, ...
    {'Requirement_speed','Requirement_speed_per_session'}, 'Requirement speed');

%% ---------- output dir ----------
outDir = fullfile(runDir,'figs','dashboard_newcohort', tag);
if ~exist(outDir,'dir'), mkdir(outDir); end

% Spaghetti plots for all metrics shown in the dashboard
for si=1:numel(metrics)
    ycol = metrics{si};
    ylab = labels{si};
    if ismember(ycol, D.Properties.VariableNames)
        plotSpaghettiByDayGroup(D, ycol, ylab, outDir);
    end
end

% Cumulative RequirementLast across days (Active vs Passive)
if ismember('Requirement_cum', D.Properties.VariableNames)
    plotRequirementCumulativeAcrossDays(D, outDir);
end

% Report missing RequirementLast by mouse/day/period for raw data checks
if ismember('RequirementLast', D.Properties.VariableNames)
    writeMissingRequirementReport(D, outDir);
end

%% ---------- main plots + stats ----------
allStats = table();
for k = 1:numel(metrics)
    ycol = metrics{k}; ylab = labels{k};

    W = perMousePeriodTable(D, ycol);     % one row per mouse×period (median across days)
    if isRequirementMetric(ycol)
        W = W(~(W.Group=="Passive" & W.Period=="During"), :);
    end
    if isempty(W), continue; end

    % Dual-panel period lines (FORCE same Y)
    fh = figure('Color','w','Position',[80 60 900 900]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    dualPanelPeriodLines(W, ylab);
    printpng(fh, fullfile(outDir, sprintf('dual_period_%s.png', safeName(ycol))));
    close(fh);

    % Change-from-Pre strip/box + p-values (no grid)
    fh = figure('Color','w','Position',[100 80 800 500]);
    changeFromPreStrip_Period(W, ylab, 'none');  % 'holm'|'fdr'|'none'
    printpng(fh, fullfile(outDir, sprintf('delta_from_pre_%s.png', safeName(ycol))));
    close(fh);

    % RM-ANOVA within group
    % 1) RAW version
    rmOneWayWithinGroupFigure_Period(W, ylab, ycol, outDir, false);
    % 2) DELTA version
    rmOneWayWithinGroupFigure_Period(W, ylab, ycol, outDir, true);

    % Stats tables (NOTE: different sets of columns => harmonize before vertcat)
    ST  = runStats_Period(W, ycol);
    STd = runDeltaStats_Period(W, ycol);

    ST.metric  = repmat(string(ycol), height(ST),1);
    ST.label   = repmat(string(ylab), height(ST),1);
    STd.metric = repmat(string(ycol), height(STd),1);
    STd.label  = repmat(string(ylab), height(STd),1);

    STrm_raw   = runRMAnovaWithinGroup_Period(W, ycol, false);
    STrm_delta = runRMAnovaWithinGroup_Period(W, ycol, true);
    STrm_raw.label   = repmat(string(ylab), height(STrm_raw),1);
    STrm_delta.label = repmat(string(ylab), height(STrm_delta),1);
    STrm_raw.variant   = repmat("RAW",   height(STrm_raw),1);
    STrm_delta.variant = repmat("DELTA", height(STrm_delta),1);

    % FIX: align columns before vertcat
    chunk = alignAndVertcatStats({ST, STd, STrm_raw, STrm_delta});
    allStats = [allStats; chunk]; %#ok<AGROW>
end

% Global multiple-comparison corrections across ALL tests for ALL metrics
if ~isempty(allStats) && ismember('p', allStats.Properties.VariableNames)
    allStats.p_FDR  = fdrBH(allStats.p);
    allStats.p_Holm = holmBonferroni(allStats.p);
    writetable(allStats, fullfile(outDir,'stats_summary.csv'));
    disp(allStats(:, intersect({'metric','label','variant','test','effect','level','p','p_FDR','p_Holm','N'}, allStats.Properties.VariableNames)));
end

% Correlation + PCA examples (per mouse-period, core behavior metrics)
core = {'lick_freq_per_min','rew_freq_per_min','lick_meanDur_s','lick_medianIEI_s','pupil_mean'};
core = core(ismember(core, D.Properties.VariableNames));
if numel(core) >= 3
    Wcore = perMousePeriodTable(D, core);
    corrPCA(Wcore, outDir);
end

fprintf('Dashboard (%s) saved to:\n  %s\n', tag, outDir);
end

%% ========================================================================
% FIX helper: harmonize stats table columns then concatenate
function out = alignAndVertcatStats(tbls)
% tbls: cell array of tables
% Ensures all have same variables (union), fills missing columns with missing,
% and orders columns consistently.

% union of variable names
allVars = strings(0,1);
for i=1:numel(tbls)
    if isempty(tbls{i}), continue; end
    v = string(tbls{i}.Properties.VariableNames);
    v = v(:);
    allVars = unique([allVars; v], 'stable');
end

% preferred front order (if present)
pref = ["metric","label","variant","test","effect","level","p","N"];
front = pref(ismember(pref, allVars));
front = front(:);
rest  = setdiff(allVars, front, 'stable');
rest = rest(:);
ordered = [front; rest];

out = table();
for i=1:numel(tbls)
    Ti = tbls{i};
    if isempty(Ti) || height(Ti)==0
        continue
    end

    % add missing vars
    missingVars = setdiff(ordered, string(Ti.Properties.VariableNames), 'stable');
    for v = missingVars'
        Ti.(v) = missingColOfType(height(Ti));
    end

    % reorder
    Ti = Ti(:, cellstr(ordered));

    out = [out; Ti]; %#ok<AGROW>
end
end

function col = missingColOfType(n)
% Use string missing; numeric ones like p,N will be converted later if needed.
% Safer: return NaN double for p/N columns? But we don’t know column name here.
% Using string missing works fine for text columns; MATLAB will error if later expects numeric.
% So we return NaN double, which is acceptable for p/N and also for string columns (they can store NaN? no).
% Better: return cellstr missing? We'll use "missing string" (""), but p/N need NaN.
% To avoid type issues, we return a cell array of missing, which makes table variable type cell.
col = repmat(missing, n, 1); %#ok<MISSL>
end

function [metrics, labels] = addMetricIfExists(metrics, labels, D, candidates, label)
if any(strcmp(metrics, label))
    return
end
for i = 1:numel(candidates)
    c = candidates{i};
    if ismember(c, D.Properties.VariableNames) && ~ismember(c, metrics)
        metrics{end+1} = c; %#ok<AGROW>
        labels{end+1}  = label; %#ok<AGROW>
        break
    end
end
end

function tf = isRequirementMetric(ycol)
tf = startsWith(string(ycol), "Requirement");
end

function perCounts = periodFiniteCounts(Tg)
if isempty(Tg)
    perCounts = table(categorical([],[],[]), zeros(0,1), 'VariableNames', {'Period','count'});
    return
end
per = categorical(Tg.Period);
cats = categories(per);
cnt = zeros(numel(cats),1);
for i=1:numel(cats)
    mask = per==cats{i};
    cnt(i) = sum(isfinite(Tg.value(mask)));
end
perCounts = table(categorical(cats, categories(per), 'Ordinal', isordinal(per)), cnt, ...
    'VariableNames', {'Period','count'});
end

function writeMissingRequirementReport(D, outDir)
T = D(:, {'mouse_key','Group','day_index','Period','RequirementLast'});
T.Group = string(T.Group);
T.Period = string(T.Period);
missMask = ~isfinite(double(T.RequirementLast));
missing = T(missMask, {'mouse_key','Group','day_index','Period'});
if ~isempty(missing)
    writetable(missing, fullfile(outDir, 'RequirementLast_missing_rows.csv'));
end

% Summary: count missing by mouse/period (compatible with older MATLAB)
if ~isempty(missing)
    [Gid, mKey, mGroup, mPeriod] = findgroups(missing.mouse_key, missing.Group, missing.Period);
    cnt = splitapply(@numel, missing.day_index, Gid);
    G = table(mKey, mGroup, mPeriod, cnt, ...
        'VariableNames', {'mouse_key','Group','Period','missing_days'});
    writetable(G, fullfile(outDir, 'RequirementLast_missing_summary.csv'));
end
end

function plotSpaghettiByDayGroup(D, ycol, ylab, outDir)
if ~ismember(ycol, D.Properties.VariableNames), return; end

T = D(:, {'mouse_key','Group','day_index', ycol});
T.Group = string(T.Group);
T = T(isfinite(T.day_index), :);
if isempty(T), return; end

fig = figure('Color','w','Position',[80 80 950 520]); hold on
cols = struct('Active',[0.9 0.3 0.2], 'Passive',[0.2 0.4 0.9]);

groups = ["Active","Passive"];
for gi=1:numel(groups)
    g = groups(gi);
    Tg = T(T.Group==g & isfinite(T.(ycol)), :);
    if isempty(Tg), continue; end
    mk = unique(Tg.mouse_key, 'stable');
    for i=1:numel(mk)
        r = Tg.mouse_key==mk(i);
        days = Tg.day_index(r);
        vals = double(Tg.(ycol)(r));
        [dU,~,idx] = unique(days);
        vU = accumarray(idx, vals, [], @(x) mean(x,'omitnan'));
        [dU,ord] = sort(dU); vU=vU(ord);
        plot(dU, vU, '-', 'Color', cols.(g), 'LineWidth', 0.8, 'HandleVisibility','off');
    end

    % group mean by day
    dAll = unique(Tg.day_index); dAll = sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(Tg.(ycol)(Tg.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, '-', 'Color', cols.(g), 'LineWidth', 2.5);
end

xlabel('Day'); ylabel(ylab, 'Interpreter','none');
title(sprintf('%s (spaghetti by mouse)', ylab), 'Interpreter','none');
legend(groups, 'Location','best');
grid off; box off
exportgraphics(fig, fullfile(outDir, sprintf('SPAGHETTI_%s.png', safeName(ycol))), 'Resolution', 220);
close(fig);
end

function plotRequirementCumulativeAcrossDays(D, outDir)
if ~ismember('Requirement_cum', D.Properties.VariableNames), return; end
fig = figure('Color','w','Position',[80 80 950 520]); hold on
groups = ["Active","Passive"];
cols = struct('Active',[0.9 0.3 0.2], 'Passive',[0.2 0.4 0.9]);
for gi=1:numel(groups)
    g = groups(gi);
    rG = string(D.Group)==g;
    if ~any(rG), continue; end
    mk = unique(D.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & D.mouse_key==mk(i);
        [d,ord] = sort(D.day_index(r));
        y = double(D.Requirement_cum(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',3,'Color',cols.(g), 'HandleVisibility','off');
    end
    dAll = unique(D.day_index(rG)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(D.Requirement_cum(rG & D.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, '-', 'LineWidth', 2.5, 'Color', cols.(g));
end
xlabel('Day'); ylabel('Cumulative Requirement');
title('Cumulative Requirement across days');
legend(groups,'Location','best');
grid off; box off
exportgraphics(fig, fullfile(outDir,'Requirement_cumulative_across_days.png'), 'Resolution', 220);
close(fig);
end

%% ========================================================================
% COHORT MAP: explicit, prevents "messing up" Active vs Passive
function cohort = buildNewCohortMap()
% Columns: Cage (string), Color (string), Sex (string 'F'/'M'), GroupType ('Active'/'Passive'), PairID (string)

rows = {
    % cage   color    sex  group     pair
    "6100", "red",    "F", "Passive", "6100_pair1"
    "6100", "orange", "F", "Passive", "6100_pair1"
    "6100", "black",  "F", "Active",  "6100_pair1"

    "0911", "red",    "F", "Active",  "0911_pair1"
    "0911", "orange", "F", "Passive", "0911_pair1"

    "0911", "black",  "F", "Passive", "0911_pair2"
    "0911", "white",  "F", "Active",  "0911_pair2"

    "0910", "black",  "M", "Active",  "0910_pair1"
    "0910", "red",    "M", "Passive", "0910_pair1"
    "0910", "orange", "M", "Passive", "0910_pair1"

    "6099", "red",    "M", "Passive", "6099_pair1"
    "6099", "orange", "M", "Active",  "6099_pair1"

    "6099", "black",  "M", "Active",  "6099_pair2"
    "6099", "white",  "M", "Passive", "6099_pair2"   % died day13 (handled naturally by missing later days)
    };

cohort = cell2table(rows, 'VariableNames', {'Cage','Color','Sex','GroupType','PairID'});
end

function [S, D] = applyCohortMapToSD(S, D, cohort)
% Parse mouse_key -> Cage/Color, then join with cohort map.

% --- apply to S ---
[Scage, Scolor] = parseCageColor(string(S.mouse_key));
S.Cage  = Scage;
S.Color = Scolor;
S = addMapCols(S, cohort);

% --- apply to D ---
[Dcage, Dcolor] = parseCageColor(string(D.mouse_key));
D.Cage  = Dcage;
D.Color = Dcolor;
D = addMapCols(D, cohort);
end

function T = addMapCols(T, cohort)
% Adds Sex, GroupType, PairID by matching Cage+Color. Unmatched remain missing.
T.Sex       = strings(height(T),1);
T.GroupType = strings(height(T),1);
T.PairID    = strings(height(T),1);

% normalize for match
cageT  = string(T.Cage);   colorT = lower(string(T.Color));
cageC  = string(cohort.Cage); colorC = lower(string(cohort.Color));

for i = 1:height(cohort)
    m = (cageT==cageC(i)) & (colorT==colorC(i));
    if any(m)
        T.Sex(m)       = string(cohort.Sex(i));
        T.GroupType(m) = string(cohort.GroupType(i));
        T.PairID(m)    = string(cohort.PairID(i));
    end
end

% convert to categoricals where useful
T.Sex       = categorical(T.Sex,       {'F','M'});
T.GroupType = categorical(T.GroupType, {'Active','Passive'});
end

function [cage, color] = parseCageColor(mouse_key)
% Robust parsing from mouse_key string.
% Extract first 4-digit cage number and a color token among {black, red, orange, white}.
mk = lower(regexprep(strtrim(mouse_key),'[_\-]+',' '));

% cage: first 4-digit block
cage = strings(numel(mk),1);
for i=1:numel(mk)
    tok = regexp(mk(i), '(\d{4})', 'tokens', 'once');
    if ~isempty(tok), cage(i) = string(tok{1}); else, cage(i) = ""; end
end

% color: word-boundary match (avoids 'forange' confusion)
colors = ["black","red","orange","white"];
color = strings(numel(mk),1);
for i=1:numel(mk)
    found = "";
    for c = 1:numel(colors)
        pat = "\<" + colors(c) + "\>";
        if ~isempty(regexp(mk(i), pat, 'once'))
            found = colors(c);
            break
        end
    end
    color(i) = found;
end
end

%% ========================================================================
% PERIODS (new)
function P = periodOfDay(d)
p = strings(size(d));
p(d>=3  & d<=5 )  = "Pre";
p(d>=6  & d<=10)  = "During";
p(d>=11 & d<=13)  = "Post";
p(d>=14 & d<=16)  = "Withdrawal";
p(d>=17 & d<=18)  = "Re-exposure";
P = categorical(p, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

%% ========================================================================
% PER-MOUSE × PERIOD table
function W = perMousePeriodTable(D, ycols)
if ischar(ycols) || isstring(ycols), ycols = cellstr(string(ycols)); end
D2 = D;

if numel(ycols)==1
    y = ycols{1};
    G = groupsummary(D2(:,{'mouse_key','Group','Period',y}), {'mouse_key','Group','Period'}, 'median', y);
    W = renamevars(G, "median_"+y, "value");
else
    key = unique(D2(:,{'mouse_key','Group','Period'}),'rows','stable');
    W = key;
    for j=1:numel(ycols)
        y = ycols{j};
        G = groupsummary(D2(:,{'mouse_key','Group','Period',y}), {'mouse_key','Group','Period'}, 'median', y);
        W = outerjoin(W, renamevars(G, "median_"+y, y), 'Keys',{'mouse_key','Group','Period'}, 'MergeKeys',true);
    end
end
end

%% ========================================================================
% PLOTTING for PERIODS
function COL = palette()
COL.passive = [0 0.45 0.74];
COL.active  = [0.85 0.33 0.10];
COL.grey    = 0.7*[1 1 1];
COL.periodOrder = {'Pre','During','Post','Withdrawal','Re-exposure'};
COL.periodMark  = {'o','^','s','d','p'};
end

function dualPanelPeriodLines(W, ylab)
COL = palette(); periods = COL.periodOrder;

ax1 = nexttile; hold(ax1,'on'); title(ax1, sprintf('Passive — %s', ylab));
drawPanelPeriod(ax1, W(W.Group=="Passive",:), periods, COL, ylab);

ax2 = nexttile; hold(ax2,'on'); title(ax2, sprintf('Active — %s', ylab));
drawPanelPeriod(ax2, W(W.Group=="Active",:), periods, COL, ylab);

% FORCE SAME Y LIMITS (Passive vs Active)
yl1 = ylim(ax1); yl2 = ylim(ax2);
yl  = [min(yl1(1),yl2(1)), max(yl1(2),yl2(2))];
ylim(ax1, yl);
ylim(ax2, yl);
end

function drawPanelPeriod(ax, Wg, periods, COL, ylab)
axes(ax);

if isempty(Wg)
    text(0.02,0.95,'N=0 mice','Units','normalized');
    set(ax,'XTick',1:numel(periods),'XTickLabel',periods);
    xlim(ax,[0.7 numel(periods)+0.3]);
    box(ax,'on');
    grid(ax,'off');
    ylabel(ax, ylab);
    return;
end

mice = unique(Wg.mouse_key,'stable');
for i=1:numel(mice)
    s = Wg(Wg.mouse_key==mice(i),:);
    [x,ord] = sort(categorical(s.Period,periods,'Ordinal',true));
    y = s.value(ord);
    plot(ax, double(x), y, '-o', 'Color',COL.grey, 'MarkerSize',3, 'LineWidth',0.9, ...
         'MarkerFaceColor',COL.grey*0.9, 'MarkerEdgeColor','none');
end

M = groupsummary(Wg, 'Period', 'mean', 'value');
E = groupsummary(Wg, 'Period', @(x) std(x,'omitnan')./sqrt(sum(isfinite(x))), 'value');

x = double(categorical(M.Period,periods,'Ordinal',true));
y = M.mean_value;
e = E.fun1_value;

fill(ax, [x; flipud(x)], [y-e; flipud(y+e)], [0 0 0], 'FaceAlpha',0.12, 'EdgeColor','none');
plot(ax, x, y, 'k-', 'LineWidth',2.0);
errorbar(ax, x, y, e, 'k','LineWidth',1.1,'CapSize',8);

set(ax,'XTick',1:numel(periods),'XTickLabel',periods);
xlim(ax,[0.7 numel(periods)+0.3]);
box(ax,'on');
grid(ax,'off');
ylabel(ax, ylab);

text(0.02,0.95,sprintf('N=%d mice', numel(mice)),'Units','normalized','FontWeight','bold');
end

function changeFromPreStrip_Period(W, ylab, adjustMethod)
if nargin<3 || isempty(adjustMethod), adjustMethod = 'holm'; end

COL     = palette();
periods = COL.periodOrder;
W = W(~isnan(W.value),:);

base = groupsummary(W(W.Period=="Pre",:), 'mouse_key', 'median', 'value');
base.Properties.VariableNames{'median_value'} = 'base';
W = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key', 'MergeKeys',true);
W.delta = W.value - W.base;

hold on
groups = {'Passive','Active'};
off    = [-0.12 0.12];
hG(1) = scatter(nan,nan,36,'o','filled','MarkerFaceColor',COL.passive,'MarkerEdgeColor','k','DisplayName','Passive');
hG(2) = scatter(nan,nan,36,'o','filled','MarkerFaceColor',COL.active, 'MarkerEdgeColor','k','DisplayName','Active');

topY = -inf(1,numel(periods));
for g=1:numel(groups)
    G = W(W.Group==groups{g},:);
    thisCol = COL.passive; if g==2, thisCol = COL.active; end
    for e=1:numel(periods)
        Ge = G(G.Period==periods{e} & isfinite(G.delta),:);
        if isempty(Ge), continue; end
        x = e + off(g);
        b = boxchart(repmat(x,height(Ge),1), Ge.delta, 'BoxWidth',0.18, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
        if isprop(b,'BoxFaceColor'), b.BoxFaceColor = thisCol; end
        if isprop(b,'BoxEdgeColor'), b.BoxEdgeColor = thisCol; end
        scatter(repmat(x,height(Ge),1), Ge.delta, 24, 'filled', 'MarkerFaceColor',thisCol, ...
            'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.75);
        topY(e) = max(topY(e), max(Ge.delta));
    end
end

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlim([0.5 numel(periods)+0.5]);
yline(0,'k:');
box on;
grid off;
ylabel(sprintf('%s – Pre', ylab));
legend(hG,'Location','best');
title('Change from Pre (per mouse)');

% Passive vs Active per period: rank-sum (skip Pre)
pRaw  = nan(1,numel(periods)); valid = false(1,numel(periods));
for e = 1:numel(periods)
    if strcmp(periods{e},'Pre'), continue; end
    Pa = W(W.Group=="Passive" & W.Period==periods{e} & isfinite(W.delta),:);
    Ac = W(W.Group=="Active"  & W.Period==periods{e} & isfinite(W.delta),:);
    x = Pa.delta; y = Ac.delta;
    if numel(x)>=2 && numel(y)>=2
        pRaw(e) = ranksum(x,y);
        valid(e)= true;
    end
end

pAdj = pRaw;
switch lower(string(adjustMethod))
    case "holm", if any(valid), pAdj(valid) = holmBonferroni(pRaw(valid)); end
    case "fdr",  if any(valid), pAdj(valid) = fdrBH(pRaw(valid));          end
    otherwise
end

yl = ylim; rngY = diff(yl); pad = 0.06*rngY;
for e = 1:numel(periods)
    if ~valid(e), continue; end
    x1 = e + off(1);  x2 = e + off(2);
    y  = max(topY(e), yl(2)-2*pad) + pad;
    plot([x1 x1 x2 x2],[y-0.5*pad y y y-0.5*pad],'k-','LineWidth',1);
    tag = upper(char(adjustMethod)); if tag=="NONE", tag="raw"; end
    txt = sprintf('%s (%s %.3g)', starStr(pAdj(e)), tag, pAdj(e));
    text(mean([x1 x2]), y+0.2*pad, txt, 'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', 'FontSize',9, 'FontWeight','bold');
    yl(2) = max(yl(2), y+0.8*pad);
end
ylim(yl);
end

%% ========================================================================
% RM-ANOVA plot (RAW or DELTA) with:
%   - no grid
%   - same Y across Passive/Active panels
%   - unique posthoc bars (no duplicated A–B / B–A)
function rmOneWayWithinGroupFigure_Period(W, ylab, ycol, outDir, useDelta)
if nargin<5, useDelta = true; end

Eall = ["Pre","During","Post","Withdrawal","Re-exposure"];
if useDelta
    E = Eall(2:end);
else
    E = Eall;
end
groups = ["Passive","Active"];
safeE  = matlab.lang.makeValidName(cellstr(E));
xMap   = containers.Map(cellstr(E), num2cell(1:numel(E)));

fh = figure('Color','w','Position',[90 90 980 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax = gobjects(1,2);

for gi = 1:numel(groups)
    G = groups(gi);
    ax(gi) = nexttile; hold(ax(gi),'on');

    if useDelta
        ttlBase = sprintf('%s — %s (Δ from Pre)', ylab, G);
    else
        ttlBase = sprintf('%s — %s (RAW)', ylab, G);
    end
    title(ax(gi), ttlBase);

    if useDelta
        base = groupsummary(W(W.Group==G & W.Period=="Pre",:), 'mouse_key', 'median', 'value');
        base.Properties.VariableNames{'median_value'}='base';
        Tg = outerjoin(W(W.Group==G & ismember(W.Period,E),:), base(:,{'mouse_key','base'}), ...
            'Keys','mouse_key','MergeKeys',true);
        Tg.value = Tg.value - Tg.base;
    else
        Tg = W(W.Group==G & ismember(W.Period,E),:);
    end

    % Drop periods with insufficient data for this group
    perCounts = periodFiniteCounts(Tg);
    keepPer = perCounts.Period(perCounts.count >= 3);
    Euse = E(ismember(E, string(keepPer)));
    if numel(Euse) < 3
        msg = sprintf('N<3 complete mice (RM-ANOVA): N=%d', 0);
        text(ax(gi), 0.02,0.95, msg, 'Units','normalized','FontWeight','bold');
        continue
    end
    safeE  = matlab.lang.makeValidName(cellstr(Euse));
    xMap   = containers.Map(cellstr(Euse), num2cell(1:numel(Euse)));

    % boxcharts
    for ei = 1:numel(Euse)
        Ge = Tg(Tg.Period==Euse(ei) & isfinite(Tg.value),:);
        if isempty(Ge), continue; end
        b = boxchart(ax(gi), repmat(ei,height(Ge),1), Ge.value, 'BoxWidth',0.25, 'MarkerStyle','none','BoxFaceAlpha',0.22);
        if G=="Passive", b.BoxFaceColor=[0 0.45 0.74]; else, b.BoxFaceColor=[0.85 0.33 0.10]; end
        scatter(ax(gi), repmat(ei,height(Ge),1), Ge.value, 26, 'k', 'filled', 'MarkerFaceAlpha',0.75);
    end
    set(ax(gi),'XTick',1:numel(Euse),'XTickLabel',Euse);
    box(ax(gi),'on');
    grid(ax(gi),'off');

    if useDelta
        yline(ax(gi), 0,'k:');
        ylabel(ax(gi), sprintf('%s – Pre', ylab));
    else
        ylabel(ax(gi), ylab);
    end

    % complete cases only
    Wg   = unstack(Tg(:,{'mouse_key','Period','value'}), 'value', 'Period');
    have = safeE(ismember(safeE, Wg.Properties.VariableNames));
    ok   = all(isfinite(Wg{:,have}),2);
    Wc   = Wg(ok, [{'mouse_key'}, have]);

    if height(Wc) < 3
        msg = sprintf('N<3 complete mice (RM-ANOVA): N=%d', height(Wc));
        text(ax(gi), 0.02,0.95, msg, 'Units','normalized','FontWeight','bold');
        continue
    end

    orig = strings(numel(have),1);
    for i = 1:numel(have)
        orig(i) = Euse(strcmp(have{i}, safeE));
    end
    within = table(categorical(orig,'Ordinal',true), 'VariableNames', {'Period'});

    rm  = fitrm(Wc, sprintf('%s ~ 1', strjoin(have,',')), 'WithinDesign', within);
    RA  = ranova(rm,'WithinModel','Period');

    termLabels = ranovaTermLabels(RA);
    idx = find(strcmpi(termLabels,'Period') | contains(lower(termLabels),'period'), 1, 'first');
    pOmni = NaN;
    if ~isempty(idx) && ismember('pValue', RA.Properties.VariableNames)
        pOmni = RA.pValue(idx);
    end

    if isfinite(pOmni)
        title(ax(gi), sprintf('%s (RM-ANOVA p=%.3g)', ttlBase, pOmni), 'FontWeight','bold');
    else
        title(ax(gi), sprintf('%s (RM-ANOVA p=n/a)', ttlBase), 'FontWeight','bold');
    end

    % posthoc
    try
        C = multcompare(rm,'Period','ComparisonType','tukey-kramer');
    catch
        C = multcompare(rm,'Period','ComparisonType','bonferroni');
    end

    vnames = string(C.Properties.VariableNames);
    L1 = 'LevelA'; L2 = 'LevelB';
    if ~ismember(L1,vnames)
        i1 = find(endsWith(vnames,'_1'),1,'first'); i2 = find(endsWith(vnames,'_2'),1,'first');
        if ~isempty(i1), L1 = char(vnames(i1)); end
        if ~isempty(i2), L2 = char(vnames(i2)); end
    end
    if ismember('pValue',vnames), pv = C.pValue;
    elseif ismember('pValueAdj',vnames), pv = C.pValueAdj;
    else, pv = NaN(height(C),1);
    end

    yl = ylim(ax(gi)); baseY = yl(2); pad = 0.07*diff(yl);

    % UNIQUE comparisons only (no duplicates A–B and B–A)
    seen = strings(0,1);
    used = zeros(1,numel(E));  % track stacking per span
    for r = 1:height(C)
        a = string(C.(L1)(r)); b = string(C.(L2)(r));
        if ~isKey(xMap, char(a)) || ~isKey(xMap, char(b)), continue; end
        if ~isfinite(pv(r)) || pv(r) >= 0.05, continue; end

        pair = sort([a b]);
        key  = pair(1) + "_" + pair(2);
        if any(seen==key), continue; end
        seen(end+1) = key; %#ok<AGROW>

        x1 = xMap(char(a)); x2 = xMap(char(b));
        span = abs(x2-x1); if span==0, span=1; end
        used(span) = used(span) + 1;
        y = baseY + pad*(used(span)-0.5 + 0.8*(span-1));

        plot(ax(gi), [x1 x1 x2 x2],[y-0.4*pad y y y-0.4*pad],'k-','LineWidth',0.9);
        text(ax(gi), mean([x1 x2]), y+0.05*pad, starStr(pv(r)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9,'FontWeight','bold');

        yl(2) = max(yl(2), y+pad);
    end
    ylim(ax(gi), yl);
end

% FORCE SAME Y across the two panels (critical)
if all(isgraphics(ax))
    yl1 = ylim(ax(1)); yl2 = ylim(ax(2));
    yl  = [min(yl1(1),yl2(1)), max(yl1(2),yl2(2))];
    ylim(ax(1), yl);
    ylim(ax(2), yl);
end

if useDelta
    fn = fullfile(outDir, sprintf('rm_oneway_%s_DELTA.png', safeName(ycol)));
else
    fn = fullfile(outDir, sprintf('rm_oneway_%s_RAW.png',   safeName(ycol)));
end

printpng(fh, fn);
close(fh);
end

%% ========================================================================
% STATS for PERIODS
function ST = runStats_Period(W, ycol)
T = W(:,{'mouse_key','Group','Period','value'});
T.mouse_key = categorical(T.mouse_key);
T.Group     = categorical(T.Group);
T.Period    = categorical(T.Period, {'Pre','During','Post','Withdrawal','Re-exposure'}, 'Ordinal',true);

rows = {};

% Within-group paired: Pre vs {During, Post, Withdrawal, Re-exposure}
G = ["Passive","Active"];
Ecmp = ["During","Post","Withdrawal","Re-exposure"];
for g = 1:numel(G)
    for e = 1:numel(Ecmp)
        A = groupsummary(T(T.Group==G(g) & T.Period=="Pre",:),    'mouse_key','median','value');
        B = groupsummary(T(T.Group==G(g) & T.Period==Ecmp(e),:),  'mouse_key','median','value');
        M = outerjoin(A(:,{'mouse_key','median_value'}), B(:,{'mouse_key','median_value'}), ...
            'Keys','mouse_key','MergeKeys',true);
        xa = M.median_value_left; xb = M.median_value_right;
        good = isfinite(xa) & isfinite(xb);
        p = NaN; N = nnz(good);
        if N >= 3, p = signrank(xa(good), xb(good)); end
        rows(end+1,:) = {'WithinNP','Group='+string(G(g)),'Pre|'+Ecmp(e), p, N}; %#ok<AGROW>
    end
end

% Between-group at each period
Eall = ["Pre","During","Post","Withdrawal","Re-exposure"];
for e = 1:numel(Eall)
    Pa = T(T.Group=="Passive" & T.Period==Eall(e) & isfinite(T.value),:);
    Ac = T(T.Group=="Active"  & T.Period==Eall(e) & isfinite(T.value),:);
    x = Pa.value; y = Ac.value;
    p = NaN; N = numel(x)+numel(y);
    if numel(x)>=2 && numel(y)>=2, p = ranksum(x,y); end
    rows(end+1,:) = {'BetweenNP','Period='+string(Eall(e)),'Passive vs Active', p, N}; %#ok<AGROW>
end

ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
ST.metric = repmat(string(ycol), height(ST),1);
end

function ST = runDeltaStats_Period(W, ycol)
base = groupsummary(W(W.Period=="Pre",:), 'mouse_key', 'median', 'value');
base.Properties.VariableNames{'median_value'} = 'base';
T = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key', 'MergeKeys',true);
T.delta = T.value - T.base;

E = ["During","Post","Withdrawal","Re-exposure"];
T = T(ismember(string(T.Period), E) & isfinite(T.delta),:);

rows = {};
G = ["Passive","Active"];

% within-group, paired signrank across period pairs
for g = 1:numel(G)
    Tg = T(T.Group==G(g),:);
    for i = 1:numel(E)-1
        for j = i+1:numel(E)
            A = outerjoin( ...
                Tg(Tg.Period==E(i), {'mouse_key','delta'}), ...
                Tg(Tg.Period==E(j), {'mouse_key','delta'}), ...
                'Keys','mouse_key','MergeKeys',true);
            xa = A.delta_left; xb = A.delta_right;
            good = isfinite(xa) & isfinite(xb); N = nnz(good); p = NaN;
            if N >= 3, p = signrank(xa(good), xb(good)); end
            rows(end+1,:) = {'DeltaWithinNP', 'Group='+string(G(g)), E(i)+"|"+E(j), p, N}; %#ok<AGROW>
        end
    end
end

% between groups at each period (ranksum)
for e = 1:numel(E)
    Pa = T(T.Group=="Passive" & T.Period==E(e),:);
    Ac = T(T.Group=="Active"  & T.Period==E(e),:);
    x = Pa.delta; y = Ac.delta; x=x(isfinite(x)); y=y(isfinite(y)); p = NaN; N = numel(x)+numel(y);
    if numel(x)>=2 && numel(y)>=2, p = ranksum(x, y); end
    rows(end+1,:) = {'DeltaBetweenNP', 'Period='+string(E(e)), 'Passive vs Active', p, N}; %#ok<AGROW>
end

% all cell-pairs (unpaired ranksum)
cells = unique(T(:,{'Group','Period'}),'rows','stable');
for i = 1:height(cells)-1
    for j = i+1:height(cells)
        A = T(T.Group==cells.Group(i) & T.Period==cells.Period(i),:);
        B = T(T.Group==cells.Group(j) & T.Period==cells.Period(j),:);
        x = A.delta; y = B.delta; x=x(isfinite(x)); y=y(isfinite(y)); p = NaN; N = numel(x)+numel(y);
        if numel(x)>=2 && numel(y)>=2, p = ranksum(x, y); end
        lab = string(cells.Group(i))+"|"+string(cells.Period(i))+" vs "+ ...
              string(cells.Group(j))+"|"+string(cells.Period(j));
        rows(end+1,:) = {'DeltaAllPairsNP', 'CellPair', lab, p, N}; %#ok<AGROW>
    end
end

ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
ST.metric = repmat(string(ycol), height(ST),1);
end

function STrm = runRMAnovaWithinGroup_Period(W, ycol, useDelta)
if nargin<3, useDelta = true; end
Eall = ["Pre","During","Post","Withdrawal","Re-exposure"];
if useDelta
    E = Eall(2:end);
else
    E = Eall;
end
rows = {};

for G = ["Passive","Active"]

    if useDelta
        base = groupsummary(W(W.Group==G & W.Period=="Pre",:), 'mouse_key', 'median', 'value');
        base.Properties.VariableNames{'median_value'} = 'base';
        Tg = outerjoin(W(W.Group==G & ismember(W.Period,E),:), base(:,{'mouse_key','base'}), ...
            'Keys','mouse_key','MergeKeys',true);
        Tg.value = Tg.value - Tg.base;
    else
        Tg = W(W.Group==G & ismember(W.Period,E),:);
    end

    % Drop periods with insufficient data for this group
    perCounts = periodFiniteCounts(Tg);
    keepPer = perCounts.Period(perCounts.count >= 3);
    Euse = E(ismember(E, string(keepPer)));
    if numel(Euse) < 3
        rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Period', NaN, 0}; %#ok<AGROW>
        continue
    end

    Wg    = unstack(Tg(:,{'mouse_key','Period','value'}), 'value', 'Period');
    safeE = matlab.lang.makeValidName(cellstr(Euse));
    have  = safeE(ismember(safeE, Wg.Properties.VariableNames));
    if numel(have) < 3
        rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Period', NaN, height(Wg)}; %#ok<AGROW>
        continue
    end
    mask = all(isfinite(Wg{:,have}),2);
    Wc   = Wg(mask, [{'mouse_key'}, have]);
    if height(Wc) < 3
        rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Period', NaN, height(Wc)}; %#ok<AGROW>
        continue
    end

    orig = strings(numel(have),1);
    for i=1:numel(have)
        orig(i) = Euse(strcmp(have{i}, safeE));
    end
    within = table(categorical(orig,'Ordinal',true), 'VariableNames', {'Period'});

    measStr = strjoin(have, ',');
    rm  = fitrm(Wc, sprintf('%s ~ 1', measStr), 'WithinDesign', within);
    RA  = ranova(rm, 'WithinModel', 'Period');

    termLabels = ranovaTermLabels(RA);
    idx = find(strcmpi(termLabels,'Period') | contains(lower(termLabels),'period'), 1, 'first');
    if isempty(idx)
        idx = find(~contains(lower(termLabels),'error'), 1, 'first');
    end
    if isempty(idx) || ~ismember('pValue', RA.Properties.VariableNames)
        pOmni = NaN;
    else
        pOmni = RA.pValue(idx);
    end
    rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Period', pOmni, height(Wc)}; %#ok<AGROW>

    try
        C = multcompare(rm, 'Period', 'ComparisonType', 'tukey-kramer');
    catch
        C = multcompare(rm, 'Period', 'ComparisonType', 'bonferroni');
    end

    L1 = 'LevelA'; L2 = 'LevelB';
    vnames = string(C.Properties.VariableNames);
    if ~ismember(L1, vnames)
        i1 = find(endsWith(vnames,'_1'), 1, 'first');
        i2 = find(endsWith(vnames,'_2'), 1, 'first');
        if ~isempty(i1) && ~isempty(i2)
            L1 = char(vnames(i1)); L2 = char(vnames(i2));
        end
    end
    if ismember('pValue', vnames), pv = C.pValue;
    elseif ismember('pValueAdj', vnames), pv = C.pValueAdj;
    else, pv = NaN(height(C),1);
    end

    % record UNIQUE only (no duplicates)
    seen = strings(0,1);
    for r = 1:height(C)
        a = string(C.(L1)(r)); b = string(C.(L2)(r));
        pair = sort([a b]);
        key  = pair(1) + "_" + pair(2);
        if any(seen==key), continue; end
        seen(end+1) = key; %#ok<AGROW>
        lvl = a + "|" + b;
        rows(end+1,:) = {'RM-ANOVA posthoc','Group='+string(G), lvl, pv(r), height(Wc)}; %#ok<AGROW>
    end
end

STrm = cell2table(rows, 'VariableNames', {'test','effect','level','p','N'});
STrm.metric = repmat(string(ycol), height(STrm), 1);
end

%% ========================================================================
% FAST SESSION->DAY METRICS (UPDATED with lick/bout freq per 10s + PR derived)
function [S, D] = fast_session_day_metrics(T, runDir)
cacheMat = fullfile(runDir, 'S_D_cache.mat');
if exist(cacheMat,'file')
    L = load(cacheMat,'S','D','cache_hash');
    if isfield(L,'cache_hash') && isequal(L.cache_hash, local_hash(T))
        % Validate cache has required derived metrics
        requiredCols = {
            'lick_freq_per_10s','lick_meanDur_s','bout_meanDur_s', ...
            'bout_freq_per_10s','bout_freq_per_min', ...
            'Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min'};
        missingCols = setdiff(requiredCols, L.D.Properties.VariableNames);
        if isempty(missingCols)
            fprintf('Loaded S/D from cache: %s\n', cacheMat);
            S = L.S; D = L.D; return;
        else
            fprintf('Cache missing %d derived columns; recomputing S/D.\n', numel(missingCols));
        end
    end
end

% bout definition
boutGap_s = 1; % edit if your standard differs

tic;
need = {'mouse_key','day_index','day_name','session_idx','Diameter_px', ...
        'Lick_TTL','Injector_TTL','CamTime_rel_s','PupilTimestamp_s', ...
        'CamTime_s','PlotTime_s_30fps','Session_Paradigm','isPassive','RequirementLast'};
V = T.Properties.VariableNames;
keepExtras = V( contains(V,'Immersion_Latency','IgnoreCase',true) | ...
                contains(V,'TST_','IgnoreCase',true)            | ...
                contains(V,'HOT_','IgnoreCase',true) );
need = unique([need, keepExtras], 'stable');
need = intersect(need, V, 'stable');
T = T(:, need);

if ~isstring(T.mouse_key),        T.mouse_key = string(T.mouse_key); end
if ~isstring(T.day_name),         T.day_name  = string(T.day_name);  end
if ismember('Session_Paradigm', T.Properties.VariableNames) && ~isstring(T.Session_Paradigm)
    T.Session_Paradigm = string(T.Session_Paradigm);
end
T.mouse_key  = categorical(T.mouse_key);
T.day_name   = categorical(T.day_name);
if ismember('Session_Paradigm', T.Properties.VariableNames)
    T.Session_Paradigm = categorical(T.Session_Paradigm);
end
if ismember('isPassive', T.Properties.VariableNames) && ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

if ismember('Lick_TTL', T.Properties.VariableNames)
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;
end
if ismember('Injector_TTL', T.Properties.VariableNames)
    T.Injector_TTL(isnan(T.Injector_TTL)) = 0; T.Injector_TTL = T.Injector_TTL > 0.5;
end

[g, keys_mouse, keys_day, keys_dayname, keys_sess] = findgroups( ...
    T.mouse_key, T.day_index, T.day_name, T.session_idx);
nG = max(g);
fprintf('Computing per-session metrics for %d sessions...\n', nG);

S = table();
S.mouse_key        = removecats(keys_mouse);
S.day_index        = double(keys_day);
S.day_name         = removecats(keys_dayname);
S.session_idx      = double(keys_sess);
S.RequirementLast  = nan(nG,1);
S.isPassive        = nan(nG,1);
S.SessionMinutes   = nan(nG,1);
S.Session_Paradigm = strings(nG,1);

vars = { ...
    'lick_n','lick_freq_per_min','lick_freq_per_10s', ...
    'lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
    'bout_n','bout_freq_per_10s','bout_freq_per_min', ...
    'bout_meanDur_s','bout_totalDur_s', ...
    'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s', ...
    'pupil_mean'};
for v = vars, S.(v{1}) = nan(nG,1); end
for j=1:numel(keepExtras), S.(keepExtras{j}) = nan(nG,1); end

tb_all = pickTimebase_fast(T);
lastTick = tic; step = max(1, floor(nG/50));
for k = 1:nG
    idx = (g == k);
    tb  = tb_all(idx);
    dur_s = finiteRange_fast(tb);
    S.SessionMinutes(k) = dur_s/60;

    if ismember('RequirementLast',T.Properties.VariableNames)
        S.RequirementLast(k) = mean(double(T.RequirementLast(idx)),'omitnan');
    end
    if ismember('isPassive',T.Properties.VariableNames)
        ip = double(T.isPassive(idx)); ip = ip(isfinite(ip));
        if ~isempty(ip), S.isPassive(k) = mode(round(ip)); end
    end
    if ismember('Session_Paradigm', T.Properties.VariableNames)
        sp = T.Session_Paradigm(idx);
        if iscategorical(sp), sp = string(mode(sp(~ismissing(sp))));
        else, sp = string(sp); sp = sp(~ismissing(sp)); if isempty(sp), sp=""; else, sp=sp(1); end
        end
        S.Session_Paradigm(k) = sp;
    end

    if ismember('Diameter_px',T.Properties.VariableNames)
        S.pupil_mean(k) = mean(double(T.Diameter_px(idx)),'omitnan');
    end

    if ismember('Lick_TTL',T.Properties.VariableNames)
        ttlL = logical(T.Lick_TTL(idx));
        [n,md,td,iei,onL,offL] = eventMetrics_fast(tb, ttlL);
        S.lick_n(k)=n;
        S.lick_meanDur_s(k)=md;
        S.lick_totalDur_s(k)=td;
        S.lick_medianIEI_s(k)=iei;

        if S.SessionMinutes(k)>0
            S.lick_freq_per_min(k)=n/S.SessionMinutes(k);
        end
        if dur_s > 0
            S.lick_freq_per_10s(k) = n / (dur_s/10);
        end

        [bn,bmd,btd] = boutMetrics_fromEvents(tb, onL, offL, boutGap_s);
        S.bout_n(k)         = bn;
        S.bout_meanDur_s(k) = bmd;
        S.bout_totalDur_s(k)= btd;

        if S.SessionMinutes(k)>0
            S.bout_freq_per_min(k) = bn/S.SessionMinutes(k);
        end
        if dur_s > 0
            S.bout_freq_per_10s(k) = bn / (dur_s/10);
        end
    end

    if ismember('Injector_TTL',T.Properties.VariableNames)
        [n,md,td,iri] = eventMetrics_fast(tb, logical(T.Injector_TTL(idx)));
        S.rew_n(k)=n; S.rew_meanDur_s(k)=md; S.rew_totalDur_s(k)=td; S.rew_medianIRI_s(k)=iri;
        if S.SessionMinutes(k)>0, S.rew_freq_per_min(k)=n/S.SessionMinutes(k); end
    end

    for j=1:numel(keepExtras)
        col = keepExtras{j};
        if ismember(col, T.Properties.VariableNames)
            S.(col)(k) = mean(double(T.(col)(idx)),'omitnan');
        end
    end

    if mod(k, step)==0 && toc(lastTick)>0.5
        fprintf('  %d/%d (%.0f%%)\n', k, nG, 100*k/nG); lastTick = tic;
    end
end
fprintf('Per-session metrics done in %.1f s.\n', toc);

fprintf('Collapsing to per-day medians...\n');
[g2, mk2, di2, dn2] = findgroups(S.mouse_key, S.day_index, S.day_name);
D = table(removecats(mk2), double(di2), removecats(dn2), ...
          'VariableNames',{'mouse_key','day_index','day_name'});

baseList  = [{'RequirementLast','isPassive','SessionMinutes'}, vars];
extraList = intersect(keepExtras, S.Properties.VariableNames, 'stable');
list = unique([baseList, extraList], 'stable');

for v = list
    D.(v{1}) = splitapply(@(x) median(x,'omitnan'), S.(v{1}), g2);
end

% --- PR-derived daily metrics ---
if ismember('RequirementLast', D.Properties.VariableNames)
    D.Requirement_cum = nan(height(D),1);
    D.Requirement_speed_per_day = nan(height(D),1);
    D.Requirement_speed_per_min = nan(height(D),1);

    mice = unique(D.mouse_key,'stable');
    for i = 1:numel(mice)
        msk = D.mouse_key==mice(i);
        [dd,ord] = sort(double(D.day_index(msk)));
        idx = find(msk); idx = idx(ord);

        r = double(D.RequirementLast(idx));
        r(~isfinite(r)) = NaN;

        % cumulative
        rcum = nan(size(r));
        acc = 0;
        for j = 1:numel(r)
            if isfinite(r(j)), acc = acc + r(j); end
            rcum(j) = acc;
        end
        D.Requirement_cum(idx) = rcum;

        % speed per day
        spd = nan(size(r));
        for j = 2:numel(r)
            if isfinite(r(j)) && isfinite(r(j-1)) && isfinite(dd(j)) && isfinite(dd(j-1)) && (dd(j)-dd(j-1))~=0
                spd(j) = (r(j)-r(j-1)) / (dd(j)-dd(j-1));
            end
        end
        D.Requirement_speed_per_day(idx) = spd;

        % speed per min proxy
        if ismember('SessionMinutes', D.Properties.VariableNames)
            mins = double(D.SessionMinutes(idx));
            good = isfinite(r) & isfinite(mins) & mins>0;
            tmp = nan(size(r));
            tmp(good) = r(good) ./ mins(good);
            D.Requirement_speed_per_min(idx) = tmp;
        end
    end
end

cache_hash = local_hash(T(:, {'mouse_key','day_index','session_idx'}));
try, save(cacheMat,'S','D','cache_hash','-v7.3'); fprintf('Cached S/D to: %s\n', cacheMat); catch, end
end

%% ========================================================================
% Bout helper
function [boutN, boutMeanDur, boutTotalDur] = boutMetrics_fromEvents(t, on, off, boutGap_s)
t = double(t(:));
boutN = 0; boutMeanDur = NaN; boutTotalDur = 0;
if isempty(on) || isempty(off) || numel(on)~=numel(off)
    return;
end

onT  = t(on);
% conservative: off uses t(off) (frame timestamp), not off+1
offT = t(off);
bad = ~isfinite(onT) | ~isfinite(offT);
onT(bad) = []; offT(bad) = [];
if isempty(onT), return; end

gaps = [Inf; diff(onT)];
isNew = gaps > boutGap_s;

bStart = find(isNew);
bEnd   = [bStart(2:end)-1; numel(onT)];

boutN = numel(bStart);
boutDur = nan(boutN,1);
for b = 1:boutN
    i1 = bStart(b);
    i2 = bEnd(b);
    boutDur(b) = offT(i2) - onT(i1);
end
boutDur = boutDur(isfinite(boutDur) & boutDur>=0);

if isempty(boutDur)
    boutMeanDur = NaN; boutTotalDur = 0;
else
    boutMeanDur  = mean(boutDur,'omitnan');
    boutTotalDur = sum(boutDur,'omitnan');
end
end

%% ========================================================================
% SMALL HELPERS
function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end

function pick = pickPctColumn(cols, preferKey)
pick = '';
if isempty(cols), return; end
pctCols = cols(contains(cols,'Pct'));
if isempty(pctCols), return; end
pref = pctCols(contains(lower(pctCols), lower(preferKey)));
if ~isempty(pref), pick = pref{1}; else, pick = pctCols{1}; end
end

function nm = colLike(T, key)
nm = '';
if istable(T)
    V = T.Properties.VariableNames;
    hit = find(strcmpi(V,key) | contains(lower(V), lower(key)), 1, 'first');
    if ~isempty(hit), nm = V{hit}; end
end
end

function s = safeName(nm), s = regexprep(nm,'[^a-zA-Z0-9]+','_'); end

function printpng(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end

function tb = pickTimebase_fast(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
tb = nan(height(T),1);
for i=1:numel(cands)
    if ismember(cands{i}, T.Properties.VariableNames)
        v = double(T.(cands{i})); if any(isfinite(v)), tb = v; return; end
    end
end
end

function r = finiteRange_fast(x)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), r = 0; else, r = max(x)-min(x); end
end

function [n, meanDur, totalDur, medianIEI, on, off] = eventMetrics_fast(t, ttl)
t   = double(t(:)); ttl = logical(ttl(:));
good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);

on = []; off = [];
if numel(t)<2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end

dt  = diff(t); md = median(dt(isfinite(dt))); if ~isfinite(md), md=1/30; end
t(2:end) = max(t(2:end), t(1:end-1)+md*0.5);

d   = diff([false; ttl; false]); on = find(d==1); off = find(d==-1)-1;
n   = numel(on);
if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; on=[]; off=[]; return; end

edges  = [t; t(end)+md];
segDur = sum(edges(off+1) - edges(on), 2, 'omitnan');
meanDur  = mean(segDur,'omitnan'); totalDur = sum(segDur,'omitnan');

if n>=2, medianIEI = median(diff(t(on)),'omitnan'); else, medianIEI = NaN; end
end

function h = local_hash(Tkeys)
try
    raw = [uint32(double(grp2idx(categorical(Tkeys.mouse_key)))), ...
           uint32(double(Tkeys.day_index)), uint32(double(Tkeys.session_idx))];
    h = uint64(sum(uint64(raw(:)).*1664525 + 1013904223));
catch, h = now;
end
end

function labs = ranovaTermLabels(tbl)
if ismember('Term', tbl.Properties.VariableNames)
    labs = string(tbl.Term);
elseif ~isempty(tbl.Properties.RowNames)
    labs = string(tbl.Properties.RowNames);
else
    v = tbl.Properties.VariableNames;
    cand = v(contains(lower(v),'term') | contains(lower(v),'effect') | contains(lower(v),'within'));
    if ~isempty(cand)
        labs = string(tbl.(cand{1}));
    else
        labs = strings(height(tbl),1);
    end
end
end

function s = starStr(p)
if ~isfinite(p), s = 'n/a'; return; end
if p < 1e-4, s = '****';
elseif p < 1e-3, s = '***';
elseif p < 1e-2, s = '**';
elseif p < 0.05, s = '*';
else, s = 'n.s.';
end
end

function q = fdrBH(p)
p = double(p); m = numel(p); [ps,idx] = sort(p); q = nan(size(p));
ranks = (1:m)'; adj = ps.*m./ranks;
for i=m-1:-1:1, adj(i)=min(adj(i),adj(i+1)); end
q(idx) = min(adj,1);
end

function p_holm = holmBonferroni(p)
p = double(p(:)); [ps,idx] = sort(p); m = numel(p);
adj = (m - (1:m)' + 1) .* ps;
for k = 2:m, adj(k) = max(adj(k), adj(k-1)); end
p_holm = zeros(size(p)); p_holm(idx) = min(adj,1);
end

%% ========================================================================
% Legacy passive classifier kept as fallback only
function HadPassive = classifyHadPassive(S, D)
mice = unique(S.mouse_key,'stable');
HadPassive = false(height(D),1);
havePar = ismember('Session_Paradigm', S.Properties.VariableNames);
for i = 1:numel(mice)
    rS  = S.mouse_key==mice(i);
    ip  = S.isPassive(rS);
    di  = S.day_index(rS);
    par = strings(nnz(rS),1);
    if havePar
        p = S.Session_Paradigm(rS);
        if iscategorical(p), p = string(p); end
        par = string(p);
    end
    flag = false;
    if any(ip==1,'all')
        flag = true;
    elseif havePar && any(contains(lower(par),'passive'),'all')
        flag = true;
    elseif any(ismember(double(di),6:8) & (ip>=1 | isnan(ip)),'all')
        flag = true;
    end
    Dmask = D.mouse_key==mice(i);
    HadPassive(Dmask) = flag;
end
end

%% ========================================================================
% CORR + PCA (grid removed)
function corrPCA(Wwide, outDir)
COL = palette();
dat = Wwide; key = dat(:,1:3); X = dat{:,4:end};
vars = dat.Properties.VariableNames(4:end);

C = corr(X,'rows','pairwise');
fh = figure('Color','w','Position',[60 60 650 550]);
imagesc(C); axis image; colorbar; caxis([-1 1]);
set(gca,'XTick',1:numel(vars),'XTickLabel',vars,'XTickLabelRotation',45);
set(gca,'YTick',1:numel(vars),'YTickLabel',vars);
title('Metric correlation (per mouse×period)');
grid off;
box on;
printpng(fh, fullfile(outDir,'corr_heatmap.png')); close(fh);

Xz = X;
for j=1:size(Xz,2)
    mu = nanmean(Xz(:,j)); sd = nanstd(Xz(:,j));
    if sd>0 && isfinite(sd)
        Xz(:,j) = (Xz(:,j)-mu)./sd;
    else
        Xz(:,j) = 0;
    end
end
Xz(~isfinite(Xz))=0;

[~,score,~,~,expl] = pca(Xz);

plotPcaScatter(score, expl, key, COL, [1 2], outDir, 'pca_scatter_pc1_pc2.png');
if size(score,2) >= 3
    plotPcaScatter(score, expl, key, COL, [2 3], outDir, 'pca_scatter_pc2_pc3.png');
end
if size(score,2) >= 4
    plotPcaScatter(score, expl, key, COL, [3 4], outDir, 'pca_scatter_pc3_pc4.png');
end
end

function plotPcaScatter(score, expl, key, COL, dims, outDir, fname)
fh = figure('Color','w','Position',[80 80 740 540]); hold on;
grp = string(key.Group); per = string(key.Period);
colors = zeros(numel(grp),3);
colors(grp=="Passive",:) = repmat(COL.passive, sum(grp=="Passive"), 1);
colors(grp=="Active", :) = repmat(COL.active,  sum(grp=="Active"),  1);

for e = 1:numel(COL.periodOrder)
    mask = (per==COL.periodOrder{e});
    if any(mask)
        scatter(score(mask,dims(1)), score(mask,dims(2)), 46, colors(mask,:), 'filled', ...
            'Marker', COL.periodMark{e}, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.85);
    end
end
hG(1) = scatter(nan,nan,60,'o','filled','MarkerFaceColor',COL.passive,'MarkerEdgeColor','k','DisplayName','Passive');
hG(2) = scatter(nan,nan,60,'o','filled','MarkerFaceColor',COL.active, 'MarkerEdgeColor','k','DisplayName','Active');
hE = gobjects(numel(COL.periodOrder),1);
for e=1:numel(COL.periodOrder)
    hE(e) = plot(nan,nan,'k','LineStyle','none','Marker',COL.periodMark{e}, 'MarkerSize',8,'DisplayName',COL.periodOrder{e});
end
lg1 = legend(hG, 'Location','northeastoutside'); title(lg1,'Group');
lg2 = legend(hE, 'Location','southeastoutside');  title(lg2,'Period');
uistack(lg1,'top'); uistack(lg2,'top');

xlabel(sprintf('PC%d (%.1f%%)', dims(1), expl(dims(1))));
ylabel(sprintf('PC%d (%.1f%%)', dims(2), expl(dims(2))));
grid off;
box on;
title(sprintf('PCA of mouse×period metric vectors (PC%d vs PC%d)', dims(1), dims(2)));
printpng(fh, fullfile(outDir, fname)); close(fh);
end
