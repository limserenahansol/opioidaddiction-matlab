function make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203()
% Uses YOUR longitudinal CSV schema (no guessing):
% Required:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1), Injector_TTL (0/1), Diameter_px (numeric), RequirementLast (numeric)
%   Time axis: CamTime_rel_s OR PlotTime_s_30fps  (either one must exist)
%
% Outputs: QC + requested analyses (licking + pupil + PR + period comparisons + optional tests)

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Session-time binning for cumulative lick curves
dt_bin_s = 60;          % 1-min bins for cumulative licking across session
bout_gap_s = 2.0;       % lick bout definition: pause > 2s splits bouts
pupil_win = [-2 2];     % seconds around event for pupil event-lock
min_event_separation_s = 0.5; % to prevent degenerate overlaps in event extraction

% Period definitions (days)
PER.pre       = 3:5;     % water PR
PER.during    = 6:10;    % morphine period (active PR; passive forced/yoked)
PER.post      = 11:13;   % morphine PR all
PER.withdraw  = 14:16;   % withdrawal PR (water)
PER.reexpo    = 17:18;   % morphine PR re-exposure

% "transition/unreliable days" you asked to optionally exclude
% (You mentioned day4, day6, day11, day14 are less reliable)
transitionDays = [4 6 11 14];

% Passive-specific period comparisons you requested
% Passive has no PR score during 6-10 in your design, but we still analyze licking/pupil.
PASSIVE_COMPARE.focus_during_days = 7:9;    % “passive focus” subset within 6–10
PASSIVE_COMPARE.focus_post_days   = 12:13;

% Event-lock pupil focus comparison you requested:
% passive mouse only: day7-9 vs day12-13; also passive vs active on those days
FOCUS_PASSIVE_DAYS_A = 7:9;
FOCUS_PASSIVE_DAYS_B = 12:13;

% If you want cumulative licking to exclude habituation days 1-2:
excludeDays_licking_global = [1 2];

%% ===================== FIND LATEST run_* =====================
if ~exist(rootTry,'dir')
    error('rootTry not found: %s', rootTry);
end
D = dir(fullfile(rootTry,'run_*'));
assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')==2, 'Expected CSV not found: %s', csvPath);

ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(runDir, ['QC_AND_REQUESTED_ANALYSES_' ts]);
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\nUsing runDir:\n  %s\n', runDir);
fprintf('Reading CSV:\n  %s\n', csvPath);
fprintf('Saving outputs to:\n  %s\n\n', outDir);

%% ===================== LOAD TABLE =====================
T = readtable(csvPath, 'VariableNamingRule','preserve');

% enforce mouse_key string
if ismember('mouse_key', T.Properties.VariableNames)
    if ~isstring(T.mouse_key), T.mouse_key = string(T.mouse_key); end
else
    error('Missing required column: mouse_key');
end

% required columns (schema-true)
REQ = {'mouse_key','day_index','session_idx','Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, REQ, 'Core longitudinal schema');

% time axis (explicitly allowed alternatives from YOUR pipeline)
timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end

% cohort metadata (YOUR new cohort)
COH = buildNewCohortTable(); % cage/color/sex/group/pair

% join cohort metadata onto rows by mouse_key
T = addCohortMeta(T, COH);

% drop days with NaN day_index
T = T(~isnan(T.day_index),:);

%% ===================== QC: BASIC SANITY =====================
% 1) Missingness heatmaps (lick, pupil, reward) by mouse x day
doMissingnessQC(T, outDir);

% 2) Session-level summaries (lick count, lick/min, IEI stats, PR score, pupil mean)
S = buildSessionSummary(T, timeVar, dt_bin_s);

% 3) Outlier QC plots (robust MAD) for key metrics
doOutlierQC(S, outDir);

% 4) Spaghetti QC: lick/min and PR score across days (active/passive)
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'RequirementLast', outDir, 'QC_spaghetti_PRscore_RequirementLast.png');
plotPupilByGroupDay(S, outDir); % day-level pupil mean (raw)

%% ===================== REQUESTED ANALYSES: LICKING =====================

% A) Cumulative licking across session time: day-by-day, Active vs Passive
plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global);

% B) Early vs Late comparisons (example like your figure): Active early(days3-5) vs late(days6-10)
plotEarlyLateCumulative(T, timeVar, outDir, dt_bin_s, 3:5, 6:10);

% C) Lick count per day spaghetti (active vs passive)
plotDailyLickCountSpaghetti(S, outDir);

% D) Lick raster (per day) – useful QC + “something novel” detection
% (creates a couple representative days if available)
plotLickRasterExamples(T, timeVar, outDir);

%% ===================== REQUESTED ANALYSES: PERIOD COMPARISONS =====================

% Make two “day sets”: include all days vs exclude transition days
DAYSETS.allDays.name  = 'ALLDAYS';
DAYSETS.allDays.mask  = true(size(S.day_index));

DAYSETS.noTransition.name = 'NO_TRANSITION';
DAYSETS.noTransition.mask = ~ismember(S.day_index, transitionDays);

% Active: Day5 vs during(6-10) vs post(11-13) vs withdraw(14-16) vs reexpo(day17)
% Passive: Day5 vs post(11-13) vs withdraw(14-16) vs reexpo(17-18)
doPeriodComparisons(S, outDir, PER, DAYSETS, transitionDays);

%% ===================== REQUESTED ANALYSES: PUPIL =====================

% 1) Day-level pupil normalization using Day3 baseline per mouse
S = addPupilBaselineNormalization(S); % adds pupil_norm_rel (delta/base)

plotPupilNormalizedDayTraces(S, outDir);

% 2) Event-locked pupil: reward-locked and lick-bout-locked
%    - baseline subtract using pre 2s
%    - split rewarded vs non-rewarded lick bouts
doEventLockedPupilAnalyses(T, timeVar, outDir, pupil_win, bout_gap_s, min_event_separation_s, ...
    FOCUS_PASSIVE_DAYS_A, FOCUS_PASSIVE_DAYS_B);

%% ===================== REQUESTED ANALYSES: TAIL IMMERSION / HOT / TST / STRAUB =====================
% These are included if your final CSV contains the columns (from your pipeline).
doOptionalPainAndBehaviorTests(T, outDir);

fprintf('\nDONE. Outputs saved in:\n  %s\n', outDir);
end

%% =====================================================================
%% =========================== HELPERS =================================
%% =====================================================================

function assertHasVars(T, vars, context)
missing = vars(~ismember(vars, T.Properties.VariableNames));
if ~isempty(missing)
    fprintf('\nSCHEMA ERROR in %s.\nMissing required columns:\n', context);
    disp(missing(:));
    error('Fix CSV schema or update required variable list in the script.');
end
end

function COH = buildNewCohortTable()
% mouse_key format expected in your pipeline is typically: "6100_red" or "0911_black" etc.
% If your pipeline uses "6100red" without underscore, we handle that by normalizing later.

% Define each animal
rows = {
% cage   color   sex  group  pair_id  role_in_pair
'6100','red',    'F', 'P',   '6100_pairA', 'Passive';
'6100','orange', 'F', 'P',   '6100_pairA', 'Passive';
'6100','black',  'F', 'A',   '6100_pairA', 'Active';

'0911','red',    'F', 'A',   '0911_pairA', 'Active';
'0911','orange', 'F', 'P',   '0911_pairA', 'Passive';

'0911','black',  'F', 'P',   '0911_pairB', 'Passive';
'0911','white',  'F', 'A',   '0911_pairB', 'Active';

'0910','red',    'M', 'P',   '0910_pairA', 'Passive';
'0910','orange', 'M', 'P',   '0910_pairA', 'Passive';
'0910','black',  'M', 'A',   '0910_pairA', 'Active';

'6099','red',    'M', 'P',   '6099_pairA', 'Passive';
'6099','orange', 'M', 'A',   '6099_pairA', 'Active';

'6099','black',  'M', 'A',   '6099_pairB', 'Active';
'6099','white',  'M', 'P',   '6099_pairB', 'Passive'; % died after day13 (handled by missingness)
};

COH = cell2table(rows, 'VariableNames', {'cage','color','sex','group','pair_id','pair_role'});

% build canonical mouse_key as cage_color
COH.mouse_key = string(COH.cage) + "_" + string(COH.color);
end

function T = addCohortMeta(T, COH)
% Normalize mouse_key to cage_color form for joining
mk = string(T.mouse_key);
mk = regexprep(mk,'\s+','');
mk = regexprep(mk,'-','_');

% if keys are like "6100red" -> convert to "6100_red"
mk2 = mk;
isNoUnd = ~contains(mk2,'_');
mk2(isNoUnd) = regexprep(mk2(isNoUnd), '^(\d{4})([A-Za-z]+)$', '$1_$2');

% lowercase color part for robust join
mk2 = lower(mk2);

% normalize COH key similarly
cohKey = lower(string(COH.mouse_key));

T.mouse_key_norm = mk2;
COH.mouse_key_norm = cohKey;

% left join (keeps all T rows)
T = outerjoin(T, COH(:,{'mouse_key_norm','cage','color','sex','group','pair_id','pair_role'}), ...
    'Keys','mouse_key_norm', 'MergeKeys', true, 'Type','left');

% Define GroupMouse categorical for plotting
% group: 'A'/'P'
g = string(T.group);
g(ismissing(g)) = "U";
T.GroupMouse = categorical(g, ["A","P","U"], {'Active','Passive','Unknown'});

end

function doMissingnessQC(T, outDir)
% Missingness per mouse x day for key signals
vars = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
vars = vars(ismember(vars, T.Properties.VariableNames));

mk = unique(T.mouse_key_norm,'stable');
days = unique(T.day_index);
days = sort(days(~isnan(days)));

M = nan(numel(mk), numel(days), numel(vars));
for i=1:numel(mk)
    for j=1:numel(days)
        r = T.mouse_key_norm==mk(i) & T.day_index==days(j);
        if ~any(r), continue; end
        for v=1:numel(vars)
            x = T.(vars{v})(r);
            M(i,j,v) = mean(isnan(x));
        end
    end
end

for v=1:numel(vars)
    fig = figure('Color','w','Position',[80 80 1100 520]);
    imagesc(days, 1:numel(mk), M(:,:,v));
    colormap(parula); colorbar;
    xlabel('Day'); ylabel('Mouse'); title(['Missing fraction: ' vars{v}]);
    set(gca,'YTick',1:numel(mk),'YTickLabel',mk, 'TickLabelInterpreter','none');
    exportgraphics(fig, fullfile(outDir, ['QC_missing_' vars{v} '.png']), 'Resolution', 200);
    close(fig);
end
end

function S = buildSessionSummary(T, timeVar, dt_bin_s)
% Builds one row per (mouse, day, session) with lick + pupil + PR metrics.

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

S = table;
S.mouse_key = splitapply(@(x)x(1), T.mouse_key_norm, G);
S.day_index = splitapply(@(x)x(1), T.day_index, G);
S.session_idx = splitapply(@(x)x(1), T.session_idx, G);

% group labels
if ismember('GroupMouse', T.Properties.VariableNames)
    S.GroupMouse = splitapply(@(x)x(1), T.GroupMouse, G);
else
    S.GroupMouse = categorical(repmat("Unknown", height(S),1));
end

% PR score (RequirementLast) – session-level scalar
S.RequirementLast = splitapply(@nanmean, T.RequirementLast, G);

% pupil mean
S.pupil_mean_px = splitapply(@nanmean, T.Diameter_px, G);

% lick times per session
S.nLicks = zeros(height(S),1);
S.licks_per_min = nan(height(S),1);
S.iei_median = nan(height(S),1);
S.iei_mean = nan(height(S),1);
S.iei_cv = nan(height(S),1);

for i=1:height(S)
    r = (T.mouse_key_norm==S.mouse_key(i) & T.day_index==S.day_index(i) & T.session_idx==S.session_idx(i));
    t = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    if isempty(t) || all(isnan(t))
        continue;
    end
    % sort by time
    [t,ord] = sort(t); lk = lk(ord);

    % rising edges = lick events
    lickTimes = detectRisingEdges(t, lk);

    S.nLicks(i) = numel(lickTimes);
    dur = max(t)-min(t);
    if dur > 1
        S.licks_per_min(i) = (numel(lickTimes) / dur) * 60;
    end
    if numel(lickTimes) >= 2
        iei = diff(lickTimes);
        S.iei_median(i) = median(iei);
        S.iei_mean(i)   = mean(iei);
        S.iei_cv(i)     = std(iei)/max(eps,mean(iei));
    end
end

end

function lickTimes = detectRisingEdges(t, x01)
% robust rising-edge detection for TTL-like signals
x01 = x01(:);
x01(~isfinite(x01)) = 0;
x01 = x01 > 0.5;
dx = diff([false; x01]);
idx = find(dx==1);
lickTimes = t(idx);
% drop non-monotonic issues
lickTimes = lickTimes(isfinite(lickTimes));
end

function doOutlierQC(S, outDir)
% Robust MAD-based outlier flagging for key metrics
metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

for k=1:numel(metrics)
    x = double(S.(metrics{k}));
    [isOut, z] = madOutliers(x, 4.5);
    fig = figure('Color','w','Position',[80 80 950 420]); hold on
    scatter(1:numel(x), x, 18, 'filled','MarkerFaceAlpha',0.5);
    scatter(find(isOut), x(isOut), 45, 'o','LineWidth',1.8);
    title(sprintf('QC outliers (MAD z>4.5): %s', metrics{k}), 'Interpreter','none');
    xlabel('Session index'); ylabel(metrics{k}, 'Interpreter','none'); grid on; box on
    exportgraphics(fig, fullfile(outDir, ['QC_outliers_' metrics{k} '.png']), 'Resolution', 200);
    close(fig);

    % save outlier list
    if any(isOut)
        Tout = S(isOut, {'mouse_key','day_index','session_idx','GroupMouse'});
        Tout.metric = repmat(string(metrics{k}), height(Tout),1);
        Tout.value  = x(isOut);
        Tout.madz   = z(isOut);
        writetable(Tout, fullfile(outDir, ['QC_outliers_' metrics{k} '.csv']));
    end
end
end

function [isOut, z] = madOutliers(x, thr)
x = x(:);
m = median(x,'omitnan');
madv = median(abs(x-m),'omitnan');
if madv < eps
    z = zeros(size(x));
else
    z = 0.6745*(x-m)/madv;
end
isOut = abs(z) > thr & isfinite(z);
end

function plotSpaghettiByGroup(S, yvar, outDir, fname)
if ~ismember(yvar, S.Properties.VariableNames), return; end
fig = figure('Color','w','Position',[80 80 1100 520]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

groups = categories(S.GroupMouse);
groups = groups(~strcmp(groups,'Unknown'));

for gi=1:numel(groups)
    nexttile; hold on
    g = groups{gi};
    rG = S.GroupMouse==g;
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & S.mouse_key==mk(i);
        [d,ord] = sort(S.day_index(r));
        y = double(S.(yvar)(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    % mean
    dAll = unique(S.day_index(rG)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.(yvar)(rG & S.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, 'k-', 'LineWidth',3);
    title(sprintf('%s: %s', yvar, g), 'Interpreter','none');
    xlabel('Day'); ylabel(yvar, 'Interpreter','none'); grid on; box on
end

exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotPupilByGroupDay(S, outDir)
if ~ismember('pupil_mean_px', S.Properties.VariableNames), return; end
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
groups = categories(S.GroupMouse);
groups = groups(~strcmp(groups,'Unknown'));
for gi=1:numel(groups)
    nexttile; hold on
    g = groups{gi};
    rG = S.GroupMouse==g;
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & S.mouse_key==mk(i);
        [d,ord] = sort(S.day_index(r));
        y = double(S.pupil_mean_px(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(S.day_index(rG)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.pupil_mean_px(rG & S.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, 'k-', 'LineWidth',3);
    title(sprintf('Pupil mean (px): %s', g));
    xlabel('Day'); ylabel('Pupil diameter (px)'); grid on; box on
end
exportgraphics(fig, fullfile(outDir,'QC_pupil_mean_px_by_group.png'), 'Resolution', 220);
close(fig);
end

function plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays)
% Produces:
% - Day-by-day mean cumulative licking across session time (Active)
% - Day-by-day mean cumulative licking across session time (Passive)
% - Cumulative licking across all days (excluding habituation)

% Filter out excluded days
R0 = ~ismember(T.day_index, excludeDays);

for gname = ["Active","Passive"]
    Rg = R0 & (string(T.GroupMouse)==gname);
    if ~any(Rg), continue; end

    days = unique(T.day_index(Rg)); days = sort(days);
    fig = figure('Color','w','Position',[80 80 1250 520]); hold on

    maxT = 0;
    curves = {};
    for di=1:numel(days)
        d = days(di);
        R = Rg & T.day_index==d;

        [tgrid, muCum] = meanCumulativeLicksByDay(T(R,:), timeVar, dt_bin_s);
        if isempty(tgrid), continue; end
        maxT = max(maxT, max(tgrid));

        plot(tgrid/60, muCum, 'LineWidth', 2);
        curves{end+1} = sprintf('Day %d', d); %#ok<AGROW>
    end
    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title(sprintf('Cumulative licking across session (%s mice)', gname));
    grid on; box on
    legend(curves, 'Location','eastoutside');
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_cumulative_across_session_%s.png', gname)), 'Resolution', 220);
    close(fig);
end

% Cumulative lick per day (spaghetti) – combined by group
% (already in S plots, but this one is “cumulative across days”)
% We derive daily totals from session summaries.
S = buildSessionSummary(T(R0,:), timeVar, dt_bin_s);
S.daily_licks = S.nLicks;
plotDailyCumulativeAcrossDays(S, outDir);

end

function [tgrid, muCum] = meanCumulativeLicksByDay(Tday, timeVar, dt_bin_s)
% Average cumulative lick curve across sessions (equal weight per session)
if isempty(Tday), tgrid=[]; muCum=[]; return; end

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(Tday(:,keys));

% build each session curve
sessionCurves = {};
tmaxAll = 0;
for gi=1:max(G)
    r = (G==gi);
    t = double(Tday.(timeVar)(r));
    lk = double(Tday.Lick_TTL(r));
    if isempty(t) || all(isnan(t)), continue; end
    [t,ord] = sort(t); lk=lk(ord);
    lickTimes = detectRisingEdges(t, lk);

    if isempty(lickTimes)
        tmaxAll = max(tmaxAll, max(t)-min(t));
        sessionCurves{end+1} = struct('dur', max(t)-min(t), 'lickTimes', lickTimes); %#ok<AGROW>
    else
        dur = max(t)-min(t);
        tmaxAll = max(tmaxAll, dur);
        sessionCurves{end+1} = struct('dur', dur, 'lickTimes', lickTimes - min(t)); %#ok<AGROW>
    end
end

if isempty(sessionCurves)
    tgrid=[]; muCum=[]; return;
end

% common grid
tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;
cumMat = nan(numel(sessionCurves), numel(tgrid));

for i=1:numel(sessionCurves)
    lt = sessionCurves{i}.lickTimes;
    if isempty(lt)
        cumMat(i,:) = 0;
        continue;
    end
    % cumulative count at each grid time
    cumMat(i,:) = arrayfun(@(tt) sum(lt <= tt), tgrid);
end

muCum = mean(cumMat, 1, 'omitnan');
end

function plotEarlyLateCumulative(T, timeVar, outDir, dt_bin_s, earlyDays, lateDays)
% Matches your “Early vs Late” style plot
for gname = ["Active","Passive"]
    Rg = (string(T.GroupMouse)==gname);

    fig = figure('Color','w','Position',[80 80 900 520]); hold on
    [tE, muE] = meanCumulativeLicksByDay(T(Rg & ismember(T.day_index, earlyDays),:), timeVar, dt_bin_s);
    [tL, muL] = meanCumulativeLicksByDay(T(Rg & ismember(T.day_index, lateDays),:), timeVar, dt_bin_s);

    if ~isempty(tE), plot(tE/60, muE, 'LineWidth', 3); end
    if ~isempty(tL), plot(tL/60, muL, 'LineWidth', 3); end
    xlabel('Session time (min)'); ylabel('Mean cumulative licks');
    title(sprintf('Early (days %s) vs Late (days %s) sessions (%s mice)', ...
        rangeStr(earlyDays), rangeStr(lateDays), gname));
    grid on; box on
    legend({sprintf('Days %s',rangeStr(earlyDays)), sprintf('Days %s',rangeStr(lateDays))}, 'Location','southeast');
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_early_vs_late_cumulative_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
% compact display
if all(diff(v)==1), s=sprintf('%d-%d', v(1), v(end));
else, s=strjoin(string(v),'-'); end
end

function plotDailyLickCountSpaghetti(S, outDir)
fig = figure('Color','w','Position',[80 80 900 520]); hold on
mk = unique(S.mouse_key,'stable');
for i=1:numel(mk)
    r = (S.mouse_key==mk(i));
    [d,ord] = sort(S.day_index(r));
    y = double(S.nLicks(r)); y=y(ord);
    plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
end
dAll = unique(S.day_index); dAll=sort(dAll);
mu = nan(size(dAll));
for j=1:numel(dAll)
    mu(j) = mean(double(S.nLicks(S.day_index==dAll(j))), 'omitnan');
end
plot(dAll, mu, 'k-', 'LineWidth', 3);
xlabel('Day index'); ylabel('Lick count per day'); title('Spaghetti: Lick count per day');
grid on; box on
exportgraphics(fig, fullfile(outDir,'LICK_spaghetti_lickcount_per_day.png'), 'Resolution', 220);
close(fig);
end

function plotDailyCumulativeAcrossDays(S, outDir)
% cumulative sum of daily totals per mouse (excluding missing days)
fig = figure('Color','w','Position',[80 80 1000 520]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

for gname = ["Active","Passive"]
    nexttile; hold on
    rG = (string(S.GroupMouse)==gname);
    if ~any(rG), continue; end
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & S.mouse_key==mk(i);
        [d,ord] = sort(S.day_index(r));
        y = double(S.daily_licks(r)); y=y(ord);
        ycum = cumsum(y,'omitnan');
        plot(d,ycum,'-o','LineWidth',1,'MarkerSize',3);
    end
    xlabel('Day'); ylabel('Cumulative licks'); title(sprintf('Cumulative licks across days (%s)', gname));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'LICK_cumulative_across_days_active_vs_passive.png'), 'Resolution', 220);
close(fig);
end

function plotLickRasterExamples(T, timeVar, outDir)
% Raster plots for a couple representative days (if present)
days = unique(T.day_index);
days = sort(days(~isnan(days)));
if isempty(days), return; end
pickDays = intersect(days, [4 9 12]); % examples similar to your screenshot
if isempty(pickDays), pickDays = days(min(2,numel(days))); end

for pd = pickDays(:)'
    fig = figure('Color','w','Position',[80 80 1100 520]); hold on
    % choose a few mice with that day
    mk = unique(T.mouse_key_norm(T.day_index==pd),'stable');
    mk = mk(1:min(10,numel(mk)));

    y = 0;
    for i=1:numel(mk)
        r = (T.mouse_key_norm==mk(i) & T.day_index==pd);
        t = double(T.(timeVar)(r));
        lk = double(T.Lick_TTL(r));
        if isempty(t), continue; end
        [t,ord] = sort(t); lk=lk(ord);
        lickTimes = detectRisingEdges(t, lk);
        y = y + 1;
        plot(lickTimes/60, y*ones(size(lickTimes)), 'k.', 'MarkerSize', 6);
    end
    xlabel('Session time (min)'); ylabel('Mouse (row)');
    title(sprintf('Lick raster (Day %d) – first %d mice', pd, numel(mk)));
    grid on; box on
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_raster_day%d.png', pd)), 'Resolution', 220);
    close(fig);
end
end

function doPeriodComparisons(S, outDir, PER, DAYSETS, transitionDays)
% Your requested period comparisons (with and without transition days)

% Build “period label” per row
S.period = strings(height(S),1);
S.period(ismember(S.day_index, PER.pre))      = "pre";
S.period(ismember(S.day_index, PER.during))   = "during";
S.period(ismember(S.day_index, PER.post))     = "post";
S.period(ismember(S.day_index, PER.withdraw)) = "withdrawal";
S.period(ismember(S.day_index, PER.reexpo))   = "reexposure";
S.period = categorical(S.period, ["pre","during","post","withdrawal","reexposure","<undef>"], ...
                                  {"pre","during","post","withdrawal","reexposure","<undef>"});

metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

% Active vs Passive period boxes
for dsName = fieldnames(DAYSETS)'
    ds = DAYSETS.(dsName{1});
    for m=1:numel(metrics)
        yvar = metrics{m};
        fig = figure('Color','w','Position',[80 80 1200 520]);
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        for gi=1:2
            gLabel = ["Active","Passive"];
            g = gLabel(gi);
            nexttile; hold on
            Rg = (string(S.GroupMouse)==g) & ds.mask & (S.period~="<undef>");

            % apply your passive rule for PR score: ignore during for passive PR comparisons
            Rp = Rg;
            if strcmp(yvar,'RequirementLast') && strcmp(g,'Passive')
                Rp = Rg & (S.period~="during");
            end

            if ~any(Rp)
                title(sprintf('%s (%s) – no data', yvar, g));
                axis off; continue;
            end

            % per-mouse mean per period (avoids overweighting multiple sessions)
            [Gmp, keyTbl] = findgroups(S.mouse_key(Rp), S.period(Rp));
            val = splitapply(@(x) mean(double(x),'omitnan'), S.(yvar)(Rp), Gmp);
            mk  = splitapply(@(x) x(1), S.mouse_key(Rp), Gmp);
            pr  = splitapply(@(x) x(1), S.period(Rp), Gmp);
            P = table(mk, pr, val, 'VariableNames', {'mouse','period','value'});

            periodsOrder = categories(S.period);
            periodsOrder = periodsOrder(~strcmp(periodsOrder,'<undef>'));
            boxchart(double(P.period), P.value);
            set(gca,'XTick',1:numel(periodsOrder),'XTickLabel',periodsOrder);
            ylabel(yvar,'Interpreter','none');
            title(sprintf('%s (%s) – %s', yvar, g, ds.name));
            grid on; box on
        end

        exportgraphics(fig, fullfile(outDir, sprintf('PERIOD_box_%s_%s.png', yvar, ds.name)), 'Resolution', 220);
        close(fig);
    end
end

% Save a human-readable “what days are excluded” note
fid = fopen(fullfile(outDir,'PERIOD_transition_days_note.txt'),'w');
fprintf(fid,'Transition/unreliable days excluded in NO_TRANSITION set: %s\n', mat2str(transitionDays));
fprintf(fid,'Active requested compare: Day5 vs (6-10) vs (11-13) vs (14-16) vs Day17\n');
fprintf(fid,'Passive requested compare: Day5 vs (11-13) vs (14-16) vs (17-18)\n');
fclose(fid);

end

function S = addPupilBaselineNormalization(S)
% baseline day3 pupil mean per mouse
S.pupil_base_day3 = nan(height(S),1);
S.pupil_norm_rel  = nan(height(S),1);

mk = unique(S.mouse_key,'stable');
for i=1:numel(mk)
    rM = (S.mouse_key==mk(i));
    base = mean(double(S.pupil_mean_px(rM & S.day_index==3)), 'omitnan');
    if ~isfinite(base) || base<=0
        continue;
    end
    S.pupil_base_day3(rM) = base;
    S.pupil_norm_rel(rM)  = (double(S.pupil_mean_px(rM)) - base) ./ base; % relative change
end
end

function plotPupilNormalizedDayTraces(S, outDir)
if ~ismember('pupil_norm_rel', S.Properties.VariableNames), return; end

fig = figure('Color','w','Position',[80 80 1100 650]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

for gi=1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);
    nexttile; hold on
    rG = (string(S.GroupMouse)==g);
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & (S.mouse_key==mk(i));
        [d,ord] = sort(S.day_index(r));
        y = double(S.pupil_norm_rel(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(S.day_index(rG)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.pupil_norm_rel(rG & S.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, 'k-', 'LineWidth', 3);
    yline(0,'k:');
    xlabel('Day'); ylabel('\Delta pupil / baseline(D3)');
    title(sprintf('Pupil normalized by Day3 baseline (%s)', g));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_daylevel_normalized_by_day3.png'), 'Resolution', 220);
close(fig);
end

function doEventLockedPupilAnalyses(T, timeVar, outDir, win_s, bout_gap_s, minSep_s, focusDaysA, focusDaysB)
% Event-locked pupil:
% 1) reward-locked using Injector_TTL rising edges
% 2) lick-bout-locked (bout start/end), split into rewarded vs non-reward bouts
% and the special comparisons for passive days 7-9 vs 12-13.

% Build per-session eventlocked averages, then aggregate across days/groups
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

% time base for extraction: assume approximately uniform; use interpolation
tAxis = win_s(1):0.1:win_s(2); % 10 Hz resampled for clean averaging

E = table; % event summaries
row = 0;

for gi=1:max(G)
    r = (G==gi);
    t = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    rw = double(T.Injector_TTL(r));
    pup = double(T.Diameter_px(r));

    if isempty(t) || all(isnan(t)) || all(isnan(pup)), continue; end
    [t,ord] = sort(t); lk=lk(ord); rw=rw(ord); pup=pup(ord);

    mk = string(T.mouse_key_norm(find(r,1,'first')));
    day = double(T.day_index(find(r,1,'first')));
    ses = double(T.session_idx(find(r,1,'first')));
    grp = string(T.GroupMouse(find(r,1,'first')));

    % ---------- reward events ----------
    rewardTimes = detectRisingEdges(t, rw);
    rewardTimes = enforceMinSeparation(rewardTimes, minSep_s);

    if ~isempty(rewardTimes)
        ER = extractEventLocked(t, pup, rewardTimes, tAxis);
        if ~isempty(ER)
            row=row+1;
            E.mouse_key(row,1) = mk;
            E.day_index(row,1) = day;
            E.session_idx(row,1) = ses;
            E.GroupMouse(row,1) = grp;
            E.eventType(row,1) = "reward";
            E.subType(row,1) = "all";
            E.trace{row,1} = ER;
        end
    end

    % ---------- lick bouts ----------
    lickTimes = detectRisingEdges(t, lk);
    lickTimes = enforceMinSeparation(lickTimes, 0.02); % keep all real licks
    if numel(lickTimes)>=1
        [boutStart, boutEnd] = makeLickBouts(lickTimes, bout_gap_s);

        % classify rewarded vs non-rewarded bout:
        % rewarded if any reward within [boutStart, boutEnd+1s]
        isRewarded = false(size(boutStart));
        for bi=1:numel(boutStart)
            isRewarded(bi) = any(rewardTimes>=boutStart(bi) & rewardTimes<=boutEnd(bi)+1.0);
        end

        % Start-locked
        EB_rs = extractEventLocked(t, pup, boutStart(isRewarded), tAxis);
        EB_ns = extractEventLocked(t, pup, boutStart(~isRewarded), tAxis);
        if ~isempty(EB_rs)
            row=row+1; E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
            E.GroupMouse(row,1)=grp; E.eventType(row,1)="lickBoutStart"; E.subType(row,1)="rewarded";
            E.trace{row,1}=EB_rs;
        end
        if ~isempty(EB_ns)
            row=row+1; E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
            E.GroupMouse(row,1)=grp; E.eventType(row,1)="lickBoutStart"; E.subType(row,1)="nonreward";
            E.trace{row,1}=EB_ns;
        end

        % End-locked
        EB_re = extractEventLocked(t, pup, boutEnd(isRewarded), tAxis);
        EB_ne = extractEventLocked(t, pup, boutEnd(~isRewarded), tAxis);
        if ~isempty(EB_re)
            row=row+1; E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
            E.GroupMouse(row,1)=grp; E.eventType(row,1)="lickBoutEnd"; E.subType(row,1)="rewarded";
            E.trace{row,1}=EB_re;
        end
        if ~isempty(EB_ne)
            row=row+1; E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
            E.GroupMouse(row,1)=grp; E.eventType(row,1)="lickBoutEnd"; E.subType(row,1)="nonreward";
            E.trace{row,1}=EB_ne;
        end
    end
end

if isempty(E)
    warning('No event-locked pupil traces could be computed.');
    return;
end

% Plot group-level means for key comparisons
% 1) Passive only: days 7-9 vs 12-13 (lickBoutStart, all bouts pooled by subtype)
plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "rewarded", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_rewarded_PASSIVE_days7-9_vs_12-13.png');
plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "nonreward", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_nonreward_PASSIVE_days7-9_vs_12-13.png');

% 2) Reward-locked: Active vs Passive, early(days3-5) vs late(days6-10)
plotRewardLockedEarlyLate(E, tAxis, outDir);

end

function times = enforceMinSeparation(times, minSep)
if isempty(times), return; end
times = sort(times(:));
keep = true(size(times));
last = -Inf;
for i=1:numel(times)
    if times(i)-last < minSep
        keep(i) = false;
    else
        last = times(i);
    end
end
times = times(keep);
end

function [boutStart, boutEnd] = makeLickBouts(lickTimes, gap_s)
lickTimes = sort(lickTimes(:));
if isempty(lickTimes)
    boutStart=[]; boutEnd=[]; return;
end
boutStart = lickTimes(1);
boutEnd   = lickTimes(1);
starts = [];
ends = [];
for i=2:numel(lickTimes)
    if lickTimes(i) - lickTimes(i-1) > gap_s
        starts(end+1,1) = boutStart; %#ok<AGROW>
        ends(end+1,1)   = boutEnd;   %#ok<AGROW>
        boutStart = lickTimes(i);
        boutEnd   = lickTimes(i);
    else
        boutEnd = lickTimes(i);
    end
end
starts(end+1,1)=boutStart;
ends(end+1,1)=boutEnd;
boutStart = starts; boutEnd = ends;
end

function trace = extractEventLocked(t, pup, eventTimes, tAxis)
% Returns baseline-subtracted average trace:
% baseline = mean of tAxis(t<0) window for each event
if isempty(eventTimes), trace=[]; return; end

% interpolate pupil to event-centered axis for each event, then average
M = nan(numel(eventTimes), numel(tAxis));
for i=1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    % baseline subtract using pre window (-2..0)
    pre = (tAxis >= min(tAxis) & tAxis < 0);
    b = mean(pi(pre), 'omitnan');
    pi = pi - b;
    M(i,:) = pi;
end
trace = mean(M, 1, 'omitnan');
end

function plotEventLockedComparison(E, tAxis, outDir, groupName, eventType, subType, daysA, daysB, fname)
rG = (string(E.GroupMouse)==groupName) & (string(E.eventType)==eventType) & (string(E.subType)==subType);

A = stackTraces(E.trace(rG & ismember(E.day_index, daysA)));
B = stackTraces(E.trace(rG & ismember(E.day_index, daysB)));

if isempty(A) || isempty(B)
    warning('Not enough data for %s %s %s comparison.', groupName, eventType, subType);
    return;
end

muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));

fig = figure('Color','w','Position',[80 80 900 520]); hold on
shaded(tAxis, muA, seA);
plot(tAxis, muA, 'LineWidth', 2.5);
shaded(tAxis, muB, seB);
plot(tAxis, muB, 'LineWidth', 2.5);
xline(0,'k-'); yline(0,'k:');
xlabel('Time from event (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('%s %s (%s): Days %s vs %s', groupName, eventType, subType, rangeStr(daysA), rangeStr(daysB)));
legend({sprintf('Days %s',rangeStr(daysA)), sprintf('Days %s',rangeStr(daysB))}, 'Location','best');
grid on; box on
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotRewardLockedEarlyLate(E, tAxis, outDir)
eventType = "reward";
subType   = "all";
earlyDays = 3:5;
lateDays  = 6:10;

fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi=1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);

    nexttile; hold on
    rG = (string(E.GroupMouse)==g) & (string(E.eventType)==eventType) & (string(E.subType)==subType);

    A = stackTraces(E.trace(rG & ismember(E.day_index, earlyDays)));
    B = stackTraces(E.trace(rG & ismember(E.day_index, lateDays)));

    if ~isempty(A)
        muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
        shaded(tAxis, muA, seA); plot(tAxis, muA, 'LineWidth', 2.5);
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shaded(tAxis, muB, seB); plot(tAxis, muB, 'LineWidth', 2.5);
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil (%s): Days %s vs %s', g, rangeStr(earlyDays), rangeStr(lateDays)));
    legend({sprintf('Days %s',rangeStr(earlyDays)), sprintf('Days %s',rangeStr(lateDays))}, 'Location','best');
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_reward_locked_early_vs_late_active_passive.png'), 'Resolution', 220);
close(fig);
end

function M = stackTraces(tracesCell)
if isempty(tracesCell), M=[]; return; end
n = numel(tracesCell{1});
M = nan(numel(tracesCell), n);
for i=1:numel(tracesCell)
    v = tracesCell{i};
    if isempty(v) || numel(v)~=n, continue; end
    M(i,:) = v(:)';
end
M = M(any(isfinite(M),2),:);
end

function shaded(x, mu, se)
x=x(:)'; mu=mu(:)'; se=se(:)';
fill([x fliplr(x)], [mu-se fliplr(mu+se)], 0.9*[1 1 1], 'EdgeColor','none', 'FaceAlpha',0.35);
end

function doOptionalPainAndBehaviorTests(T, outDir)
% Uses these ONLY if present in your final CSV.
optionalVars = {
    'Immersion_Latency_s'
    'HOT_Frames_Non_moving'
    'TST_Frames_Non_moving'
    % Straub tail: your pipeline merges a STRAUB summary, column name depends on that summary CSV.
    % We do NOT guess names; instead we only look for any column containing "STRAUB" that is numeric.
};

% plot Immersion latency by day active vs passive if present
if ismember('Immersion_Latency_s', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'Immersion_Latency_s', outDir, 'TEST_tailImmersion_latency.png');
end
if ismember('HOT_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'HOT_Frames_Non_moving', outDir, 'TEST_hotPlate_nonmoving.png');
end
if ismember('TST_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'TST_Frames_Non_moving', outDir, 'TEST_TST_nonmoving.png');
end

% Straub: strictly schema-safe: only if a numeric column exists that contains "straub"
vn = T.Properties.VariableNames;
isStraub = contains(lower(vn),'straub');
straubCols = vn(isStraub);
straubCols = straubCols(cellfun(@(c) isnumeric(T.(c)), straubCols));

if ~isempty(straubCols)
    % plot the first one (or you can edit to choose)
    plotScalarByDayGroup(T, straubCols{1}, outDir, 'TEST_straub_metric.png');
else
    % not an error; just informative
    fid = fopen(fullfile(outDir,'TEST_straub_note.txt'),'w');
    fprintf(fid,'No numeric Straub-related column found in ALL_mice_longitudinal.csv.\n');
    fprintf(fid,'If your STRAUB summary is merged, confirm the output column name and add it explicitly.\n');
    fclose(fid);
end
end

function plotScalarByDayGroup(T, varName, outDir, fname)
% Aggregate per mouse per day for scalar var; plot group means
assert(ismember(varName, T.Properties.VariableNames), 'Missing %s', varName);

% reduce to one value per mouse-day (mean)
G = findgroups(T.mouse_key_norm, T.day_index, T.GroupMouse);
mk = splitapply(@(x)x(1), T.mouse_key_norm, G);
dy = splitapply(@(x)x(1), T.day_index, G);
gp = splitapply(@(x)x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), T.(varName), G);

S = table(mk, dy, gp, val, 'VariableNames', {'mouse','day','group','value'});

fig = figure('Color','w','Position',[80 80 950 420]); hold on
groups = categories(S.group);
groups = groups(~strcmp(groups,'Unknown'));

for gi=1:numel(groups)
    g = groups{gi};
    rG = (S.group==g);
    days = unique(S.day(rG)); days=sort(days);
    mu = nan(size(days)); se = nan(size(days));
    for j=1:numel(days)
        x = S.value(rG & S.day==days(j));
        mu(j) = mean(x,'omitnan');
        se(j) = std(x,0,'omitnan')/sqrt(max(1,sum(isfinite(x))));
    end
    errorbar(days, mu, se, '-o', 'LineWidth', 2, 'MarkerSize',6);
end
xlabel('Day'); ylabel(varName,'Interpreter','none');
title(sprintf('%s: Active vs Passive (mouse-day means)', varName), 'Interpreter','none');
grid on; box on
legend(groups,'Location','best');
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end
