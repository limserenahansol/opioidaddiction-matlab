function make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v5()

% Uses YOUR longitudinal CSV schema (no guessing).
%
% Required columns:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1 numeric or logical), Injector_TTL (0/1), Diameter_px (numeric), RequirementLast (numeric)
%   Time axis: CamTime_rel_s OR PlotTime_s_30fps (either one must exist)
%
% Requested updates in v5:
%   - Fix legend mismatch in event-locked plots (line colors match legend)
%   - Fix Straub plotting explicitly using:
%       STRAUB_Frames_Non_moving, STRAUB_Pct_Non_moving
%   - Add 0911_red breakdown (which day/session has missing Lick_TTL and how much)
%   - Add bad-time check (time corruption vs TTL corruption)
%   - Fix categorical blank-category error
%   - Fix schema error printing (no mat2str on string arrays)
%   - Fix string/cell issues that caused brace-indexing / cell2mat failures

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

dt_bin_s = 60;          % 1-min bins
bout_gap_s = 2.0;       % lick bout definition
pupil_win = [-2 2];
min_event_separation_s = 0.5;

% Period definitions
PER.pre       = 3:5;
PER.during    = 6:10;
PER.post      = 11:13;
PER.withdraw  = 14:16;
PER.reexpo    = 17:18;

transitionDays = [4 6 11 14];

FOCUS_PASSIVE_DAYS_A = 7:9;
FOCUS_PASSIVE_DAYS_B = 12:13;

excludeDays_licking_global = [1 2];

% Target mouse for debugging
DEBUG_MOUSE = "0911_red";   % normalized key expected like "0911_red"

%% ===================== FIND LATEST run_* =====================
assert(exist(rootTry,'dir')==7, 'rootTry not found: %s', rootTry);

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

% Enforce mouse_key as string
assertHasVars(T, {'mouse_key','day_index','session_idx'}, 'Top-level');
T.mouse_key = string(T.mouse_key);

% Time axis: must exist (explicit allowed alternatives)
timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end

% Required core signals (schema-true)
REQ = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, [{'mouse_key','day_index','session_idx'}, REQ], 'Core longitudinal schema');

% Normalize mouse key and attach cohort metadata
COH = buildNewCohortTable();
T = addCohortMeta(T, COH);

% Drop rows with NaN day_index
T = T(~isnan(double(T.day_index)),:);

%% ===================== BAD TIME CHECK (before filtering) =====================
reportBadTime(T, timeVar, outDir, DEBUG_MOUSE);

%% ===================== REMOVE NON-FINITE TIME ROWS (but log them) =====================
tcol = double(T.(timeVar));
badTimeRows = ~isfinite(tcol);
if any(badTimeRows)
    fprintf('[WARN] Removing %d rows with non-finite %s\n', nnz(badTimeRows), timeVar);
    Tw = T(badTimeRows, {'mouse_key_norm','day_index','session_idx'});
    Tw.timeVar = repmat(string(timeVar), height(Tw),1);
    writetable(Tw, fullfile(outDir, sprintf('WARN_removed_nonfinite_%s_rows.csv', timeVar)));
    T = T(~badTimeRows,:);
end

%% ===================== QC: MISSINGNESS =====================
doMissingnessQC(T, outDir, {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'});

%% ===================== SESSION SUMMARY =====================
S = buildSessionSummary(T, timeVar);

%% ===================== OUTLIER QC =====================
doOutlierQC(S, outDir);

%% ===================== SPAGHETTI QC =====================
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'RequirementLast', outDir, 'QC_spaghetti_PRscore_RequirementLast.png');
plotPupilByGroupDay(S, outDir);

%% ===================== REQUESTED ANALYSES: LICKING =====================
plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global);
plotEarlyLateCumulative(T, timeVar, outDir, dt_bin_s, 3:5, 6:10);
plotDailyLickCountSpaghetti(S, outDir);
plotLickRasterExamples(T, timeVar, outDir);

%% ===================== PERIOD COMPARISONS =====================
DAYSETS.allDays.name  = 'ALLDAYS';
DAYSETS.allDays.mask  = true(height(S),1);

DAYSETS.noTransition.name = 'NO_TRANSITION';
DAYSETS.noTransition.mask = ~ismember(double(S.day_index), transitionDays);

doPeriodComparisons(S, outDir, PER, DAYSETS, transitionDays);

%% ===================== PUPIL =====================
S = addPupilBaselineNormalization(S);
plotPupilNormalizedDayTraces(S, outDir);

doEventLockedPupilAnalyses(T, timeVar, outDir, pupil_win, bout_gap_s, min_event_separation_s, ...
    FOCUS_PASSIVE_DAYS_A, FOCUS_PASSIVE_DAYS_B);

%% ===================== STRAUB + HOTPLATE + TST + TAIL IMMERSION =====================
doOptionalPainAndBehaviorTests_v5(T, outDir);

%% ===================== DEBUG: 0911_red missing lick TTL breakdown =====================
debugMissingLickTTL(T, timeVar, outDir, DEBUG_MOUSE);

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
    % robust printing (no mat2str on strings)
    missStr = strjoin(string(missing), ', ');
    error('Fix CSV schema or update required variable list. Missing: %s', missStr);
end
end

function COH = buildNewCohortTable()
rows = {
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
'6099','white',  'M', 'P',   '6099_pairB', 'Passive';
};
COH = cell2table(rows, 'VariableNames', {'cage','color','sex','group','pair_id','pair_role'});
COH.mouse_key = lower(string(COH.cage) + "_" + string(COH.color));
end

function T = addCohortMeta(T, COH)
mk = string(T.mouse_key);
mk = regexprep(mk,'\s+','');
mk = regexprep(mk,'-','_');

% if keys are like "6100red" -> "6100_red"
mk2 = mk;
isNoUnd = ~contains(mk2,'_');
mk2(isNoUnd) = regexprep(mk2(isNoUnd), '^(\d{4})([A-Za-z]+)$', '$1_$2');
mk2 = lower(mk2);

T.mouse_key_norm = mk2;
COH.mouse_key_norm = COH.mouse_key;

T = outerjoin(T, COH(:,{'mouse_key_norm','cage','color','sex','group','pair_id','pair_role'}), ...
    'Keys','mouse_key_norm', 'MergeKeys', true, 'Type','left');

g = string(T.group);
g(ismissing(g)) = "U";
T.GroupMouse = categorical(g, ["A","P","U"], {'Active','Passive','Unknown'});
end

function doMissingnessQC(T, outDir, vars)
mk = unique(string(T.mouse_key_norm),'stable');
days = unique(double(T.day_index));
days = sort(days(isfinite(days)));

M = nan(numel(mk), numel(days), numel(vars));
for i=1:numel(mk)
    for j=1:numel(days)
        r = (string(T.mouse_key_norm)==mk(i)) & (double(T.day_index)==days(j));
        if ~any(r), continue; end
        for v=1:numel(vars)
            x = T.(vars{v})(r);
            M(i,j,v) = missingFraction(x);
        end
    end
end

for v=1:numel(vars)
    fig = figure('Color','w','Position',[80 80 1100 520]);
    imagesc(days, 1:numel(mk), M(:,:,v));
    colormap(parula); colorbar;
    xlabel('Day'); ylabel('Mouse');
    title(['Missing fraction: ' vars{v}], 'Interpreter','none');
    set(gca,'YTick',1:numel(mk),'YTickLabel',mk, 'TickLabelInterpreter','none');
    exportgraphics(fig, fullfile(outDir, ['QC_missing_' vars{v} '.png']), 'Resolution', 200);
    close(fig);
end
end

function f = missingFraction(x)
% fraction of entries that are missing/NaN
if isnumeric(x) || islogical(x)
    f = mean(~isfinite(double(x)));
elseif isstring(x)
    f = mean(ismissing(x));
elseif iscell(x)
    xs = string(x);
    f = mean(ismissing(xs));
else
    % fallback
    try
        f = mean(~isfinite(double(x)));
    catch
        f = NaN;
    end
end
end

function S = buildSessionSummary(T, timeVar)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

S = table;
S.mouse_key = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
S.day_index = splitapply(@(x) double(x(1)), double(T.day_index), G);
S.session_idx = splitapply(@(x) double(x(1)), double(T.session_idx), G);

if ismember('GroupMouse', T.Properties.VariableNames)
    S.GroupMouse = splitapply(@(x) x(1), T.GroupMouse, G);
else
    S.GroupMouse = categorical(repmat("Unknown", height(S),1));
end

S.RequirementLast = splitapply(@(x) mean(double(x),'omitnan'), T.RequirementLast, G);
S.pupil_mean_px   = splitapply(@(x) mean(double(x),'omitnan'), T.Diameter_px, G);

S.nLicks = zeros(height(S),1);
S.licks_per_min = nan(height(S),1);
S.iei_median = nan(height(S),1);

for i=1:height(S)
    r = (string(T.mouse_key_norm)==S.mouse_key(i) & double(T.day_index)==S.day_index(i) & double(T.session_idx)==S.session_idx(i));
    t  = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));

    if isempty(t) || all(~isfinite(t)), continue; end

    ok = isfinite(t);
    t = t(ok); lk = lk(ok);

    [t,ord] = sort(t); lk = lk(ord);
    lickTimes = detectRisingEdges(t, lk);

    S.nLicks(i) = numel(lickTimes);
    dur = max(t) - min(t);
    if dur > 1
        S.licks_per_min(i) = (numel(lickTimes)/dur) * 60;
    end
    if numel(lickTimes) >= 2
        iei = diff(lickTimes);
        S.iei_median(i) = median(iei);
    end
end
end

function lickTimes = detectRisingEdges(t, x01)
x01 = x01(:);
x01(~isfinite(x01)) = 0;
x01 = x01 > 0.5;
dx = diff([false; x01]);
idx = find(dx==1);
lickTimes = t(idx);
lickTimes = lickTimes(isfinite(lickTimes));
end

function doOutlierQC(S, outDir)
metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

for k=1:numel(metrics)
    x = double(S.(metrics{k}));
    [isOut, z] = madOutliers(x, 4.5);

    fig = figure('Color','w','Position',[80 80 950 420]); hold on
    scatter(1:numel(x), x, 18, 'filled','MarkerFaceAlpha',0.5);
    scatter(find(isOut), x(isOut), 45, 'o','LineWidth',1.8);
    title(sprintf('QC outliers (MAD z>4.5): %s', metrics{k}), 'Interpreter','none');
    xlabel('Session index'); ylabel(metrics{k}, 'Interpreter','none');
    grid on; box on
    exportgraphics(fig, fullfile(outDir, ['QC_outliers_' metrics{k} '.png']), 'Resolution', 200);
    close(fig);

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
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

groups = categories(S.GroupMouse);
groups = groups(~strcmp(groups,'Unknown'));

for gi=1:numel(groups)
    nexttile; hold on
    g = groups{gi};
    rG = S.GroupMouse==g;
    mk = unique(string(S.mouse_key(rG)),'stable');
    for i=1:numel(mk)
        r = rG & (string(S.mouse_key)==mk(i));
        [d,ord] = sort(double(S.day_index(r)));
        y = double(S.(yvar)(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(double(S.day_index(rG))); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.(yvar)(rG & double(S.day_index)==dAll(j))), 'omitnan');
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
    mk = unique(string(S.mouse_key(rG)),'stable');
    for i=1:numel(mk)
        r = rG & (string(S.mouse_key)==mk(i));
        [d,ord] = sort(double(S.day_index(r)));
        y = double(S.pupil_mean_px(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(double(S.day_index(rG))); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.pupil_mean_px(rG & double(S.day_index)==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, 'k-', 'LineWidth',3);
    title(sprintf('Pupil mean (px): %s', g));
    xlabel('Day'); ylabel('Pupil diameter (px)'); grid on; box on
end
exportgraphics(fig, fullfile(outDir,'QC_pupil_mean_px_by_group.png'), 'Resolution', 220);
close(fig);
end

function plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays)
R0 = ~ismember(double(T.day_index), excludeDays);

for gname = ["Active","Passive"]
    Rg = R0 & (string(T.GroupMouse)==gname);
    if ~any(Rg), continue; end

    days = unique(double(T.day_index(Rg))); days = sort(days);
    fig = figure('Color','w','Position',[80 80 1250 520]); hold on

    for di=1:numel(days)
        d = days(di);
        R = Rg & (double(T.day_index)==d);

        [tgrid, muCum] = meanCumulativeLicksByDay(T(R,:), timeVar, dt_bin_s);
        if isempty(tgrid), continue; end
        plot(tgrid/60, muCum, 'LineWidth', 2, 'DisplayName', sprintf('Day %d', d));
    end
    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title(sprintf('Cumulative licking across session (%s mice)', gname));
    grid on; box on
    legend('show','Location','eastoutside');
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_cumulative_across_session_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function [tgrid, muCum] = meanCumulativeLicksByDay(Tday, timeVar, dt_bin_s)
if isempty(Tday), tgrid=[]; muCum=[]; return; end

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(Tday(:,keys));

sessionCurves = {};
tmaxAll = 0;
for gi=1:max(G)
    r = (G==gi);
    t = double(Tday.(timeVar)(r));
    lk = double(Tday.Lick_TTL(r));

    ok = isfinite(t);
    t = t(ok); lk = lk(ok);
    if isempty(t), continue; end

    [t,ord] = sort(t); lk=lk(ord);
    lickTimes = detectRisingEdges(t, lk);

    dur = max(t)-min(t);
    tmaxAll = max(tmaxAll, dur);
    if isempty(lickTimes)
        sessionCurves{end+1} = struct('dur', dur, 'lickTimes', []); %#ok<AGROW>
    else
        sessionCurves{end+1} = struct('dur', dur, 'lickTimes', lickTimes - min(t)); %#ok<AGROW>
    end
end

if isempty(sessionCurves)
    tgrid=[]; muCum=[]; return;
end

tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;
cumMat = nan(numel(sessionCurves), numel(tgrid));

for i=1:numel(sessionCurves)
    lt = sessionCurves{i}.lickTimes;
    if isempty(lt)
        cumMat(i,:) = 0;
    else
        cumMat(i,:) = arrayfun(@(tt) sum(lt <= tt), tgrid);
    end
end

muCum = mean(cumMat, 1, 'omitnan');
end

function plotEarlyLateCumulative(T, timeVar, outDir, dt_bin_s, earlyDays, lateDays)
for gname = ["Active","Passive"]
    Rg = (string(T.GroupMouse)==gname);

    fig = figure('Color','w','Position',[80 80 900 520]); hold on
    [tE, muE] = meanCumulativeLicksByDay(T(Rg & ismember(double(T.day_index), earlyDays),:), timeVar, dt_bin_s);
    [tL, muL] = meanCumulativeLicksByDay(T(Rg & ismember(double(T.day_index), lateDays),:), timeVar, dt_bin_s);

    if ~isempty(tE)
        plot(tE/60, muE, 'LineWidth', 3, 'DisplayName', sprintf('Days %s', rangeStr(earlyDays)));
    end
    if ~isempty(tL)
        plot(tL/60, muL, 'LineWidth', 3, 'DisplayName', sprintf('Days %s', rangeStr(lateDays)));
    end
    xlabel('Session time (min)'); ylabel('Mean cumulative licks');
    title(sprintf('Early vs Late cumulative (%s mice)', gname));
    grid on; box on
    legend('show','Location','southeast');
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_early_vs_late_cumulative_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1), s=sprintf('%d-%d', v(1), v(end));
else, s=strjoin(string(v),'-'); end
end

function plotDailyLickCountSpaghetti(S, outDir)
fig = figure('Color','w','Position',[80 80 900 520]); hold on
mk = unique(string(S.mouse_key),'stable');
for i=1:numel(mk)
    r = (string(S.mouse_key)==mk(i));
    [d,ord] = sort(double(S.day_index(r)));
    y = double(S.nLicks(r)); y=y(ord);
    plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
end
dAll = unique(double(S.day_index)); dAll=sort(dAll);
mu = nan(size(dAll));
for j=1:numel(dAll)
    mu(j) = mean(double(S.nLicks(double(S.day_index)==dAll(j))), 'omitnan');
end
plot(dAll, mu, 'k-', 'LineWidth', 3);
xlabel('Day index'); ylabel('Lick count per day'); title('Spaghetti: Lick count per day');
grid on; box on
exportgraphics(fig, fullfile(outDir,'LICK_spaghetti_lickcount_per_day.png'), 'Resolution', 220);
close(fig);
end

function plotLickRasterExamples(T, timeVar, outDir)
days = unique(double(T.day_index));
days = sort(days(isfinite(days)));
if isempty(days), return; end
pickDays = intersect(days, [4 9 12]);
if isempty(pickDays), pickDays = days(min(2,numel(days))); end

for pd = pickDays(:)'
    fig = figure('Color','w','Position',[80 80 1100 520]); hold on
    mk = unique(string(T.mouse_key_norm(double(T.day_index)==pd)),'stable');
    mk = mk(1:min(10,numel(mk)));

    y = 0;
    for i=1:numel(mk)
        r = (string(T.mouse_key_norm)==mk(i) & double(T.day_index)==pd);
        t = double(T.(timeVar)(r));
        lk = double(T.Lick_TTL(r));

        ok = isfinite(t);
        t = t(ok); lk = lk(ok);
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
S.period = strings(height(S),1);
S.period(ismember(double(S.day_index), PER.pre))      = "pre";
S.period(ismember(double(S.day_index), PER.during))   = "during";
S.period(ismember(double(S.day_index), PER.post))     = "post";
S.period(ismember(double(S.day_index), PER.withdraw)) = "withdrawal";
S.period(ismember(double(S.day_index), PER.reexpo))   = "reexposure";
S.period(S.period=="") = "<undef>";

% IMPORTANT FIX: no blank category names
S.period = categorical(S.period, ...
    ["pre","during","post","withdrawal","reexposure","<undef>"], ...
    ["pre","during","post","withdrawal","reexposure","<undef>"]);

metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

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

            Rp = Rg;
            if strcmp(yvar,'RequirementLast') && strcmp(g,'Passive')
                Rp = Rg & (S.period~="during");
            end

            if ~any(Rp)
                title(sprintf('%s (%s) – no data', yvar, g));
                axis off; continue;
            end

            [Gmp,~] = findgroups(string(S.mouse_key(Rp)), S.period(Rp));
            val = splitapply(@(x) mean(double(x),'omitnan'), S.(yvar)(Rp), Gmp);
            pr  = splitapply(@(x) x(1), S.period(Rp), Gmp);

            boxchart(double(pr), val);
            periodsOrder = categories(S.period);
            periodsOrder = periodsOrder(~strcmp(periodsOrder,'<undef>'));
            set(gca,'XTick',1:numel(periodsOrder),'XTickLabel',periodsOrder);
            ylabel(yvar,'Interpreter','none');
            title(sprintf('%s (%s) – %s', yvar, g, ds.name));
            grid on; box on
        end

        exportgraphics(fig, fullfile(outDir, sprintf('PERIOD_box_%s_%s.png', yvar, ds.name)), 'Resolution', 220);
        close(fig);
    end
end

fid = fopen(fullfile(outDir,'PERIOD_transition_days_note.txt'),'w');
fprintf(fid,'Transition/unreliable days excluded in NO_TRANSITION set: %s\n', mat2str(transitionDays));
fprintf(fid,'Active compare: Day5 vs (6-10) vs (11-13) vs (14-16) vs (17)\n');
fprintf(fid,'Passive compare: Day5 vs (11-13) vs (14-16) vs (17-18)\n');
fclose(fid);
end

function S = addPupilBaselineNormalization(S)
S.pupil_base_day3 = nan(height(S),1);
S.pupil_norm_rel  = nan(height(S),1);

mk = unique(string(S.mouse_key),'stable');
for i=1:numel(mk)
    rM = (string(S.mouse_key)==mk(i));
    base = mean(double(S.pupil_mean_px(rM & double(S.day_index)==3)), 'omitnan');
    if ~isfinite(base) || base<=0, continue; end
    S.pupil_base_day3(rM) = base;
    S.pupil_norm_rel(rM)  = (double(S.pupil_mean_px(rM)) - base) ./ base;
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
    mk = unique(string(S.mouse_key(rG)),'stable');
    for i=1:numel(mk)
        r = rG & (string(S.mouse_key)==mk(i));
        [d,ord] = sort(double(S.day_index(r)));
        y = double(S.pupil_norm_rel(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(double(S.day_index(rG))); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.pupil_norm_rel(rG & double(S.day_index)==dAll(j))), 'omitnan');
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
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

tAxis = win_s(1):0.1:win_s(2);

E = table; row = 0;

for gi=1:max(G)
    r = (G==gi);
    t = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    rw = double(T.Injector_TTL(r));
    pup = double(T.Diameter_px(r));

    ok = isfinite(t) & isfinite(pup);
    t=t(ok); lk=lk(ok); rw=rw(ok); pup=pup(ok);
    if isempty(t), continue; end

    [t,ord] = sort(t); lk=lk(ord); rw=rw(ord); pup=pup(ord);

    mk  = string(T.mouse_key_norm(find(r,1,'first')));
    day = double(T.day_index(find(r,1,'first')));
    ses = double(T.session_idx(find(r,1,'first')));
    grp = string(T.GroupMouse(find(r,1,'first')));

    rewardTimes = enforceMinSeparation(detectRisingEdges(t, rw), minSep_s);

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

    lickTimes = enforceMinSeparation(detectRisingEdges(t, lk), 0.02);
    if ~isempty(lickTimes)
        [boutStart, boutEnd] = makeLickBouts(lickTimes, bout_gap_s);

        isRewarded = false(size(boutStart));
        for bi=1:numel(boutStart)
            isRewarded(bi) = any(rewardTimes>=boutStart(bi) & rewardTimes<=boutEnd(bi)+1.0);
        end

        addTrace("lickBoutStart","rewarded",  boutStart(isRewarded));
        addTrace("lickBoutStart","nonreward", boutStart(~isRewarded));
        addTrace("lickBoutEnd","rewarded",    boutEnd(isRewarded));
        addTrace("lickBoutEnd","nonreward",   boutEnd(~isRewarded));
    end

    function addTrace(ev, sub, times)
        tr = extractEventLocked(t, pup, times, tAxis);
        if isempty(tr), return; end
        row=row+1;
        E.mouse_key(row,1)=mk;
        E.day_index(row,1)=day;
        E.session_idx(row,1)=ses;
        E.GroupMouse(row,1)=grp;
        E.eventType(row,1)=string(ev);
        E.subType(row,1)=string(sub);
        E.trace{row,1}=tr;
    end
end

if isempty(E)
    warning('No event-locked pupil traces could be computed.');
    return;
end

plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "rewarded", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_rewarded_PASSIVE_days7-9_vs_12-13.png');
plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "nonreward", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_nonreward_PASSIVE_days7-9_vs_12-13.png');

plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir);
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
starts = lickTimes(1);
ends   = lickTimes(1);
for i=2:numel(lickTimes)
    if lickTimes(i) - lickTimes(i-1) > gap_s
        starts(end+1,1) = lickTimes(i); %#ok<AGROW>
        ends(end+1,1)   = lickTimes(i); %#ok<AGROW>
    else
        ends(end) = lickTimes(i);
    end
end
boutStart = starts;
boutEnd   = ends;
end

function trace = extractEventLocked(t, pup, eventTimes, tAxis)
if isempty(eventTimes), trace=[]; return; end
M = nan(numel(eventTimes), numel(tAxis));
pre = (tAxis < 0);

for i=1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    b = mean(pi(pre), 'omitnan');
    pi = pi - b;
    M(i,:) = pi;
end
trace = mean(M, 1, 'omitnan');
end

function plotEventLockedComparison(E, tAxis, outDir, groupName, eventType, subType, daysA, daysB, fname)
rG = (string(E.GroupMouse)==groupName) & (string(E.eventType)==eventType) & (string(E.subType)==subType);

A = stackTraces(E.trace(rG & ismember(double(E.day_index), daysA)));
B = stackTraces(E.trace(rG & ismember(double(E.day_index), daysB)));

if isempty(A) || isempty(B)
    warning('Not enough data for %s %s %s comparison.', groupName, eventType, subType);
    return;
end

muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));

% fixed colors + legend binding to line handles
cA = [0.8500 0.3250 0.0980]; % orange-ish
cB = [0.0000 0.4470 0.7410]; % blue

fig = figure('Color','w','Position',[80 80 900 520]); hold on
shaded(tAxis, muA, seA);
hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', cA, 'DisplayName', sprintf('Days %s',rangeStr(daysA)));
shaded(tAxis, muB, seB);
hB = plot(tAxis, muB, 'LineWidth', 2.5, 'Color', cB, 'DisplayName', sprintf('Days %s',rangeStr(daysB)));

xline(0,'k-'); yline(0,'k:');
xlabel('Time from event (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('%s %s (%s)', groupName, eventType, subType));
legend([hA hB],'Location','best');
grid on; box on
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir)
eventType = "reward";
subType   = "all";
earlyDays = 3:5;
lateDays  = 6:10;

% fixed colors (match your uploaded figure expectation)
cEarly = [0.8500 0.3250 0.0980]; % orange
cLate  = [0.0000 0.4470 0.7410]; % blue

fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi=1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);

    nexttile; hold on
    rG = (string(E.GroupMouse)==g) & (string(E.eventType)==eventType) & (string(E.subType)==subType);

    A = stackTraces(E.trace(rG & ismember(double(E.day_index), earlyDays)));
    B = stackTraces(E.trace(rG & ismember(double(E.day_index), lateDays)));

    hA=[]; hB=[];
    if ~isempty(A)
        muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
        shaded(tAxis, muA, seA);
        hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', cEarly, 'DisplayName', sprintf('Days %s',rangeStr(earlyDays)));
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shaded(tAxis, muB, seB);
        hB = plot(tAxis, muB, 'LineWidth', 2.5, 'Color', cLate, 'DisplayName', sprintf('Days %s',rangeStr(lateDays)));
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil (%s)', g));
    if ~isempty(hA) && ~isempty(hB)
        legend([hA hB],'Location','best');
    elseif ~isempty(hA)
        legend(hA,'Location','best');
    elseif ~isempty(hB)
        legend(hB,'Location','best');
    end
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

function doOptionalPainAndBehaviorTests_v5(T, outDir)
% Explicitly include Straub metrics if present
% Tail immersion / hotplate / TST kept same pattern

if ismember('Immersion_Latency_s', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'Immersion_Latency_s', outDir, 'TEST_tailImmersion_latency.png');
end
if ismember('HOT_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'HOT_Frames_Non_moving', outDir, 'TEST_hotPlate_nonmoving_frames.png');
end
if ismember('TST_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'TST_Frames_Non_moving', outDir, 'TEST_TST_nonmoving_frames.png');
end

% FIX: explicitly plot the two Straub metrics you requested
hasStraubFrames = ismember('STRAUB_Frames_Non_moving', T.Properties.VariableNames);
hasStraubPct    = ismember('STRAUB_Pct_Non_moving',    T.Properties.VariableNames);

if hasStraubFrames
    plotScalarByDayGroup(T, 'STRAUB_Frames_Non_moving', outDir, 'TEST_STRAUB_nonmoving_frames.png');
end
if hasStraubPct
    plotScalarByDayGroup(T, 'STRAUB_Pct_Non_moving', outDir, 'TEST_STRAUB_nonmoving_pct.png');
end

if ~hasStraubFrames && ~hasStraubPct
    vn = string(T.Properties.VariableNames);
    straubLike = vn(contains(lower(vn),'straub'));
    fid = fopen(fullfile(outDir,'TEST_straub_note.txt'),'w');
    fprintf(fid,'No STRAUB_Frames_Non_moving or STRAUB_Pct_Non_moving in ALL_mice_longitudinal.csv.\n');
    fprintf(fid,'STRAUB* columns present (if any):\n');
    for i=1:numel(straubLike)
        fprintf(fid,'  %s\n', straubLike(i));
    end
    fclose(fid);
end
end

function plotScalarByDayGroup(T, varName, outDir, fname)
assert(ismember(varName, T.Properties.VariableNames), 'Missing %s', varName);

G = findgroups(string(T.mouse_key_norm), double(T.day_index), T.GroupMouse);
mk = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
dy = splitapply(@(x) double(x(1)), double(T.day_index), G);
gp = splitapply(@(x) x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), double(T.(varName)), G);

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
    errorbar(days, mu, se, '-o', 'LineWidth', 2, 'MarkerSize',6, 'DisplayName', g);
end
xlabel('Day'); ylabel(varName,'Interpreter','none');
title(sprintf('%s: Active vs Passive (mouse-day means)', varName), 'Interpreter','none');
grid on; box on
legend('show','Location','best');
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function debugMissingLickTTL(T, timeVar, outDir, mouseKeyNorm)
mk = lower(string(mouseKeyNorm));
rM = (lower(string(T.mouse_key_norm)) == mk);
if ~any(rM)
    warning('debugMissingLickTTL: mouse %s not found in mouse_key_norm.', mk);
    return;
end

t = T(rM,:);
t.day_index = double(t.day_index);
t.session_idx = double(t.session_idx);

G = findgroups(t.day_index, t.session_idx);
day = splitapply(@(x) x(1), t.day_index, G);
ses = splitapply(@(x) x(1), t.session_idx, G);

lick = t.Lick_TTL;
time = t.(timeVar);

nRows      = splitapply(@numel, lick, G);
nLickNaN   = splitapply(@(x) nnz(~isfinite(double(x))), lick, G);
fracLickNaN= nLickNaN ./ max(1,nRows);

nTimeBad   = splitapply(@(x) nnz(~isfinite(double(x))), time, G);
fracTimeBad= nTimeBad ./ max(1,nRows);

% monotonicity check (only on finite time)
nonMono = splitapply(@(x) any(diff(double(x(isfinite(double(x)))))<0), time, G);

R = table(day, ses, nRows, nLickNaN, fracLickNaN, nTimeBad, fracTimeBad, nonMono);
R = sortrows(R, {'fracLickNaN','fracTimeBad','nonMono'},{'descend','descend','descend'});

writetable(R, fullfile(outDir, sprintf('DEBUG_%s_missing_LickTTL_by_day_session.csv', mk)));

% Print top suspicious to console
fprintf('\n[DEBUG] %s: top day/session with missing Lick_TTL or bad time\n', mk);
disp(R(1:min(20,height(R)),:));

% Also save a quick plot of missing fraction across days
fig = figure('Color','w','Position',[80 80 900 350]); hold on
uDays = unique(R.day); uDays=sort(uDays);
y = nan(size(uDays));
for i=1:numel(uDays)
    y(i) = max(R.fracLickNaN(R.day==uDays(i)));
end
plot(uDays, y, '-o','LineWidth',2);
xlabel('Day'); ylabel('Max frac NaN in Lick\_TTL (across sessions)');
title(sprintf('DEBUG %s: Missing Lick_TTL by day (max across sessions)', mk), 'Interpreter','none');
grid on; box on
exportgraphics(fig, fullfile(outDir, sprintf('DEBUG_%s_missing_LickTTL_by_day.png', mk)), 'Resolution', 220);
close(fig);
end

function reportBadTime(T, timeVar, outDir, debugMouse)
% Reports sessions with:
%   - non-finite time
%   - non-monotonic time
%   - huge time jumps (optional heuristic)

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

mk  = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
day = splitapply(@(x) double(x(1)), double(T.day_index), G);
ses = splitapply(@(x) double(x(1)), double(T.session_idx), G);

tt  = splitapply(@(x) {double(x)}, double(T.(timeVar)), G);

nRows    = cellfun(@numel, tt);
nBad     = cellfun(@(x) nnz(~isfinite(x)), tt);

nonMono  = cellfun(@(x) any(diff(x(isfinite(x)))<0), tt);

% jump heuristic: 99th percentile of dt
p99dt = cellfun(@(x) prctile(diff(x(isfinite(x))), 99), tt);
p99dt(~isfinite(p99dt)) = NaN;

R = table(mk, day, ses, nRows, nBad, nBad./max(1,nRows), nonMono, p99dt, ...
    'VariableNames', {'mouse','day','session','nRows','nNonFinite','fracNonFinite','nonMonotonic','p99_dt'});

% Save full report
writetable(R, fullfile(outDir, sprintf('QC_badTime_report_%s.csv', timeVar)));

% Focused report for debug mouse
dm = lower(string(debugMouse));
rD = lower(R.mouse) == dm;
if any(rD)
    writetable(R(rD,:), fullfile(outDir, sprintf('QC_badTime_DEBUG_%s_%s.csv', dm, timeVar)));
end
end
