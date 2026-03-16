function make_longitudinal_QC_and_requested_analysesv4()
% Uses YOUR longitudinal CSV schema (no guessing beyond explicit alternatives):
% Required columns:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1), Injector_TTL (0/1), Diameter_px (numeric), RequirementLast (numeric)
% Time axis (explicit alternatives only):
%   CamTime_rel_s OR PlotTime_s_30fps  (either one must exist)
%
% Optional tests (ONLY if these exact columns exist):
%   Immersion_Latency_s
%   HOT_Frames_Non_moving
%   TST_Frames_Non_moving
%   STRAUB_Frames_Non_moving
%   STRAUB_Pct_Non_moving
%
% Adds:
%   - Legend fix for event-locked pupil plots (shading no longer hijacks legend)
%   - Categorical fix (no blank categories)
%   - Robust assertHasVars (no mat2str crash)
%   - Robust string handling for mouse_key comparisons
%   - 0911_red missing lick breakdown report
%   - badT time-axis corruption check (distinguish time-broken vs TTL-broken)

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Session-time binning for cumulative lick curves
dt_bin_s = 60;                % 1-min bins
bout_gap_s = 2.0;             % lick bout: pause > 2s splits bouts
pupil_win = [-2 2];           % seconds around event
min_event_separation_s = 0.5; % avoid degenerate overlaps

% Period definitions (days)
PER.pre       = 3:5;
PER.during    = 6:10;
PER.post      = 11:13;
PER.withdraw  = 14:16;
PER.reexpo    = 17:18;

% Transition/unreliable days
transitionDays = [4 6 11 14];

% Passive-specific comparisons you requested
FOCUS_PASSIVE_DAYS_A = 7:9;
FOCUS_PASSIVE_DAYS_B = 12:13;

% Exclude habituation days for some licking plots
excludeDays_licking_global = [1 2];

% Target debug mouse
DEBUG_MOUSE = "0911_red"; % your request
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
assertHasVars(T, {'mouse_key'}, 'Top-level');
T.mouse_key = toStringCol(T.mouse_key);

% required columns (schema-true)
REQ = {'mouse_key','day_index','session_idx','Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, REQ, 'Core longitudinal schema');

% enforce core numeric columns robustly
T.day_index        = toDoubleCol(T.day_index);
T.session_idx      = toDoubleCol(T.session_idx);
T.Lick_TTL         = toDoubleCol(T.Lick_TTL);
T.Injector_TTL     = toDoubleCol(T.Injector_TTL);
T.Diameter_px      = toDoubleCol(T.Diameter_px);
T.RequirementLast  = toDoubleCol(T.RequirementLast);

% time axis (explicitly allowed alternatives)
timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end
T.(timeVar) = toDoubleCol(T.(timeVar));

% cohort metadata (YOUR new cohort)
COH = buildNewCohortTable();

% join cohort metadata onto rows by mouse_key
T = addCohortMeta(T, COH);

% drop days with NaN day_index
T = T(isfinite(T.day_index),:);

%% ===================== TIME-AXIS INTEGRITY CHECK (badT) =====================
badT = checkTimeAxisIntegrity(T, timeVar);
if ~isempty(badT)
    writetable(badT, fullfile(outDir, 'QC_badTimeAxis_sessions.csv'));
    fprintf('[WARN] Found %d sessions with time-axis integrity problems. Saved QC_badTimeAxis_sessions.csv\n', height(badT));
end

%% ===================== QC: BASIC SANITY =====================
doMissingnessQC(T, outDir);

% 0911_red breakdown you requested (lick TTL missing vs time corruption)
debugMissingLickForMouse(T, timeVar, outDir, DEBUG_MOUSE);

% Session-level summaries
S = buildSessionSummary(T, timeVar);

doOutlierQC(S, outDir);

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
DAYSETS.noTransition.mask = ~ismember(S.day_index, transitionDays);

doPeriodComparisons(S, outDir, PER, DAYSETS, transitionDays);

%% ===================== PUPIL =====================
S = addPupilBaselineNormalization(S);
plotPupilNormalizedDayTraces(S, outDir);

doEventLockedPupilAnalyses(T, timeVar, outDir, pupil_win, bout_gap_s, min_event_separation_s, ...
    FOCUS_PASSIVE_DAYS_A, FOCUS_PASSIVE_DAYS_B);

%% ===================== TAIL IMMERSION / HOT / TST / STRAUB =====================
doOptionalPainAndBehaviorTests_explicit(T, outDir);

fprintf('\nDONE. Outputs saved in:\n  %s\n', outDir);
end

%% =====================================================================
%% =========================== HELPERS =================================
%% =====================================================================

function assertHasVars(T, vars, context)
missing = vars(~ismember(vars, T.Properties.VariableNames));
if ~isempty(missing)
    fprintf('\nSCHEMA ERROR in %s.\nMissing required columns:\n', context);
    for i=1:numel(missing)
        fprintf('  - %s\n', missing{i});
    end
    error('Fix CSV schema or update required variable list in the script.');
end
end

function s = toStringCol(x)
% robust string conversion
if isstring(x)
    s = x;
elseif iscellstr(x) || iscell(x)
    s = string(x);
elseif ischar(x)
    s = string(cellstr(x));
else
    s = string(x);
end
s = strip(s);
end

function d = toDoubleCol(x)
% robust numeric conversion (logical, numeric, string, cell)
if isnumeric(x)
    d = double(x);
elseif islogical(x)
    d = double(x);
elseif isstring(x) || ischar(x)
    d = str2double(string(x));
elseif iscell(x)
    % cell can contain numeric/logical/char
    try
        d = nan(size(x));
        for i=1:numel(x)
            xi = x{i};
            if isnumeric(xi) || islogical(xi)
                d(i) = double(xi);
            else
                d(i) = str2double(string(xi));
            end
        end
    catch
        d = str2double(string(x));
    end
else
    d = str2double(string(x));
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
COH.mouse_key = string(COH.cage) + "_" + string(COH.color);
end

function T = addCohortMeta(T, COH)
mk = string(T.mouse_key);
mk = regexprep(mk,'\s+','');
mk = regexprep(mk,'-','_');

mk2 = mk;
isNoUnd = ~contains(mk2,'_');
mk2(isNoUnd) = regexprep(mk2(isNoUnd), '^(\d{4})([A-Za-z]+)$', '$1_$2');
mk2 = lower(mk2);

cohKey = lower(string(COH.mouse_key));

T.mouse_key_norm = mk2;
COH.mouse_key_norm = cohKey;

T = outerjoin(T, COH(:,{'mouse_key_norm','cage','color','sex','group','pair_id','pair_role'}), ...
    'Keys','mouse_key_norm', 'MergeKeys', true, 'Type','left');

g = string(T.group);
g(ismissing(g)) = "U";
T.GroupMouse = categorical(g, ["A","P","U"], {'Active','Passive','Unknown'});
end

function badT = checkTimeAxisIntegrity(T, timeVar)
% returns one row per bad session (mouse/day/session) if time axis is corrupted
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

mk = splitapply(@(x) x(1), T.mouse_key_norm, G);
dy = splitapply(@(x) x(1), T.day_index, G);
ss = splitapply(@(x) x(1), T.session_idx, G);

bad_nonfinite = false(max(G),1);
bad_nonmono   = false(max(G),1);
bad_dur       = false(max(G),1);
dur_s         = nan(max(G),1);
n_rows        = nan(max(G),1);

for i=1:max(G)
    r = (G==i);
    t = double(T.(timeVar)(r));
    n_rows(i) = numel(t);

    bad_nonfinite(i) = any(~isfinite(t));
    t2 = t(isfinite(t));
    if numel(t2) >= 2
        bad_nonmono(i) = any(diff(t2) < 0);
        dur_s(i) = max(t2) - min(t2);
        bad_dur(i) = ~(dur_s(i) > 0);
    else
        bad_nonmono(i) = true;
        dur_s(i) = NaN;
        bad_dur(i) = true;
    end
end

keep = bad_nonfinite | bad_nonmono | bad_dur;
badT = table(string(mk(keep)), dy(keep), ss(keep), n_rows(keep), dur_s(keep), ...
    bad_nonfinite(keep), bad_nonmono(keep), bad_dur(keep), ...
    'VariableNames', {'mouse_key','day_index','session_idx','n_rows','duration_s', ...
                      'has_nonfinite_time','has_nonmonotonic_time','bad_duration'});
end

function doMissingnessQC(T, outDir)
vars = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
vars = vars(ismember(vars, T.Properties.VariableNames));

mk = unique(T.mouse_key_norm,'stable');
days = unique(T.day_index);
days = sort(days(isfinite(days)));

M = nan(numel(mk), numel(days), numel(vars));
for i=1:numel(mk)
    for j=1:numel(days)
        r = (T.mouse_key_norm==mk(i) & T.day_index==days(j));
        if ~any(r), continue; end
        for v=1:numel(vars)
            x = T.(vars{v})(r);
            M(i,j,v) = mean(~isfinite(x)); % treat NaN/Inf as missing
        end
    end
end

for v=1:numel(vars)
    fig = figure('Color','w','Position',[80 80 1100 520]);
    imagesc(days, 1:numel(mk), M(:,:,v));
    colormap(parula); colorbar;
    xlabel('Day'); ylabel('Mouse'); title(['Missing fraction: ' vars{v}], 'Interpreter','none');
    set(gca,'YTick',1:numel(mk),'YTickLabel',mk, 'TickLabelInterpreter','none');
    exportgraphics(fig, fullfile(outDir, ['QC_missing_' vars{v} '.png']), 'Resolution', 200);
    close(fig);
end
end

function debugMissingLickForMouse(T, timeVar, outDir, mouseKey)
% Prints which day/session has missing Lick_TTL and whether time is corrupted.
mouseKey = lower(string(mouseKey));
mouseKey = regexprep(mouseKey,'-','_');
if ~contains(mouseKey,'_')
    mouseKey = regexprep(mouseKey, '^(\d{4})([A-Za-z]+)$', '$1_$2');
end

rM = (T.mouse_key_norm == mouseKey);
if ~any(rM)
    fprintf('[DEBUG] Mouse %s not found in mouse_key_norm.\n', mouseKey);
    return;
end

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(rM,keys));

dy = splitapply(@(x) x(1), T.day_index(rM), G);
ss = splitapply(@(x) x(1), T.session_idx(rM), G);

fracLickMissing = nan(max(G),1);
fracTimeBad     = nan(max(G),1);
nrows           = nan(max(G),1);

for i=1:max(G)
    rr = (G==i);
    lk = T.Lick_TTL(rM);
    lk = lk(rr);
    t  = T.(timeVar)(rM);
    t  = t(rr);

    nrows(i) = numel(lk);
    fracLickMissing(i) = mean(~isfinite(lk));
    fracTimeBad(i)     = mean(~isfinite(t));
end

DBG = table(repmat(mouseKey,max(G),1), dy, ss, nrows, fracLickMissing, fracTimeBad, ...
    'VariableNames', {'mouse_key','day_index','session_idx','n_rows','frac_LickTTL_missing','frac_time_missing'});
DBG = sortrows(DBG, {'day_index','session_idx'});

writetable(DBG, fullfile(outDir, sprintf('DEBUG_missingLick_%s.csv', mouseKey)));
fprintf('[DEBUG] Saved missing-lick breakdown for %s to DEBUG_missingLick_%s.csv\n', mouseKey, mouseKey);

% Print the worst offenders to console
DBG2 = DBG(DGBool(DGBool(fracLickMissing),:),:); %#ok<NASGU>
[~,ix] = sort(DBG.frac_LickTTL_missing, 'descend');
topN = min(15, height(DBG));
fprintf('\n[DEBUG] Top %d sessions by Lick_TTL missing (mouse %s):\n', topN, mouseKey);
disp(DBG(ix(1:topN),:));
end

function tf = DGBool(x)
tf = isfinite(x) & x>0;
end

function S = buildSessionSummary(T, timeVar)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

S = table;
S.mouse_key   = string(splitapply(@(x)x(1), T.mouse_key_norm, G));
S.day_index   = splitapply(@(x)x(1), T.day_index, G);
S.session_idx = splitapply(@(x)x(1), T.session_idx, G);

if ismember('GroupMouse', T.Properties.VariableNames)
    S.GroupMouse = splitapply(@(x)x(1), T.GroupMouse, G);
else
    S.GroupMouse = categorical(repmat("Unknown", height(S),1));
end

S.RequirementLast = splitapply(@nanmean, T.RequirementLast, G);
S.pupil_mean_px   = splitapply(@nanmean, T.Diameter_px, G);

S.nLicks = zeros(height(S),1);
S.licks_per_min = nan(height(S),1);
S.iei_median = nan(height(S),1);

for i=1:height(S)
    r = (T.mouse_key_norm==S.mouse_key(i) & T.day_index==S.day_index(i) & T.session_idx==S.session_idx(i));
    t = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    if isempty(t) || all(~isfinite(t)), continue; end

    ok = isfinite(t);
    t = t(ok); lk = lk(ok);

    [t,ord] = sort(t); lk = lk(ord);

    lickTimes = detectRisingEdges(t, lk);
    S.nLicks(i) = numel(lickTimes);

    dur = max(t)-min(t);
    if isfinite(dur) && dur > 1
        S.licks_per_min(i) = (numel(lickTimes) / dur) * 60;
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
    xlabel('Session index'); ylabel(metrics{k}, 'Interpreter','none'); grid on; box on
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
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & S.mouse_key==mk(i);
        [d,ord] = sort(S.day_index(r));
        y = double(S.(yvar)(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
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
        r = rG & (S.mouse_key==mk(i));
        [d,ord] = sort(S.day_index(r));
        y = double(S.pupil_mean_px(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end
    dAll = unique(S.day_index(rG)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(S.pupil_mean_px(rG & S.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, 'k-', 'LineWidth', 3);
    title(sprintf('Pupil mean (px): %s', g));
    xlabel('Day'); ylabel('Pupil diameter (px)'); grid on; box on
end
exportgraphics(fig, fullfile(outDir,'QC_pupil_mean_px_by_group.png'), 'Resolution', 220);
close(fig);
end

function plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays)
R0 = ~ismember(T.day_index, excludeDays);

for gname = ["Active","Passive"]
    Rg = R0 & (string(T.GroupMouse)==gname);
    if ~any(Rg), continue; end

    days = unique(T.day_index(Rg)); days = sort(days);
    fig = figure('Color','w','Position',[80 80 1250 520]); hold on

    labels = strings(0);
    for di=1:numel(days)
        d = days(di);
        R = Rg & T.day_index==d;
        [tgrid, muCum] = meanCumulativeLicksByDay(T(R,:), timeVar, dt_bin_s);
        if isempty(tgrid), continue; end
        plot(tgrid/60, muCum, 'LineWidth', 2);
        labels(end+1) = "Day " + string(d); %#ok<AGROW>
    end

    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title(sprintf('Cumulative licking across session (%s mice)', gname));
    grid on; box on
    if ~isempty(labels), legend(labels, 'Location','eastoutside'); end
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

    dur = max(t) - min(t);
    if ~isfinite(dur) || dur<=0, continue; end
    tmaxAll = max(tmaxAll, dur);

    sessionCurves{end+1} = struct('dur', dur, 'lickTimes', lickTimes - min(t)); %#ok<AGROW>
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

    [tE, muE] = meanCumulativeLicksByDay(T(Rg & ismember(T.day_index, earlyDays),:), timeVar, dt_bin_s);
    [tL, muL] = meanCumulativeLicksByDay(T(Rg & ismember(T.day_index, lateDays),:), timeVar, dt_bin_s);

    h = gobjects(0);
    labels = strings(0);

    if ~isempty(tE)
        h(end+1) = plot(tE/60, muE, 'LineWidth', 3); %#ok<AGROW>
        labels(end+1) = "Days " + rangeStr(earlyDays); %#ok<AGROW>
    end
    if ~isempty(tL)
        h(end+1) = plot(tL/60, muL, 'LineWidth', 3); %#ok<AGROW>
        labels(end+1) = "Days " + rangeStr(lateDays); %#ok<AGROW>
    end

    xlabel('Session time (min)'); ylabel('Mean cumulative licks');
    title(sprintf('Early (days %s) vs Late (days %s) sessions (%s mice)', ...
        rangeStr(earlyDays), rangeStr(lateDays), gname));
    grid on; box on
    if ~isempty(h), legend(h, labels, 'Location','southeast'); end
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_early_vs_late_cumulative_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1), s=sprintf('%d-%d', v(1), v(end));
else, s=strjoin(string(v),','); end
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

function plotLickRasterExamples(T, timeVar, outDir)
days = unique(T.day_index);
days = sort(days(isfinite(days)));
if isempty(days), return; end
pickDays = intersect(days, [4 9 12]);
if isempty(pickDays), pickDays = days(min(2,numel(days))); end

for pd = pickDays(:)'
    fig = figure('Color','w','Position',[80 80 1100 520]); hold on
    mk = unique(T.mouse_key_norm(T.day_index==pd),'stable');
    mk = mk(1:min(10,numel(mk)));

    y = 0;
    for i=1:numel(mk)
        r = (T.mouse_key_norm==mk(i) & T.day_index==pd);
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
S.period = repmat("<undef>", height(S), 1);
S.period(ismember(S.day_index, PER.pre))      = "pre";
S.period(ismember(S.day_index, PER.during))   = "during";
S.period(ismember(S.day_index, PER.post))     = "post";
S.period(ismember(S.day_index, PER.withdraw)) = "withdrawal";
S.period(ismember(S.day_index, PER.reexpo))   = "reexposure";

% FIX: categorical must not contain blank category names
S.period = categorical(S.period, ["pre","during","post","withdrawal","reexposure","<undef>"], ...
                                  {"pre","during","post","withdrawal","reexposure","<undef>"});

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

            [Gmp,~] = findgroups(S.mouse_key(Rp), S.period(Rp));
            val = splitapply(@(x) mean(double(x),'omitnan'), S.(yvar)(Rp), Gmp);
            pr  = splitapply(@(x) x(1), S.period(Rp), Gmp);
            P = table(pr, val, 'VariableNames', {'period','value'});

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

fid = fopen(fullfile(outDir,'PERIOD_transition_days_note.txt'),'w');
fprintf(fid,'Transition/unreliable days excluded in NO_TRANSITION set: %s\n', mat2str(transitionDays));
fprintf(fid,'Active requested compare: Day5 vs (6-10) vs (11-13) vs (14-16) vs Day17\n');
fprintf(fid,'Passive requested compare: Day5 vs (11-13) vs (14-16) vs (17-18)\n');
fclose(fid);
end

function S = addPupilBaselineNormalization(S)
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
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));
tAxis = win_s(1):0.1:win_s(2);

E = table; row = 0;

for gi=1:max(G)
    r = (G==gi);
    t = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    rw = double(T.Injector_TTL(r));
    pup= double(T.Diameter_px(r));

    ok = isfinite(t) & isfinite(pup);
    t=t(ok); lk=lk(ok); rw=rw(ok); pup=pup(ok);
    if numel(t)<3, continue; end

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
            E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
            E.GroupMouse(row,1)=grp; E.eventType(row,1)="reward"; E.subType(row,1)="all";
            E.trace{row,1}=ER;
        end
    end

    lickTimes = enforceMinSeparation(detectRisingEdges(t, lk), 0.02);
    if ~isempty(lickTimes)
        [boutStart, boutEnd] = makeLickBouts(lickTimes, bout_gap_s);

        isRewarded = false(size(boutStart));
        for bi=1:numel(boutStart)
            isRewarded(bi) = any(rewardTimes>=boutStart(bi) & rewardTimes<=boutEnd(bi)+1.0);
        end

        addTrace("lickBoutStart","rewarded", extractEventLocked(t,pup,boutStart(isRewarded), tAxis));
        addTrace("lickBoutStart","nonreward",extractEventLocked(t,pup,boutStart(~isRewarded),tAxis));
        addTrace("lickBoutEnd","rewarded",   extractEventLocked(t,pup,boutEnd(isRewarded),  tAxis));
        addTrace("lickBoutEnd","nonreward",  extractEventLocked(t,pup,boutEnd(~isRewarded), tAxis));
    end

    function addTrace(ev, sub, tr)
        if isempty(tr), return; end
        row=row+1;
        E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.session_idx(row,1)=ses;
        E.GroupMouse(row,1)=grp; E.eventType(row,1)=ev; E.subType(row,1)=sub;
        E.trace{row,1}=tr;
    end
end

if isempty(E)
    warning('No event-locked pupil traces could be computed.');
    return;
end

% Passive days 7-9 vs 12-13
plotEventLockedComparison_fixedLegend(E, tAxis, outDir, "Passive", "lickBoutStart", "rewarded", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_rewarded_PASSIVE_days7-9_vs_12-13.png');
plotEventLockedComparison_fixedLegend(E, tAxis, outDir, "Passive", "lickBoutStart", "nonreward", focusDaysA, focusDaysB, ...
    'PUPIL_lickBoutStart_nonreward_PASSIVE_days7-9_vs_12-13.png');

% Reward-locked early vs late (legend fixed)
plotRewardLockedEarlyLate_fixedLegend(E, tAxis, outDir);
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
if isempty(eventTimes), trace=[]; return; end
M = nan(numel(eventTimes), numel(tAxis));
for i=1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    pre = (tAxis < 0);
    b = mean(pi(pre), 'omitnan');
    pi = pi - b;
    M(i,:) = pi;
end
trace = mean(M, 1, 'omitnan');
end

function plotEventLockedComparison_fixedLegend(E, tAxis, outDir, groupName, eventType, subType, daysA, daysB, fname)
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

% shaded areas should not appear in legend
shaded_noLegend(tAxis, muA, seA);
hA = plot(tAxis, muA, 'LineWidth', 2.5);

shaded_noLegend(tAxis, muB, seB);
hB = plot(tAxis, muB, 'LineWidth', 2.5);

xline(0,'k-'); yline(0,'k:');
xlabel('Time from event (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('%s %s (%s): Days %s vs %s', groupName, eventType, subType, rangeStr(daysA), rangeStr(daysB)));
legend([hA hB], {sprintf('Days %s',rangeStr(daysA)), sprintf('Days %s',rangeStr(daysB))}, 'Location','best');
grid on; box on
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotRewardLockedEarlyLate_fixedLegend(E, tAxis, outDir)
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

    h = gobjects(0);
    labels = strings(0);

    if ~isempty(A)
        muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
        shaded_noLegend(tAxis, muA, seA);
        h(end+1) = plot(tAxis, muA, 'LineWidth', 2.5); %#ok<AGROW>
        labels(end+1) = "Days " + rangeStr(earlyDays); %#ok<AGROW>
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shaded_noLegend(tAxis, muB, seB);
        h(end+1) = plot(tAxis, muB, 'LineWidth', 2.5); %#ok<AGROW>
        labels(end+1) = "Days " + rangeStr(lateDays); %#ok<AGROW>
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil (%s): %s vs %s', g, rangeStr(earlyDays), rangeStr(lateDays)));
    if ~isempty(h), legend(h, labels, 'Location','best'); end
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

function shaded_noLegend(x, mu, se)
x=x(:)'; mu=mu(:)'; se=se(:)';
hFill = fill([x fliplr(x)], [mu-se fliplr(mu+se)], 0.9*[1 1 1], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
set(hFill, 'HandleVisibility','off'); % FIX: legend should not use shading
end

function doOptionalPainAndBehaviorTests_explicit(T, outDir)
% Explicit columns only (as you requested)
plotIfExists('Immersion_Latency_s',      'TEST_tailImmersion_latency.png');
plotIfExists('HOT_Frames_Non_moving',    'TEST_hotPlate_nonmoving.png');
plotIfExists('TST_Frames_Non_moving',    'TEST_TST_nonmoving.png');

% STRAUB explicitly (your request)
plotIfExists('STRAUB_Frames_Non_moving', 'TEST_straub_frames_nonmoving.png');
plotIfExists('STRAUB_Pct_Non_moving',    'TEST_straub_pct_nonmoving.png');

    function plotIfExists(varName, fname)
        if ismember(varName, T.Properties.VariableNames)
            plotScalarByDayGroup(T, varName, outDir, fname);
        else
            fid = fopen(fullfile(outDir, ['NOTE_missing_' varName '.txt']),'w');
            fprintf(fid,'Column not found in ALL_mice_longitudinal.csv: %s\n', varName);
            fclose(fid);
        end
    end
end

function plotScalarByDayGroup(T, varName, outDir, fname)
G = findgroups(T.mouse_key_norm, T.day_index, T.GroupMouse);
mk = splitapply(@(x)x(1), T.mouse_key_norm, G);
dy = splitapply(@(x)x(1), T.day_index, G);
gp = splitapply(@(x)x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), T.(varName), G);
S = table(string(mk), dy, gp, val, 'VariableNames', {'mouse','day','group','value'});

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
