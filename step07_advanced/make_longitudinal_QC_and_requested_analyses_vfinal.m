function make_longitudinal_QC_and_requested_analyses_vfinal()
% make_longitudinal_QC_and_requested_analyses_NEWCOHORT_vFINAL
%
% STRICT schema (based on your pipeline output ALL_mice_longitudinal.csv).
% NO guessing for CORE columns.
%
% Required core columns:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1), Injector_TTL (0/1), Diameter_px (numeric), RequirementLast (numeric)
% Time axis:
%   CamTime_rel_s OR PlotTime_s_30fps  (either one must exist)
%
% Outputs:
%   QC + requested analyses:
%     - lick freq, lick count, IEI (lick period), cumulative licks within session
%     - cumulative licks across days (exclude day1-2)
%     - period comparisons (with/without transition days)
%     - pupil normalization by day3 baseline
%     - event-locked pupil: reward-locked + lick-bout-locked (rewarded vs nonreward)
%     - optional: Tail immersion / TST / Hot plate / Straub if present (prefix-based, no name invention)

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Cumulative licking within session
dt_bin_s   = 60;     % 1-min bins
bout_gap_s = 2.0;    % lick bout: pause > 2s splits bouts

% Pupil event-lock window
pupil_win = [-2 2];          % seconds
resample_dt = 0.05;          % 20 Hz for event-locked traces
min_reward_separation_s = 0.5; % drop rewards closer than this (rare TTL glitches)

% Period definitions
PER.pre       = 3:5;
PER.during    = 6:10;
PER.post      = 11:13;
PER.withdraw  = 14:16;
PER.reexpo    = 17:18;

% Transition / unreliable days (optionally exclude)
transitionDays = [4 6 11 14];

% Passive focus comparison days you requested
FOCUS_PASSIVE_DAYS_A = 7:9;     % passive during subset
FOCUS_PASSIVE_DAYS_B = 12:13;   % post subset

% Exclude habituation from cumulative-across-days
excludeDays_global = [1 2];

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
if exist(outDir,'dir')~=7, mkdir(outDir); end

fprintf('\nUsing runDir:\n  %s\n', runDir);
fprintf('Reading CSV:\n  %s\n', csvPath);
fprintf('Saving outputs to:\n  %s\n\n', outDir);

%% ===================== LOAD TABLE =====================
T = readtable(csvPath, 'VariableNamingRule','preserve');

% Enforce string mouse_key
assertHasVars(T, {'mouse_key','day_index','session_idx'}, 'Top-level');
T.mouse_key = string(T.mouse_key);

% Required core columns (STRICT)
REQ = {'mouse_key','day_index','session_idx','Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, REQ, 'Core longitudinal schema');

% Time axis (one must exist)
timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end

% Remove non-finite time rows (explains your warning)
badT = ~isfinite(double(T.(timeVar)));
if any(badT)
    fprintf('[WARN] Removing %d rows with non-finite %s\n', sum(badT), timeVar);
    T = T(~badT,:);
end

% Cohort metadata (your NEW cohort)
COH = buildNewCohortTable();
T = addCohortMeta_strict(T, COH);

% Drop NaN day rows (safety)
T = T(isfinite(double(T.day_index)),:);

%% ===================== QC 1: missingness heatmaps =====================
doMissingnessQC(T, outDir);

%% ===================== Build session summary =====================
S = buildSessionSummary(T, timeVar);

%% ===================== QC 2: outliers =====================
doOutlierQC(S, outDir);

%% ===================== QC 3: spaghetti (mouse-level) =====================
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'nLicks',        outDir, 'QC_spaghetti_nLicks.png');
plotSpaghettiByGroup(S, 'iei_median',    outDir, 'QC_spaghetti_iei_median.png');
plotSpaghettiByGroup(S, 'RequirementLast', outDir, 'QC_spaghetti_PRscore_RequirementLast.png');
plotPupilByGroupDay(S, outDir);

%% ===================== REQUESTED LICKING ANALYSES =====================
% A) Within-session cumulative licking curves: day-by-day, Active vs Passive
plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays_global);

% B) Cumulative licking across days (exclude day1-2): Active vs Passive
plotCumulativeAcrossDays(S, outDir, excludeDays_global);

% C) Lick freq + lick count + lick period summary plots (by day, group)
plotDailyMetricGroupMeans(S, outDir, excludeDays_global);

%% ===================== REQUESTED PERIOD COMPARISONS =====================
DAYSETS.allDays.name = 'ALLDAYS';
DAYSETS.allDays.mask = true(height(S),1);

DAYSETS.noTransition.name = 'NO_TRANSITION';
DAYSETS.noTransition.mask = ~ismember(double(S.day_index), transitionDays);

doPeriodComparisons_strict(S, outDir, PER, DAYSETS, transitionDays);

%% ===================== REQUESTED PUPIL ANALYSES =====================
% 1) Normalize pupil by Day3 baseline (per mouse)
S = addPupilBaselineNormalization(S);

plotPupilNormalizedDayTraces(S, outDir);

% 2) Event-locked pupil: reward-locked + lick-bout-locked (rewarded vs nonreward)
doEventLockedPupilAnalyses_strict( ...
    T, timeVar, outDir, pupil_win, resample_dt, bout_gap_s, ...
    min_reward_separation_s, FOCUS_PASSIVE_DAYS_A, FOCUS_PASSIVE_DAYS_B);

%% ===================== OPTIONAL: Immersion / HOT / TST / STRAUB =====================
doOptionalPainAndBehaviorTests_prefixStrict(T, outDir);

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
    msg = strjoin(string(missing), ', ');
    error('Fix CSV schema or update required variable list in the script. Missing: %s', msg);
end
end

function COH = buildNewCohortTable()
% Canonical mouse_key in this script: "6100_red" etc.
rows = {
'6100','red',    'F','P','6100_pairA','Passive';
'6100','orange', 'F','P','6100_pairA','Passive';
'6100','black',  'F','A','6100_pairA','Active';

'0911','red',    'F','A','0911_pairA','Active';
'0911','orange', 'F','P','0911_pairA','Passive';

'0911','black',  'F','P','0911_pairB','Passive';
'0911','white',  'F','A','0911_pairB','Active';

'0910','red',    'M','P','0910_pairA','Passive';
'0910','orange', 'M','P','0910_pairA','Passive';
'0910','black',  'M','A','0910_pairA','Active';

'6099','red',    'M','P','6099_pairA','Passive';
'6099','orange', 'M','A','6099_pairA','Active';

'6099','black',  'M','A','6099_pairB','Active';
'6099','white',  'M','P','6099_pairB','Passive'; % died after day13 => naturally missing later days
};
COH = cell2table(rows, 'VariableNames', {'cage','color','sex','group','pair_id','pair_role'});
COH.mouse_key_norm = lower(string(COH.cage) + "_" + string(COH.color));
end

function T = addCohortMeta_strict(T, COH)
% Normalize incoming mouse_key to cage_color (strict transform, not guessing columns)
mk = lower(strrep(strrep(strtrim(T.mouse_key),'-','_'),' ',''));

% Handle "6100red" -> "6100_red" (this is not schema guessing; it’s ID normalization)
noUnd = ~contains(mk,'_');
mk(noUnd) = regexprep(mk(noUnd), '^(\d{4})([a-zA-Z]+)$', '$1_$2');

T.mouse_key_norm = mk;

% Join (left join)
J = COH(:, {'mouse_key_norm','cage','color','sex','group','pair_id','pair_role'});
T = outerjoin(T, J, 'Keys','mouse_key_norm', 'MergeKeys', true, 'Type','left');

% Group label
g = string(T.group);
g(ismissing(g) | g=="") = "U";
T.GroupMouse = categorical(g, ["A","P","U"], {'Active','Passive','Unknown'});
end

function doMissingnessQC(T, outDir)
vars = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
mk = unique(T.mouse_key_norm,'stable');
days = unique(double(T.day_index)); days = sort(days(isfinite(days)));

M = nan(numel(mk), numel(days), numel(vars));
for i=1:numel(mk)
    for j=1:numel(days)
        r = (T.mouse_key_norm==mk(i)) & (double(T.day_index)==days(j));
        if ~any(r), continue; end
        for v=1:numel(vars)
            x = T.(vars{v})(r);
            M(i,j,v) = mean(~isfinite(double(x)));
        end
    end
end

for v=1:numel(vars)
    fig = figure('Color','w','Position',[80 80 1200 520]);
    imagesc(days, 1:numel(mk), M(:,:,v));
    colormap(parula); colorbar;
    xlabel('Day'); ylabel('Mouse');
    title(['Missing fraction: ' vars{v}], 'Interpreter','none');
    set(gca,'YTick',1:numel(mk),'YTickLabel',mk,'TickLabelInterpreter','none');
    exportgraphics(fig, fullfile(outDir, ['QC_missing_' vars{v} '.png']), 'Resolution', 200);
    close(fig);
end
end

function S = buildSessionSummary(T, timeVar)
% One row per (mouse, day, session)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

S = table;
S.mouse_key    = splitapply(@(x) string(x(1)), T.mouse_key_norm, G);
S.day_index    = splitapply(@(x) double(x(1)), T.day_index,     G);
S.session_idx  = splitapply(@(x) double(x(1)), T.session_idx,   G);
S.GroupMouse   = splitapply(@(x) x(1),         T.GroupMouse,    G);

% PR score
S.RequirementLast = splitapply(@(x) mean(double(x),'omitnan'), T.RequirementLast, G);

% Pupil mean
S.pupil_mean_px   = splitapply(@(x) mean(double(x),'omitnan'), T.Diameter_px, G);

% Lick metrics
n = height(S);
S.nLicks      = zeros(n,1);
S.licks_per_min = nan(n,1);
S.iei_median  = nan(n,1);
S.iei_mean    = nan(n,1);
S.iei_cv      = nan(n,1);

for i=1:n
    r = (T.mouse_key_norm==S.mouse_key(i)) & (double(T.day_index)==S.day_index(i)) & (double(T.session_idx)==S.session_idx(i));
    t  = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));
    if isempty(t) || all(~isfinite(t)), continue; end

    [t,ord] = sort(t); lk = lk(ord);
    lickTimes = detectRisingEdges(t, lk);

    S.nLicks(i) = numel(lickTimes);
    dur = max(t)-min(t);
    if isfinite(dur) && dur>1
        S.licks_per_min(i) = (numel(lickTimes)/dur)*60;
    end
    if numel(lickTimes)>=2
        iei = diff(lickTimes);
        S.iei_median(i) = median(iei,'omitnan');
        S.iei_mean(i)   = mean(iei,'omitnan');
        S.iei_cv(i)     = std(iei,0,'omitnan')/max(eps,mean(iei,'omitnan'));
    end
end
end

function times = detectRisingEdges(t, x01)
x01 = x01(:);
x01(~isfinite(x01)) = 0;
x01 = x01 > 0.5;
dx = diff([false; x01]);
idx = find(dx==1);
times = t(idx);
times = times(isfinite(times));
end

function doOutlierQC(S, outDir)
metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
for k=1:numel(metrics)
    yvar = metrics{k};
    if ~ismember(yvar, S.Properties.VariableNames), continue; end
    x = double(S.(yvar));
    [isOut, z] = madOutliers(x, 4.5);

    fig = figure('Color','w','Position',[80 80 950 420]); hold on
    scatter(1:numel(x), x, 18, 'filled','MarkerFaceAlpha',0.5);
    scatter(find(isOut), x(isOut), 55, 'o','LineWidth',1.8);
    title(sprintf('QC outliers (MAD z>4.5): %s', yvar), 'Interpreter','none');
    xlabel('Session #'); ylabel(yvar, 'Interpreter','none'); grid on; box on
    exportgraphics(fig, fullfile(outDir, ['QC_outliers_' yvar '.png']), 'Resolution', 200);
    close(fig);

    if any(isOut)
        Tout = S(isOut, {'mouse_key','day_index','session_idx','GroupMouse'});
        Tout.metric = repmat(string(yvar), height(Tout),1);
        Tout.value  = x(isOut);
        Tout.madz   = z(isOut);
        writetable(Tout, fullfile(outDir, ['QC_outliers_' yvar '.csv']));
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
fig = figure('Color','w','Position',[80 80 1150 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

groups = {'Active','Passive'};
for gi=1:2
    g = groups{gi};
    nexttile; hold on
    rG = string(S.GroupMouse)==g;
    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & (S.mouse_key==mk(i));
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
    xlabel('Day'); ylabel(yvar,'Interpreter','none'); grid on; box on
end

exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotPupilByGroupDay(S, outDir)
if ~ismember('pupil_mean_px', S.Properties.VariableNames), return; end
fig = figure('Color','w','Position',[80 80 1150 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
groups = {'Active','Passive'};
for gi=1:2
    g = groups{gi};
    nexttile; hold on
    rG = string(S.GroupMouse)==g;
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
    plot(dAll, mu, 'k-', 'LineWidth',3);
    title(sprintf('Pupil mean (px): %s', g));
    xlabel('Day'); ylabel('Pupil diameter (px)'); grid on; box on
end
exportgraphics(fig, fullfile(outDir,'QC_pupil_mean_px_by_group.png'), 'Resolution', 220);
close(fig);
end

function plotCumulativeLickingAcrossSession(T, timeVar, outDir, dt_bin_s, excludeDays)
for gname = ["Active","Passive"]
    Rg = (string(T.GroupMouse)==gname) & ~ismember(double(T.day_index), excludeDays);
    if ~any(Rg), continue; end
    days = unique(double(T.day_index(Rg))); days = sort(days);

    fig = figure('Color','w','Position',[80 80 1200 520]); hold on
    leg = strings(0,1);

    for di=1:numel(days)
        d = days(di);
        Td = T(Rg & double(T.day_index)==d,:);
        [tgrid, muCum] = meanCumulativeLicks_sessions(Td, timeVar, dt_bin_s);
        if isempty(tgrid), continue; end
        plot(tgrid/60, muCum, 'LineWidth', 2);
        leg(end+1) = "Day " + d; %#ok<AGROW>
    end

    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title(sprintf('Within-session cumulative licks (%s)', gname));
    grid on; box on
    if ~isempty(leg), legend(leg,'Location','eastoutside'); end
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_cum_withinSession_dayByDay_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function [tgrid, muCum] = meanCumulativeLicks_sessions(Tsub, timeVar, dt_bin_s)
if isempty(Tsub), tgrid=[]; muCum=[]; return; end
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(Tsub(:,keys));

sess = {};
tmaxAll = 0;

for gi=1:max(G)
    r = (G==gi);
    t  = double(Tsub.(timeVar)(r));
    lk = double(Tsub.Lick_TTL(r));
    if isempty(t) || all(~isfinite(t)), continue; end
    [t,ord] = sort(t); lk=lk(ord);
    lickTimes = detectRisingEdges(t, lk);

    t0 = min(t);
    dur = max(t)-t0;
    if ~isfinite(dur), continue; end
    tmaxAll = max(tmaxAll, dur);

    sess{end+1} = lickTimes - t0; %#ok<AGROW>
end

if isempty(sess), tgrid=[]; muCum=[]; return; end
tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;

M = nan(numel(sess), numel(tgrid));
for i=1:numel(sess)
    lt = sess{i};
    if isempty(lt)
        M(i,:) = 0;
    else
        M(i,:) = arrayfun(@(tt) sum(lt<=tt), tgrid);
    end
end
muCum = mean(M,1,'omitnan');
end

function plotCumulativeAcrossDays(S, outDir, excludeDays)
fig = figure('Color','w','Position',[80 80 1100 560]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

for gname = ["Active","Passive"]
    nexttile; hold on
    rG = (string(S.GroupMouse)==gname) & ~ismember(double(S.day_index), excludeDays);
    if ~any(rG), title(sprintf('%s: no data', gname)); axis off; continue; end

    mk = unique(S.mouse_key(rG),'stable');
    for i=1:numel(mk)
        r = rG & (S.mouse_key==mk(i));
        [d,ord] = sort(S.day_index(r));
        y = double(S.nLicks(r)); y=y(ord);
        plot(d, cumsum(y,'omitnan'), '-o', 'LineWidth',1, 'MarkerSize',3);
    end
    xlabel('Day'); ylabel('Cumulative licks'); title(sprintf('Cumulative licks across days (%s)', gname));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'LICK_cumulative_acrossDays_active_vs_passive.png'), 'Resolution', 220);
close(fig);
end

function plotDailyMetricGroupMeans(S, outDir, excludeDays)
metrics = {'licks_per_min','nLicks','iei_median'};
for m=1:numel(metrics)
    yvar = metrics{m};
    if ~ismember(yvar, S.Properties.VariableNames), continue; end

    fig = figure('Color','w','Position',[80 80 1000 420]); hold on
    for gname = ["Active","Passive"]
        rG = (string(S.GroupMouse)==gname) & ~ismember(double(S.day_index), excludeDays);
        if ~any(rG), continue; end
        days = unique(S.day_index(rG)); days = sort(days);
        mu = nan(size(days)); se = nan(size(days));
        for j=1:numel(days)
            x = double(S.(yvar)(rG & S.day_index==days(j)));
            mu(j) = mean(x,'omitnan');
            se(j) = std(x,0,'omitnan')/sqrt(max(1,sum(isfinite(x))));
        end
        errorbar(days, mu, se, '-o', 'LineWidth',2, 'MarkerSize',6);
    end
    xlabel('Day'); ylabel(yvar,'Interpreter','none');
    title(sprintf('Daily %s (mean±SEM): Active vs Passive', yvar), 'Interpreter','none');
    grid on; box on; legend({'Active','Passive'},'Location','best');
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_daily_%s_meanSEM.png', yvar)), 'Resolution', 220);
    close(fig);
end
end

function doPeriodComparisons_strict(S, outDir, PER, DAYSETS, transitionDays)
% Assign period label (no blank category)
period = strings(height(S),1);
period(ismember(S.day_index, PER.pre))       = "pre";
period(ismember(S.day_index, PER.during))    = "during";
period(ismember(S.day_index, PER.post))      = "post";
period(ismember(S.day_index, PER.withdraw))  = "withdrawal";
period(ismember(S.day_index, PER.reexpo))    = "reexposure";
period(period=="") = "<undef>";
S.period = categorical(period, ...
    ["pre","during","post","withdrawal","reexposure","<undef>"], ...
    ["pre","during","post","withdrawal","reexposure","<undef>"]);

metrics = {'licks_per_min','nLicks','iei_median','pupil_mean_px','RequirementLast'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

groups = ["Active","Passive"];
for dsName = fieldnames(DAYSETS)'
    ds = DAYSETS.(dsName{1});
    for mi=1:numel(metrics)
        yvar = metrics{mi};

        fig = figure('Color','w','Position',[80 80 1200 520]);
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        for gi=1:2
            g = groups(gi);
            nexttile; hold on
            Rg = (string(S.GroupMouse)==g) & ds.mask & (S.period~="<undef>");

            % Passive PR rule: exclude "during" for RequirementLast (PR score not meaningful there)
            if yvar=="RequirementLast" && g=="Passive"
                Rg = Rg & (S.period~="during");
            end

            if ~any(Rg)
                title(sprintf('%s (%s) – no data', yvar, g), 'Interpreter','none');
                axis off; continue;
            end

            % mouse-level mean per period
            [Gmp, ~] = findgroups(S.mouse_key(Rg), S.period(Rg));
            val = splitapply(@(x) mean(double(x),'omitnan'), S.(yvar)(Rg), Gmp);
            pr  = splitapply(@(x) x(1), S.period(Rg), Gmp);

            boxchart(double(pr), val);
            cats = categories(S.period); cats = cats(~strcmp(cats,'<undef>'));
            set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
            ylabel(yvar,'Interpreter','none');
            title(sprintf('%s (%s) – %s', yvar, g, ds.name), 'Interpreter','none');
            grid on; box on
        end

        exportgraphics(fig, fullfile(outDir, sprintf('PERIOD_box_%s_%s.png', yvar, ds.name)), 'Resolution', 220);
        close(fig);
    end
end

fid = fopen(fullfile(outDir,'PERIOD_transition_days_note.txt'),'w');
fprintf(fid,'Transition/unreliable days excluded in NO_TRANSITION: %s\n', mat2str(transitionDays));
fprintf(fid,'Active compare: Day5 vs (6-10) vs (11-13) vs (14-16) vs Day17\n');
fprintf(fid,'Passive compare: Day5 vs (11-13) vs (14-16) vs (17-18)\n');
fclose(fid);
end

function S = addPupilBaselineNormalization(S)
S.pupil_base_day3 = nan(height(S),1);
S.pupil_norm_rel  = nan(height(S),1);

mk = unique(S.mouse_key,'stable');
for i=1:numel(mk)
    rM = (S.mouse_key==mk(i));
    base = mean(double(S.pupil_mean_px(rM & S.day_index==3)), 'omitnan');
    if ~isfinite(base) || base<=0, continue; end
    S.pupil_base_day3(rM) = base;
    S.pupil_norm_rel(rM)  = (double(S.pupil_mean_px(rM)) - base) ./ base;
end
end

function plotPupilNormalizedDayTraces(S, outDir)
if ~ismember('pupil_norm_rel', S.Properties.VariableNames), return; end

fig = figure('Color','w','Position',[80 80 1150 640]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
groups = {'Active','Passive'};

for gi=1:2
    g = groups{gi};
    nexttile; hold on
    rG = string(S.GroupMouse)==g;
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
    plot(dAll, mu, 'k-', 'LineWidth',3);
    yline(0,'k:');
    xlabel('Day'); ylabel('\Delta pupil / baseline(D3)');
    title(sprintf('Pupil normalized by Day3 baseline (%s)', g));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_dayLevel_normalized_by_day3.png'), 'Resolution', 220);
close(fig);
end

function doEventLockedPupilAnalyses_strict(T, timeVar, outDir, win_s, dt, bout_gap_s, min_reward_sep, daysA, daysB)
tAxis = win_s(1):dt:win_s(2);

% group per session
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

E = table;
row = 0;

for gi=1:max(G)
    r = (G==gi);

    t   = double(T.(timeVar)(r));
    lk  = double(T.Lick_TTL(r));
    rw  = double(T.Injector_TTL(r));
    pup = double(T.Diameter_px(r));

    if isempty(t) || all(~isfinite(t)) || all(~isfinite(pup)), continue; end
    [t,ord] = sort(t); lk=lk(ord); rw=rw(ord); pup=pup(ord);

    mk  = string(T.mouse_key_norm(find(r,1,'first')));
    day = double(T.day_index(find(r,1,'first')));
    grp = string(T.GroupMouse(find(r,1,'first')));

    rewardTimes = detectRisingEdges(t, rw);
    rewardTimes = enforceMinSeparation(rewardTimes, min_reward_sep);

    lickTimes = detectRisingEdges(t, lk);
    lickTimes = enforceMinSeparation(lickTimes, 0.02);

    % Reward-locked
    if ~isempty(rewardTimes)
        tr = extractEventLockedTrace(t, pup, rewardTimes, tAxis);
        if ~isempty(tr)
            row=row+1;
            E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.GroupMouse(row,1)=grp;
            E.eventType(row,1)="reward"; E.subType(row,1)="all"; E.trace{row,1}=tr;
        end
    end

    % Lick-bout locked
    if ~isempty(lickTimes)
        [boutStart, boutEnd] = makeLickBouts(lickTimes, bout_gap_s);
        isRewarded = false(size(boutStart));
        for bi=1:numel(boutStart)
            isRewarded(bi) = any(rewardTimes>=boutStart(bi) & rewardTimes<=boutEnd(bi)+1.0);
        end

        tr_rs = extractEventLockedTrace(t, pup, boutStart(isRewarded), tAxis);
        tr_ns = extractEventLockedTrace(t, pup, boutStart(~isRewarded), tAxis);

        if ~isempty(tr_rs)
            row=row+1;
            E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.GroupMouse(row,1)=grp;
            E.eventType(row,1)="lickBoutStart"; E.subType(row,1)="rewarded"; E.trace{row,1}=tr_rs;
        end
        if ~isempty(tr_ns)
            row=row+1;
            E.mouse_key(row,1)=mk; E.day_index(row,1)=day; E.GroupMouse(row,1)=grp;
            E.eventType(row,1)="lickBoutStart"; E.subType(row,1)="nonreward"; E.trace{row,1}=tr_ns;
        end
    end
end

if isempty(E)
    warning('No event-locked pupil traces computed.');
    return;
end

% Passive-only focus: days 7-9 vs 12-13 for lickBoutStart rewarded/nonreward
plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "rewarded", daysA, daysB, ...
    'PUPIL_lickBoutStart_rewarded_PASSIVE_days7-9_vs_12-13.png');
plotEventLockedComparison(E, tAxis, outDir, "Passive", "lickBoutStart", "nonreward", daysA, daysB, ...
    'PUPIL_lickBoutStart_nonreward_PASSIVE_days7-9_vs_12-13.png');

% Reward-locked early vs during (example overview)
plotRewardLockedOverview(E, tAxis, outDir);
end

function times = enforceMinSeparation(times, minSep)
if isempty(times), return; end
times = sort(times(:));
keep = true(size(times));
last = -Inf;
for i=1:numel(times)
    if times(i)-last < minSep
        keep(i)=false;
    else
        last = times(i);
    end
end
times = times(keep);
end

function [boutStart, boutEnd] = makeLickBouts(lickTimes, gap_s)
lickTimes = sort(lickTimes(:));
if isempty(lickTimes), boutStart=[]; boutEnd=[]; return; end
bs = lickTimes(1); be = lickTimes(1);
starts = []; ends = [];
for i=2:numel(lickTimes)
    if lickTimes(i)-lickTimes(i-1) > gap_s
        starts(end+1,1)=bs; ends(end+1,1)=be; %#ok<AGROW>
        bs = lickTimes(i); be = lickTimes(i);
    else
        be = lickTimes(i);
    end
end
starts(end+1,1)=bs; ends(end+1,1)=be;
boutStart = starts; boutEnd = ends;
end

function tr = extractEventLockedTrace(t, pup, eventTimes, tAxis)
if isempty(eventTimes), tr=[]; return; end
M = nan(numel(eventTimes), numel(tAxis));
preMask = (tAxis < 0);

for i=1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    b  = mean(pi(preMask), 'omitnan');
    pi = pi - b;
    M(i,:) = pi;
end
tr = mean(M,1,'omitnan');
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

fig = figure('Color','w','Position',[80 80 900 520]); hold on
shaded(tAxis, muA, seA); plot(tAxis, muA, 'LineWidth',2.5);
shaded(tAxis, muB, seB); plot(tAxis, muB, 'LineWidth',2.5);
xline(0,'k-'); yline(0,'k:');
xlabel('Time from event (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('%s %s (%s): Days %s vs %s', groupName, eventType, subType, rangeStr(daysA), rangeStr(daysB)));
legend({sprintf('Days %s',rangeStr(daysA)), sprintf('Days %s',rangeStr(daysB))}, 'Location','best');
grid on; box on
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotRewardLockedOverview(E, tAxis, outDir)
fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

groups = ["Active","Passive"];
for gi=1:2
    g = groups(gi);
    nexttile; hold on
    r = (string(E.GroupMouse)==g) & (string(E.eventType)=="reward");
    M = stackTraces(E.trace(r));
    if ~isempty(M)
        mu = mean(M,1,'omitnan'); se = std(M,0,1,'omitnan')/sqrt(size(M,1));
        shaded(tAxis, mu, se); plot(tAxis, mu, 'LineWidth',2.5);
    end
    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil');
    title(sprintf('Reward-locked pupil (all days): %s', g));
    grid on; box on
end
exportgraphics(fig, fullfile(outDir,'PUPIL_rewardLocked_overview_active_passive.png'), 'Resolution', 220);
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

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1), s=sprintf('%d-%d', v(1), v(end));
else, s=strjoin(string(v),'-'); end
end

function doOptionalPainAndBehaviorTests_prefixStrict(T, outDir)
% Do not invent names.
% Only summarize if columns exist.

% Tail immersion (your pipeline mentions immersion latency)
if ismember('Immersion_Latency_s', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'Immersion_Latency_s', outDir, 'TEST_tailImmersion_latency.png');
end

% Manual scoring summaries often create multiple columns with these prefixes
prefixes = ["TST_","HOT_","STRAUB_"];
found = strings(0,1);
for p = prefixes
    cols = string(T.Properties.VariableNames);
    hit = cols(startsWith(cols, p));
    hit = hit(arrayfun(@(c) isnumeric(T.(c)), cellstr(hit)));
    if ~isempty(hit)
        found = [found; hit(:)]; %#ok<AGROW>
        % Plot first few columns as daily mean overview
        for k=1:min(4,numel(hit))
            plotScalarByDayGroup(T, hit(k), outDir, sprintf('TEST_%s.png', hit(k)));
        end
    end
end

fid = fopen(fullfile(outDir,'OPTIONAL_tests_columns_found.txt'),'w');
fprintf(fid,'Columns found & used (prefix-based):\n');
if isempty(found)
    fprintf(fid,'  None found (this is OK if not merged into ALL_mice_longitudinal.csv yet).\n');
else
    for i=1:numel(found), fprintf(fid,'  %s\n', found(i)); end
end
fclose(fid);
end

function plotScalarByDayGroup(T, varName, outDir, fname)
varName = string(varName);
assert(ismember(varName, string(T.Properties.VariableNames)), 'Missing %s', varName);

% Reduce to mouse-day mean
G = findgroups(T.mouse_key_norm, double(T.day_index), T.GroupMouse);
mk  = splitapply(@(x) string(x(1)), T.mouse_key_norm, G);
dy  = splitapply(@(x) double(x(1)), T.day_index, G);
gp  = splitapply(@(x) x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), T.(varName), G);

S = table(mk, dy, gp, val, 'VariableNames', {'mouse','day','group','value'});

fig = figure('Color','w','Position',[80 80 950 420]); hold on
for gname = ["Active","Passive"]
    rG = string(S.group)==gname;
    if ~any(rG), continue; end
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
legend({'Active','Passive'},'Location','best');
exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end
