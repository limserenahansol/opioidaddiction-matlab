function make_longitudinal_QC_and_requested_analyses_NEWCOHORT_FINAL()
% NEW MASTER SCRIPT (FINAL):
% QC + Licking metrics + Period comparisons (with/without transition days)
% + Pupil day-normalization (baseline Day3) + Event-locked pupil (reward + lick-bout)
% + Straub / TST / Tail immersion / Hot plate summary plots + robust stats
%
% INPUT:
%   rootTry/run_*/ALL_mice_longitudinal.csv
%
% REQUIRED columns (minimum):
%   mouse_key, day_index, session_idx
%   Lick_TTL (0/1), Injector_TTL (0/1), Diameter_px (numeric)
%   plus a time vector: CamTime_rel_s OR PlotTime_s_30fps
%
% OPTIONAL columns for assays (any subset OK):
%   STRAUB_*, TST_*, Immersion_Latency_s, HOT_*
%
% OUTPUT:
%   runDir/QC_AND_REQUESTED_ANALYSES_FINAL_yyyymmdd_HHMMSS/

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Your task phases
PER.pre       = 3:5;      % water PR
PER.during    = 6:10;     % passive training morphine (active has PR, passive forced reward)
PER.post      = 11:13;    % morphine PR all
PER.withdraw  = 14:16;    % withdrawal (water)
PER.reexpo    = 17:18;    % re-exposure morphine

% Transition days you said are less reliable (edit if needed)
transitionDays = [6 11 14];

% Exclusions for licking plots
excludeDays_licking_global = [1 2];

% Lick bout definition
bout_gap_s = 2.0;               % if no lick for >2s => new bout
min_event_separation_s = 0.50;  % for reward events de-dup

% Pupil event window
pupil_win = [-2 2];
dt_event  = 0.1;

% Cumulative licking binning
dt_bin_s = 60;

% Passive “during” handling:
% true  = keep passive during (it will show forced reward based effects)
% false = remove passive during from phase stats (often what you want)
INCLUDE_PASSIVE_DURING_IN_PHASE_STATS = true;

%% ===================== FIND LATEST run_* =====================
assert(exist(rootTry,'dir')==7, 'rootTry not found: %s', rootTry);
D = dir(fullfile(rootTry,'run_*'));  assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')==2, 'Expected CSV not found: %s', csvPath);

ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(runDir, ['QC_AND_REQUESTED_ANALYSES_FINAL_' ts]);
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\n[INFO] runDir:\n  %s\n', runDir);
fprintf('[INFO] CSV:\n  %s\n', csvPath);
fprintf('[INFO] outDir:\n  %s\n\n', outDir);

%% ===================== LOAD TABLE =====================
T = readtable(csvPath, 'VariableNamingRule','preserve');
assertHasVars(T, {'mouse_key','day_index','session_idx'}, 'Top-level');

T.mouse_key = string(T.mouse_key);
T.day_index = double(T.day_index);
T.session_idx = double(T.session_idx);

timeVar = pickTimeVar(T);
assertHasVars(T, {'Lick_TTL','Injector_TTL','Diameter_px'}, 'Core schema');

% Cohort mapping (your new cohort)
COH = buildNewCohortTable_NEW();
T = addCohortMeta(T, COH);

% Remove rows with NaN day or non-finite time
T = T(isfinite(T.day_index),:);
badTime = ~isfinite(double(T.(timeVar)));
if any(badTime)
    Tw = T(badTime, {'mouse_key_norm','day_index','session_idx'});
    writetable(Tw, fullfile(outDir, sprintf('WARN_removed_nonfinite_%s_rows.csv', timeVar)));
    T = T(~badTime,:);
end

%% ===================== QC 1: MISSINGNESS MAPS =====================
doMissingnessQC(T, outDir, {'Lick_TTL','Injector_TTL','Diameter_px'});

%% ===================== BUILD SESSION SUMMARY (lick + pupil day means) =====================
S = buildSessionSummary_FULL(T, timeVar, bout_gap_s);

% Save a compact summary to inspect quickly
writetable(S, fullfile(outDir, 'SESSION_SUMMARY_FULL.csv'));

%% ===================== QC 2: OUTLIER FLAGS (MAD) =====================
doOutlierQC(S, outDir);

%% ===================== QC 3: SPAGHETTI QUICKLOOKS =====================
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'nLicks',        outDir, 'QC_spaghetti_nLicks.png');
plotSpaghettiByGroup(S, 'bout_count',    outDir, 'QC_spaghetti_bout_count.png');
plotSpaghettiByGroup(S, 'iei_median',    outDir, 'QC_spaghetti_iei_median.png');
plotSpaghettiByGroup(S, 'pupil_mean_px', outDir, 'QC_spaghetti_pupil_mean_px.png');

%% ===================== LICKING: DAY-BY-DAY CUMULATIVE CURVES =====================
% Like your example: each day is a different colored curve (per group)
plotCumulativeLickingAcrossSession_byDay(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global, "Active");
plotCumulativeLickingAcrossSession_byDay(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global, "Passive");

% Like your example: passive=blue, active=red, averaged across (all days except 1-2)
plotCumulativeLicking_GroupOverlay(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global);

%% ===================== PHASE / PERIOD COMPARISONS (WITH + WITHOUT TRANSITION DAYS) =====================
DAYSETS.allDays.name  = 'ALLDAYS';
DAYSETS.allDays.keepDay = @(d) true(size(d));

DAYSETS.noTransition.name = 'NO_TRANSITION';
DAYSETS.noTransition.keepDay = @(d) ~ismember(d, transitionDays);

doPhaseComparisons_withStats_FINAL(S, outDir, PER, DAYSETS, INCLUDE_PASSIVE_DURING_IN_PHASE_STATS);

%% ===================== PUPIL: DAY-NORMALIZED (BASELINE = DAY3) =====================
S = addPupilBaselineNormalization_DAY3(S);
writetable(S, fullfile(outDir, 'SESSION_SUMMARY_WITH_PUPIL_NORMALIZATION.csv'));

plotPupilNormalizedDayTraces(S, outDir);                       % active vs passive by day
plotPupilMorphineVsWater_DayCompare(S, outDir, [4 5], [12 13]);% your requested comparison style

%% ===================== PUPIL: EVENT-LOCKED (REWARD + LICK-BOUT; REWARDED vs UNREWARDED) =====================
tAxis = pupil_win(1):dt_event:pupil_win(2);
E = buildEventLockedTable_FULL(T, timeVar, tAxis, bout_gap_s, min_event_separation_s);

if ~isempty(E)
    writetable(E(:,{'mouse_key','day_index','session_idx','GroupMouse','eventType','subType'}), ...
        fullfile(outDir,'EVENTLOCK_index.csv'));
    plotRewardLockedEarlyLate(E, tAxis, outDir, PER);         % your “reward-locked pupil” style
    plotLickBoutLocked_RewardedVsUnrewarded(E, tAxis, outDir, PER);
else
    fid=fopen(fullfile(outDir,'PUPIL_EVENTLOCK_NOTE.txt'),'w');
    fprintf(fid,'No event-locked pupil traces computed. Check Injector_TTL/Lick_TTL and Diameter_px availability.\n');
    fclose(fid);
end

%% ===================== STRAUB / TST / TAIL IMMERSION / HOT PLATE =====================
doAssayPlotsAndStats(T, outDir);

fprintf('\nDONE.\nOutputs saved in:\n  %s\n', outDir);
end

%% =====================================================================
%% =============================== HELPERS ==============================
%% =====================================================================

function timeVar = pickTimeVar(T)
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end
end

function assertHasVars(T, vars, context)
missing = vars(~ismember(vars, T.Properties.VariableNames));
if ~isempty(missing)
    fprintf('\nSCHEMA ERROR in %s.\nMissing required columns:\n', context);
    disp(missing(:));
    error('Missing required columns: %s', strjoin(string(missing), ', '));
end
end

function COH = buildNewCohortTable_NEW()
% Your mapping (cage + color) -> sex/group/pair
rows = {
'6100','red',    'F','P','6100_pairA';
'6100','orange', 'F','P','6100_pairA';
'6100','black',  'F','A','6100_pairA';

'0911','red',    'F','A','0911_pairA';
'0911','orange', 'F','P','0911_pairA';

'0911','black',  'F','P','0911_pairB';
'0911','white',  'F','A','0911_pairB';

'0910','red',    'M','P','0910_pairA';
'0910','orange', 'M','P','0910_pairA';
'0910','black',  'M','A','0910_pairA';

'6099','red',    'M','P','6099_pairA';
'6099','orange', 'M','A','6099_pairA';

'6099','black',  'M','A','6099_pairB';
'6099','white',  'M','P','6099_pairB';
};
COH = cell2table(rows, 'VariableNames', {'cage','color','sex','ap','pair_id'});
COH.mouse_key_norm = lower(string(COH.cage) + "_" + string(COH.color));
COH.GroupMouse = categorical(COH.ap, ["A","P"], {'Active','Passive'});
end

function T = addCohortMeta(T, COH)
mk = string(T.mouse_key);
mk = regexprep(mk,'\s+','');
mk = regexprep(mk,'-','_');
mk = lower(mk);
% If someone stored "6100red" style, convert to "6100_red"
isNoUnd = ~contains(mk,'_');
mk(isNoUnd) = regexprep(mk(isNoUnd), '^(\d{4})([a-z]+)$', '$1_$2');
T.mouse_key_norm = mk;

T = outerjoin(T, COH(:,{'mouse_key_norm','cage','color','sex','pair_id','GroupMouse'}), ...
    'Keys','mouse_key_norm','MergeKeys',true,'Type','left');

% Unknown group fallback
g = string(T.GroupMouse);
g(ismissing(g)) = "Unknown";
T.GroupMouse = categorical(g, ["Active","Passive","Unknown"]);
end

function doMissingnessQC(T, outDir, vars)
mk = unique(string(T.mouse_key_norm),'stable');
days = unique(double(T.day_index)); days = sort(days(isfinite(days)));

for v=1:numel(vars)
    M = nan(numel(mk), numel(days));
    for i=1:numel(mk)
        for j=1:numel(days)
            r = (string(T.mouse_key_norm)==mk(i)) & (double(T.day_index)==days(j));
            if ~any(r), continue; end
            x = T.(vars{v})(r);
            M(i,j) = mean(~isfinite(double(x)));
        end
    end
    fig = figure('Color','w','Position',[80 80 1100 520]);
    imagesc(days, 1:numel(mk), M); colormap(parula); colorbar;
    xlabel('Day'); ylabel('Mouse');
    title(['Missing fraction: ' vars{v}], 'Interpreter','none');
    set(gca,'YTick',1:numel(mk),'YTickLabel',mk,'TickLabelInterpreter','none');
    exportgraphics(fig, fullfile(outDir, ['QC_missing_' vars{v} '.png']), 'Resolution', 200);
    close(fig);
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
boutStart = lickTimes(1); boutEnd = lickTimes(1);
for i=2:numel(lickTimes)
    if lickTimes(i)-lickTimes(i-1) > gap_s
        boutStart(end+1,1)=lickTimes(i); %#ok<AGROW>
        boutEnd(end+1,1)=lickTimes(i); %#ok<AGROW>
    else
        boutEnd(end)=lickTimes(i);
    end
end
end

function S = buildSessionSummary_FULL(T, timeVar, bout_gap_s)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

S = table;
S.mouse_key = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
S.day_index = splitapply(@(x) double(x(1)), double(T.day_index), G);
S.session_idx = splitapply(@(x) double(x(1)), double(T.session_idx), G);
S.GroupMouse = splitapply(@(x) x(1), T.GroupMouse, G);

S.pupil_mean_px = splitapply(@(x) mean(double(x),'omitnan'), T.Diameter_px, G);

% Lick-derived per-session metrics
n = height(S);
S.nLicks = nan(n,1);
S.licks_per_min = nan(n,1);
S.iei_median = nan(n,1);
S.bout_count = nan(n,1);
S.bout_mean_dur_s = nan(n,1);
S.total_lick_duration_s = nan(n,1);

for i=1:n
    r = (string(T.mouse_key_norm)==S.mouse_key(i)) & (double(T.day_index)==S.day_index(i)) & (double(T.session_idx)==S.session_idx(i));
    t  = double(T.(timeVar)(r));
    lk = double(T.Lick_TTL(r));

    ok = isfinite(t);
    t=t(ok); lk=lk(ok);
    if isempty(t), continue; end
    [t,ord]=sort(t); lk=lk(ord);

    lickTimes = detectRisingEdges(t, lk);
    S.nLicks(i) = numel(lickTimes);

    dur = max(t)-min(t);
    if isfinite(dur) && dur>1
        S.licks_per_min(i) = (numel(lickTimes)/dur)*60;
    end

    if numel(lickTimes)>=2
        S.iei_median(i) = median(diff(lickTimes));
    end

    [bs,be] = makeLickBouts(lickTimes, bout_gap_s);
    S.bout_count(i) = numel(bs);
    if ~isempty(bs)
        bd = be-bs;
        S.bout_mean_dur_s(i) = mean(bd,'omitnan');
        S.total_lick_duration_s(i) = sum(bd,'omitnan');
    else
        S.bout_mean_dur_s(i) = 0;
        S.total_lick_duration_s(i) = 0;
    end
end
end

function doOutlierQC(S, outDir)
metrics = {'licks_per_min','nLicks','iei_median','bout_count','bout_mean_dur_s','total_lick_duration_s','pupil_mean_px'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

for k=1:numel(metrics)
    x = double(S.(metrics{k}));
    [isOut, z] = madOutliers(x, 4.5);

    fig = figure('Color','w','Position',[80 80 950 420]); hold on
    scatter(1:numel(x), x, 18, 'filled','MarkerFaceAlpha',0.5);
    scatter(find(isOut), x(isOut), 55, 'o','LineWidth',1.8);
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
    xlabel('Day'); ylabel(yvar,'Interpreter','none');
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotCumulativeLickingAcrossSession_byDay(T, timeVar, outDir, dt_bin_s, excludeDays, groupName)
R = (string(T.GroupMouse)==string(groupName)) & ~ismember(double(T.day_index), excludeDays);
if ~any(R), return; end

days = unique(double(T.day_index(R))); days=sort(days);
fig = figure('Color','w','Position',[80 80 1150 520]); hold on
for di=1:numel(days)
    Td = T(R & double(T.day_index)==days(di),:);
    [tgrid, mu] = meanCumulativeLicks(Td, timeVar, dt_bin_s);
    if isempty(tgrid), continue; end
    plot(tgrid/60, mu, 'LineWidth', 1.8, 'DisplayName', sprintf('Day %d', days(di)));
end
xlabel('Session time (min)');
ylabel('Mean cumulative licks');
title(sprintf('Cumulative licking across session (%s mice) - day-by-day', groupName));
grid on; box on
legend('show','Location','eastoutside');
exportgraphics(fig, fullfile(outDir, sprintf('LICK_cumulative_byDay_%s.png', groupName)), 'Resolution', 220);
close(fig);
end

function plotCumulativeLicking_GroupOverlay(T, timeVar, outDir, dt_bin_s, excludeDays)
T0 = T(~ismember(double(T.day_index), excludeDays), :);

fig = figure('Color','w','Position',[80 80 1050 520]); hold on
for gname = ["Passive","Active"]
    R = (string(T0.GroupMouse)==gname);
    if ~any(R), continue; end
    [tgrid, mu, sem] = meanCumulativeLicks_SEM(T0(R,:), timeVar, dt_bin_s);
    if isempty(tgrid), continue; end
    if gname=="Passive"
        shadedSEM(tgrid/60, mu, sem, [0 0.4470 0.7410]); % blue-ish
        plot(tgrid/60, mu, 'LineWidth', 3, 'Color', [0 0.4470 0.7410], 'DisplayName','Passive');
    else
        shadedSEM(tgrid/60, mu, sem, [0.8500 0.3250 0.0980]); % red-ish
        plot(tgrid/60, mu, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980], 'DisplayName','Active');
    end
end
xlabel('Session time (min)');
ylabel('Mean cumulative licks');
title('Cumulative licking (days excluding 1-2): Passive (blue) vs Active (red)');
grid on; box on
legend('show','Location','best');
exportgraphics(fig, fullfile(outDir, 'LICK_cumulative_passive_vs_active_overlay.png'), 'Resolution', 220);
close(fig);
end

function [tgrid, mu] = meanCumulativeLicks(Tsub, timeVar, dt_bin_s)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(Tsub(:,keys));
nG = max(G);

curves = {};
tmaxAll = 0;
for gi=1:nG
    r = (G==gi);
    t = double(Tsub.(timeVar)(r));
    lk = double(Tsub.Lick_TTL(r));
    ok = isfinite(t);
    t=t(ok); lk=lk(ok);
    if isempty(t), continue; end
    [t,ord]=sort(t); lk=lk(ord);

    lt = detectRisingEdges(t, lk);
    t0 = min(t);
    dur = max(t)-t0;
    tmaxAll = max(tmaxAll, dur);
    lt = lt - t0;
    curves{end+1} = lt; %#ok<AGROW>
end

if isempty(curves), tgrid=[]; mu=[]; return; end
tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;

cumMat = nan(numel(curves), numel(tgrid));
for i=1:numel(curves)
    lt = curves{i};
    if isempty(lt)
        cumMat(i,:) = 0;
    else
        cumMat(i,:) = arrayfun(@(tt) sum(lt<=tt), tgrid);
    end
end
mu = mean(cumMat,1,'omitnan');
end

function [tgrid, mu, sem] = meanCumulativeLicks_SEM(Tsub, timeVar, dt_bin_s)
[tgrid, mu] = meanCumulativeLicks(Tsub, timeVar, dt_bin_s);
if isempty(tgrid), sem=[]; return; end

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(Tsub(:,keys));
nG = max(G);

curves = {};
tmaxAll = 0;
for gi=1:nG
    r=(G==gi);
    t = double(Tsub.(timeVar)(r));
    lk = double(Tsub.Lick_TTL(r));
    ok = isfinite(t);
    t=t(ok); lk=lk(ok);
    if isempty(t), continue; end
    [t,ord]=sort(t); lk=lk(ord);
    lt = detectRisingEdges(t, lk);
    t0=min(t);
    dur=max(t)-t0;
    tmaxAll=max(tmaxAll,dur);
    curves{end+1}=lt-t0; %#ok<AGROW>
end
tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;

cumMat = nan(numel(curves), numel(tgrid));
for i=1:numel(curves)
    lt=curves{i};
    if isempty(lt), cumMat(i,:)=0;
    else, cumMat(i,:)=arrayfun(@(tt) sum(lt<=tt), tgrid);
    end
end

nEff = sum(any(isfinite(cumMat),2));
if nEff<=1
    sem = zeros(size(mu));
else
    sem = std(cumMat,0,1,'omitnan') ./ sqrt(nEff);
end
end

function shadedSEM(x, mu, sem, colorRGB)
x = x(:)'; mu=mu(:)'; sem=sem(:)';
fill([x fliplr(x)], [mu-sem fliplr(mu+sem)], colorRGB, ...
    'EdgeColor','none','FaceAlpha',0.20);
end

function doPhaseComparisons_withStats_FINAL(S, outDir, PER, DAYSETS, includePassiveDuring)
S.phase = assignPhaseFromDay(S.day_index, PER);

% Metrics to compare across phases
metrics = {'licks_per_min','nLicks','iei_median','bout_count','bout_mean_dur_s','total_lick_duration_s','pupil_mean_px','pupil_norm_rel'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

for dsName = fieldnames(DAYSETS)'
    ds = DAYSETS.(dsName{1});

    for m=1:numel(metrics)
        yvar = metrics{m};

        % Save stats table: Active vs Passive per phase
        phasesOrder = categories(S.phase);
        phasesOrder = phasesOrder(~strcmp(phasesOrder,'<undef>'));

        POUT = table;
        row=0;

        fig = figure('Color','w','Position',[80 80 1200 520]);
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        A_val=[]; A_ph=categorical([]);
        P_val=[]; P_ph=categorical([]);

        for gi=1:2
            gLabel = ["Active","Passive"];
            g = gLabel(gi);
            nexttile; hold on

            rG = (string(S.GroupMouse)==g) ...
                 & ds.keepDay(S.day_index) ...
                 & (S.phase~="<undef>");

            % Optional: drop Passive during if requested
            if g=="Passive" && ~includePassiveDuring
                rG = rG & (S.phase~="during");
            end

            if ~any(rG)
                title(sprintf('%s (%s) – no data', yvar, g));
                axis off; continue;
            end

            % Collapse to mouse x phase means
            [Gmp,~] = findgroups(string(S.mouse_key(rG)), S.phase(rG));
            val = splitapply(@(x) mean(double(x),'omitnan'), double(S.(yvar)(rG)), Gmp);
            ph  = splitapply(@(x) x(1), S.phase(rG), Gmp);

            % Plot
            boxchart(double(ph), val);
            set(gca,'XTick',1:numel(phasesOrder),'XTickLabel',phasesOrder);
            ylabel(yvar,'Interpreter','none');
            title(sprintf('%s (%s) – %s', yvar, g, ds.name));
            grid on; box on

            if gi==1
                A_val = val; A_ph = ph;
            else
                P_val = val; P_ph = ph;
            end
        end

        % Stats: ranksum per phase (robust; fixes your crash)
        for pi=1:numel(phasesOrder)
            phName = phasesOrder{pi};
            a = A_val(A_ph==phName);
            p = P_val(P_ph==phName);
            pval = safeRanksum(a,p);

            row=row+1;
            POUT.phase(row,1) = string(phName);
            POUT.metric(row,1)= string(yvar);
            POUT.dataset(row,1)=string(ds.name);
            POUT.p_ranksum_active_vs_passive(row,1)=pval;
            POUT.n_active(row,1)=sum(isfinite(a));
            POUT.n_passive(row,1)=sum(isfinite(p));
        end

        writetable(POUT, fullfile(outDir, sprintf('STATS_phase_%s_%s_active_vs_passive.csv', yvar, ds.name)));
        exportgraphics(fig, fullfile(outDir, sprintf('PHASE_box_%s_%s.png', yvar, ds.name)), 'Resolution', 220);
        close(fig);
    end
end
end

function p = safeRanksum(a,b)
a = a(:); b=b(:);
a = a(isfinite(a)); b=b(isfinite(b));
if isempty(a) || isempty(b)
    p = NaN;
    return;
end
p = ranksum(a,b);
end

function ph = assignPhaseFromDay(day_index, PER)
phase = strings(numel(day_index),1);
d = double(day_index);

phase(ismember(d, PER.pre))      = "pre";
phase(ismember(d, PER.during))   = "during";
phase(ismember(d, PER.post))     = "post";
phase(ismember(d, PER.withdraw)) = "withdrawal";
phase(ismember(d, PER.reexpo))   = "reexposure";
phase(phase=="") = "<undef>";

ph = categorical(phase, ...
    ["pre","during","post","withdrawal","reexposure","<undef>"], ...
    ["pre","during","post","withdrawal","reexposure","<undef>"]);
end

function S = addPupilBaselineNormalization_DAY3(S)
% Baseline = Day3 mean pupil per mouse (your request)
S.pupil_base_day3 = nan(height(S),1);
S.pupil_norm_rel  = nan(height(S),1);

mk = unique(string(S.mouse_key),'stable');
for i=1:numel(mk)
    rM = (string(S.mouse_key)==mk(i));
    base = mean(double(S.pupil_mean_px(rM & S.day_index==3)), 'omitnan');
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
        [d,ord]=sort(double(S.day_index(r)));
        y = double(S.pupil_norm_rel(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',1,'MarkerSize',4);
    end

    dAll = unique(double(S.day_index(rG))); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j)=mean(double(S.pupil_norm_rel(rG & double(S.day_index)==dAll(j))),'omitnan');
    end
    plot(dAll, mu, 'k-','LineWidth',3);
    yline(0,'k:');
    xlabel('Day'); ylabel('\Delta pupil / baseline(D3)');
    title(sprintf('Pupil normalized by Day3 baseline (%s)', g));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_daylevel_normalized_by_day3.png'), 'Resolution', 220);
close(fig);
end

function plotPupilMorphineVsWater_DayCompare(S, outDir, daysA, daysB)
% Example: water (days 4-5) vs morphine PR (days 12-13) using normalized pupil
if ~ismember('pupil_norm_rel', S.Properties.VariableNames), return; end

fig = figure('Color','w','Position',[80 80 900 420]);
tiledlayout(1,1,'TileSpacing','compact','Padding','compact'); nexttile; hold on

for gname = ["Active","Passive"]
    rG = (string(S.GroupMouse)==gname);
    A = meanByMouseOverDays(S, 'pupil_norm_rel', rG, daysA);
    B = meanByMouseOverDays(S, 'pupil_norm_rel', rG, daysB);

    % paired-ish scatter (not forcing pairing)
    x1 = repmat(1, numel(A), 1);
    x2 = repmat(2, numel(B), 1);
    scatter(x1 + 0.05*(gname=="Active"), A, 40, 'filled');
    scatter(x2 + 0.05*(gname=="Active"), B, 40, 'filled');

    % summary
    muA=mean(A,'omitnan'); seA=std(A,0,'omitnan')/sqrt(max(1,sum(isfinite(A))));
    muB=mean(B,'omitnan'); seB=std(B,0,'omitnan')/sqrt(max(1,sum(isfinite(B))));
    errorbar([1 2] + 0.05*(gname=="Active"), [muA muB], [seA seB], '-o', 'LineWidth',2, 'DisplayName', char(gname));
end

set(gca,'XTick',[1 2], 'XTickLabel',{sprintf('Days %s',rangeStr(daysA)), sprintf('Days %s',rangeStr(daysB))});
ylabel('Pupil (normalized to Day3)');
title('Pupil comparison (water vs morphine PR windows)');
grid on; box on
legend('show','Location','best');

exportgraphics(fig, fullfile(outDir, sprintf('PUPIL_compare_days_%s_vs_%s.png', rangeStr(daysA), rangeStr(daysB))), 'Resolution', 220);
close(fig);
end

function V = meanByMouseOverDays(S, varName, rowMask, days)
mk = unique(string(S.mouse_key(rowMask)),'stable');
V = nan(numel(mk),1);
for i=1:numel(mk)
    r = rowMask & string(S.mouse_key)==mk(i) & ismember(double(S.day_index), days);
    V(i) = mean(double(S.(varName)(r)),'omitnan');
end
V = V(isfinite(V));
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1), s=sprintf('%d-%d',v(1),v(end));
else, s=strjoin(string(v),',');
end
end

function E = buildEventLockedTable_FULL(T, timeVar, tAxis, bout_gap_s, minSep_s)
% Builds event-locked pupil traces:
% - reward events (Injector_TTL rising edges)
% - lick-bout events: bout start, classified rewarded vs unrewarded
%
% Normalization: baseline = mean pupil in (-2..0) window relative to event
% For lick bouts: baseline is 2s before FIRST lick, and window is relative to BOUT START
% (You can extend to include end-of-bout logic later if you want; this already matches your requested baseline idea.)

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));
nG = max(G);

rows = {};
for gi=1:nG
    r = (G==gi);
    t   = double(T.(timeVar)(r));
    lk  = double(T.Lick_TTL(r));
    rw  = double(T.Injector_TTL(r));
    pup = double(T.Diameter_px(r));

    ok = isfinite(t) & isfinite(pup);
    t=t(ok); lk=lk(ok); rw=rw(ok); pup=pup(ok);
    if isempty(t), continue; end
    [t,ord]=sort(t); lk=lk(ord); rw=rw(ord); pup=pup(ord);

    mk  = string(T.mouse_key_norm(find(r,1,'first')));
    day = double(T.day_index(find(r,1,'first')));
    ses = double(T.session_idx(find(r,1,'first')));
    grp = string(T.GroupMouse(find(r,1,'first')));

    rewardTimes = enforceMinSeparation(detectRisingEdges(t, rw), minSep_s);

    % Reward-locked
    if ~isempty(rewardTimes)
        tr = extractEventLocked(t, pup, rewardTimes, tAxis);
        if ~isempty(tr)
            rows(end+1,:) = {mk, day, ses, grp, "reward", "all", tr}; %#ok<AGROW>
        end
    end

    % Lick bouts (start)
    lickTimes = enforceMinSeparation(detectRisingEdges(t, lk), 0.02);
    if ~isempty(lickTimes)
        [bs,be] = makeLickBouts(lickTimes, bout_gap_s);

        isRewarded = false(size(bs));
        for bi=1:numel(bs)
            % if any reward occurs during bout or within 1s after bout end, call it rewarded
            isRewarded(bi) = any(rewardTimes >= bs(bi) & rewardTimes <= (be(bi)+1.0));
        end

        trR = extractEventLocked(t, pup, bs(isRewarded), tAxis);
        if ~isempty(trR)
            rows(end+1,:) = {mk, day, ses, grp, "lickBoutStart", "rewarded", trR}; %#ok<AGROW>
        end
        trU = extractEventLocked(t, pup, bs(~isRewarded), tAxis);
        if ~isempty(trU)
            rows(end+1,:) = {mk, day, ses, grp, "lickBoutStart", "unrewarded", trU}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    E = table;
else
    E = cell2table(rows, 'VariableNames', {'mouse_key','day_index','session_idx','GroupMouse','eventType','subType','trace'});
end
end

function trace = extractEventLocked(t, pup, eventTimes, tAxis)
if isempty(eventTimes), trace=[]; return; end
M = nan(numel(eventTimes), numel(tAxis));
pre = (tAxis < 0);
for i=1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    b  = mean(pi(pre), 'omitnan');  % baseline: 2s pre
    pi = pi - b;
    M(i,:) = pi;
end
% session-level mean trace
trace = mean(M, 1, 'omitnan');
end

function plotRewardLockedEarlyLate(E, tAxis, outDir, PER)
% “Reward-locked pupil” plot like your example:
% compare early(pre) vs during (or pre vs during) separately for active and passive
earlyDays = PER.pre;     % 3-5
lateDays  = PER.during;  % 6-10

fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi=1:2
    gLabel = ["Passive","Active"];  % match your typical figure layout (passive left)
    g = gLabel(gi);
    nexttile; hold on

    rG = (string(E.GroupMouse)==g) & (string(E.eventType)=="reward") & (string(E.subType)=="all");

    A = stackTraces(E.trace(rG & ismember(double(E.day_index), earlyDays)));
    B = stackTraces(E.trace(rG & ismember(double(E.day_index), lateDays)));

    h = gobjects(0);

    if ~isempty(A)
        muA=mean(A,1,'omitnan'); seA=std(A,0,1,'omitnan')/sqrt(size(A,1));
        shadedSEM(tAxis, muA, seA, [0.2 0.6 1.0]);
        h(end+1)=plot(tAxis, muA, 'LineWidth',2.5, 'DisplayName',sprintf('%s Early (days %s)', g, rangeStr(earlyDays)));
    end
    if ~isempty(B)
        muB=mean(B,1,'omitnan'); seB=std(B,0,1,'omitnan')/sqrt(size(B,1));
        shadedSEM(tAxis, muB, seB, [0.1 0.3 0.8]);
        h(end+1)=plot(tAxis, muB, '--', 'LineWidth',2.5, 'DisplayName',sprintf('%s Late (days %s)', g, rangeStr(lateDays)));
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil (%s only)', g));
    grid on; box on
    if ~isempty(h), legend(h,'Location','best'); end
end

exportgraphics(fig, fullfile(outDir,'PUPIL_reward_locked_pre_vs_during.png'), 'Resolution', 220);
close(fig);
end

function plotLickBoutLocked_RewardedVsUnrewarded(E, tAxis, outDir, PER)
% Focus days you care about for passive forced period vs later morphine PR:
% Example plot: Passive days 7-9 vs Passive days 12-13 (both for rewarded/unrewarded)
passiveFocusA = intersect(PER.during, [7 8 9]);
passiveFocusB = PER.post; % 11-13

fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for col=1:2
    if col==1
        daysUse = passiveFocusA;
        ttl = sprintf('Passive days %s (during)', rangeStr(daysUse));
    else
        daysUse = passiveFocusB;
        ttl = sprintf('Passive days %s (post)', rangeStr(daysUse));
    end

    nexttile; hold on
    rBase = (string(E.GroupMouse)=="Passive") & (string(E.eventType)=="lickBoutStart") & ismember(double(E.day_index), daysUse);

    R = stackTraces(E.trace(rBase & string(E.subType)=="rewarded"));
    U = stackTraces(E.trace(rBase & string(E.subType)=="unrewarded"));

    h = gobjects(0);
    if ~isempty(R)
        mu=mean(R,1,'omitnan'); se=std(R,0,1,'omitnan')/sqrt(size(R,1));
        shadedSEM(tAxis, mu, se, [0.2 0.8 0.2]);
        h(end+1)=plot(tAxis, mu, 'LineWidth',2.5, 'DisplayName','Rewarded bout');
    end
    if ~isempty(U)
        mu=mean(U,1,'omitnan'); se=std(U,0,1,'omitnan')/sqrt(size(U,1));
        shadedSEM(tAxis, mu, se, [0.7 0.7 0.7]);
        h(end+1)=plot(tAxis, mu, '--', 'LineWidth',2.5, 'DisplayName','Unrewarded bout');
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from lick-bout start (s)');
    ylabel('\Delta pupil (baseline-subtracted)');
    title(ttl);
    grid on; box on
    if ~isempty(h), legend(h,'Location','best'); end
end

exportgraphics(fig, fullfile(outDir,'PUPIL_lickbout_locked_rewarded_vs_unrewarded_PASSIVE.png'), 'Resolution', 220);
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

function doAssayPlotsAndStats(T, outDir)
% Auto-detect relevant columns and produce Active vs Passive day-level plots
cands = string(T.Properties.VariableNames);

assays = {};
% Straub candidates
assays = [assays, cands(contains(lower(cands),'straub'))];
% TST candidates
assays = [assays, cands(contains(lower(cands),'tst'))];
% Hot plate candidates
assays = [assays, cands(contains(lower(cands),'hot'))];
% Tail immersion latency
if any(cands=="Immersion_Latency_s"), assays=[assays, "Immersion_Latency_s"]; end

assays = unique(assays,'stable');
assays = assays(~ismember(assays, ["mouse_key","mouse_key_norm","day_index","session_idx"]));

if isempty(assays)
    fid=fopen(fullfile(outDir,'ASSAY_NOTE.txt'),'w');
    fprintf(fid,'No assay columns detected (STRAUB/TST/HOT/Immersion_Latency_s).\n');
    fclose(fid);
    return;
end

for i=1:numel(assays)
    varName = assays(i);
    if ~isnumeric(T.(varName)) && ~islogical(T.(varName))
        continue;
    end
    plotScalarByDayGroup(T, varName, outDir);
end
end

function plotScalarByDayGroup(T, varName, outDir)
% Collapses to mouse-day means, then plots mean±SEM by group across days.
G = findgroups(string(T.mouse_key_norm), double(T.day_index), T.GroupMouse);
mk = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
dy = splitapply(@(x) double(x(1)), double(T.day_index), G);
gp = splitapply(@(x) x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), double(T.(varName)), G);

S2 = table(mk, dy, gp, val, 'VariableNames', {'mouse','day','group','value'});
S2 = S2(S2.group~="Unknown",:);

fig = figure('Color','w','Position',[80 80 1050 420]); hold on
groups = categories(S2.group);

for gi=1:numel(groups)
    g = groups{gi};
    rG = (S2.group==g);
    days = unique(S2.day(rG)); days=sort(days);
    mu = nan(size(days)); se = nan(size(days));
    for j=1:numel(days)
        x = S2.value(rG & S2.day==days(j));
        x = x(isfinite(x));
        mu(j)=mean(x,'omitnan');
        se(j)=std(x,0,'omitnan')/sqrt(max(1,numel(x)));
    end
    errorbar(days, mu, se, '-o', 'LineWidth', 2, 'MarkerSize',6, 'DisplayName', g);
end

xlabel('Day'); ylabel(varName,'Interpreter','none');
title(sprintf('%s: Active vs Passive (mouse-day means)', varName), 'Interpreter','none');
grid on; box on
legend('show','Location','best');

exportgraphics(fig, fullfile(outDir, sprintf('ASSAY_%s_byDay_ActiveVsPassive.png', varName)), 'Resolution', 220);
close(fig);

% Basic per-day stats (Active vs Passive) using safeRanksum
daysAll = unique(S2.day); daysAll=sort(daysAll);
ST = table;
row=0;
for j=1:numel(daysAll)
    d = daysAll(j);
    a = S2.value(S2.day==d & S2.group=="Active");
    p = S2.value(S2.day==d & S2.group=="Passive");
    row=row+1;
    ST.day(row,1)=d;
    ST.p_ranksum_active_vs_passive(row,1)=safeRanksum(a,p);
    ST.n_active(row,1)=sum(isfinite(a));
    ST.n_passive(row,1)=sum(isfinite(p));
end
writetable(ST, fullfile(outDir, sprintf('ASSAY_STATS_%s_active_vs_passive_byDay.csv', varName)));
end
