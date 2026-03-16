function make_longitudinal_QC_v6()
% make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v6
%
% Reads:
%   <rootTry>/run_*/ALL_mice_longitudinal.csv   (latest run_* folder)
%
% REQUIRED columns:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1 numeric or logical), Injector_TTL (0/1),
%   Diameter_px (numeric), RequirementLast (numeric)
%   Time axis: CamTime_rel_s OR PlotTime_s_30fps  (one must exist)
%
% Requested updates implemented:
%   - Fix legend mismatch in event-locked plots
%   - Fix Straub plotting using STRAUB_Frames_Non_moving and STRAUB_Pct_Non_moving
%   - Add 0911_red breakdown (missing Lick_TTL and time corruption)
%   - Add bad-time report (time corruption vs TTL corruption)
%   - Cumulative licking with SEM, by phases:
%       pre=4-5, during=7-10, post=12-13, withdrawal=15-16, reexposure=17-18
%   - Phase boxplots + p-values (ttest2 + ranksum) saved to CSV

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

dt_bin_s = 60;                % 1-min bins for cumulative licking curves
bout_gap_s = 2.0;             % lick bout definition
pupil_win = [-2 2];
min_event_separation_s = 0.5; % seconds

% Phase day-sets (YOUR REQUEST)
PH.pre       = 4:5;
PH.during    = 7:10;
PH.post      = 12:13;
PH.withdraw  = 15:16;
PH.reexpo    = 17:18;

% For reward-locked pupil comparison figure (requested earlier)
earlyDays_pupil = 3:5;
lateDays_pupil  = 6:10;

% Days to exclude globally for licking metrics (if desired)
excludeDays_licking_global = [1 2];

% Debug target
DEBUG_MOUSE = "0911_red";

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

assertHasVars(T, {'mouse_key','day_index','session_idx'}, 'Top-level');
T.mouse_key = string(T.mouse_key);

% Time axis choice
timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end

REQ = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, [{'mouse_key','day_index','session_idx'}, REQ], 'Core longitudinal schema');

% Normalize mouse key + attach cohort meta
COH = buildNewCohortTable();
T = addCohortMeta(T, COH);

% Remove NaN day rows
T = T(~isnan(double(T.day_index)),:);

%% ===================== BAD TIME REPORT (before filtering) =====================
reportBadTime(T, timeVar, outDir, DEBUG_MOUSE);

%% ===================== REMOVE NON-FINITE TIME ROWS (log) =====================
tcol = double(T.(timeVar));
badTimeRows = ~isfinite(tcol);
if any(badTimeRows)
    fprintf('[WARN] Removing %d rows with non-finite %s\n', nnz(badTimeRows), timeVar);
    Tw = T(badTimeRows, {'mouse_key_norm','day_index','session_idx'});
    Tw.timeVar = repmat(string(timeVar), height(Tw),1);
    writetable(Tw, fullfile(outDir, sprintf('WARN_removed_nonfinite_%s_rows.csv', timeVar)));
    T = T(~badTimeRows,:);
end

%% ===================== QC: MISSINGNESS MAPS =====================
doMissingnessQC(T, outDir, {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'});

%% ===================== SESSION SUMMARY (per mouse-day-session) =====================
S = buildSessionSummary(T, timeVar);

%% ===================== BASIC QC PLOTS =====================
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'RequirementLast', outDir, 'QC_spaghetti_PRscore_RequirementLast.png');

%% ===================== CUMULATIVE LICKING with SEM by PHASE =====================
plotCumulativeLickingByPhaseWithSEM(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global, PH);

%% ===================== PHASE COMPARISONS (BOXPLOTS + STATS) =====================
phaseBoxplotsAndStats(S, outDir, PH, 'licks_per_min');

%% ===================== EVENT-LOCKED PUPIL (legend fixed) =====================
doEventLockedPupilAnalyses(T, timeVar, outDir, pupil_win, bout_gap_s, min_event_separation_s, ...
    earlyDays_pupil, lateDays_pupil);

%% ===================== STRAUB + HOTPLATE + TST + TAIL IMMERSION =====================
doOptionalPainAndBehaviorTests_v6(T, outDir);

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
COH.mouse_key_norm = COH.mouse_key;
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
if isnumeric(x) || islogical(x)
    f = mean(~isfinite(double(x)));
elseif isstring(x)
    f = mean(ismissing(x));
elseif iscell(x)
    xs = string(x);
    f = mean(ismissing(xs));
else
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
S.GroupMouse = splitapply(@(x) x(1), T.GroupMouse, G);

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

function plotCumulativeLickingByPhaseWithSEM(T, timeVar, outDir, dt_bin_s, excludeDays, PH)
% Mean cumulative lick curves (mean ± SEM across mice) per PHASE, separately for Active/Passive.

phaseNames = fieldnames(PH);
for gname = ["Active","Passive"]
    for pi=1:numel(phaseNames)
        phName = phaseNames{pi};
        daysPH = PH.(phName);

        R = (string(T.GroupMouse)==gname) & ismember(double(T.day_index), daysPH) & ~ismember(double(T.day_index), excludeDays);
        if ~any(R), continue; end

        [tgrid, mu, sem, nMice] = meanCumLicksAcrossMice(T(R,:), timeVar, dt_bin_s);

        if isempty(tgrid), continue; end
        fig = figure('Color','w','Position',[80 80 900 520]); hold on
        shadedSEM(tgrid/60, mu, sem);
        plot(tgrid/60, mu, 'LineWidth', 2.5);

        xlabel('Session time (min)');
        ylabel('Cumulative licks');
        title(sprintf('Cumulative licking (mean \\pm SEM): %s – %s (n=%d mice)', gname, phName, nMice), 'Interpreter','none');
        grid on; box on
        exportgraphics(fig, fullfile(outDir, sprintf('LICK_cumulative_%s_%s_meanSEM.png', gname, phName)), 'Resolution', 220);
        close(fig);
    end
end
end

function [tgrid, mu, sem, nMice] = meanCumLicksAcrossMice(Tsub, timeVar, dt_bin_s)
% For each mouse: build a mean cumulative curve across its sessions in Tsub.
% Then compute mean ± SEM across mice.

mice = unique(string(Tsub.mouse_key_norm),'stable');
mouseCurves = [];
tgrid = [];

for i=1:numel(mice)
    rM = (string(Tsub.mouse_key_norm)==mice(i));
    Tm = Tsub(rM,:);

    [tgrid_i, muCum_sessionMean] = meanCumulativeLicksAcrossSessions(Tm, timeVar, dt_bin_s);
    if isempty(tgrid_i), continue; end

    if isempty(tgrid)
        tgrid = tgrid_i;
        mouseCurves = muCum_sessionMean(:)';
    else
        % align by interpolation to common grid if needed
        if numel(tgrid_i) ~= numel(tgrid) || any(tgrid_i(:) ~= tgrid(:))
            muCum_sessionMean = interp1(tgrid_i, muCum_sessionMean, tgrid, 'linear', nan);
        end
        mouseCurves(end+1,:) = muCum_sessionMean(:)'; %#ok<AGROW>
    end
end

if isempty(tgrid) || isempty(mouseCurves)
    tgrid=[]; mu=[]; sem=[]; nMice=0; return;
end

mouseCurves = mouseCurves(any(isfinite(mouseCurves),2),:);
nMice = size(mouseCurves,1);

mu  = mean(mouseCurves, 1, 'omitnan');
sem = std(mouseCurves, 0, 1, 'omitnan') ./ sqrt(max(1,nMice));
end

function [tgrid, muCum] = meanCumulativeLicksAcrossSessions(Tm, timeVar, dt_bin_s)
% Compute session-level cumulative curves, then average across sessions (within a mouse)

keys = {'day_index','session_idx'};
G = findgroups(Tm(:,keys));

sessionCurves = {};
tmaxAll = 0;

for gi=1:max(G)
    r = (G==gi);
    t = double(Tm.(timeVar)(r));
    lk = double(Tm.Lick_TTL(r));

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

function shadedSEM(x, mu, sem)
x=x(:)'; mu=mu(:)'; sem=sem(:)';
fill([x fliplr(x)], [mu-sem fliplr(mu+sem)], 0.9*[1 1 1], 'EdgeColor','none', 'FaceAlpha',0.35);
end

function phaseBoxplotsAndStats(S, outDir, PH, metricName)
% Boxplots by phase for Active and Passive, and p-values Active vs Passive per phase.
if ~ismember(metricName, S.Properties.VariableNames)
    warning('phaseBoxplotsAndStats: missing metric %s', metricName);
    return;
end

% Assign phase label by day
phase = strings(height(S),1);
phase(ismember(double(S.day_index), PH.pre))      = "pre";
phase(ismember(double(S.day_index), PH.during))   = "during";
phase(ismember(double(S.day_index), PH.post))     = "post";
phase(ismember(double(S.day_index), PH.withdraw)) = "withdrawal";
phase(ismember(double(S.day_index), PH.reexpo))   = "reexposure";
phase(phase=="") = "<undef>";

S.phase = categorical(phase, ["pre","during","post","withdrawal","reexposure","<undef>"]);

% Per-mouse mean within phase (mouse-day means already in S; we average across days inside phase)
R = S(S.phase~="<undef>" & S.GroupMouse~="Unknown", :);

G = findgroups(string(R.mouse_key), R.GroupMouse, R.phase);
mouse = splitapply(@(x) string(x(1)), string(R.mouse_key), G);
grp   = splitapply(@(x) x(1), R.GroupMouse, G);
ph    = splitapply(@(x) x(1), R.phase, G);
val   = splitapply(@(x) mean(double(x),'omitnan'), R.(metricName), G);

M = table(mouse, grp, ph, val, 'VariableNames', {'mouse','group','phase','value'});

% Plot: two panels Active vs Passive boxplots
fig = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi=1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);
    nexttile; hold on
    rG = (string(M.group)==g);
    if ~any(rG), axis off; title(sprintf('%s: no data', g)); continue; end

    phasesOrder = categories(M.phase);
    phasesOrder = phasesOrder(~strcmp(phasesOrder,"<undef>"));

    % boxchart needs numeric x; map categorical to 1..K
    x = double(M.phase(rG));
    y = M.value(rG);
    boxchart(x, y);
    set(gca,'XTick',1:numel(phasesOrder),'XTickLabel',phasesOrder);
    ylabel(metricName,'Interpreter','none');
    title(sprintf('%s – mouse means by phase', g), 'Interpreter','none');
    grid on; box on
end

exportgraphics(fig, fullfile(outDir, sprintf('PHASE_box_%s_mouseMeans.png', metricName)), 'Resolution', 220);
close(fig);

% Stats: Active vs Passive per phase (mouse-level)
phasesOrder = categories(M.phase);
phasesOrder = phasesOrder(~strcmp(phasesOrder,"<undef>"));

statsRows = [];
for i=1:numel(phasesOrder)
    phName = phasesOrder{i};
    a = M.value(M.group=="Active"  & string(M.phase)==phName);
    p = M.value(M.group=="Passive" & string(M.phase)==phName);

    a = a(isfinite(a)); p = p(isfinite(p));
    if numel(a)<2 || numel(p)<2
        t_p = NaN; r_p = NaN;
    else
        try
            [~,t_p] = ttest2(a,p);
        catch
            t_p = NaN;
        end
        try
            r_p = ranksum(a,p);
        catch
            r_p = NaN;
        end
    end

    statsRows = [statsRows; {phName, numel(a), numel(p), mean(a,'omitnan'), mean(p,'omitnan'), t_p, r_p}]; %#ok<AGROW>
end

ST = cell2table(statsRows, 'VariableNames', ...
    {'phase','nActive','nPassive','meanActive','meanPassive','p_ttest2','p_ranksum'});

writetable(ST, fullfile(outDir, sprintf('STATS_phase_Active_vs_Passive_%s.csv', metricName)));
end

function doEventLockedPupilAnalyses(T, timeVar, outDir, win_s, bout_gap_s, minSep_s, earlyDays, lateDays)
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

        addRow("lickBoutStart","rewarded",  boutStart(isRewarded));
        addRow("lickBoutStart","nonreward", boutStart(~isRewarded));
        addRow("lickBoutEnd","rewarded",    boutEnd(isRewarded));
        addRow("lickBoutEnd","nonreward",   boutEnd(~isRewarded));
    end

    function addRow(ev, sub, times)
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

plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir, earlyDays, lateDays);
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
if isempty(lickTimes), boutStart=[]; boutEnd=[]; return; end
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

function plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir, earlyDays, lateDays)
eventType = "reward";
subType   = "all";

% Orange = early; Blue = late
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
        shadedSEM(tAxis, muA, seA);
        hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', cEarly, ...
            'DisplayName', sprintf('Days %d-%d', earlyDays(1), earlyDays(end)));
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shadedSEM(tAxis, muB, seB);
        hB = plot(tAxis, muB, 'LineWidth', 2.5, 'Color', cLate, ...
            'DisplayName', sprintf('Days %d-%d', lateDays(1), lateDays(end)));
    end

    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil (%s)', g));

    % Legend bound to line handles only (no confusion)
    if ~isempty(hA) && ~isempty(hB)
        legend([hA hB],'Location','best');
    elseif ~isempty(hA)
        legend(hA,'Location','best');
    elseif ~isempty(hB)
        legend(hB,'Location','best');
    end
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_reward_locked_early_vs_late_active_passive_FIXEDLEGEND.png'), 'Resolution', 220);
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

function doOptionalPainAndBehaviorTests_v6(T, outDir)
% Hotplate / TST / Tail immersion kept; Straub explicitly fixed.

if ismember('Immersion_Latency_s', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'Immersion_Latency_s', outDir, 'TEST_tailImmersion_latency.png');
end
if ismember('HOT_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'HOT_Frames_Non_moving', outDir, 'TEST_hotPlate_nonmoving_frames.png');
end
if ismember('TST_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'TST_Frames_Non_moving', outDir, 'TEST_TST_nonmoving_frames.png');
end

% STRAUB FIX
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

% non-monotonic time flag
nonMono = splitapply(@(x) any(diff(double(x(isfinite(double(x)))))<0), time, G);

% dt jump heuristic: p99 dt
p99dt = splitapply(@(x) prctile(diff(double(x(isfinite(double(x))))),99), time, G);
p99dt(~isfinite(p99dt)) = NaN;

R = table(day, ses, nRows, nLickNaN, fracLickNaN, nTimeBad, fracTimeBad, nonMono, p99dt);
R = sortrows(R, {'fracLickNaN','fracTimeBad','nonMono','p99dt'},{'descend','descend','descend','descend'});

writetable(R, fullfile(outDir, sprintf('DEBUG_%s_missing_LickTTL_by_day_session.csv', mk)));

fprintf('\n[DEBUG] %s: top day/session with missing Lick_TTL or bad time\n', mk);
disp(R(1:min(20,height(R)),:));

% day-level plot: max frac missing lick per day
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
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));

mk  = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
day = splitapply(@(x) double(x(1)), double(T.day_index), G);
ses = splitapply(@(x) double(x(1)), double(T.session_idx), G);

tt  = splitapply(@(x) {double(x)}, double(T.(timeVar)), G);

nRows   = cellfun(@numel, tt);
nBad    = cellfun(@(x) nnz(~isfinite(x)), tt);
fracBad = nBad ./ max(1,nRows);

nonMono = cellfun(@(x) any(diff(x(isfinite(x)))<0), tt);

% p99 dt as jump indicator
p99dt = cellfun(@(x) localP99dt(x), tt);
p99dt(~isfinite(p99dt)) = NaN;

R = table(mk, day, ses, nRows, nBad, fracBad, nonMono, p99dt, ...
    'VariableNames', {'mouse','day','session','nRows','nNonFinite','fracNonFinite','nonMonotonic','p99_dt'});

writetable(R, fullfile(outDir, sprintf('QC_badTime_report_%s.csv', timeVar)));

dm = lower(string(debugMouse));
rD = lower(R.mouse) == dm;
if any(rD)
    writetable(R(rD,:), fullfile(outDir, sprintf('QC_badTime_DEBUG_%s_%s.csv', dm, timeVar)));
end
end

function p = localP99dt(x)
x = x(isfinite(x));
if numel(x) < 3, p = NaN; return; end
dt = diff(x);
dt = dt(isfinite(dt));
if isempty(dt), p = NaN; return; end
p = prctile(dt,99);
end
