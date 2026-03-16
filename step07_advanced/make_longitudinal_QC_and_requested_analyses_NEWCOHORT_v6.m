function make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v6()
% make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v6
% Full revised v6: fixes function nesting errors + adds SEM + phase windows + stats + Straub + debug.
%
% Required columns in ALL_mice_longitudinal.csv:
%   mouse_key (string/cellstr), day_index (numeric), session_idx (numeric)
%   Lick_TTL (0/1 numeric or logical), Injector_TTL (0/1), Diameter_px (numeric), RequirementLast (numeric)
%   Time axis: CamTime_rel_s OR PlotTime_s_30fps (either one must exist)

%% ===================== USER SETTINGS =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

dt_bin_s = 60;                 % bin for cumulative curve
pupil_win = [-2 2];
bout_gap_s = 2.0;
min_event_separation_s = 0.5;

% You requested these phase windows (use these everywhere)
PER.pre       = 4:5;
PER.during    = 7:10;
PER.post      = 12:13;
PER.withdraw  = 15:16;
PER.reexpo    = 17:18;

% If you still want to define "transition" days, keep it separate.
transitionDays = [6 11 14];    % optional, used only for "NO_TRANSITION"

excludeDays_licking_global = [1 2];  % optional exclude

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

timeVar = '';
if ismember('CamTime_rel_s', T.Properties.VariableNames)
    timeVar = 'CamTime_rel_s';
elseif ismember('PlotTime_s_30fps', T.Properties.VariableNames)
    timeVar = 'PlotTime_s_30fps';
else
    error('Missing time column: need CamTime_rel_s OR PlotTime_s_30fps');
end

REQ = {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'};
assertHasVars(T, [{'mouse_key','day_index','session_idx'}, REQ], 'Core schema');

% cohort meta (NEW cohort mapping)
COH = buildNewCohortTable();
T = addCohortMeta(T, COH);

% Drop NaN day_index
T = T(~isnan(double(T.day_index)),:);

%% ===================== BAD TIME CHECK (before filtering) =====================
reportBadTime(T, timeVar, outDir, DEBUG_MOUSE);

%% ===================== REMOVE NON-FINITE TIME ROWS (log them) =====================
tcol = double(T.(timeVar));
badTimeRows = ~isfinite(tcol);
if any(badTimeRows)
    fprintf('[WARN] Removing %d rows with non-finite %s\n', nnz(badTimeRows), timeVar);
    Tw = T(badTimeRows, {'mouse_key_norm','day_index','session_idx'});
    Tw.timeVar = repmat(string(timeVar), height(Tw),1);
    writetable(Tw, fullfile(outDir, sprintf('WARN_removed_nonfinite_%s_rows.csv', timeVar)));
    T = T(~badTimeRows,:);
end

%% ===================== QC: MISSINGNESS HEATMAP =====================
doMissingnessQC(T, outDir, {'Lick_TTL','Injector_TTL','Diameter_px','RequirementLast'});

%% ===================== SESSION SUMMARY =====================
S = buildSessionSummary(T, timeVar);

%% ===================== OUTLIER QC =====================
doOutlierQC(S, outDir);

%% ===================== BASIC SPAGHETTI QC =====================
plotSpaghettiByGroup(S, 'licks_per_min', outDir, 'QC_spaghetti_licks_per_min.png');
plotSpaghettiByGroup(S, 'RequirementLast', outDir, 'QC_spaghetti_PRscore_RequirementLast.png');
plotPupilByGroupDay(S, outDir);

%% ===================== CUMULATIVE LICKING (with SEM) =====================
% 1) across session by phase windows (Active and Passive separately)
plotCumulativeLickingByPhases_SEM(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global, PER);

% (optional) also keep the older day-by-day version if you want:
% plotCumulativeLickingAcrossSession_SEM_byDay(T, timeVar, outDir, dt_bin_s, excludeDays_licking_global);

%% ===================== PERIOD/PHASE COMPARISONS + STATS =====================
DAYSETS.allDays.name  = 'ALLDAYS';
DAYSETS.allDays.mask  = true(height(S),1);

DAYSETS.noTransition.name = 'NO_TRANSITION';
DAYSETS.noTransition.mask = ~ismember(double(S.day_index), transitionDays);

doPhaseComparisons_withStats(S, outDir, PER, DAYSETS);

%% ===================== PUPIL =====================
S = addPupilBaselineNormalization(S);
plotPupilNormalizedDayTraces(S, outDir);

doEventLockedPupilAnalyses_v6(T, timeVar, outDir, pupil_win, bout_gap_s, min_event_separation_s);

%% ===================== STRAUB + HOTPLATE + TST + TAIL IMMERSION =====================
doOptionalPainAndBehaviorTests_v6(T, outDir);

%% ===================== DEBUG: missing lick TTL breakdown =====================
debugMissingLickTTL(T, timeVar, outDir, DEBUG_MOUSE);

fprintf('\nDONE. Outputs saved in:\n  %s\n', outDir);
end

%% =====================================================================
%% =========================== LOCAL FUNCTIONS ==========================
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

function lickTimes = detectRisingEdges(t, x01)
x01 = x01(:);
x01(~isfinite(x01)) = 0;
x01 = x01 > 0.5;
dx = diff([false; x01]);
idx = find(dx==1);
lickTimes = t(idx);
lickTimes = lickTimes(isfinite(lickTimes));
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

    ok = isfinite(t);
    t = t(ok); lk = lk(ok);
    if isempty(t), continue; end

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

function plotCumulativeLickingByPhases_SEM(T, timeVar, outDir, dt_bin_s, excludeDays, PER)
% Cumulative licking curves with SEM for each phase window, per group.
phaseNames = ["pre","during","post","withdrawal","reexposure"];
phaseDays = {PER.pre, PER.during, PER.post, PER.withdraw, PER.reexpo};

T0 = T(~ismember(double(T.day_index), excludeDays), :);

for gname = ["Active","Passive"]
    Rg = (string(T0.GroupMouse)==gname);
    if ~any(Rg), continue; end

    fig = figure('Color','w','Position',[80 80 1200 650]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for pi=1:numel(phaseNames)
        nexttile; hold on
        daysThis = phaseDays{pi};
        R = Rg & ismember(double(T0.day_index), daysThis);
        if ~any(R)
            title(sprintf('%s (%s): no data', phaseNames(pi), gname));
            axis off;
            continue;
        end

        [tgrid, mu, sem] = meanCumulativeLicks_SEM(T0(R,:), timeVar, dt_bin_s);
        if isempty(tgrid)
            title(sprintf('%s (%s): no curves', phaseNames(pi), gname));
            axis off;
            continue;
        end

        shadedSEM(tgrid/60, mu, sem);
        plot(tgrid/60, mu, 'LineWidth', 2.5);

        xlabel('Session time (min)');
        ylabel('Cumulative licks');
        title(sprintf('%s days %s (%s)', phaseNames(pi), rangeStr(daysThis), gname), 'Interpreter','none');
        grid on; box on
    end

    exportgraphics(fig, fullfile(outDir, sprintf('LICK_cumulative_byPhase_SEM_%s.png', gname)), 'Resolution', 220);
    close(fig);
end
end

function [tgrid, mu, sem] = meanCumulativeLicks_SEM(Tsub, timeVar, dt_bin_s)
% Build per-session cumulative lick curves aligned to session start, then average + SEM.
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
    t = t(ok); lk = lk(ok);
    if isempty(t), continue; end

    [t,ord] = sort(t); lk = lk(ord);
    lt = detectRisingEdges(t, lk);

    t0 = min(t);
    dur = max(t) - t0;
    tmaxAll = max(tmaxAll, dur);

    lt = lt - t0; % align
    curves{end+1} = lt; %#ok<AGROW>
end

if isempty(curves)
    tgrid=[]; mu=[]; sem=[]; return;
end

tgrid = 0:dt_bin_s:ceil(tmaxAll/dt_bin_s)*dt_bin_s;
cumMat = nan(numel(curves), numel(tgrid));

for i=1:numel(curves)
    lt = curves{i};
    if isempty(lt)
        cumMat(i,:) = 0;
    else
        cumMat(i,:) = arrayfun(@(tt) sum(lt <= tt), tgrid);
    end
end

mu = mean(cumMat, 1, 'omitnan');
nEff = sum(any(isfinite(cumMat),2));
if nEff <= 1
    sem = zeros(size(mu));
else
    sem = std(cumMat, 0, 1, 'omitnan') ./ sqrt(nEff);
end
end

function shadedSEM(x, mu, sem)
x = x(:)'; mu = mu(:)'; sem = sem(:)';
fill([x fliplr(x)], [mu-sem fliplr(mu+sem)], 0.9*[1 1 1], ...
    'EdgeColor','none', 'FaceAlpha',0.35);
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1)
    s = sprintf('%d-%d', v(1), v(end));
else
    s = strjoin(string(v), ',');
end
end

function doPhaseComparisons_withStats(S, outDir, PER, DAYSETS)
% Creates phase labels using your windows and then:
% - boxplots (mouse means per phase)
% - stats: Active vs Passive within each phase using ranksum

% Assign phase
phase = strings(height(S),1);
phase(ismember(double(S.day_index), PER.pre))      = "pre";
phase(ismember(double(S.day_index), PER.during))   = "during";
phase(ismember(double(S.day_index), PER.post))     = "post";
phase(ismember(double(S.day_index), PER.withdraw)) = "withdrawal";
phase(ismember(double(S.day_index), PER.reexpo))   = "reexposure";
phase(phase=="") = "<undef>";

S.phase = categorical(phase, ...
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

        POUT = table; prow = 0;

        for gi=1:2
            gLabel = ["Active","Passive"];
            g = gLabel(gi);
            nexttile; hold on

            Rg = (string(S.GroupMouse)==g) & ds.mask & (S.phase~="<undef>");
            if ~any(Rg)
                title(sprintf('%s (%s) – no data', yvar, g));
                axis off; continue;
            end

            % collapse to mouse x phase means
            [Gmp,~] = findgroups(string(S.mouse_key(Rg)), S.phase(Rg));
            val = splitapply(@(x) mean(double(x),'omitnan'), S.(yvar)(Rg), Gmp);
            ph  = splitapply(@(x) x(1), S.phase(Rg), Gmp);

            boxchart(double(ph), val);
            phasesOrder = categories(S.phase);
            phasesOrder = phasesOrder(~strcmp(phasesOrder,'<undef>'));
            set(gca,'XTick',1:numel(phasesOrder),'XTickLabel',phasesOrder);
            ylabel(yvar,'Interpreter','none');
            title(sprintf('%s (%s) – %s', yvar, g, ds.name));
            grid on; box on

            % Store for between-group test
            if gi==1
                A_val = val; A_ph = ph;
            else
                P_val = val; P_ph = ph;
            end
        end

        % Between-group p-values per phase (ranksum)
        phasesOrder = categories(S.phase);
        phasesOrder = phasesOrder(~strcmp(phasesOrder,'<undef>'));

        for pi=1:numel(phasesOrder)
            phName = phasesOrder{pi};
            a = A_val(A_ph==phName);
            p = P_val(P_ph==phName);
            if numel(a)>=2 && numel(p)>=2
                pval = ranksum(a,p);
            else
                pval = NaN;
            end
            prow = prow + 1;
            POUT.phase(prow,1) = string(phName);
            POUT.metric(prow,1) = string(yvar);
            POUT.dataset(prow,1) = string(ds.name);
            POUT.p_ranksum_active_vs_passive(prow,1) = pval;
            POUT.n_active(prow,1) = numel(a);
            POUT.n_passive(prow,1) = numel(p);
        end

        writetable(POUT, fullfile(outDir, sprintf('STATS_phase_%s_%s_active_vs_passive.csv', yvar, ds.name)));

        exportgraphics(fig, fullfile(outDir, sprintf('PHASE_box_%s_%s.png', yvar, ds.name)), 'Resolution', 220);
        close(fig);
    end
end
end

function S = addPupilBaselineNormalization(S)
S.pupil_base_day4 = nan(height(S),1);
S.pupil_norm_rel  = nan(height(S),1);

mk = unique(string(S.mouse_key),'stable');
for i=1:numel(mk)
    rM = (string(S.mouse_key)==mk(i));
    base = mean(double(S.pupil_mean_px(rM & ismember(double(S.day_index),4))), 'omitnan');
    if ~isfinite(base) || base<=0, continue; end
    S.pupil_base_day4(rM) = base;
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
    xlabel('Day'); ylabel('\Delta pupil / baseline(D4)');
    title(sprintf('Pupil normalized by Day4 baseline (%s)', g));
    grid on; box on
end

exportgraphics(fig, fullfile(outDir,'PUPIL_daylevel_normalized_by_day4.png'), 'Resolution', 220);
close(fig);
end

function doEventLockedPupilAnalyses_v6(T, timeVar, outDir, win_s, bout_gap_s, minSep_s)
% Minimal: fixes legend by binding to handles. Reward-locked early vs late example kept.
tAxis = win_s(1):0.1:win_s(2);
E = buildEventLockedTable(T, timeVar, tAxis, bout_gap_s, minSep_s);

if isempty(E)
    warning('No event-locked pupil traces computed.');
    return;
end

plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir);
end

function E = buildEventLockedTable(T, timeVar, tAxis, bout_gap_s, minSep_s)
keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));
nG = max(G);

rows = {};
for gi=1:nG
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
        tr = extractEventLocked(t, pup, rewardTimes, tAxis);
        if ~isempty(tr)
            rows(end+1,:) = {mk, day, ses, grp, "reward", "all", tr}; %#ok<AGROW>
        end
    end

    lickTimes = enforceMinSeparation(detectRisingEdges(t, lk), 0.02);
    if ~isempty(lickTimes)
        [boutStart, boutEnd] = makeLickBouts(lickTimes, bout_gap_s);
        isRewarded = false(size(boutStart));
        for bi=1:numel(boutStart)
            isRewarded(bi) = any(rewardTimes>=boutStart(bi) & rewardTimes<=boutEnd(bi)+1.0);
        end

        tr1 = extractEventLocked(t, pup, boutStart(isRewarded), tAxis);
        if ~isempty(tr1)
            rows(end+1,:) = {mk, day, ses, grp, "lickBoutStart", "rewarded", tr1}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    E = table;
else
    E = cell2table(rows, 'VariableNames', {'mouse_key','day_index','session_idx','GroupMouse','eventType','subType','trace'});
end
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

function plotRewardLockedEarlyLate_FIXEDLEGEND(E, tAxis, outDir)
eventType = "reward";
subType   = "all";
earlyDays = 4:5;      % match your new phase definition
lateDays  = 7:10;

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
        hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', cEarly, 'DisplayName', sprintf('Days %s',rangeStr(earlyDays)));
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shadedSEM(tAxis, muB, seB);
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

function doOptionalPainAndBehaviorTests_v6(T, outDir)
if ismember('Immersion_Latency_s', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'Immersion_Latency_s', outDir, 'TEST_tailImmersion_latency.png');
end
if ismember('HOT_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'HOT_Frames_Non_moving', outDir, 'TEST_hotPlate_nonmoving_frames.png');
end
if ismember('TST_Frames_Non_moving', T.Properties.VariableNames)
    plotScalarByDayGroup(T, 'TST_Frames_Non_moving', outDir, 'TEST_TST_nonmoving_frames.png');
end

% Requested: Straub explicit metrics
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
G = findgroups(string(T.mouse_key_norm), double(T.day_index), T.GroupMouse);
mk = splitapply(@(x) string(x(1)), string(T.mouse_key_norm), G);
dy = splitapply(@(x) double(x(1)), double(T.day_index), G);
gp = splitapply(@(x) x(1), T.GroupMouse, G);
val = splitapply(@(x) mean(double(x),'omitnan'), double(T.(varName)), G);

S2 = table(mk, dy, gp, val, 'VariableNames', {'mouse','day','group','value'});

fig = figure('Color','w','Position',[80 80 950 420]); hold on
groups = categories(S2.group);
groups = groups(~strcmp(groups,'Unknown'));

for gi=1:numel(groups)
    g = groups{gi};
    rG = (S2.group==g);
    days = unique(S2.day(rG)); days=sort(days);
    mu = nan(size(days)); se = nan(size(days));
    for j=1:numel(days)
        x = S2.value(rG & S2.day==days(j));
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

nRows       = splitapply(@numel, lick, G);
nLickNaN    = splitapply(@(x) nnz(~isfinite(double(x))), lick, G);
fracLickNaN = nLickNaN ./ max(1,nRows);

nTimeBad    = splitapply(@(x) nnz(~isfinite(double(x))), time, G);
fracTimeBad = nTimeBad ./ max(1,nRows);

nonMono = splitapply(@(x) any(diff(double(x(isfinite(double(x)))))<0), time, G);

R = table(day, ses, nRows, nLickNaN, fracLickNaN, nTimeBad, fracTimeBad, nonMono);
R = sortrows(R, {'fracLickNaN','fracTimeBad','nonMono'},{'descend','descend','descend'});

writetable(R, fullfile(outDir, sprintf('DEBUG_%s_missing_LickTTL_by_day_session.csv', mk)));

fprintf('\n[DEBUG] %s: top day/session with missing Lick_TTL or bad time\n', mk);
disp(R(1:min(30,height(R)),:));

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

nRows = cellfun(@numel, tt);
nBad  = cellfun(@(x) nnz(~isfinite(x)), tt);
nonMono = cellfun(@(x) any(diff(x(isfinite(x)))<0), tt);

p99dt = cellfun(@(x) safeP99dt(x), tt);

R = table(mk, day, ses, nRows, nBad, nBad./max(1,nRows), nonMono, p99dt, ...
    'VariableNames', {'mouse','day','session','nRows','nNonFinite','fracNonFinite','nonMonotonic','p99_dt'});

writetable(R, fullfile(outDir, sprintf('QC_badTime_report_%s.csv', timeVar)));

dm = lower(string(debugMouse));
rD = lower(R.mouse) == dm;
if any(rD)
    writetable(R(rD,:), fullfile(outDir, sprintf('QC_badTime_DEBUG_%s_%s.csv', dm, timeVar)));
end
end

function v = safeP99dt(x)
x = x(isfinite(x));
if numel(x) < 3
    v = NaN;
    return;
end
d = diff(x);
d = d(isfinite(d));
if isempty(d)
    v = NaN;
else
    v = prctile(d,99);
end
end
