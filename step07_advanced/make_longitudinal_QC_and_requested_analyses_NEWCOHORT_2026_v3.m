function make_longitudinal_QC_and_requested_analyses_NEWCOHORT_2026_v3()
% make_longitudinal_QC_and_requested_analyses_NEWCOHORT_2026_v3
%
% Schema-STRICT (no guessing):
%   Reads ALL_mice_longitudinal.csv from the latest run_* folder
%   Requires columns (frame-level):
%     mouse_key, day_index, session_idx, Frame, Lick_TTL, Injector_TTL, Diameter_px
%   Optional (if present): RequirementLast, TrialRequirement, Trial, Session_Paradigm
%
% Uses explicit cohort + explicit active/passive pairing (NEWCOHORT).
% Produces:
%   - QC summaries (coverage, missingness, outlier flags)
%   - Lick metrics: lick count, lick freq, lick IEI/period, cumulative lick curves
%   - Period comparisons (with and without transition days)
%   - Pupil day trends normalized to Day3 baseline per mouse
%   - Event-locked pupil: lick-bout locked & reward-locked; rewarded vs unrewarded
%   - Tail immersion / TST / Hotplate / Straub if columns exist in a session-level table
%
% NOTES:
%   - Timebase from Frame/30 (consistent with your existing script).
%   - Lick bouts: consecutive licks with inter-lick interval <= 2.0 s
%   - Rewarded lick-bout: reward edge within [bout_start, bout_end + 1.0 s]
%
% Output:
%   Creates outDir under figs/lick_patterns_MASTER/plots_matlab_requested_v3/

%% -------------------- USER SETTINGS --------------------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Session time binning for cumulative curves
binMin = 1;           % minutes
maxPlotMin = 20;      % x-axis cap (still handles longer sessions, but plots up to this)

% Pupil event-lock windows
preS  = 2.0;
postS = 2.0;
boutGapS = 2.0;       % licks within 2s are same bout
rewardWinAfterBoutS = 1.0;  % reward within 1s after bout end => "rewarded bout"

% Transition days you called unreliable (you explicitly listed these)
transitionDays = [4 6 11 14 17];

% Passive focus comparisons you asked for
passiveFocus_days_during = [7 8 9];
passiveFocus_days_post   = [12 13];

% Assay day pairs (example you showed day5 vs day9)
assayDayA = 5;
assayDayB = 9;

rng(1);

%% -------------------- LOCATE LATEST RUN --------------------
if ~exist(rootTry,'dir')
    error('rootTry not found: %s', rootTry);
end

D = dir(fullfile(rootTry,'run_*'));
assert(~isempty(D), 'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);

dataDir = fullfile(runDir,'figs','lick_patterns_MASTER');
assert(exist(dataDir,'dir')==7, 'Expected folder not found: %s', dataDir);

ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(dataDir,'plots_matlab_requested_v3', ts);
if ~exist(outDir,'dir'); mkdir(outDir); end

%% -------------------- LOAD FRAME-LEVEL TABLE --------------------
longCSV = fullfile(dataDir,'ALL_mice_longitudinal.csv');
assert(exist(longCSV,'file')==2, 'Missing file: %s', longCSV);

T = readtable(longCSV);

% REQUIRED schema (no guessing)
REQ = {'mouse_key','day_index','session_idx','Frame','Lick_TTL','Injector_TTL','Diameter_px'};
requireVars(T, REQ, 'ALL_mice_longitudinal.csv');

% Timebase
T.tb_s = T.Frame ./ 30;  % consistent with your prior code

% Ensure sane ordering
T = sortrows(T, {'mouse_key','day_index','session_idx','tb_s'});

%% -------------------- EXPLICIT COHORT + PAIRING (NO GUESSING) --------------------
Coh = buildExplicitCohortTable_NEWCOHORT();
% Coh has: mouse_key, cage, color, sex, group, pair_id

% Join and validate
T = outerjoin(T, Coh, 'Keys','mouse_key', 'MergeKeys',true);
if any(ismissing(T.group))
    missingKeys = unique(T.mouse_key(ismissing(T.group)));
    error('These mouse_key were found in CSV but NOT in explicit cohort list:\n%s', strjoin(cellstr(missingKeys), ', '));
end

isActive  = strcmpi(T.group,'active');
isPassive = strcmpi(T.group,'passive');

%% -------------------- QC: COVERAGE + BASIC INTEGRITY --------------------
% Per (mouse,day,session): coverage + key signal counts
G = findgroups(T.mouse_key, T.day_index, T.session_idx);

S = table();
S.mouse_key   = splitapply(@unique, T.mouse_key, G);
S.day_index   = splitapply(@unique, T.day_index, G);
S.session_idx = splitapply(@unique, T.session_idx, G);

S.group = splitapply(@unique, T.group, G);
S.sex   = splitapply(@unique, T.sex,   G);
S.pair_id = splitapply(@unique, T.pair_id, G);

S.nFrames = splitapply(@numel, T.Frame, G);

% Licks as rising edges
S.nLicks = splitapply(@(x) countEdges(x), T.Lick_TTL, G);
S.nRewards = splitapply(@(x) countEdges(x), T.Injector_TTL, G);

S.duration_min = splitapply(@(x) (max(x)-min(x))/60, T.tb_s, G);

S.pupil_mean = splitapply(@(x) mean(x,'omitnan'), T.Diameter_px, G);
S.pupil_nanFrac = splitapply(@(x) mean(isnan(x)), T.Diameter_px, G);

% Lick IEI / lick frequency
S.lick_freq_per_min = S.nLicks ./ max(S.duration_min, 1e-6);

S.lick_IEI_mean_s = splitapply(@(ttl,tb) mean(diff(tb(edgeTimes(ttl,tb))),'omitnan'), T.Lick_TTL, T.tb_s, G);
S.lick_IEI_cv     = splitapply(@(ttl,tb) safeCV(diff(tb(edgeTimes(ttl,tb)))), T.Lick_TTL, T.tb_s, G);

% PR score if present (RequirementLast preferred)
S.PR = nan(height(S),1);
if ismember('RequirementLast', T.Properties.VariableNames)
    S.PR = splitapply(@(x) max(x,[],'omitnan'), T.RequirementLast, G);
elseif ismember('TrialRequirement', T.Properties.VariableNames)
    S.PR = splitapply(@(x) max(x,[],'omitnan'), T.TrialRequirement, G);
end

writetable(S, fullfile(outDir, sprintf('QC_session_summary_%s.csv', ts)));

% QC flags
QC = S(:, {'mouse_key','day_index','session_idx','group','nFrames','duration_min','nLicks','nRewards','pupil_mean','pupil_nanFrac','lick_freq_per_min','lick_IEI_mean_s','lick_IEI_cv','PR'});
QC.flag_short = QC.duration_min < 5;
QC.flag_noPupil = QC.pupil_nanFrac > 0.50;
QC.flag_noLicks = QC.nLicks == 0;

writetable(QC, fullfile(outDir, sprintf('QC_flags_%s.csv', ts)));

%% -------------------- PERIOD DEFINITIONS (EXPLICIT, AS YOU REQUESTED) --------------------
% Your timeline:
% day1-2 habituation (ignore for licking analyses)
% day3-5 water PR
% day6-10 passive training morphine (active has PR; passive forced replay)
% day11-13 morphine PR for all
% day14-16 withdrawal (water)
% day17-18 re-exposure morphine

Periods = buildPeriods(transitionDays);

% Two modes: include transition days vs exclude them
modes = {'includeTransitions','excludeTransitions'};
includeTransitionDaysVec = [true, false];

%% -------------------- PLOTS: SPAGHETTI (LICK COUNT PER DAY) --------------------
plotSpaghettiMetricPerDay(S, outDir, ts, 'nLicks', 'Lick count per day', true);

%% -------------------- PLOTS: CUMULATIVE LICKING CURVES (DAY-BY-DAY) --------------------
% Cumulative curves: for each session, bin licks in 1-min bins, then cumulative sum.
% Then average across mice within day/group.

plotCumulativeByDay(T, outDir, ts, maxPlotMin, binMin, 'active');
plotCumulativeByDay(T, outDir, ts, maxPlotMin, binMin, 'passive');

%% -------------------- PERIOD COMPARISONS (ACTIVE/PASSIVE DIFFERENT DAY SETS) --------------------
for mi = 1:numel(modes)
    includeTransitions = includeTransitionDaysVec(mi);
    modeName = modes{mi};

    % Active comparison (your request):
    % Active: Day5 vs Day6-10 vs Day11-13 vs Day14-16 vs Day17
    activeDays = selectDaysForGroup('active', Periods, includeTransitions);

    % Passive comparison (your request):
    % Passive: Day5 vs Day11-13 vs Day14-16 vs Day17-18
    passiveDays = selectDaysForGroup('passive', Periods, includeTransitions);

    plotPeriodCumulativeCurves(T, outDir, ts, maxPlotMin, binMin, activeDays,  sprintf('ACTIVE period cumulative (%s)', modeName),  sprintf('ACTIVE_period_cumulative_%s_%s.png', modeName, ts));
    plotPeriodCumulativeCurves(T, outDir, ts, maxPlotMin, binMin, passiveDays, sprintf('PASSIVE period cumulative (%s)', modeName), sprintf('PASSIVE_period_cumulative_%s_%s.png', modeName, ts));

    % Period-level box/scatter summaries for lick_freq, lick_IEI_mean, PR (if present)
    periodSummaryAndStats(S, outDir, ts, activeDays,  'active',  modeName);
    periodSummaryAndStats(S, outDir, ts, passiveDays, 'passive', modeName);
end

%% -------------------- PUPIL: DAY TREND + NORMALIZED TO DAY3 BASELINE --------------------
plotPupilDayTrendNormalized(S, outDir, ts);

%% -------------------- EVENT-LOCKED PUPIL: LICK-BOUT LOCKED + REWARD LOCKED --------------------
% We compute:
% - lick-bout locked pupil, separated by rewarded vs non-rewarded bouts
% - reward-locked pupil
% Then we do special passive focus comparison (days 7-8-9 vs 12-13) as requested.

eventLockedPupil_All(T, outDir, ts, preS, postS, boutGapS, rewardWinAfterBoutS);

eventLockedPupil_PassiveFocus(T, outDir, ts, preS, postS, boutGapS, rewardWinAfterBoutS, passiveFocus_days_during, passiveFocus_days_post);

%% -------------------- ASSAYS: TAIL IMMERSION / TST / HOTPLATE / STRAUB (IF PRESENT) --------------------
% These are typically session-level; we will look for a session-level file if you have it.
% If your assay variables are already inside S (because ALL_mice_longitudinal had them), this will still work.
%
% If the variables do not exist, we skip with a clear message.
assayCompareIfAvailable(S, outDir, ts, assayDayA, assayDayB);

disp('DONE. Outputs in:');
disp(outDir);

end

%% ========================= HELPERS =========================

function requireVars(T, vars, context)
missing = vars(~ismember(vars, T.Properties.VariableNames));
if ~isempty(missing)
    error('%s missing required columns:\n%s', context, strjoin(missing, ', '));
end
end

function n = countEdges(x)
% counts rising edges from 0->1 (robust for numeric TTL)
x = x(:);
x = double(x > 0);
n = sum(diff([0; x]) == 1);
end

function t = edgeTimes(ttl, tb)
% return event times (s) for rising edges
ttl = double(ttl(:) > 0);
ix = find(diff([0; ttl]) == 1);
t  = tb(ix);
end

function c = safeCV(x)
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    c = NaN;
else
    mu = mean(x);
    if mu == 0, c = NaN; else, c = std(x)/mu; end
end
end

function Coh = buildExplicitCohortTable_NEWCOHORT()
% EXPLICIT list from your message (NO GUESSING)
% mouse_key format is "cage_color" based on your pipeline.
%
% You said:
% 6100 orange, red (passive) vs 6100 black (active)
% 0911 red (active) vs 0911 orange (passive)
% 0911 black (passive) vs 0911 white (active)
% 0910 black (active) vs 0910 orange, red (passive)
% 6099 red (passive) vs 6099 orange (active)
% 6099 black (active) vs 6099 white (passive) BUT 6099 white died day1-13 (partial)
%
% Also: f/m sex.

rows = {
    % mouse_key       cage   color    sex    group     pair_id
    '6100_red',       '6100','red',   'f',  'passive', 'pair1';
    '6100_orange',    '6100','orange','f',  'passive', 'pair1';
    '6100_black',     '6100','black', 'f',  'active',  'pair1';

    '0911_red',       '0911','red',   'f',  'active',  'pair2';
    '0911_orange',    '0911','orange','f',  'passive', 'pair2';

    '0911_black',     '0911','black', 'f',  'passive', 'pair3';
    '0911_white',     '0911','white', 'f',  'active',  'pair3';

    '0910_red',       '0910','red',   'm',  'passive', 'pair4';
    '0910_orange',    '0910','orange','m',  'passive', 'pair4';
    '0910_black',     '0910','black', 'm',  'active',  'pair4';

    '6099_red',       '6099','red',   'm',  'passive', 'pair5';
    '6099_orange',    '6099','orange','m',  'active',  'pair5';

    '6099_black',     '6099','black', 'm',  'active',  'pair6';
    '6099_white',     '6099','white', 'm',  'passive', 'pair6'; % partial OK
};

Coh = cell2table(rows, 'VariableNames', {'mouse_key','cage','color','sex','group','pair_id'});
end

function Periods = buildPeriods(transitionDays)
Periods = struct();

Periods.habituation = 1:2;

Periods.pre      = 3:5;        % water PR
Periods.during   = 6:10;       % passive training morphine (active PR, passive forced)
Periods.post     = 11:13;      % morphine PR all
Periods.withdraw = 14:16;      % water withdrawal PR
Periods.reexp    = 17:18;      % morphine PR

Periods.transitionDays = transitionDays;

% For your “compare with/without unreliable days” mode:
% We will drop these from each period when excludeTransitions mode is used.
end

function daySets = selectDaysForGroup(groupName, Periods, includeTransitions)
% Returns a struct of named day sets for your requested period comparisons
% (Active and Passive are different).
%
% Active request:
%   Day5 vs Day6-10 vs Day11-13 vs Day14-16 vs Day17
% Passive request:
%   Day5 vs Day11-13 vs Day14-16 vs Day17-18

drop = Periods.transitionDays;
if includeTransitions
    drop = [];
end

if strcmpi(groupName,'active')
    daySets = struct();
    daySets.label = {'Day5','During(6-10)','Post(11-13)','Withdrawal(14-16)','Reexp(17)'};
    daySets.days  = {setdiff(5,drop), setdiff(6:10,drop), setdiff(11:13,drop), setdiff(14:16,drop), setdiff(17,drop)};
elseif strcmpi(groupName,'passive')
    daySets = struct();
    daySets.label = {'Day5','Post(11-13)','Withdrawal(14-16)','Reexp(17-18)'};
    daySets.days  = {setdiff(5,drop), setdiff(11:13,drop), setdiff(14:16,drop), setdiff(17:18,drop)};
else
    error('Unknown groupName: %s', groupName);
end

% Remove empty sets (if a day was dropped leaving empty)
keep = ~cellfun(@isempty, daySets.days);
daySets.label = daySets.label(keep);
daySets.days  = daySets.days(keep);
end

function plotSpaghettiMetricPerDay(S, outDir, ts, varName, titleStr, splitByGroup)
% Spaghetti plot of a session-level metric by day (averaged per mouse/day first)
if ~ismember(varName, S.Properties.VariableNames)
    warning('Skipping %s spaghetti: column not found.', varName);
    return;
end

% average across sessions within day for each mouse
G = findgroups(S.mouse_key, S.day_index);
M = table();
M.mouse_key = splitapply(@unique, S.mouse_key, G);
M.day_index = splitapply(@unique, S.day_index, G);
M.group     = splitapply(@unique, S.group, G);
M.val       = splitapply(@(x) mean(x,'omitnan'), S.(varName), G);

% drop habituation days 1-2 for licking-related plots (your request)
M = M(M.day_index >= 3, :);

if splitByGroup
    groups = {'active','passive'};
else
    groups = {'all'};
end

for gi = 1:numel(groups)
    if strcmp(groups{gi},'all')
        R = M;
        subTitle = 'All mice';
        fileTag  = 'ALL';
    else
        R = M(strcmpi(M.group,groups{gi}),:);
        subTitle = upper(groups{gi});
        fileTag  = upper(groups{gi});
    end

    if isempty(R), continue; end

    f = figure('Color','w','Position',[80 80 900 450]);
    hold on;

    mice = unique(R.mouse_key);
    for i = 1:numel(mice)
        rr = R(strcmp(R.mouse_key,mice{i}),:);
        [d,ix] = sort(rr.day_index);
        plot(d, rr.val(ix), '-o', 'LineWidth',1.2, 'MarkerSize',4);
    end

    % thick mean across mice per day
    Gd = findgroups(R.day_index);
    mu = splitapply(@(x) mean(x,'omitnan'), R.val, Gd);
    dd = splitapply(@unique, R.day_index, Gd);
    [dd,ix] = sort(dd);
    mu = mu(ix);
    plot(dd, mu, 'k-', 'LineWidth',3);

    xlabel('Day index');
    ylabel(varName);
    title(sprintf('%s (%s)', titleStr, subTitle));
    grid on;

    outP = fullfile(outDir, sprintf('Spaghetti_%s_%s_%s.png', varName, fileTag, ts));
    exportgraphics(f, outP, 'Resolution', 300);
    close(f);
end

end

function plotCumulativeByDay(T, outDir, ts, maxPlotMin, binMin, groupName)
% Day-by-day cumulative lick curves, like your example figure.
R = T(strcmpi(T.group,groupName), :);
if isempty(R)
    warning('No rows for group=%s', groupName);
    return;
end

% ignore day1-2
R = R(R.day_index >= 3, :);

% build per-session cumulative curve
G = findgroups(R.mouse_key, R.day_index, R.session_idx);
sessKey = splitapply(@(a,b,c) sprintf('%s_d%d_s%d', a{1}, b(1), c(1)), R.mouse_key, R.day_index, R.session_idx, G);
% precompute curves
edgesMin = 0:binMin:maxPlotMin;
nb = numel(edgesMin);

U = table();
U.sessKey = unique(sessKey);
U.mouse_key = splitapply(@unique, R.mouse_key, G);
U.day_index = splitapply(@unique, R.day_index, G);
U.group = splitapply(@unique, R.group, G);

U.cum = cell(height(U),1);

for i = 1:height(U)
    sk = U.sessKey{i};
    idx = strcmp(sessKey, sk);
    tb = R.tb_s(idx);
    ttl = R.Lick_TTL(idx);
    et = edgeTimes(ttl, tb);

    % bin into minutes
    etMin = et/60;
    h = histcounts(etMin, [edgesMin, edgesMin(end)+binMin]); %#ok<HISTC>
    c = cumsum(h(:));
    U.cum{i} = c(1:nb).';
end

% average by day across sessions
days = unique(U.day_index);
f = figure('Color','w','Position',[70 70 1000 550]);
hold on;

for d = reshape(sort(days),1,[])
    rr = U(U.day_index==d,:);
    if isempty(rr), continue; end
    C = cell2mat(rr.cum);
    mu = mean(C,1,'omitnan');
    plot(edgesMin, mu, 'LineWidth',2);
    text(edgesMin(end), mu(end), sprintf('Day %d', d), 'FontSize',9);
end

xlabel('Session time (min)');
ylabel('Mean cumulative licks');
title(sprintf('Cumulative licking across session (%s mice)', upper(groupName)));
grid on;

outP = fullfile(outDir, sprintf('Cumulative_byDay_%s_%s.png', upper(groupName), ts));
exportgraphics(f, outP, 'Resolution', 300);
close(f);

end

function plotPeriodCumulativeCurves(T, outDir, ts, maxPlotMin, binMin, daySets, figTitle, fileName)
% Plots mean cumulative curves for each period label
R = T(T.day_index >= 3, :);

edgesMin = 0:binMin:maxPlotMin;
nb = numel(edgesMin);

f = figure('Color','w','Position',[80 80 1000 550]);
hold on;

for k = 1:numel(daySets.label)
    days = daySets.days{k};
    rr = R(ismember(R.day_index, days) & strcmpi(R.group, inferGroupFromDaySets(daySets)), :);
    if isempty(rr), continue; end

    % per-session curves
    G = findgroups(rr.mouse_key, rr.day_index, rr.session_idx);
    cumCell = splitapply(@(ttl,tb) {cumCurve(ttl,tb,edgesMin,binMin,nb)}, rr.Lick_TTL, rr.tb_s, G);
    C = cell2mat(cumCell);
    mu = mean(C,1,'omitnan');
    se = std(C,0,1,'omitnan') ./ sqrt(size(C,1));

    plot(edgesMin, mu, 'LineWidth',2);
    % sparse errorbars
    ii = 1:max(1,round(numel(edgesMin)/10)):numel(edgesMin);
    errorbar(edgesMin(ii), mu(ii), se(ii), 'LineStyle','none');

end

xlabel('Session time (min)');
ylabel('Mean cumulative licks');
title(figTitle);
grid on;
legend(daySets.label,'Location','southeast');

outP = fullfile(outDir, fileName);
exportgraphics(f, outP, 'Resolution', 300);
close(f);

end

function grp = inferGroupFromDaySets(daySets)
% daySets created for either active or passive
% We infer by label count/pattern; safer: check first label list
L = daySets.label;
if any(contains(L,'During','IgnoreCase',true)) || any(contains(L,'Reexp(17)','IgnoreCase',true))
    grp = 'active';
else
    grp = 'passive';
end
end

function c = cumCurve(ttl, tb, edgesMin, binMin, nb)
et = edgeTimes(ttl, tb);
etMin = et/60;
h = histcounts(etMin, [edgesMin, edgesMin(end)+binMin]);
c = cumsum(h(:));
c = c(1:nb).';
end

function periodSummaryAndStats(S, outDir, ts, daySets, groupName, modeName)
% For each mouse, compute mean metric per period day set, then plot + write stats
metrics = {'nLicks','lick_freq_per_min','lick_IEI_mean_s','PR','pupil_mean'};
metrics = metrics(ismember(metrics, S.Properties.VariableNames));

R = S(strcmpi(S.group,groupName), :);
R = R(R.day_index >= 3, :);

if isempty(R) || isempty(metrics)
    return;
end

% average across sessions within (mouse,day)
G = findgroups(R.mouse_key, R.day_index);
D = table();
D.mouse_key = splitapply(@unique, R.mouse_key, G);
D.day_index = splitapply(@unique, R.day_index, G);
for m = 1:numel(metrics)
    D.(metrics{m}) = splitapply(@(x) mean(x,'omitnan'), R.(metrics{m}), G);
end

mice = unique(D.mouse_key);

for m = 1:numel(metrics)
    M = nan(numel(mice), numel(daySets.label));
    for i = 1:numel(mice)
        for k = 1:numel(daySets.label)
            days = daySets.days{k};
            rr = D(strcmp(D.mouse_key,mice{i}) & ismember(D.day_index,days), :);
            M(i,k) = mean(rr.(metrics{m}), 'omitnan');
        end
    end

    % Plot
    f = figure('Color','w','Position',[90 90 1100 420]);
    hold on;

    % scatter each mouse
    for i = 1:size(M,1)
        plot(1:size(M,2), M(i,:), '-o', 'LineWidth',1.2, 'MarkerSize',4);
    end

    mu = mean(M,1,'omitnan');
    se = std(M,0,1,'omitnan') ./ sqrt(sum(isfinite(M),1));
    errorbar(1:numel(mu), mu, se, 'k-', 'LineWidth',3);

    xticks(1:numel(daySets.label));
    xticklabels(daySets.label);
    xtickangle(20);
    ylabel(metrics{m});
    title(sprintf('%s %s (%s) — %s', upper(groupName), metrics{m}, modeName, 'per-mouse period means'));
    grid on;

    outP = fullfile(outDir, sprintf('PeriodSummary_%s_%s_%s_%s.png', upper(groupName), metrics{m}, modeName, ts));
    exportgraphics(f, outP, 'Resolution', 300);
    close(f);

    % Save table
    Tout = array2table(M, 'VariableNames', matlab.lang.makeValidName(daySets.label));
    Tout.mouse_key = mice;
    Tout = movevars(Tout,'mouse_key','Before',1);
    writetable(Tout, fullfile(outDir, sprintf('PeriodSummary_%s_%s_%s_%s.csv', upper(groupName), metrics{m}, modeName, ts)));

end
end

function plotPupilDayTrendNormalized(S, outDir, ts)
% Raw and Day3-normalized pupil diameter per mouse, split active/passive
if ~ismember('pupil_mean', S.Properties.VariableNames)
    warning('No pupil_mean in session summary; skipping pupil day plots.');
    return;
end

% average across sessions within (mouse,day)
G = findgroups(S.mouse_key, S.day_index);
D = table();
D.mouse_key = splitapply(@unique, S.mouse_key, G);
D.day_index = splitapply(@unique, S.day_index, G);
D.group     = splitapply(@unique, S.group, G);
D.pupil_mean = splitapply(@(x) mean(x,'omitnan'), S.pupil_mean, G);

% baseline day3 per mouse
mice = unique(D.mouse_key);
base = nan(numel(mice),1);
for i = 1:numel(mice)
    rr = D(strcmp(D.mouse_key,mice{i}) & D.day_index==3, :);
    base(i) = mean(rr.pupil_mean,'omitnan');
end
D.base_day3 = nan(height(D),1);
for i = 1:numel(mice)
    D.base_day3(strcmp(D.mouse_key,mice{i})) = base(i);
end
D.pupil_norm = (D.pupil_mean - D.base_day3) ./ D.base_day3;  % delta/baseline

% Plot active and passive
for grp = {'active','passive'}
    R = D(strcmpi(D.group, grp{1}), :);
    if isempty(R), continue; end

    f = figure('Color','w','Position',[80 80 950 500]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % raw
    nexttile; hold on;
    miceG = unique(R.mouse_key);
    for i = 1:numel(miceG)
        rr = R(strcmp(R.mouse_key,miceG{i}), :);
        [d,ix] = sort(rr.day_index);
        plot(d, rr.pupil_mean(ix), '-o', 'LineWidth',1.2);
    end
    Gd = findgroups(R.day_index);
    mu = splitapply(@(x) mean(x,'omitnan'), R.pupil_mean, Gd);
    dd = splitapply(@unique, R.day_index, Gd);
    [dd,ix] = sort(dd); mu = mu(ix);
    plot(dd, mu, 'k-', 'LineWidth',3);
    xlabel('Day'); ylabel('Pupil diameter (px)');
    title(sprintf('Pupil day trend (RAW) — %s', upper(grp{1})));
    grid on;

    % normalized
    nexttile; hold on;
    for i = 1:numel(miceG)
        rr = R(strcmp(R.mouse_key,miceG{i}), :);
        [d,ix] = sort(rr.day_index);
        plot(d, rr.pupil_norm(ix), '-o', 'LineWidth',1.2);
    end
    mu = splitapply(@(x) mean(x,'omitnan'), R.pupil_norm, Gd);
    [dd,ix] = sort(dd); mu = mu(ix);
    plot(dd, mu, 'k-', 'LineWidth',3);
    xlabel('Day'); ylabel('\Delta Pupil / Day3 baseline');
    title(sprintf('Pupil day trend (Day3-normalized) — %s', upper(grp{1})));
    grid on;

    outP = fullfile(outDir, sprintf('Pupil_DayTrend_%s_%s.png', upper(grp{1}), ts));
    exportgraphics(f, outP, 'Resolution', 300);
    close(f);
end

end

function eventLockedPupil_All(T, outDir, ts, preS, postS, boutGapS, rewardWinAfterBoutS)
% Event-locked pupil (lick-bout locked + reward locked) for active vs passive

% We will compute mean traces across all eligible sessions per group:
groups = {'active','passive'};

for gi = 1:numel(groups)
    R = T(strcmpi(T.group,groups{gi}), :);
    if isempty(R), continue; end

    % ignore day1-2
    R = R(R.day_index >= 3, :);

    [tvec, mu_bout_rew, se_bout_rew, mu_bout_norew, se_bout_norew, mu_rew, se_rew] = ...
        computeEventLockedTraces(R, preS, postS, boutGapS, rewardWinAfterBoutS);

    % Plot
    f = figure('Color','w','Position',[90 90 1100 420]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % lick-bout locked
    nexttile; hold on;
    plot(tvec, mu_bout_rew, 'LineWidth',2);
    plot(tvec, mu_bout_norew, 'LineWidth',2);
    xlabel('Time from bout start (s)');
    ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Lick-bout locked pupil — %s', upper(groups{gi})));
    grid on;
    legend({'Rewarded bouts','Non-reward bouts'}, 'Location','best');

    % reward-locked
    nexttile; hold on;
    plot(tvec, mu_rew, 'LineWidth',2);
    xlabel('Time from reward (s)');
    ylabel('\Delta pupil (baseline-subtracted)');
    title(sprintf('Reward-locked pupil — %s', upper(groups{gi})));
    grid on;

    outP = fullfile(outDir, sprintf('EventLockedPupil_ALL_%s_%s.png', upper(groups{gi}), ts));
    exportgraphics(f, outP, 'Resolution', 300);
    close(f);

    % Save traces
    Tout = table(tvec(:), mu_bout_rew(:), se_bout_rew(:), mu_bout_norew(:), se_bout_norew(:), mu_rew(:), se_rew(:), ...
        'VariableNames', {'t_s','mu_bout_rewarded','se_bout_rewarded','mu_bout_nonreward','se_bout_nonreward','mu_reward','se_reward'});
    writetable(Tout, fullfile(outDir, sprintf('EventLockedPupil_ALL_%s_%s.csv', upper(groups{gi}), ts)));

end

end

function eventLockedPupil_PassiveFocus(T, outDir, ts, preS, postS, boutGapS, rewardWinAfterBoutS, daysA, daysB)
% Passive focus: compare days 7-8-9 vs 12-13 (your request)
R = T(strcmpi(T.group,'passive') & T.day_index >= 3, :);
if isempty(R), return; end

RA = R(ismember(R.day_index, daysA), :);
RB = R(ismember(R.day_index, daysB), :);

if isempty(RA) || isempty(RB)
    warning('Passive focus event-locked pupil: missing data for requested days.');
    return;
end

[tvec, muA_bout_rew, ~, muA_bout_norew, ~, muA_rew, ~] = computeEventLockedTraces(RA, preS, postS, boutGapS, rewardWinAfterBoutS);
[~,   muB_bout_rew, ~, muB_bout_norew, ~, muB_rew, ~] = computeEventLockedTraces(RB, preS, postS, boutGapS, rewardWinAfterBoutS);

f = figure('Color','w','Position',[90 90 1100 420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(tvec, muA_bout_rew,  'LineWidth',2);
plot(tvec, muA_bout_norew,'LineWidth',2);
plot(tvec, muB_bout_rew,  '--', 'LineWidth',2);
plot(tvec, muB_bout_norew,'--', 'LineWidth',2);
xlabel('Time from bout start (s)');
ylabel('\Delta pupil');
title(sprintf('PASSIVE lick-bout locked: days [%s] vs [%s]', num2str(daysA), num2str(daysB)));
grid on;
legend({'A rewarded','A nonreward','B rewarded','B nonreward'}, 'Location','best');

nexttile; hold on;
plot(tvec, muA_rew, 'LineWidth',2);
plot(tvec, muB_rew, '--', 'LineWidth',2);
xlabel('Time from reward (s)');
ylabel('\Delta pupil');
title('PASSIVE reward-locked comparison');
grid on;
legend({'DaysA','DaysB'}, 'Location','best');

outP = fullfile(outDir, sprintf('EventLockedPupil_PASSIVE_FOCUS_%s.png', ts));
exportgraphics(f, outP, 'Resolution', 300);
close(f);

end

function [tvec, muBoutRew, seBoutRew, muBoutNo, seBoutNo, muRew, seRew] = ...
    computeEventLockedTraces(R, preS, postS, boutGapS, rewardWinAfterBoutS)
% Compute lick-bout locked traces (rewarded vs nonreward) and reward-locked traces.
% Aggregates across sessions.

% We resample pupil onto a fixed time grid for averaging
dt = median(diff(R.tb_s));
if ~isfinite(dt) || dt <= 0
    dt = 1/30;
end
tvec = (-preS:dt:postS)';

boutTraces_rew = [];
boutTraces_no  = [];
rewTraces      = [];

G = findgroups(R.mouse_key, R.day_index, R.session_idx);
uG = unique(G);

for gi = 1:numel(uG)
    idx = (G == uG(gi));
    tb   = R.tb_s(idx);
    pup  = R.Diameter_px(idx);
    lick = R.Lick_TTL(idx);
    rew  = R.Injector_TTL(idx);

    if numel(tb) < 10
        continue;
    end

    % Event times
    lickT = edgeTimes(lick, tb);
    rewT  = edgeTimes(rew,  tb);

    % Reward-locked traces (all rewards)
    for r = 1:numel(rewT)
        tr = extractDeltaTrace(tb, pup, rewT(r), tvec);
        if ~any(isnan(tr)), rewTraces(end+1,:) = tr; end %#ok<AGROW>
    end

    % Lick bouts
    if numel(lickT) < 2
        continue;
    end
    [boutStart, boutEnd] = lickBouts(lickT, boutGapS);

    for b = 1:numel(boutStart)
        % classify rewarded bout
        isRewBout = any(rewT >= boutStart(b) & rewT <= (boutEnd(b) + rewardWinAfterBoutS));

        tr = extractDeltaTrace(tb, pup, boutStart(b), tvec);
        if any(isnan(tr)), continue; end

        if isRewBout
            boutTraces_rew(end+1,:) = tr; %#ok<AGROW>
        else
            boutTraces_no(end+1,:)  = tr; %#ok<AGROW>
        end
    end

end

% Means/SEMs
[muBoutRew, seBoutRew] = meanSEM(boutTraces_rew);
[muBoutNo,  seBoutNo ] = meanSEM(boutTraces_no);
[muRew,     seRew    ] = meanSEM(rewTraces);

end

function tr = extractDeltaTrace(tb, pup, t0, tvec)
% Extract pupil around t0 and baseline-subtract using [-pre,0) window
tq = t0 + tvec;
tr = interp1(tb, pup, tq, 'linear', NaN);

% baseline = mean in [-pre, 0)
baseIdx = (tvec < 0);
base = mean(tr(baseIdx), 'omitnan');
tr = tr - base;
end

function [boutStart, boutEnd] = lickBouts(lickT, gapS)
% lick bout = licks separated by <= gapS
lickT = lickT(:);
d = diff(lickT);
cut = find(d > gapS);

if isempty(cut)
    boutStart = lickT(1);
    boutEnd   = lickT(end);
    return;
end

startIdx = [1; cut+1];
endIdx   = [cut; numel(lickT)];

boutStart = lickT(startIdx);
boutEnd   = lickT(endIdx);
end

function [mu, se] = meanSEM(X)
if isempty(X)
    mu = nan(1,0);
    se = nan(1,0);
    return;
end
mu = mean(X,1,'omitnan');
se = std(X,0,1,'omitnan') ./ sqrt(size(X,1));
end

function assayCompareIfAvailable(S, outDir, ts, dayA, dayB)
% This runs only if relevant assay columns exist in S.
% We check for common names seen in your pipeline.
cand = {'Immersion_Latency_s','HOT_Frames_Non_moving','TST_Frames_Non_moving','STRAUB_Score','STRAUB_Frames_Non_moving'};
present = cand(ismember(cand, S.Properties.VariableNames));

if isempty(present)
    disp('Assay columns not found in session summary; skipping Tail/TST/HOT/Straub.');
    return;
end

% average across sessions within (mouse,day)
G = findgroups(S.mouse_key, S.day_index);
D = table();
D.mouse_key = splitapply(@unique, S.mouse_key, G);
D.day_index = splitapply(@unique, S.day_index, G);
D.group     = splitapply(@unique, S.group, G);
for i = 1:numel(present)
    D.(present{i}) = splitapply(@(x) mean(x,'omitnan'), S.(present{i}), G);
end

for i = 1:numel(present)
    var = present{i};

    % Active vs Passive per day (ranksum + BH-FDR)
    days = unique(D.day_index);
    pvals = nan(numel(days),1);
    for di = 1:numel(days)
        d = days(di);
        a = D.(var)(D.day_index==d & strcmpi(D.group,'active'));
        p = D.(var)(D.day_index==d & strcmpi(D.group,'passive'));
        if numel(a)>=2 && numel(p)>=2
            pvals(di) = ranksum(a,p);
        end
    end
    q = bhFDR(pvals);

    Tout = table(days, pvals, q, 'VariableNames', {'day','p_ranksum','q_BH'});
    writetable(Tout, fullfile(outDir, sprintf('Assay_%s_ranksum_byDay_%s.csv', var, ts)));

    % 4-group ANOVA style view for dayA/dayB (Active/Passive x dayA/dayB)
    xa = D.(var)(D.day_index==dayA & strcmpi(D.group,'active'));
    xp = D.(var)(D.day_index==dayA & strcmpi(D.group,'passive'));
    ya = D.(var)(D.day_index==dayB & strcmpi(D.group,'active'));
    yp = D.(var)(D.day_index==dayB & strcmpi(D.group,'passive'));

    if numel(xa)>=2 && numel(xp)>=2 && numel(ya)>=2 && numel(yp)>=2
        X = [xa; xp; ya; yp];
        Gg = [repmat({'DayA_active'},numel(xa),1);
              repmat({'DayA_passive'},numel(xp),1);
              repmat({'DayB_active'},numel(ya),1);
              repmat({'DayB_passive'},numel(yp),1)];
        pAn = anova1(X, Gg, 'off');

        f = figure('Color','w','Position',[90 90 900 450]);
        boxplot(X, Gg);
        ylabel(var);
        title(sprintf('%s: 4-group ANOVA p=%.4g', var, pAn));
        grid on;

        outP = fullfile(outDir, sprintf('Assay_%s_box_ANOVA_day%d_day%d_%s.png', var, dayA, dayB, ts));
        exportgraphics(f, outP, 'Resolution', 300);
        close(f);
    end

end

end

function q = bhFDR(p)
% Benjamini-Hochberg FDR correction
p = p(:);
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv), return; end
[sp,ix] = sort(pv);
m = numel(sp);
qq = sp .* (m ./ (1:m)');
qq = min(cummin(flipud(qq)), 1);
qq = flipud(qq);
tmp = nan(size(pv));
tmp(ix) = qq;
q(valid) = tmp;
end
