function arranged_plots_all_ids_new()
% arranged_plots_all_ids_new
% Fully self-contained, robust version that:
% 1) Finds the latest run_* folder
% 2) Locates per_session_features.csv anywhere under that run
% 3) Locates a trial-level CSV (auto-detect from per_trial*.csv candidates)
% 4) Handles naming differences in lick/min and lick-rate columns
% 5) Adds your NEW cohort mapping (sex, active/passive, PairID)
% 6) Adds your NEW experimental periods:
%       Pre:        D3–5
%       During:     D6–10
%       Post:       D11–13
%       Withdrawal: D14–16
%       Reexposure: D17–18
%    And supports excluding “less reliable” transition days (default: [4 6 11 14 17])
% 7) Produces:
%    - Progress curves + ID grids
%    - Across-day curves + ID grids
%    - Active vs Passive period summaries (including vs excluding transition days)
%    - Pairwise (Active vs Passive) comparisons per period (including vs excluding)
%
% Output folder:
%   <where per_session_features.csv lives>\plots_matlab\

%% ----- USER OPTIONS -----
NBINS = 100; rng(1);

% Less reliable days you mentioned (transition days)
TRANSITION_DAYS = [4, 6, 11, 14, 17];

% Analysis day window (you said day1-2 habituation; main is 3-18)
MIN_DAY_FOR_ANALYSIS = 3;
MAX_DAY_FOR_ANALYSIS = 18;

% Root folder holding run_* folders
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

%% ----- LOCATE LATEST RUN -----
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up = 1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry = p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*'));
assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);

%% ----- FIND SESSION CSV (must exist) -----
sessPath = findOneFile(runDir, 'per_session_features.csv');
assert(~isempty(sessPath), 'Could not find per_session_features.csv under %s', runDir);

%% ----- FIND TRIAL CSV (auto-detect) -----
trialPath = findTrialCSV(runDir);
assert(~isempty(trialPath), 'Could not find a trial CSV under %s (looked for per_trial*.csv)', runDir);

dataDir = fileparts(sessPath);  % where session CSV lives
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(dataDir,'plots_matlab');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\nUsing runDir:\n  %s\n', runDir);
fprintf('Session CSV:\n  %s\n', sessPath);
fprintf('Trial CSV:\n  %s\n', trialPath);
fprintf('Output dir:\n  %s\n\n', outDir);

%% ----- READ TABLES -----
Tses   = readtable(sessPath);
Ttrial = readtable(trialPath);

keys = {'mouse_key','day_index','session_idx'};
assert(all(ismember(keys, Tses.Properties.VariableNames)), 'per_session_features missing keys');
assert(all(ismember(keys, Ttrial.Properties.VariableNames)), 'trial CSV missing keys: mouse_key/day_index/session_idx');

% Robust metric names
lickPM_col = pickVar(Tses, {'lick_per_min','licks_per_min','lick_per_minute','licks_per_minute'});
assert(~isempty(lickPM_col), ...
    'Cannot find lick_per_min-like column in per_session_features. Columns are:\n%s', ...
    strjoin(Tses.Properties.VariableNames, ', '));

% Trial lick rate column might differ; if missing, try to compute it
lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_rate_Hz','lick_hz','rate_hz'});
if isempty(lickRate_col)
    Ttrial = tryComputeLickRate(Ttrial);
    lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_rate_Hz','lick_hz','rate_hz'});
end
assert(~isempty(lickRate_col), ...
    'Cannot find or compute lick-rate column in trial CSV. Columns are:\n%s', ...
    strjoin(Ttrial.Properties.VariableNames, ', '));

% Ensure Trial column
assert(ismember('Trial', Ttrial.Properties.VariableNames), ...
    'Trial CSV is missing column "Trial". Columns are:\n%s', ...
    strjoin(Ttrial.Properties.VariableNames, ', '));

%% ----- ATTACH COHORT + PERIOD BINS -----
cohort = buildCohortTable();   % your new cohort mapping
Tses   = attachCohortAndPeriods(Tses, cohort, TRANSITION_DAYS);
Ttrial = attachCohortAndPeriods(Ttrial, cohort, TRANSITION_DAYS);

% Restrict to analysis day window
TsesA   = Tses(Tses.day_index>=MIN_DAY_FOR_ANALYSIS & Tses.day_index<=MAX_DAY_FOR_ANALYSIS, :);
TtrialA = Ttrial(Ttrial.day_index>=MIN_DAY_FOR_ANALYSIS & Ttrial.day_index<=MAX_DAY_FOR_ANALYSIS, :);

%% ----- CLUSTER LABELS IF MISSING -----
if ~ismember('Cluster', TsesA.Properties.VariableNames)
    feat = {'lick_per_min','licks_per_min','iei_cv','cv2_median','rhythm_index', ...
            'burst_fraction','bout_rate_per_min','iei_median','iei_mean'};
    use  = feat(ismember(feat, TsesA.Properties.VariableNames));
    assert(~isempty(use), 'No usable features found to compute clusters in per_session_features.');

    X = TsesA{:,use};
    X = zscoreNaN(X);
    X = fillmissing(X,'constant',0);

    TsesA.Cluster = kmeans(X,3,'Replicates',50,'MaxIter',1000,'Display','off');
end

% Propagate clusters back to the full Tses using keys
if ~ismember('Cluster', Tses.Properties.VariableNames)
    Tses.Cluster = nan(height(Tses),1);
end
TkeyA = TsesA(:,keys);   TkeyA.keyStr   = keyString(TkeyA);
TkeyAll = Tses(:,keys);  TkeyAll.keyStr = keyString(TkeyAll);
[tf,loc] = ismember(TkeyAll.keyStr, TkeyA.keyStr);
Tses.Cluster(tf) = TsesA.Cluster(loc(tf));

validClusters = unique(Tses.Cluster(isfinite(Tses.Cluster)));
validClusters = validClusters(:)';
assert(~isempty(validClusters), 'No valid Cluster labels found.');
K = max(validClusters);

C = lines(K);
mk = {'o','s','^','d','>','<'};

%% ==== Mouse-level cluster assignment (majority) ====
TsesCl = Tses(isfinite(Tses.Cluster) & Tses.day_index>=MIN_DAY_FOR_ANALYSIS & Tses.day_index<=MAX_DAY_FOR_ANALYSIS, :);

Gm = findgroups(TsesCl.mouse_key);
mouse_ids = splitapply(@unique, TsesCl.mouse_key, Gm);

counts = zeros(numel(mouse_ids), K);
for k = 1:K
    counts(:,k) = splitapply(@(x)sum(x==k), TsesCl.Cluster, Gm);
end
[bestCount, bestK] = max(counts, [], 2);
totalCount = sum(counts,2);
purity = bestCount ./ max(1,totalCount);

MouseAssign = table(mouse_ids, totalCount, bestK, purity, ...
    'VariableNames', {'mouse_key','n_sessions','majority_cluster','purity'});

writetable(MouseAssign, fullfile(outDir, ['mouse_majority_clusters_' ts '.csv']));

%% ==== Figure: Across-day using mouse-majority clusters (including vs excluding transition days) ====
Tmap = MouseAssign(:,{'mouse_key','majority_cluster'});
Tmap.Properties.VariableNames{2} = 'ClusterMajority';
TsesMaj = outerjoin(Tses, Tmap, 'Keys','mouse_key', 'MergeKeys',true);

Dsum_inc   = daySummaryByCluster(TsesMaj, lickPM_col, 'ClusterMajority', ...
    MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS, false, TRANSITION_DAYS);
Dsum_clean = daySummaryByCluster(TsesMaj, lickPM_col, 'ClusterMajority', ...
    MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS, true,  TRANSITION_DAYS);

plotAcrossDay(Dsum_inc,   K, C, mk, outDir, ['across_day_majority_includingTransitions_' ts '.png'], ...
    'Across day (mouse-level majority clusters) — including transition days');

plotAcrossDay(Dsum_clean, K, C, mk, outDir, ['across_day_majority_EXCLUDINGTransitions_' ts '.png'], ...
    'Across day (mouse-level majority clusters) — EXCLUDING transition days');

%% ----- Build normalized progress summary (equal weight per session) -----
Ktrl = join(TtrialA, TsesCl(:,[keys {'Cluster'}]), 'Keys', keys);

Gsess = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx);
nT    = splitapply(@max, Ktrl.Trial, Gsess);
prog  = (Ktrl.Trial - 0.5) ./ nT(Gsess);

binEdges = linspace(0,1,NBINS+1);
[~,bin] = histc(prog, binEdges); %#ok<HISTC>
bin(bin<1)=1; bin(bin>NBINS)=NBINS;

% per-session per-bin mean lick rate
Gsb = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx, bin);
Ssb = table;
Ssb.mouse_key   = splitapply(@unique, Ktrl.mouse_key,   Gsb);
Ssb.day_index   = splitapply(@unique, Ktrl.day_index,   Gsb);
Ssb.session_idx = splitapply(@unique, Ktrl.session_idx, Gsb);
Ssb.Cluster     = splitapply(@unique, Ktrl.Cluster,     Gsb);
Ssb.bin         = splitapply(@unique, bin,              Gsb);
Ssb.meanLR      = splitapply(@mean,  Ktrl.(lickRate_col),    Gsb);

% cluster mean±SEM across sessions
Gcb = findgroups(Ssb.Cluster, Ssb.bin);
Csum = table;
Csum.Cluster = splitapply(@unique, Ssb.Cluster, Gcb);
Csum.bin     = splitapply(@unique, Ssb.bin,     Gcb);
Csum.nSess   = splitapply(@numel,  Ssb.meanLR,  Gcb);
Csum.mu      = splitapply(@mean,   Ssb.meanLR,  Gcb);
Csum.sem     = splitapply(@(x)nanstd(x,0)/sqrt(sum(isfinite(x))), Ssb.meanLR, Gcb);
xPct = (Csum.bin - 0.5) / NBINS * 100;

%% ===== Figure A: Normalized progress (left) + ID grids (right) =====
plotProgressWithIDs(Csum, xPct, Ssb, K, C, mk, outDir, ['progress_meanSEM_IDGRID_' ts '.png'], NBINS);

%% ===== Figure B: Across-day mean±SEM (session-level cluster) + ID grids =====
plotAcrossDayWithIDs(TsesCl, lickPM_col, K, C, mk, outDir, ['across_day_meanSEM_IDGRID_' ts '.png']);

%% ===== Period summaries: Active vs Passive (including vs excluding transition days) =====
plotPeriodActivePassive(Tses, lickPM_col, 'period', outDir, ...
    ['period_active_passive_includingTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

plotPeriodActivePassive(Tses, lickPM_col, 'period_clean', outDir, ...
    ['period_active_passive_EXCLUDINGTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

%% ===== Pairwise by PairID per period (including vs excluding transition days) =====
plotPairwiseByPeriod(Tses, lickPM_col, 'period', outDir, ...
    ['pairwise_active_passive_byPeriod_includingTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

plotPairwiseByPeriod(Tses, lickPM_col, 'period_clean', outDir, ...
    ['pairwise_active_passive_byPeriod_EXCLUDINGTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

%% ----- Save a simple cluster ID list report (session-level clusters) -----
fid = fopen(fullfile(outDir, ['mouse_ids_by_cluster_' ts '.txt']),'w');
for k=1:K
    ids = sort(unique(TsesCl.mouse_key(TsesCl.Cluster==k)));
    fprintf(fid,'Cluster %d (%d mice):\n', k, numel(ids));
    fprintf(fid,'  %s\n\n', strjoin(ids, ', '));
end
fclose(fid);

writetable(groupsummary(TsesCl(:,{'mouse_key','Cluster'}),'Cluster'), ...
    fullfile(outDir, ['mouse_counts_per_cluster_' ts '.csv']));

disp('All done.');
end

%% ====================== FILE FINDERS ======================

function p = findOneFile(rootDir, fileName)
% Recursively find fileName under rootDir; if multiple, pick most recently modified.
dd = dir(fullfile(rootDir, '**', fileName));
if isempty(dd)
    p = '';
else
    [~,ix] = max([dd.datenum]);
    p = fullfile(dd(ix).folder, dd(ix).name);
end
end

function trialPath = findTrialCSV(runDir)
% Find a reasonable trial-level CSV automatically.
% Strategy:
% 1) Prefer exact known names if present
% 2) Else pick among per_trial*.csv, preferring ones containing "rate" then "feature"
preferred = { ...
    'per_trial_rates.csv', ...
    'per_trial_features.csv', ...
    'per_trial_features_rates.csv', ...
    'per_trial_lick_rates.csv' ...
};

trialPath = '';
for i=1:numel(preferred)
    p = findOneFile(runDir, preferred{i});
    if ~isempty(p), trialPath = p; return; end
end

% Otherwise, scan per_trial*.csv
dd = dir(fullfile(runDir, '**', 'per_trial*.csv'));
if isempty(dd), return; end

names = lower(string({dd.name}));
score = zeros(numel(dd),1);

% scoring: prefer "rate", then "lick", then "feature"
score(contains(names,'rate'))    = score(contains(names,'rate'))    + 3;
score(contains(names,'lick'))    = score(contains(names,'lick'))    + 2;
score(contains(names,'feature')) = score(contains(names,'feature')) + 1;

% tie-breaker: more recent datenum
[~,ix] = max(score + 1e-6 * [dd.datenum]');
trialPath = fullfile(dd(ix).folder, dd(ix).name);
end

%% ====================== COLUMN HELPERS ======================

function col = pickVar(T, candidates)
col = '';
for i=1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        col = candidates{i};
        return
    end
end
end

function T = tryComputeLickRate(T)
% Attempt to compute a lick-rate column if absent.
% Looks for combinations like:
%   licks / duration_s
%   n_licks / trial_dur_s
% If found, adds T.lick_rate (Hz)

countCol = pickVar(T, {'n_licks','licks','lick_count','nlicks','lick_total'});
durCol   = pickVar(T, {'duration_s','dur_s','trial_duration_s','trial_dur_s','dt_s','trial_len_s'});
if isempty(countCol) || isempty(durCol)
    return
end

c = double(T.(countCol));
d = double(T.(durCol));
d(d<=0 | ~isfinite(d)) = NaN;

T.lick_rate = c ./ d;
end

function Xz = zscoreNaN(X)
mu = nanmean(X,1);
sd = nanstd(X,0,1);
sd(sd==0 | ~isfinite(sd)) = 1;
Xz = (X - mu) ./ sd;
end

function ks = keyString(Tk)
mk = string(Tk.mouse_key);
d  = string(Tk.day_index);
s  = string(Tk.session_idx);
ks = normalizeMouseKey(mk) + "__D" + d + "__S" + s;
end

function s = normalizeMouseKey(sIn)
% Makes mouse_key robust across formats:
% - lower
% - spaces/hyphens -> underscore
% - handles "0911red" -> "0911_red"
s = string(sIn);
s = strtrim(lower(s));
s = regexprep(s,'[\-]+','_');
s = regexprep(s,'\s+','_');

s2 = s;
for i=1:numel(s)
    tok = regexp(s(i), '^(\d{4})([a-z]+)$', 'tokens','once');
    if ~isempty(tok)
        s2(i) = tok{1} + "_" + tok{2};
    end
end
s = s2;
end

%% ====================== COHORT + PERIODS ======================

function cohort = buildCohortTable()
% User-provided cohort mapping
% Format: mouse_key, Sex, ActPass, PairID

rows = {
'6100_red',      'F', 'Passive', 1
'6100_orange',   'F', 'Passive', 1
'6100_black',    'F', 'Active',  1

'0911_red',      'F', 'Active',  2
'0911_orange',   'F', 'Passive', 2

'0911_white',    'F', 'Active',  3
'0911_black',    'F', 'Passive', 3

'0910_black',    'M', 'Active',  4
'0910_orange',   'M', 'Passive', 4
'0910_red',      'M', 'Passive', 4

'6099_orange',   'M', 'Active',  5
'6099_red',      'M', 'Passive', 5

'6099_black',    'M', 'Active',  6
'6099_white',    'M', 'Passive', 6   % died day1-13: handled naturally by missing days
};

mouse_key = string(rows(:,1));
Sex      = categorical(string(rows(:,2)), {'F','M'});
ActPass  = categorical(string(rows(:,3)), {'Active','Passive'});
PairID   = cell2mat(rows(:,4));

Cage  = strings(size(mouse_key));
Color = strings(size(mouse_key));
for i=1:numel(mouse_key)
    parts = split(mouse_key(i), "_");
    Cage(i)  = parts(1);
    Color(i) = parts(2);
end

cohort = table(mouse_key, Cage, Color, Sex, ActPass, PairID);
end

function T = attachCohortAndPeriods(T, cohort, transitionDays)
if ismember('mouse_key', T.Properties.VariableNames) && ~isstring(T.mouse_key)
    T.mouse_key = string(T.mouse_key);
end

mk = normalizeMouseKey(T.mouse_key);
ck = normalizeMouseKey(cohort.mouse_key);
[tf,loc] = ismember(mk, ck);

% default placeholders
T.Cage    = repmat("", height(T),1);
T.Color   = repmat("", height(T),1);
T.Sex     = categorical(repmat("U",height(T),1), {'F','M','U'});
T.ActPass = categorical(repmat("U",height(T),1), {'Active','Passive','U'});
T.PairID  = nan(height(T),1);

% fill from cohort where matched
T.Cage(tf)    = cohort.Cage(loc(tf));
T.Color(tf)   = cohort.Color(loc(tf));
T.Sex(tf)     = cohort.Sex(loc(tf));
T.ActPass(tf) = cohort.ActPass(loc(tf));
T.PairID(tf)  = cohort.PairID(loc(tf));

% periods
d = double(T.day_index);
T.period       = dayToPeriod(d, false, transitionDays);
T.period_clean = dayToPeriod(d, true,  transitionDays);
end

function cats = periodCats(withExclude)
base = ["Pre_D3-5","During_D6-10","Post_D11-13","Withdrawal_D14-16","Reexposure_D17-18","<undef>"];
if withExclude
    cats = [base(1:5), "<exclude>", base(6)];
else
    cats = base;
end
cats = cellstr(cats);
end

function P = dayToPeriod(day_index, excludeTransitions, transitionDays)
lab = strings(size(day_index)); lab(:) = "<undef>";
d = double(day_index);

lab(d>=3  & d<=5 )  = "Pre_D3-5";
lab(d>=6  & d<=10 ) = "During_D6-10";
lab(d>=11 & d<=13 ) = "Post_D11-13";
lab(d>=14 & d<=16 ) = "Withdrawal_D14-16";
lab(d>=17 & d<=18 ) = "Reexposure_D17-18";

if excludeTransitions
    lab(ismember(d, transitionDays)) = "<exclude>";
    P = categorical(lab, periodCats(true));
else
    P = categorical(lab, periodCats(false));
end
end

%% ====================== SUMMARIES ======================

function Dsum = daySummaryByCluster(T, lickCol, clusterCol, dayMin, dayMax, excludeTransitions, transitionDays)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if excludeTransitions
    T2 = T2(~ismember(double(T2.day_index), transitionDays), :);
end

cl = T2.(clusterCol);
if iscategorical(cl)
    clNum = double(cl);
elseif iscell(cl) || isstring(cl)
    clNum = double(categorical(string(cl)));
else
    clNum = double(cl);
end

keep = isfinite(clNum) & isfinite(double(T2.(lickCol))) & isfinite(double(T2.day_index));
T2 = T2(keep,:);
clNum = clNum(keep);

Gd = findgroups(T2.day_index, clNum);

Dsum = table;
Dsum.day_index = splitapply(@unique, T2.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, clNum,        Gd);
Dsum.mu        = splitapply(@mean,  T2.(lickCol),  Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), T2.(lickCol), Gd);
end

%% ====================== PLOTS ======================

function plotAcrossDay(Dsum, K, C, mk, outDir, fname, ttl)
f = figure('Color','w','Position',[60 60 900 520]); hold on
for k=1:K
    r = (Dsum.Cluster==k);
    if ~any(r), continue; end
    [d,ord] = sort(Dsum.day_index(r));
    mu = Dsum.mu(r);  mu = mu(ord);
    se = Dsum.sem(r); se = se(ord);

    plot(d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel('Day'); ylabel('Licks / min'); title(ttl);
grid on; legend('off');
exportgraphics(f, fullfile(outDir, fname), 'Resolution', 300);
close(f);
end

function plotProgressWithIDs(Csum, xPct, Ssb, K, C, mk, outDir, fname, NBINS)
fA = figure('Color','w','Position',[60 60 1400 600]);
tlo = tiledlayout(fA,1,2,'TileSpacing','compact','Padding','compact');

% Left: curves
ax1 = nexttile(tlo,1); hold(ax1,'on');
for k=1:K
    r = (Csum.Cluster==k);
    if ~any(r), continue; end
    [xx,ord] = sort(xPct(r));
    mu = Csum.mu(r);  mu = mu(ord);
    se = Csum.sem(r); se = se(ord);

    plot(ax1, xx, mu, '-', 'Color', C(k,:), 'LineWidth', 1.6, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerSize',3, 'MarkerFaceColor',C(k,:));

    ii = 1:max(1,round(numel(xx)/50)):numel(xx);
    eb = errorbar(ax1, xx(ii), mu(ii), se(ii), 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel(ax1,'Session progress (%)'); ylabel(ax1,'Lick rate (Hz)');
title(ax1, sprintf('Clusters: mean \\pm SEM across normalized session progress (%d bins)', NBINS));
grid(ax1,'on'); legend(ax1,'off');

% Right: IDs per cluster
axIDs = nexttile(tlo,2); axis(axIDs,'off'); title(axIDs,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;
for k=1:K
    ids = sort(unique(Ssb.mouse_key(Ssb.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);
    ncol = numel(txtBlocks);
    x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end

exportgraphics(fA, fullfile(outDir,fname), 'Resolution', 300);
close(fA);
end

function plotAcrossDayWithIDs(TsesCl, lickPM_col, K, C, mk, outDir, fname)
% aggregate by day x cluster
Gd = findgroups(TsesCl.day_index, TsesCl.Cluster);
Dsum = table;
Dsum.day_index = splitapply(@unique, TsesCl.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, TsesCl.Cluster,   Gd);
Dsum.mu        = splitapply(@mean,  TsesCl.(lickPM_col), Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), TsesCl.(lickPM_col), Gd);

fB = figure('Color','w','Position',[60 60 1400 600]);
t2 = tiledlayout(fB,1,2,'TileSpacing','compact','Padding','compact');

ax2 = nexttile(t2,1); hold(ax2,'on');
for k=1:K
    r = Dsum.Cluster==k;
    if ~any(r), continue; end
    [d,ord] = sort(Dsum.day_index(r));
    mu = Dsum.mu(r);  mu = mu(ord);
    se = Dsum.sem(r); se = se(ord);

    plot(ax2, d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(ax2, d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel(ax2,'Day'); ylabel(ax2,'Licks / min');
title(ax2,'Clusters: mean \pm SEM across day (session-level clusters)');
grid(ax2,'on'); legend(ax2,'off');

axIDs2 = nexttile(t2,2); axis(axIDs2,'off'); title(axIDs2,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;
for k=1:K
    ids = sort(unique(TsesCl.mouse_key(TsesCl.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);
    ncol = numel(txtBlocks); x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs2, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end

exportgraphics(fB, fullfile(outDir,fname), 'Resolution', 300);
close(fB);
end

function plotPeriodActivePassive(T, lickCol, periodVar, outDir, fname, dayMin, dayMax)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if ~ismember(periodVar, T2.Properties.VariableNames), return; end

P = T2.(periodVar);
keep = (P ~= "<undef>");
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P ~= "<exclude>");
end
T2 = T2(keep,:); P = T2.(periodVar);

% Mouse-level period average (each mouse contributes equally)
G = findgroups(T2.mouse_key, P, T2.ActPass);
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.period    = splitapply(@unique, P, G);
S.ActPass   = splitapply(@unique, T2.ActPass, G);
S.mu_mouse  = splitapply(@mean,  T2.(lickCol), G);

% Aggregate across mice
G2 = findgroups(S.period, S.ActPass);
Agg = table;
Agg.period  = splitapply(@unique, S.period, G2);
Agg.ActPass = splitapply(@unique, S.ActPass, G2);
Agg.mu      = splitapply(@mean,  S.mu_mouse, G2);
Agg.sem     = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), S.mu_mouse, G2);

periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(Agg.period)));
if isempty(periods), return; end

f = figure('Color','w','Position',[80 80 980 420]);
ax = axes(f); hold(ax,'on');

x = 1:numel(periods);
xA = x - 0.12; xP = x + 0.12;

muA = nan(1,numel(periods)); seA = muA;
muP = nan(1,numel(periods)); seP = muP;

for i=1:numel(periods)
    rrA = (Agg.ActPass=="Active")  & (Agg.period==periods{i});
    rrP = (Agg.ActPass=="Passive") & (Agg.period==periods{i});
    if any(rrA), muA(i)=Agg.mu(rrA); seA(i)=Agg.sem(rrA); end
    if any(rrP), muP(i)=Agg.mu(rrP); seP(i)=Agg.sem(rrP); end
end

plot(ax, xA, muA, '-o','LineWidth',1.8);
errorbar(ax, xA, muA, seA, 'LineStyle','none');

plot(ax, xP, muP, '-o','LineWidth',1.8);
errorbar(ax, xP, muP, seP, 'LineStyle','none');

set(ax,'XTick',x,'XTickLabel',periods); xtickangle(ax,25);
ylabel(ax,'Licks / min');
title(ax, sprintf('Period summary (Active vs Passive) — %s', periodVar), 'Interpreter','none');
grid(ax,'on');
legend(ax, {'Active','Active SEM','Passive','Passive SEM'}, 'Location','best');

exportgraphics(f, fullfile(outDir,fname), 'Resolution', 300);
close(f);
end

function plotPairwiseByPeriod(T, lickCol, periodVar, outDir, fname, dayMin, dayMax)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if ~ismember(periodVar, T2.Properties.VariableNames), return; end

P = T2.(periodVar);
keep = (P ~= "<undef>");
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P ~= "<exclude>");
end
T2 = T2(keep,:);

% Mouse-level period mean (so each mouse contributes one value per period)
G = findgroups(T2.mouse_key, T2.PairID, T2.ActPass, T2.(periodVar));
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.PairID    = splitapply(@unique, T2.PairID,    G);
S.ActPass   = splitapply(@unique, T2.ActPass,   G);
S.period    = splitapply(@unique, T2.(periodVar), G);
S.mu_mouse  = splitapply(@mean,  T2.(lickCol),  G);

periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(S.period)));
if isempty(periods), return; end

nP = numel(periods);
f = figure('Color','w','Position',[80 80 320*nP 380]);
tlo = tiledlayout(f,1,nP,'TileSpacing','compact','Padding','compact');

for ip = 1:nP
    ax = nexttile(tlo, ip); hold(ax,'on');
    pr = periods{ip};
    Sp = S(S.period==pr & isfinite(S.mu_mouse) & isfinite(S.PairID), :);

    pairIDs = unique(Sp.PairID);

    % line per passive mouse to its active partner (supports 1 active vs 2 passive)
    for i=1:numel(pairIDs)
        pid = pairIDs(i);
        A = Sp(Sp.PairID==pid & Sp.ActPass=="Active", :);
        Pm = Sp(Sp.PairID==pid & Sp.ActPass=="Passive", :);
        if isempty(A) || isempty(Pm), continue; end
        aVal = A.mu_mouse(1);
        for j=1:height(Pm)
            pVal = Pm.mu_mouse(j);
            plot(ax, [1 2], [aVal pVal], '-', 'LineWidth',1.2);
        end
    end

    aPts = Sp.mu_mouse(Sp.ActPass=="Active");
    pPts = Sp.mu_mouse(Sp.ActPass=="Passive");
    scatter(ax, 1 + 0.06*(rand(numel(aPts),1)-0.5), aPts, 28, 'filled', 'MarkerFaceAlpha',0.6);
    scatter(ax, 2 + 0.06*(rand(numel(pPts),1)-0.5), pPts, 28, 'filled', 'MarkerFaceAlpha',0.6);

    xlim(ax,[0.6 2.4]);
    set(ax,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
    ylabel(ax,'Licks / min');
    title(ax, pr, 'Interpreter','none');
    grid(ax,'on'); box(ax,'on');
end

title(tlo, sprintf('Pairwise Active vs Passive (by PairID) — %s', periodVar), 'Interpreter','none');
exportgraphics(f, fullfile(outDir,fname), 'Resolution', 300);
close(f);
end

%% ====================== ID GRID HELPER ======================

function blocks = idColumns(ids, maxRows)
% Split a long list of IDs into columns of ~maxRows lines.
ids = string(ids);
n = numel(ids);
ncol = max(1, ceil(n/maxRows));
nrow = ceil(n/ncol);

blocks = cell(1,ncol);
for c = 1:ncol
    i1 = (c-1)*nrow + 1;
    i2 = min(c*nrow, n);
    blocks{c} = strjoin(ids(i1:i2), newline);
end
end
function arranged_plots_all_ids_new()
% arranged_plots_all_ids_new
% Fully self-contained, robust version that:
% 1) Finds the latest run_* folder
% 2) Locates per_session_features.csv anywhere under that run
% 3) Locates a trial-level CSV (auto-detect from per_trial*.csv candidates)
% 4) Handles naming differences in lick/min and lick-rate columns
% 5) Adds your NEW cohort mapping (sex, active/passive, PairID)
% 6) Adds your NEW experimental periods:
%       Pre:        D3–5
%       During:     D6–10
%       Post:       D11–13
%       Withdrawal: D14–16
%       Reexposure: D17–18
%    And supports excluding “less reliable” transition days (default: [4 6 11 14 17])
% 7) Produces:
%    - Progress curves + ID grids
%    - Across-day curves + ID grids
%    - Active vs Passive period summaries (including vs excluding transition days)
%    - Pairwise (Active vs Passive) comparisons per period (including vs excluding)
%
% Output folder:
%   <where per_session_features.csv lives>\plots_matlab\

%% ----- USER OPTIONS -----
NBINS = 100; rng(1);

% Less reliable days you mentioned (transition days)
TRANSITION_DAYS = [4, 6, 11, 14, 17];

% Analysis day window (you said day1-2 habituation; main is 3-18)
MIN_DAY_FOR_ANALYSIS = 3;
MAX_DAY_FOR_ANALYSIS = 18;

% Root folder holding run_* folders
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

%% ----- LOCATE LATEST RUN -----
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up = 1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry = p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*'));
assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);

%% ----- FIND SESSION CSV (must exist) -----
sessPath = findOneFile(runDir, 'per_session_features.csv');
assert(~isempty(sessPath), 'Could not find per_session_features.csv under %s', runDir);

%% ----- FIND TRIAL CSV (auto-detect) -----
trialPath = findTrialCSV(runDir);
assert(~isempty(trialPath), 'Could not find a trial CSV under %s (looked for per_trial*.csv)', runDir);

dataDir = fileparts(sessPath);  % where session CSV lives
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(dataDir,'plots_matlab');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\nUsing runDir:\n  %s\n', runDir);
fprintf('Session CSV:\n  %s\n', sessPath);
fprintf('Trial CSV:\n  %s\n', trialPath);
fprintf('Output dir:\n  %s\n\n', outDir);

%% ----- READ TABLES -----
Tses   = readtable(sessPath);
Ttrial = readtable(trialPath);

keys = {'mouse_key','day_index','session_idx'};
assert(all(ismember(keys, Tses.Properties.VariableNames)), 'per_session_features missing keys');
assert(all(ismember(keys, Ttrial.Properties.VariableNames)), 'trial CSV missing keys: mouse_key/day_index/session_idx');

% Robust metric names
lickPM_col = pickVar(Tses, {'lick_per_min','licks_per_min','lick_per_minute','licks_per_minute'});
assert(~isempty(lickPM_col), ...
    'Cannot find lick_per_min-like column in per_session_features. Columns are:\n%s', ...
    strjoin(Tses.Properties.VariableNames, ', '));

% Trial lick rate column might differ; if missing, try to compute it
lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_rate_Hz','lick_hz','rate_hz'});
if isempty(lickRate_col)
    Ttrial = tryComputeLickRate(Ttrial);
    lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_rate_Hz','lick_hz','rate_hz'});
end
assert(~isempty(lickRate_col), ...
    'Cannot find or compute lick-rate column in trial CSV. Columns are:\n%s', ...
    strjoin(Ttrial.Properties.VariableNames, ', '));

% Ensure Trial column
assert(ismember('Trial', Ttrial.Properties.VariableNames), ...
    'Trial CSV is missing column "Trial". Columns are:\n%s', ...
    strjoin(Ttrial.Properties.VariableNames, ', '));

%% ----- ATTACH COHORT + PERIOD BINS -----
cohort = buildCohortTable();   % your new cohort mapping
Tses   = attachCohortAndPeriods(Tses, cohort, TRANSITION_DAYS);
Ttrial = attachCohortAndPeriods(Ttrial, cohort, TRANSITION_DAYS);

% Restrict to analysis day window
TsesA   = Tses(Tses.day_index>=MIN_DAY_FOR_ANALYSIS & Tses.day_index<=MAX_DAY_FOR_ANALYSIS, :);
TtrialA = Ttrial(Ttrial.day_index>=MIN_DAY_FOR_ANALYSIS & Ttrial.day_index<=MAX_DAY_FOR_ANALYSIS, :);

%% ----- CLUSTER LABELS IF MISSING -----
if ~ismember('Cluster', TsesA.Properties.VariableNames)
    feat = {'lick_per_min','licks_per_min','iei_cv','cv2_median','rhythm_index', ...
            'burst_fraction','bout_rate_per_min','iei_median','iei_mean'};
    use  = feat(ismember(feat, TsesA.Properties.VariableNames));
    assert(~isempty(use), 'No usable features found to compute clusters in per_session_features.');

    X = TsesA{:,use};
    X = zscoreNaN(X);
    X = fillmissing(X,'constant',0);

    TsesA.Cluster = kmeans(X,3,'Replicates',50,'MaxIter',1000,'Display','off');
end

% Propagate clusters back to the full Tses using keys
if ~ismember('Cluster', Tses.Properties.VariableNames)
    Tses.Cluster = nan(height(Tses),1);
end
TkeyA = TsesA(:,keys);   TkeyA.keyStr   = keyString(TkeyA);
TkeyAll = Tses(:,keys);  TkeyAll.keyStr = keyString(TkeyAll);
[tf,loc] = ismember(TkeyAll.keyStr, TkeyA.keyStr);
Tses.Cluster(tf) = TsesA.Cluster(loc(tf));

validClusters = unique(Tses.Cluster(isfinite(Tses.Cluster)));
validClusters = validClusters(:)';
assert(~isempty(validClusters), 'No valid Cluster labels found.');
K = max(validClusters);

C = lines(K);
mk = {'o','s','^','d','>','<'};

%% ==== Mouse-level cluster assignment (majority) ====
TsesCl = Tses(isfinite(Tses.Cluster) & Tses.day_index>=MIN_DAY_FOR_ANALYSIS & Tses.day_index<=MAX_DAY_FOR_ANALYSIS, :);

Gm = findgroups(TsesCl.mouse_key);
mouse_ids = splitapply(@unique, TsesCl.mouse_key, Gm);

counts = zeros(numel(mouse_ids), K);
for k = 1:K
    counts(:,k) = splitapply(@(x)sum(x==k), TsesCl.Cluster, Gm);
end
[bestCount, bestK] = max(counts, [], 2);
totalCount = sum(counts,2);
purity = bestCount ./ max(1,totalCount);

MouseAssign = table(mouse_ids, totalCount, bestK, purity, ...
    'VariableNames', {'mouse_key','n_sessions','majority_cluster','purity'});

writetable(MouseAssign, fullfile(outDir, ['mouse_majority_clusters_' ts '.csv']));

%% ==== Figure: Across-day using mouse-majority clusters (including vs excluding transition days) ====
Tmap = MouseAssign(:,{'mouse_key','majority_cluster'});
Tmap.Properties.VariableNames{2} = 'ClusterMajority';
TsesMaj = outerjoin(Tses, Tmap, 'Keys','mouse_key', 'MergeKeys',true);

Dsum_inc   = daySummaryByCluster(TsesMaj, lickPM_col, 'ClusterMajority', ...
    MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS, false, TRANSITION_DAYS);
Dsum_clean = daySummaryByCluster(TsesMaj, lickPM_col, 'ClusterMajority', ...
    MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS, true,  TRANSITION_DAYS);

plotAcrossDay(Dsum_inc,   K, C, mk, outDir, ['across_day_majority_includingTransitions_' ts '.png'], ...
    'Across day (mouse-level majority clusters) — including transition days');

plotAcrossDay(Dsum_clean, K, C, mk, outDir, ['across_day_majority_EXCLUDINGTransitions_' ts '.png'], ...
    'Across day (mouse-level majority clusters) — EXCLUDING transition days');

%% ----- Build normalized progress summary (equal weight per session) -----
Ktrl = join(TtrialA, TsesCl(:,[keys {'Cluster'}]), 'Keys', keys);

Gsess = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx);
nT    = splitapply(@max, Ktrl.Trial, Gsess);
prog  = (Ktrl.Trial - 0.5) ./ nT(Gsess);

binEdges = linspace(0,1,NBINS+1);
[~,bin] = histc(prog, binEdges); %#ok<HISTC>
bin(bin<1)=1; bin(bin>NBINS)=NBINS;

% per-session per-bin mean lick rate
Gsb = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx, bin);
Ssb = table;
Ssb.mouse_key   = splitapply(@unique, Ktrl.mouse_key,   Gsb);
Ssb.day_index   = splitapply(@unique, Ktrl.day_index,   Gsb);
Ssb.session_idx = splitapply(@unique, Ktrl.session_idx, Gsb);
Ssb.Cluster     = splitapply(@unique, Ktrl.Cluster,     Gsb);
Ssb.bin         = splitapply(@unique, bin,              Gsb);
Ssb.meanLR      = splitapply(@mean,  Ktrl.(lickRate_col),    Gsb);

% cluster mean±SEM across sessions
Gcb = findgroups(Ssb.Cluster, Ssb.bin);
Csum = table;
Csum.Cluster = splitapply(@unique, Ssb.Cluster, Gcb);
Csum.bin     = splitapply(@unique, Ssb.bin,     Gcb);
Csum.nSess   = splitapply(@numel,  Ssb.meanLR,  Gcb);
Csum.mu      = splitapply(@mean,   Ssb.meanLR,  Gcb);
Csum.sem     = splitapply(@(x)nanstd(x,0)/sqrt(sum(isfinite(x))), Ssb.meanLR, Gcb);
xPct = (Csum.bin - 0.5) / NBINS * 100;

%% ===== Figure A: Normalized progress (left) + ID grids (right) =====
plotProgressWithIDs(Csum, xPct, Ssb, K, C, mk, outDir, ['progress_meanSEM_IDGRID_' ts '.png'], NBINS);

%% ===== Figure B: Across-day mean±SEM (session-level cluster) + ID grids =====
plotAcrossDayWithIDs(TsesCl, lickPM_col, K, C, mk, outDir, ['across_day_meanSEM_IDGRID_' ts '.png']);

%% ===== Period summaries: Active vs Passive (including vs excluding transition days) =====
plotPeriodActivePassive(Tses, lickPM_col, 'period', outDir, ...
    ['period_active_passive_includingTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

plotPeriodActivePassive(Tses, lickPM_col, 'period_clean', outDir, ...
    ['period_active_passive_EXCLUDINGTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

%% ===== Pairwise by PairID per period (including vs excluding transition days) =====
plotPairwiseByPeriod(Tses, lickPM_col, 'period', outDir, ...
    ['pairwise_active_passive_byPeriod_includingTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

plotPairwiseByPeriod(Tses, lickPM_col, 'period_clean', outDir, ...
    ['pairwise_active_passive_byPeriod_EXCLUDINGTransitions_' ts '.png'], MIN_DAY_FOR_ANALYSIS, MAX_DAY_FOR_ANALYSIS);

%% ----- Save a simple cluster ID list report (session-level clusters) -----
fid = fopen(fullfile(outDir, ['mouse_ids_by_cluster_' ts '.txt']),'w');
for k=1:K
    ids = sort(unique(TsesCl.mouse_key(TsesCl.Cluster==k)));
    fprintf(fid,'Cluster %d (%d mice):\n', k, numel(ids));
    fprintf(fid,'  %s\n\n', strjoin(ids, ', '));
end
fclose(fid);

writetable(groupsummary(TsesCl(:,{'mouse_key','Cluster'}),'Cluster'), ...
    fullfile(outDir, ['mouse_counts_per_cluster_' ts '.csv']));

disp('All done.');
end

%% ====================== FILE FINDERS ======================

function p = findOneFile(rootDir, fileName)
% Recursively find fileName under rootDir; if multiple, pick most recently modified.
dd = dir(fullfile(rootDir, '**', fileName));
if isempty(dd)
    p = '';
else
    [~,ix] = max([dd.datenum]);
    p = fullfile(dd(ix).folder, dd(ix).name);
end
end

function trialPath = findTrialCSV(runDir)
% Find a reasonable trial-level CSV automatically.
% Strategy:
% 1) Prefer exact known names if present
% 2) Else pick among per_trial*.csv, preferring ones containing "rate" then "feature"
preferred = { ...
    'per_trial_rates.csv', ...
    'per_trial_features.csv', ...
    'per_trial_features_rates.csv', ...
    'per_trial_lick_rates.csv' ...
};

trialPath = '';
for i=1:numel(preferred)
    p = findOneFile(runDir, preferred{i});
    if ~isempty(p), trialPath = p; return; end
end

% Otherwise, scan per_trial*.csv
dd = dir(fullfile(runDir, '**', 'per_trial*.csv'));
if isempty(dd), return; end

names = lower(string({dd.name}));
score = zeros(numel(dd),1);

% scoring: prefer "rate", then "lick", then "feature"
score(contains(names,'rate'))    = score(contains(names,'rate'))    + 3;
score(contains(names,'lick'))    = score(contains(names,'lick'))    + 2;
score(contains(names,'feature')) = score(contains(names,'feature')) + 1;

% tie-breaker: more recent datenum
[~,ix] = max(score + 1e-6 * [dd.datenum]');
trialPath = fullfile(dd(ix).folder, dd(ix).name);
end

%% ====================== COLUMN HELPERS ======================

function col = pickVar(T, candidates)
col = '';
for i=1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        col = candidates{i};
        return
    end
end
end

function T = tryComputeLickRate(T)
% Attempt to compute a lick-rate column if absent.
% Looks for combinations like:
%   licks / duration_s
%   n_licks / trial_dur_s
% If found, adds T.lick_rate (Hz)

countCol = pickVar(T, {'n_licks','licks','lick_count','nlicks','lick_total'});
durCol   = pickVar(T, {'duration_s','dur_s','trial_duration_s','trial_dur_s','dt_s','trial_len_s'});
if isempty(countCol) || isempty(durCol)
    return
end

c = double(T.(countCol));
d = double(T.(durCol));
d(d<=0 | ~isfinite(d)) = NaN;

T.lick_rate = c ./ d;
end

function Xz = zscoreNaN(X)
mu = nanmean(X,1);
sd = nanstd(X,0,1);
sd(sd==0 | ~isfinite(sd)) = 1;
Xz = (X - mu) ./ sd;
end

function ks = keyString(Tk)
mk = string(Tk.mouse_key);
d  = string(Tk.day_index);
s  = string(Tk.session_idx);
ks = normalizeMouseKey(mk) + "__D" + d + "__S" + s;
end

function s = normalizeMouseKey(sIn)
% Makes mouse_key robust across formats:
% - lower
% - spaces/hyphens -> underscore
% - handles "0911red" -> "0911_red"
s = string(sIn);
s = strtrim(lower(s));
s = regexprep(s,'[\-]+','_');
s = regexprep(s,'\s+','_');

s2 = s;
for i=1:numel(s)
    tok = regexp(s(i), '^(\d{4})([a-z]+)$', 'tokens','once');
    if ~isempty(tok)
        s2(i) = tok{1} + "_" + tok{2};
    end
end
s = s2;
end

%% ====================== COHORT + PERIODS ======================

function cohort = buildCohortTable()
% User-provided cohort mapping
% Format: mouse_key, Sex, ActPass, PairID

rows = {
'6100_red',      'F', 'Passive', 1
'6100_orange',   'F', 'Passive', 1
'6100_black',    'F', 'Active',  1

'0911_red',      'F', 'Active',  2
'0911_orange',   'F', 'Passive', 2

'0911_white',    'F', 'Active',  3
'0911_black',    'F', 'Passive', 3

'0910_black',    'M', 'Active',  4
'0910_orange',   'M', 'Passive', 4
'0910_red',      'M', 'Passive', 4

'6099_orange',   'M', 'Active',  5
'6099_red',      'M', 'Passive', 5

'6099_black',    'M', 'Active',  6
'6099_white',    'M', 'Passive', 6   % died day1-13: handled naturally by missing days
};

mouse_key = string(rows(:,1));
Sex      = categorical(string(rows(:,2)), {'F','M'});
ActPass  = categorical(string(rows(:,3)), {'Active','Passive'});
PairID   = cell2mat(rows(:,4));

Cage  = strings(size(mouse_key));
Color = strings(size(mouse_key));
for i=1:numel(mouse_key)
    parts = split(mouse_key(i), "_");
    Cage(i)  = parts(1);
    Color(i) = parts(2);
end

cohort = table(mouse_key, Cage, Color, Sex, ActPass, PairID);
end

function T = attachCohortAndPeriods(T, cohort, transitionDays)
if ismember('mouse_key', T.Properties.VariableNames) && ~isstring(T.mouse_key)
    T.mouse_key = string(T.mouse_key);
end

mk = normalizeMouseKey(T.mouse_key);
ck = normalizeMouseKey(cohort.mouse_key);
[tf,loc] = ismember(mk, ck);

% default placeholders
T.Cage    = repmat("", height(T),1);
T.Color   = repmat("", height(T),1);
T.Sex     = categorical(repmat("U",height(T),1), {'F','M','U'});
T.ActPass = categorical(repmat("U",height(T),1), {'Active','Passive','U'});
T.PairID  = nan(height(T),1);

% fill from cohort where matched
T.Cage(tf)    = cohort.Cage(loc(tf));
T.Color(tf)   = cohort.Color(loc(tf));
T.Sex(tf)     = cohort.Sex(loc(tf));
T.ActPass(tf) = cohort.ActPass(loc(tf));
T.PairID(tf)  = cohort.PairID(loc(tf));

% periods
d = double(T.day_index);
T.period       = dayToPeriod(d, false, transitionDays);
T.period_clean = dayToPeriod(d, true,  transitionDays);
end

function cats = periodCats(withExclude)
base = ["Pre_D3-5","During_D6-10","Post_D11-13","Withdrawal_D14-16","Reexposure_D17-18","<undef>"];
if withExclude
    cats = [base(1:5), "<exclude>", base(6)];
else
    cats = base;
end
cats = cellstr(cats);
end

function P = dayToPeriod(day_index, excludeTransitions, transitionDays)
lab = strings(size(day_index)); lab(:) = "<undef>";
d = double(day_index);

lab(d>=3  & d<=5 )  = "Pre_D3-5";
lab(d>=6  & d<=10 ) = "During_D6-10";
lab(d>=11 & d<=13 ) = "Post_D11-13";
lab(d>=14 & d<=16 ) = "Withdrawal_D14-16";
lab(d>=17 & d<=18 ) = "Reexposure_D17-18";

if excludeTransitions
    lab(ismember(d, transitionDays)) = "<exclude>";
    P = categorical(lab, periodCats(true));
else
    P = categorical(lab, periodCats(false));
end
end

%% ====================== SUMMARIES ======================

function Dsum = daySummaryByCluster(T, lickCol, clusterCol, dayMin, dayMax, excludeTransitions, transitionDays)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if excludeTransitions
    T2 = T2(~ismember(double(T2.day_index), transitionDays), :);
end

cl = T2.(clusterCol);
if iscategorical(cl)
    clNum = double(cl);
elseif iscell(cl) || isstring(cl)
    clNum = double(categorical(string(cl)));
else
    clNum = double(cl);
end

keep = isfinite(clNum) & isfinite(double(T2.(lickCol))) & isfinite(double(T2.day_index));
T2 = T2(keep,:);
clNum = clNum(keep);

Gd = findgroups(T2.day_index, clNum);

Dsum = table;
Dsum.day_index = splitapply(@unique, T2.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, clNum,        Gd);
Dsum.mu        = splitapply(@mean,  T2.(lickCol),  Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), T2.(lickCol), Gd);
end

%% ====================== PLOTS ======================

function plotAcrossDay(Dsum, K, C, mk, outDir, fname, ttl)
f = figure('Color','w','Position',[60 60 900 520]); hold on
for k=1:K
    r = (Dsum.Cluster==k);
    if ~any(r), continue; end
    [d,ord] = sort(Dsum.day_index(r));
    mu = Dsum.mu(r);  mu = mu(ord);
    se = Dsum.sem(r); se = se(ord);

    plot(d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel('Day'); ylabel('Licks / min'); title(ttl);
grid on; legend('off');
exportgraphics(f, fullfile(outDir, fname), 'Resolution', 300);
close(f);
end

function plotProgressWithIDs(Csum, xPct, Ssb, K, C, mk, outDir, fname, NBINS)
fA = figure('Color','w','Position',[60 60 1400 600]);
tlo = tiledlayout(fA,1,2,'TileSpacing','compact','Padding','compact');

% Left: curves
ax1 = nexttile(tlo,1); hold(ax1,'on');
for k=1:K
    r = (Csum.Cluster==k);
    if ~any(r), continue; end
    [xx,ord] = sort(xPct(r));
    mu = Csum.mu(r);  mu = mu(ord);
    se = Csum.sem(r); se = se(ord);

    plot(ax1, xx, mu, '-', 'Color', C(k,:), 'LineWidth', 1.6, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerSize',3, 'MarkerFaceColor',C(k,:));

    ii = 1:max(1,round(numel(xx)/50)):numel(xx);
    eb = errorbar(ax1, xx(ii), mu(ii), se(ii), 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel(ax1,'Session progress (%)'); ylabel(ax1,'Lick rate (Hz)');
title(ax1, sprintf('Clusters: mean \\pm SEM across normalized session progress (%d bins)', NBINS));
grid(ax1,'on'); legend(ax1,'off');

% Right: IDs per cluster
axIDs = nexttile(tlo,2); axis(axIDs,'off'); title(axIDs,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;
for k=1:K
    ids = sort(unique(Ssb.mouse_key(Ssb.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);
    ncol = numel(txtBlocks);
    x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end

exportgraphics(fA, fullfile(outDir,fname), 'Resolution', 300);
close(fA);
end

function plotAcrossDayWithIDs(TsesCl, lickPM_col, K, C, mk, outDir, fname)
% aggregate by day x cluster
Gd = findgroups(TsesCl.day_index, TsesCl.Cluster);
Dsum = table;
Dsum.day_index = splitapply(@unique, TsesCl.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, TsesCl.Cluster,   Gd);
Dsum.mu        = splitapply(@mean,  TsesCl.(lickPM_col), Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), TsesCl.(lickPM_col), Gd);

fB = figure('Color','w','Position',[60 60 1400 600]);
t2 = tiledlayout(fB,1,2,'TileSpacing','compact','Padding','compact');

ax2 = nexttile(t2,1); hold(ax2,'on');
for k=1:K
    r = Dsum.Cluster==k;
    if ~any(r), continue; end
    [d,ord] = sort(Dsum.day_index(r));
    mu = Dsum.mu(r);  mu = mu(ord);
    se = Dsum.sem(r); se = se(ord);

    plot(ax2, d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(ax2, d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off');
end
xlabel(ax2,'Day'); ylabel(ax2,'Licks / min');
title(ax2,'Clusters: mean \pm SEM across day (session-level clusters)');
grid(ax2,'on'); legend(ax2,'off');

axIDs2 = nexttile(t2,2); axis(axIDs2,'off'); title(axIDs2,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;
for k=1:K
    ids = sort(unique(TsesCl.mouse_key(TsesCl.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);
    ncol = numel(txtBlocks); x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs2, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end

exportgraphics(fB, fullfile(outDir,fname), 'Resolution', 300);
close(fB);
end

function plotPeriodActivePassive(T, lickCol, periodVar, outDir, fname, dayMin, dayMax)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if ~ismember(periodVar, T2.Properties.VariableNames), return; end

P = T2.(periodVar);
keep = (P ~= "<undef>");
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P ~= "<exclude>");
end
T2 = T2(keep,:); P = T2.(periodVar);

% Mouse-level period average (each mouse contributes equally)
G = findgroups(T2.mouse_key, P, T2.ActPass);
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.period    = splitapply(@unique, P, G);
S.ActPass   = splitapply(@unique, T2.ActPass, G);
S.mu_mouse  = splitapply(@mean,  T2.(lickCol), G);

% Aggregate across mice
G2 = findgroups(S.period, S.ActPass);
Agg = table;
Agg.period  = splitapply(@unique, S.period, G2);
Agg.ActPass = splitapply(@unique, S.ActPass, G2);
Agg.mu      = splitapply(@mean,  S.mu_mouse, G2);
Agg.sem     = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), S.mu_mouse, G2);

periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(Agg.period)));
if isempty(periods), return; end

f = figure('Color','w','Position',[80 80 980 420]);
ax = axes(f); hold(ax,'on');

x = 1:numel(periods);
xA = x - 0.12; xP = x + 0.12;

muA = nan(1,numel(periods)); seA = muA;
muP = nan(1,numel(periods)); seP = muP;

for i=1:numel(periods)
    rrA = (Agg.ActPass=="Active")  & (Agg.period==periods{i});
    rrP = (Agg.ActPass=="Passive") & (Agg.period==periods{i});
    if any(rrA), muA(i)=Agg.mu(rrA); seA(i)=Agg.sem(rrA); end
    if any(rrP), muP(i)=Agg.mu(rrP); seP(i)=Agg.sem(rrP); end
end

plot(ax, xA, muA, '-o','LineWidth',1.8);
errorbar(ax, xA, muA, seA, 'LineStyle','none');

plot(ax, xP, muP, '-o','LineWidth',1.8);
errorbar(ax, xP, muP, seP, 'LineStyle','none');

set(ax,'XTick',x,'XTickLabel',periods); xtickangle(ax,25);
ylabel(ax,'Licks / min');
title(ax, sprintf('Period summary (Active vs Passive) — %s', periodVar), 'Interpreter','none');
grid(ax,'on');
legend(ax, {'Active','Active SEM','Passive','Passive SEM'}, 'Location','best');

exportgraphics(f, fullfile(outDir,fname), 'Resolution', 300);
close(f);
end

function plotPairwiseByPeriod(T, lickCol, periodVar, outDir, fname, dayMin, dayMax)
T2 = T(T.day_index>=dayMin & T.day_index<=dayMax, :);
if ~ismember(periodVar, T2.Properties.VariableNames), return; end

P = T2.(periodVar);
keep = (P ~= "<undef>");
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P ~= "<exclude>");
end
T2 = T2(keep,:);

% Mouse-level period mean (so each mouse contributes one value per period)
G = findgroups(T2.mouse_key, T2.PairID, T2.ActPass, T2.(periodVar));
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.PairID    = splitapply(@unique, T2.PairID,    G);
S.ActPass   = splitapply(@unique, T2.ActPass,   G);
S.period    = splitapply(@unique, T2.(periodVar), G);
S.mu_mouse  = splitapply(@mean,  T2.(lickCol),  G);

periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(S.period)));
if isempty(periods), return; end

nP = numel(periods);
f = figure('Color','w','Position',[80 80 320*nP 380]);
tlo = tiledlayout(f,1,nP,'TileSpacing','compact','Padding','compact');

for ip = 1:nP
    ax = nexttile(tlo, ip); hold(ax,'on');
    pr = periods{ip};
    Sp = S(S.period==pr & isfinite(S.mu_mouse) & isfinite(S.PairID), :);

    pairIDs = unique(Sp.PairID);

    % line per passive mouse to its active partner (supports 1 active vs 2 passive)
    for i=1:numel(pairIDs)
        pid = pairIDs(i);
        A = Sp(Sp.PairID==pid & Sp.ActPass=="Active", :);
        Pm = Sp(Sp.PairID==pid & Sp.ActPass=="Passive", :);
        if isempty(A) || isempty(Pm), continue; end
        aVal = A.mu_mouse(1);
        for j=1:height(Pm)
            pVal = Pm.mu_mouse(j);
            plot(ax, [1 2], [aVal pVal], '-', 'LineWidth',1.2);
        end
    end

    aPts = Sp.mu_mouse(Sp.ActPass=="Active");
    pPts = Sp.mu_mouse(Sp.ActPass=="Passive");
    scatter(ax, 1 + 0.06*(rand(numel(aPts),1)-0.5), aPts, 28, 'filled', 'MarkerFaceAlpha',0.6);
    scatter(ax, 2 + 0.06*(rand(numel(pPts),1)-0.5), pPts, 28, 'filled', 'MarkerFaceAlpha',0.6);

    xlim(ax,[0.6 2.4]);
    set(ax,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
    ylabel(ax,'Licks / min');
    title(ax, pr, 'Interpreter','none');
    grid(ax,'on'); box(ax,'on');
end

title(tlo, sprintf('Pairwise Active vs Passive (by PairID) — %s', periodVar), 'Interpreter','none');
exportgraphics(f, fullfile(outDir,fname), 'Resolution', 300);
close(f);
end

%% ====================== ID GRID HELPER ======================

function blocks = idColumns(ids, maxRows)
% Split a long list of IDs into columns of ~maxRows lines.
ids = string(ids);
n = numel(ids);
ncol = max(1, ceil(n/maxRows));
nrow = ceil(n/ncol);

blocks = cell(1,ncol);
for c = 1:ncol
    i1 = (c-1)*nrow + 1;
    i2 = min(c*nrow, n);
    blocks{c} = strjoin(ids(i1:i2), newline);
end
end
