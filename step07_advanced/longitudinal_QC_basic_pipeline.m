function longitudinal_QC_basic_pipeline()
% longitudinal_QC_basic_pipeline
% 목적:
%   “분석 들어가기 전에” 데이터 품질/이상치/처리 문제를 빠르게 점검하는 QC + 기본 분석 파이프라인.
%   데이터가 커서 원본을 못 주는 상황을 가정하고, longitudinal_outputs/run_*/ 아래의
%   summary CSV/feature CSV들을 자동으로 찾아서 최소한의 sanity check를 수행합니다.
%
% 핵심 기능:
%   0) 최신 run_* 자동 탐색
%   1) per_session_features.csv / per_trial*.csv 자동 탐색(파일명이 달라도 최대한 찾음)
%   2) 사용자 cohort/Pair/Active-Passive/성별 매핑 적용
%   3) 사용자 period 정의 적용:
%        Pre: D3-5
%        During: D6-10
%        Post: D11-13
%        Withdrawal: D14-16
%        Reexposure: D17-18
%      + transition day 제외 버전(period_clean)도 같이 생성
%   4) QC/이상치 체크:
%      - 누락 일수/세션 수 히트맵
%      - day별 lick/min 분포(Active vs Passive)
%      - 누적 lick(일별 누적/전체 누적)
%      - period 비교(transition day 포함/제외)
%      - Pair-wise(Active vs Passive) 연결 plot
%      - outlier 후보(robust z-score) 리스트 출력
%   5) pupil/event 관련:
%      - (선택) pupil_day_summary.csv 같은 요약이 있으면 mouse별 baseline(D3) 정규화 후 day/period 비교
%      - (선택) event-aligned pupil이 있으면(파일 존재할 때만) reward vs nonreward, lick-bout 기준 분석까지 “틀” 제공
%
% 주의:
%   - Pupil event-aligned 원본(프레임타임/이벤트타임)이 어떤 이름/형식인지 사용자 환경마다 달라서
%     이 스크립트는 “파일이 있으면 자동 처리”, 없으면 스킵합니다.
%   - lick “reward-causing vs nonreward”는 trial table에 reward flag 컬럼이 있어야 자동으로 됩니다.
%
% 사용법:
%   1) 이 파일명을 longitudinal_QC_basic_pipeline.m 로 저장
%   2) MATLAB에서 실행: longitudinal_QC_basic_pipeline
%
% 출력:
%   <latest run에서 session CSV가 있는 폴더>\QC_basic\
%     - figs/*.png
%     - tables/*.csv
%     - QC_report.txt
%
% 작성: ChatGPT (사용자 cohort/기간 규칙 반영)

%% ================= USER SETTINGS =================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% 분석에서 제외할 habituation days
ANALYSIS_DAY_MIN = 3;
ANALYSIS_DAY_MAX = 18;

% transition day (사용자 언급)
TRANSITION_DAYS = [4 6 11 14 17];

% lick bout 정의 (연속 lick 사이의 gap이 이 값보다 짧으면 같은 bout)
BOUT_GAP_SEC = 2.0;

% pupil event window (sec)
PUPIL_PRE_SEC  = 2.0;
PUPIL_POST_SEC = 2.0;

% outlier 탐지(robust z-score threshold)
OUTLIER_Z = 3.5;

% normalized progress curve bins
NBINS = 100;

rng(1);

%% ================= LOCATE LATEST RUN =================
runDir = findLatestRun(rootTry);
assert(~isempty(runDir), 'No run_* found under %s', rootTry);

% find session/trial CSVs (robust)
sessPath  = findOneFile(runDir, 'per_session_features.csv');
trialPath = findTrialCSV(runDir);

assert(~isempty(sessPath), 'Could not find per_session_features.csv under %s', runDir);
assert(~isempty(trialPath), 'Could not find per_trial*.csv under %s', runDir);

dataDir = fileparts(sessPath);
outRoot = fullfile(dataDir,'QC_basic');
figDir  = fullfile(outRoot,'figs');
tabDir  = fullfile(outRoot,'tables');
if ~exist(figDir,'dir'), mkdir(figDir); end
if ~exist(tabDir,'dir'), mkdir(tabDir); end

ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
reportPath = fullfile(outRoot, ['QC_report_' ts '.txt']);
logf = fopen(reportPath,'w');

fprintf('RUN: %s\n', runDir);
fprintf('SESSION: %s\n', sessPath);
fprintf('TRIAL: %s\n', trialPath);
fprintf('OUT: %s\n\n', outRoot);

fprintf(logf,'RUN: %s\nSESSION: %s\nTRIAL: %s\nOUT: %s\n\n', runDir, sessPath, trialPath, outRoot);

%% ================= LOAD TABLES =================
Tses   = readtable(sessPath);
Ttrial = readtable(trialPath);

keys = {'mouse_key','day_index','session_idx'};
assert(all(ismember(keys, Tses.Properties.VariableNames)), 'per_session_features missing keys');
assert(all(ismember(keys, Ttrial.Properties.VariableNames)), 'trial CSV missing keys: mouse_key/day_index/session_idx');

% normalize mouse_key format
Tses.mouse_key   = normalizeMouseKey(Tses.mouse_key);
Ttrial.mouse_key = normalizeMouseKey(Ttrial.mouse_key);

% choose lick/min column
lickPM_col = pickVar(Tses, {'lick_per_min','licks_per_min','lick_per_minute','licks_per_minute'});
assert(~isempty(lickPM_col), 'Cannot find licks/min column in per_session_features.csv');

% choose trial lick rate column (Hz)
lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_hz','rate_hz'});
if isempty(lickRate_col)
    Ttrial = tryComputeLickRate(Ttrial);
    lickRate_col = pickVar(Ttrial, {'lick_rate','lickrate','lick_rate_hz','lick_hz','rate_hz'});
end
assert(~isempty(lickRate_col), 'Cannot find or compute trial lick-rate column in trial CSV.');
assert(ismember('Trial', Ttrial.Properties.VariableNames), 'Trial CSV missing "Trial" column.');

%% ================= APPLY COHORT + PERIODS =================
cohort = buildCohortTable(); % includes PairID, Sex, ActPass
Tses   = attachCohortAndPeriods(Tses, cohort, TRANSITION_DAYS);
Ttrial = attachCohortAndPeriods(Ttrial, cohort, TRANSITION_DAYS);

% Restrict to day window for most QC plots (but keep full in memory)
TsesA   = Tses(Tses.day_index>=ANALYSIS_DAY_MIN & Tses.day_index<=ANALYSIS_DAY_MAX, :);
TtrialA = Ttrial(Ttrial.day_index>=ANALYSIS_DAY_MIN & Ttrial.day_index<=ANALYSIS_DAY_MAX, :);

% Save cohort merge coverage
mergeCoverage = table;
mergeCoverage.mouse_key = unique(Tses.mouse_key,'stable');
[tf,~] = ismember(mergeCoverage.mouse_key, cohort.mouse_key);
mergeCoverage.in_cohort_map = tf;
writetable(mergeCoverage, fullfile(tabDir, ['cohort_merge_coverage_' ts '.csv']));
fprintf(logf, 'Cohort map coverage: %d/%d mice matched\n\n', sum(tf), numel(tf));

%% ================= BASIC QC 1: MISSINGNESS HEATMAP =================
% sessions per day per mouse
H = sessionsHeatmap(TsesA);
fig = plotHeatmap(H);
savePNG(fig, fullfile(figDir, ['QC_missing_sessions_heatmap_' ts '.png']));

fprintf(logf,'[QC] Missingness heatmap saved.\n');

%% ================= BASIC QC 2: LICKS/MIN DISTRIBUTIONS BY DAY =================
% Active vs Passive per day (box+jitter)
fig = plotDayDistributions(TsesA, lickPM_col);
savePNG(fig, fullfile(figDir, ['QC_licksPerMin_byDay_ActPass_' ts '.png']));
fprintf(logf,'[QC] Licks/min distributions by day saved.\n');

%% ================= BASIC QC 3: CUMULATIVE LICKS (ACTIVE vs PASSIVE) =================
% Use per_session licks/min as proxy for "amount" (if session duration available use true licks)
% If duration exists, compute licks = licks/min * minutes
Tcum = computeDailyCumulative(TsesA, lickPM_col);
writetable(Tcum, fullfile(tabDir, ['daily_cumulative_' ts '.csv']));

fig = plotCumulative(Tcum);
savePNG(fig, fullfile(figDir, ['QC_cumulative_Act_vs_Passive_' ts '.png']));
fprintf(logf,'[QC] Cumulative plot saved.\n');

%% ================= BASIC QC 4: PERIOD COMPARISON (WITH/WITHOUT TRANSITION DAYS) =================
% Requested examples:
% Active: D5 vs D6-10 vs D11-13 vs D14-16 vs D17
% Passive: D5 vs D11-13 vs D14-16 vs D17-18
% We do:
%   - period means (mouse-level average so each mouse weight equal)
%   - two versions: period vs period_clean
%   - additionally "custom bins" for your requested comparison
%
[PerTbl_incl, PerTbl_excl] = periodSummaryActPass(TsesA, lickPM_col);
writetable(PerTbl_incl, fullfile(tabDir, ['period_summary_includingTransitions_' ts '.csv']));
writetable(PerTbl_excl, fullfile(tabDir, ['period_summary_EXCLUDINGTransitions_' ts '.csv']));

fig = plotPeriodSummary(PerTbl_incl, 'Including transition days');
savePNG(fig, fullfile(figDir, ['QC_period_summary_includingTransitions_' ts '.png']));
fig = plotPeriodSummary(PerTbl_excl, 'EXCLUDING transition days');
savePNG(fig, fullfile(figDir, ['QC_period_summary_EXCLUDINGTransitions_' ts '.png']));

% Custom bins per your request (Active vs Passive different bins)
[Cust_incl, Cust_excl] = customPeriodSummaryRequested(TsesA, lickPM_col, TRANSITION_DAYS);
writetable(Cust_incl, fullfile(tabDir, ['custom_period_requested_includingTransitions_' ts '.csv']));
writetable(Cust_excl, fullfile(tabDir, ['custom_period_requested_EXCLUDINGTransitions_' ts '.csv']));

fig = plotCustomPeriodSummary(Cust_incl, 'Including transition days');
savePNG(fig, fullfile(figDir, ['QC_custom_period_requested_includingTransitions_' ts '.png']));
fig = plotCustomPeriodSummary(Cust_excl, 'EXCLUDING transition days');
savePNG(fig, fullfile(figDir, ['QC_custom_period_requested_EXCLUDINGTransitions_' ts '.png']));

fprintf(logf,'[QC] Period comparisons saved.\n');

%% ================= BASIC QC 5: PAIRWISE (ACTIVE vs PASSIVE) BY PERIOD =================
fig = plotPairwiseByPeriod(TsesA, lickPM_col, 'period');
savePNG(fig, fullfile(figDir, ['QC_pairwise_byPeriod_includingTransitions_' ts '.png']));
fig = plotPairwiseByPeriod(TsesA, lickPM_col, 'period_clean');
savePNG(fig, fullfile(figDir, ['QC_pairwise_byPeriod_EXCLUDINGTransitions_' ts '.png']));
fprintf(logf,'[QC] Pairwise A vs P plots saved.\n');

%% ================= BASIC QC 6: OUTLIER CANDIDATES =================
% robust z-score on licks/min within each (day x ActPass) group
Out = findOutliersRobust(TsesA, lickPM_col, OUTLIER_Z);
writetable(Out, fullfile(tabDir, ['outlier_candidates_' ts '.csv']));
fprintf(logf,'[QC] Outlier candidates written: %d rows (z>|%.1f|)\n', height(Out), OUTLIER_Z);

%% ================= OPTIONAL QC 7: TRIAL-LEVEL NORMALIZED PROGRESS CURVES =================
% quick check: are there weird sessions/flat sessions?
fig = plotNormalizedProgressCurves(TtrialA, TsesA, keys, lickRate_col, NBINS);
savePNG(fig, fullfile(figDir, ['QC_normalized_progress_curves_' ts '.png']));
fprintf(logf,'[QC] Normalized progress curves saved.\n');

%% ================= OPTIONAL: PUPIL SUMMARY (IF EXISTS) =================
% 1) day-level pupil comparison morphine vs water:
%    - normalize pupil per mouse by baseline day3 (or day3 mean)
%    - compare day4-5 vs day12-13 (your request)
%
pupilDayPath = findAnyPupilDaySummary(runDir);
if ~isempty(pupilDayPath)
    fprintf('Found pupil day summary: %s\n', pupilDayPath);
    fprintf(logf,'Found pupil day summary: %s\n', pupilDayPath);

    Pup = readtable(pupilDayPath);
    Pup = sanitizePupilDayTable(Pup);
    Pup.mouse_key = normalizeMouseKey(Pup.mouse_key);
    Pup = attachCohortAndPeriods(Pup, cohort, TRANSITION_DAYS);

    PupN = normalizePupilByBaselineDay(Pup, 3); % baseline day3
    writetable(PupN, fullfile(tabDir, ['pupil_day_normalized_' ts '.csv']));

    fig = plotPupilDayCompare(PupN);
    savePNG(fig, fullfile(figDir, ['QC_pupil_day_compare_D45_vs_D1213_normByD3_' ts '.png']));

    fprintf(logf,'[QC] Pupil day-level normalized comparisons saved.\n');
else
    fprintf('No pupil day summary found. (Skipping pupil day QC)\n');
    fprintf(logf,'No pupil day summary found. (Skipping pupil day QC)\n');
end

%% ================= OPTIONAL: EVENT-ALIGNED PUPIL (IF EXISTS) =================
% This section is “best effort”:
%   - If you have an event-aligned pupil CSV already produced by your pipeline, we will use it.
%   - If not found, we skip.
%
pupilEventPath = findAnyEventAlignedPupil(runDir);
if ~isempty(pupilEventPath)
    fprintf('Found event-aligned pupil file: %s\n', pupilEventPath);
    fprintf(logf,'Found event-aligned pupil file: %s\n', pupilEventPath);

    E = readtable(pupilEventPath);
    % expected minimal columns (adaptable):
    % mouse_key, day_index, event_time, pupil, event_type(optional), reward_flag(optional), lick_time(optional)
    % We will not crash if extra columns exist; we will attempt to detect key names.
    try
        [BoutTbl, fig1, fig2] = eventAlignedPupilBoutQC(E, cohort, ...
            PUPIL_PRE_SEC, PUPIL_POST_SEC, BOUT_GAP_SEC, TRANSITION_DAYS);
        writetable(BoutTbl, fullfile(tabDir, ['pupil_event_bout_summary_' ts '.csv']));
        savePNG(fig1, fullfile(figDir, ['QC_pupil_event_aligned_reward_vs_nonreward_' ts '.png']));
        savePNG(fig2, fullfile(figDir, ['QC_pupil_event_aligned_period_focus_passiveD6_10_' ts '.png']));
        fprintf(logf,'[QC] Event-aligned pupil bout QC saved.\n');
    catch ME
        fprintf('Event-aligned pupil QC failed: %s\n', ME.message);
        fprintf(logf,'Event-aligned pupil QC failed: %s\n', ME.message);
    end
else
    fprintf('No event-aligned pupil file found. (Skipping event-aligned pupil QC)\n');
    fprintf(logf,'No event-aligned pupil file found. (Skipping event-aligned pupil QC)\n');
end

%% ================= REPORT: REMINDERS FOR NON-LICK BEHAVIOR DAYS =================
fprintf(logf,'\nNotes:\n');
fprintf(logf,'- Days 11-13: you said “not licking” but have other measures (pupil, tail immersion, TST, hot plate).\n');
fprintf(logf,'  This pipeline only auto-QC licks/pupil if pupil files exist. Add your assay CSVs into run folder to QC them similarly.\n');
fprintf(logf,'- Withdrawal: Day14 high licking, Day15-16 more reliable. Use period_clean excluding day14 if desired.\n');
fprintf(logf,'- Passive mice during Day6-10: reward is forced; compare pupil reward-aligned responses vs Post/other days.\n');

fclose(logf);

fprintf('\nDONE.\nQC outputs in:\n  %s\nReport:\n  %s\n', outRoot, reportPath);

end

%% =====================================================================
%% ============================= HELPERS ===============================
%% =====================================================================

function runDir = findLatestRun(rootTry)
runDir = '';
if ~exist(rootTry,'dir')
    here=pwd; cand=here;
    for up=1:5
        p=fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry=p; break; end
        cand=fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*'));
if isempty(D), return; end
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);
end

function p = findOneFile(rootDir, fileName)
dd = dir(fullfile(rootDir, '**', fileName));
if isempty(dd), p = ''; return; end
[~,ix] = max([dd.datenum]);
p = fullfile(dd(ix).folder, dd(ix).name);
end

function trialPath = findTrialCSV(runDir)
% robust search for trial-level file
preferred = {'per_trial_rates.csv','per_trial_features.csv','per_trial_lick_rates.csv'};
trialPath = '';
for i=1:numel(preferred)
    p = findOneFile(runDir, preferred{i});
    if ~isempty(p), trialPath = p; return; end
end
dd = dir(fullfile(runDir,'**','per_trial*.csv'));
if isempty(dd), return; end
names = lower(string({dd.name}));
score = zeros(numel(dd),1);
score(contains(names,'rate'))    = score(contains(names,'rate')) + 3;
score(contains(names,'lick'))    = score(contains(names,'lick')) + 2;
score(contains(names,'feature')) = score(contains(names,'feature')) + 1;
[~,ix] = max(score + 1e-6*[dd.datenum]');
trialPath = fullfile(dd(ix).folder, dd(ix).name);
end

function col = pickVar(T, candidates)
col = '';
for i=1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        col = candidates{i};
        return
    end
end
end

function s = normalizeMouseKey(sIn)
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

function T = tryComputeLickRate(T)
countCol = pickVar(T, {'n_licks','licks','lick_count','nlicks','lick_total'});
durCol   = pickVar(T, {'duration_s','dur_s','trial_duration_s','trial_dur_s','dt_s','trial_len_s'});
if isempty(countCol) || isempty(durCol), return; end
c = double(T.(countCol));
d = double(T.(durCol));
d(d<=0 | ~isfinite(d)) = NaN;
T.lick_rate = c ./ d;
end

function cohort = buildCohortTable()
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
'6099_white',    'M', 'Passive', 6   % died day1-13 (missing days naturally handled)
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
if ismember('mouse_key', T.Properties.VariableNames)
    T.mouse_key = normalizeMouseKey(T.mouse_key);
else
    error('Table missing mouse_key');
end

mk = T.mouse_key;
[tf,loc] = ismember(mk, cohort.mouse_key);

T.Cage    = repmat("", height(T),1);
T.Color   = repmat("", height(T),1);
T.Sex     = categorical(repmat("U",height(T),1), {'F','M','U'});
T.ActPass = categorical(repmat("U",height(T),1), {'Active','Passive','U'});
T.PairID  = nan(height(T),1);

T.Cage(tf)    = cohort.Cage(loc(tf));
T.Color(tf)   = cohort.Color(loc(tf));
T.Sex(tf)     = cohort.Sex(loc(tf));
T.ActPass(tf) = cohort.ActPass(loc(tf));
T.PairID(tf)  = cohort.PairID(loc(tf));

d = double(T.day_index);
T.period       = dayToPeriod(d, false, transitionDays);
T.period_clean = dayToPeriod(d, true,  transitionDays);
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
    P = categorical(lab, {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18','<exclude>','<undef>'});
else
    P = categorical(lab, {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18','<undef>'});
end
end

function H = sessionsHeatmap(TsesA)
mice = unique(TsesA.mouse_key,'stable');
days = unique(TsesA.day_index);
days = sort(days);
M = zeros(numel(mice), numel(days));
for i=1:numel(mice)
    for j=1:numel(days)
        M(i,j) = sum(TsesA.mouse_key==mice(i) & TsesA.day_index==days(j));
    end
end
H = struct('mice',mice,'days',days,'M',M,'ActPass',[]);
end

function fig = plotHeatmap(H)
fig = figure('Color','w','Position',[80 80 1200 520]);
imagesc(H.M);
colormap(parula);
colorbar;
set(gca,'XTick',1:numel(H.days),'XTickLabel',H.days);
xlabel('Day'); ylabel('Mouse');
yticks(1:numel(H.mice));
yticklabels(H.mice);
title('QC: #sessions per mouse per day (missingness heatmap)');
set(gca,'TickLabelInterpreter','none');
end

function fig = plotDayDistributions(TsesA, lickPM_col)
fig = figure('Color','w','Position',[80 80 1200 520]);
hold on;
days = unique(TsesA.day_index); days=sort(days);
x = []; y = []; g = [];
for i=1:numel(days)
    d = days(i);
    r = TsesA.day_index==d & isfinite(double(TsesA.(lickPM_col)));
    if ~any(r), continue; end
    xi = i + 0.25*(rand(sum(r),1)-0.5);
    yi = double(TsesA.(lickPM_col)(r));
    gi = string(TsesA.ActPass(r));
    scatter(xi, yi, 20, 'filled','MarkerFaceAlpha',0.35);
    x = [x; xi]; y = [y; yi]; g = [g; gi];
end
set(gca,'XTick',1:numel(days),'XTickLabel',days);
xlabel('Day'); ylabel('Licks / min');
title('QC: Licks/min by day (jittered; color not enforced)');
grid on; box on;
% annotate medians by ActPass
for i=1:numel(days)
    d=days(i);
    for ap = ["Active","Passive"]
        rr = TsesA.day_index==d & TsesA.ActPass==ap & isfinite(double(TsesA.(lickPM_col)));
        if any(rr)
            med = median(double(TsesA.(lickPM_col)(rr)));
            text(i, med, sprintf('%s med=%.1f', ap, med), 'FontSize',7, 'Interpreter','none');
        end
    end
end
end

function Tcum = computeDailyCumulative(TsesA, lickPM_col)
% if duration exists use it; else proxy with licks/min
durMinCol = pickVar(TsesA, {'session_minutes','duration_min','session_duration_min','minutes'});
if ~isempty(durMinCol)
    licks = double(TsesA.(lickPM_col)) .* double(TsesA.(durMinCol));
else
    licks = double(TsesA.(lickPM_col)); % proxy
end

G = findgroups(TsesA.mouse_key, TsesA.day_index, TsesA.ActPass);
S = table;
S.mouse_key = splitapply(@unique, TsesA.mouse_key, G);
S.day_index = splitapply(@unique, TsesA.day_index, G);
S.ActPass   = splitapply(@unique, TsesA.ActPass, G);
S.licks_est = splitapply(@nansum, licks, G);

% daily mean across mice for each ActPass
G2 = findgroups(S.day_index, S.ActPass);
Tcum = table;
Tcum.day_index = splitapply(@unique, S.day_index, G2);
Tcum.ActPass   = splitapply(@unique, S.ActPass, G2);
Tcum.mu        = splitapply(@mean,  S.licks_est, G2);
Tcum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), S.licks_est, G2);

% cumulative across days
acts = categories(Tcum.ActPass);
Tcum.cum = nan(height(Tcum),1);
for a=1:numel(acts)
    rr = Tcum.ActPass==acts{a};
    [d,ord]=sort(double(Tcum.day_index(rr)));
    cumv = cumsum(Tcum.mu(rr));
    tmp = nan(sum(rr),1); tmp(ord)=cumv;
    Tcum.cum(rr)=tmp;
end
end

function fig = plotCumulative(Tcum)
fig = figure('Color','w','Position',[80 80 900 520]);
hold on;
acts = categories(Tcum.ActPass);
for a=1:numel(acts)
    rr = Tcum.ActPass==acts{a};
    [d,ord]=sort(double(Tcum.day_index(rr)));
    plot(d, Tcum.cum(rr), '-o', 'LineWidth',1.8);
end
xlabel('Day'); ylabel('Cumulative (proxy) licks');
title('Cumulative licks (days>=3) Active vs Passive');
grid on; box on;
legend(acts,'Location','best');
end

function [Per_incl, Per_excl] = periodSummaryActPass(TsesA, lickPM_col)
Per_incl = periodSummaryCore(TsesA, lickPM_col, 'period');
Per_excl = periodSummaryCore(TsesA, lickPM_col, 'period_clean');
end

function Per = periodSummaryCore(T, lickPM_col, periodVar)
P = T.(periodVar);
keep = P~="<undef>";
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P~="<exclude>");
end
T2 = T(keep,:);
P2 = T2.(periodVar);

% mouse-level average first (equal weight per mouse)
G = findgroups(T2.mouse_key, P2, T2.ActPass);
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.period    = splitapply(@unique, P2, G);
S.ActPass   = splitapply(@unique, T2.ActPass, G);
S.mu_mouse  = splitapply(@mean,  T2.(lickPM_col), G);

G2 = findgroups(S.period, S.ActPass);
Per = table;
Per.period  = splitapply(@unique, S.period, G2);
Per.ActPass = splitapply(@unique, S.ActPass, G2);
Per.mu      = splitapply(@mean,  S.mu_mouse, G2);
Per.sem     = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), S.mu_mouse, G2);
Per.n_mice  = splitapply(@(x) numel(unique(x)), S.mouse_key, G2);
end

function fig = plotPeriodSummary(Per, ttl)
fig = figure('Color','w','Position',[80 80 980 420]);
hold on;
periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(Per.period)));
x = 1:numel(periods);
for ap = ["Active","Passive"]
    mu = nan(1,numel(periods)); se = mu;
    for i=1:numel(periods)
        rr = Per.ActPass==ap & Per.period==periods{i};
        if any(rr)
            mu(i)=Per.mu(rr);
            se(i)=Per.sem(rr);
        end
    end
    plot(x, mu, '-o', 'LineWidth',1.8);
    errorbar(x, mu, se, 'LineStyle','none');
end
set(gca,'XTick',x,'XTickLabel',periods); xtickangle(25);
ylabel('Licks / min'); title(['Period summary: ' ttl], 'Interpreter','none');
grid on; box on; legend({'Active','A sem','Passive','P sem'},'Location','best');
end

function [Cust_incl, Cust_excl] = customPeriodSummaryRequested(TsesA, lickPM_col, transitionDays)
Cust_incl = customPeriodCore(TsesA, lickPM_col, false, transitionDays);
Cust_excl = customPeriodCore(TsesA, lickPM_col, true,  transitionDays);
end

function Cust = customPeriodCore(T, lickPM_col, excludeTransitions, transitionDays)
T2 = T;
if excludeTransitions
    T2 = T2(~ismember(double(T2.day_index), transitionDays), :);
end

% define custom bins per ActPass
% Active: D5, D6-10, D11-13, D14-16, D17 (and D18 if exists -> merge with D17-18? user said day17 for active)
% Passive: D5, D11-13, D14-16, D17-18
bin = strings(height(T2),1); bin(:) = "<undef>";
d = double(T2.day_index);

isA = (T2.ActPass=="Active");
isP = (T2.ActPass=="Passive");

bin(isA & d==5)           = "A_D5";
bin(isA & d>=6 & d<=10)   = "A_D6-10";
bin(isA & d>=11 & d<=13)  = "A_D11-13";
bin(isA & d>=14 & d<=16)  = "A_D14-16";
bin(isA & d==17)          = "A_D17";
bin(isA & d==18)          = "A_D17-18";

bin(isP & d==5)           = "P_D5";
bin(isP & d>=11 & d<=13)  = "P_D11-13";
bin(isP & d>=14 & d<=16)  = "P_D14-16";
bin(isP & d>=17 & d<=18)  = "P_D17-18";

T2.custom_bin = categorical(bin);

keep = T2.custom_bin~="<undef>" & isfinite(double(T2.(lickPM_col)));
T2 = T2(keep,:);

% mouse-level average first
G = findgroups(T2.mouse_key, T2.custom_bin);
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.custom_bin= splitapply(@unique, T2.custom_bin, G);
S.mu_mouse  = splitapply(@mean,  T2.(lickPM_col), G);

G2 = findgroups(S.custom_bin);
Cust = table;
Cust.custom_bin = splitapply(@unique, S.custom_bin, G2);
Cust.mu         = splitapply(@mean,  S.mu_mouse, G2);
Cust.sem        = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), S.mu_mouse, G2);
Cust.n_mice     = splitapply(@(x) numel(unique(x)), S.mouse_key, G2);

Cust = sortrows(Cust, 'custom_bin');
end

function fig = plotCustomPeriodSummary(Cust, ttl)
fig = figure('Color','w','Position',[80 80 980 420]);
hold on;
x = 1:height(Cust);
plot(x, Cust.mu, '-o','LineWidth',1.8);
errorbar(x, Cust.mu, Cust.sem, 'LineStyle','none');
set(gca,'XTick',x,'XTickLabel',string(Cust.custom_bin)); xtickangle(25);
ylabel('Licks / min'); title(['Custom requested bins: ' ttl], 'Interpreter','none');
grid on; box on;
end

function fig = plotPairwiseByPeriod(T, lickPM_col, periodVar)
% mouse-level period mean then connect Passive->Active within PairID
P = T.(periodVar);
keep = P~="<undef>";
if any(strcmp("<exclude>", categories(P)))
    keep = keep & (P~="<exclude>");
end
T2 = T(keep & isfinite(double(T.(lickPM_col))) & isfinite(T.PairID), :);

G = findgroups(T2.mouse_key, T2.PairID, T2.ActPass, T2.(periodVar));
S = table;
S.mouse_key = splitapply(@unique, T2.mouse_key, G);
S.PairID    = splitapply(@unique, T2.PairID,    G);
S.ActPass   = splitapply(@unique, T2.ActPass,   G);
S.period    = splitapply(@unique, T2.(periodVar), G);
S.mu_mouse  = splitapply(@mean,  T2.(lickPM_col),  G);

periodOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
periods = periodOrder(ismember(periodOrder, categories(S.period)));
nP = numel(periods);

fig = figure('Color','w','Position',[80 80 320*max(1,nP) 380]);
tlo = tiledlayout(fig,1,max(1,nP),'TileSpacing','compact','Padding','compact');

for ip=1:nP
    ax = nexttile(tlo, ip); hold(ax,'on'); box(ax,'on');
    pr = periods{ip};
    Sp = S(S.period==pr,:);
    pairIDs = unique(Sp.PairID);

    for i=1:numel(pairIDs)
        pid = pairIDs(i);
        A = Sp(Sp.PairID==pid & Sp.ActPass=="Active", :);
        Pm= Sp(Sp.PairID==pid & Sp.ActPass=="Passive",:);
        if isempty(A) || isempty(Pm), continue; end
        aVal = A.mu_mouse(1);
        for j=1:height(Pm)
            pVal = Pm.mu_mouse(j);
            plot(ax, [1 2], [aVal pVal], '-', 'LineWidth',1.2);
        end
    end

    aPts = Sp.mu_mouse(Sp.ActPass=="Active");
    pPts = Sp.mu_mouse(Sp.ActPass=="Passive");
    scatter(ax, 1 + 0.06*(rand(numel(aPts),1)-0.5), aPts, 28, 'filled','MarkerFaceAlpha',0.6);
    scatter(ax, 2 + 0.06*(rand(numel(pPts),1)-0.5), pPts, 28, 'filled','MarkerFaceAlpha',0.6);

    xlim(ax,[0.6 2.4]);
    set(ax,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
    ylabel(ax,'Licks / min'); title(ax, pr, 'Interpreter','none');
    grid(ax,'on');
end
title(tlo, sprintf('Pairwise A vs P (by PairID) — %s', periodVar), 'Interpreter','none');
end

function Out = findOutliersRobust(T, lickPM_col, zthr)
% robust z within (day x ActPass)
keep = isfinite(double(T.(lickPM_col)));
T2 = T(keep,:);
val = double(T2.(lickPM_col));

G = findgroups(T2.day_index, T2.ActPass);
med = splitapply(@median, val, G);
mad = splitapply(@(x) median(abs(x-median(x))) + eps, val, G); % robust scale
rz  = (val - med(G)) ./ (1.4826*mad(G)); % robust z

rr = abs(rz) > zthr;
Out = T2(rr, {'mouse_key','day_index','session_idx','ActPass','PairID',lickPM_col});
Out.robust_z = rz(rr);
Out = sortrows(Out, {'day_index','ActPass','robust_z'}, {'ascend','ascend','descend'});
end

function fig = plotNormalizedProgressCurves(TtrialA, TsesA, keys, lickRate_col, NBINS)
% join trial with session cluster? Here we just plot Active vs Passive average curve (QC)
Ktrl = join(TtrialA, TsesA(:,[keys {'ActPass'}]), 'Keys', keys);

% compute normalized progress bins per session
Gsess = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx);
nT = splitapply(@max, Ktrl.Trial, Gsess);
prog = (Ktrl.Trial - 0.5) ./ nT(Gsess);

binEdges = linspace(0,1,NBINS+1);
[~,bin] = histc(prog, binEdges); %#ok<HISTC>
bin(bin<1)=1; bin(bin>NBINS)=NBINS;

% per-session per-bin
Gsb = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx, bin, Ktrl.ActPass);
Ssb = table;
Ssb.ActPass = splitapply(@unique, Ktrl.ActPass, Gsb);
Ssb.bin     = splitapply(@unique, bin, Gsb);
Ssb.meanLR  = splitapply(@mean,  Ktrl.(lickRate_col), Gsb);

% mean across sessions
Gcb = findgroups(Ssb.ActPass, Ssb.bin);
Csum = table;
Csum.ActPass = splitapply(@unique, Ssb.ActPass, Gcb);
Csum.bin     = splitapply(@unique, Ssb.bin,     Gcb);
Csum.mu      = splitapply(@mean,   Ssb.meanLR,  Gcb);
Csum.sem     = splitapply(@(x) nanstd(x,0)/sqrt(sum(isfinite(x))), Ssb.meanLR, Gcb);
xPct = (Csum.bin - 0.5) / NBINS * 100;

fig = figure('Color','w','Position',[80 80 900 420]);
hold on; box on;
acts = categories(Csum.ActPass);
for i=1:numel(acts)
    rr = Csum.ActPass==acts{i};
    [x,ord]=sort(xPct(rr));
    mu=Csum.mu(rr); mu=mu(ord);
    se=Csum.sem(rr); se=se(ord);
    plot(x,mu,'LineWidth',2);
    ii=1:max(1,round(numel(x)/30)):numel(x);
    errorbar(x(ii),mu(ii),se(ii),'LineStyle','none');
end
xlabel('Session progress (%)'); ylabel('Lick rate (Hz)');
title('QC: normalized session progress curves (Active vs Passive)');
grid on; legend(acts,'Location','best');
end

function savePNG(fig, pathOut)
set(fig,'PaperPositionMode','auto');
try
    exportgraphics(fig, pathOut, 'Resolution', 200);
catch
    print(fig, pathOut, '-dpng','-r200');
end
close(fig);
end

%% ================= PUPIL DAY SUMMARY (OPTIONAL) =================

function p = findAnyPupilDaySummary(runDir)
% tries common names; adjust if your pipeline uses different names
cands = { ...
    'pupil_day_summary.csv', ...
    'per_day_pupil.csv', ...
    'pupil_summary_by_day.csv', ...
    'pupil_features_per_day.csv' ...
};
p = '';
for i=1:numel(cands)
    pp = findOneFile(runDir, cands{i});
    if ~isempty(pp), p = pp; return; end
end
% broader search
dd = dir(fullfile(runDir,'**','*pupil*day*.csv'));
if isempty(dd), return; end
[~,ix]=max([dd.datenum]);
p = fullfile(dd(ix).folder, dd(ix).name);
end

function Pup = sanitizePupilDayTable(Pup)
% Expect at minimum: mouse_key, day_index, pupil_metric (mean/median)
% Try to find pupil column automatically and rename to pupil_mean
Pup.mouse_key = normalizeMouseKey(Pup.mouse_key);
if ~ismember('day_index', Pup.Properties.VariableNames)
    dcol = pickVar(Pup, {'day','Day','dayIdx','day_index'});
    assert(~isempty(dcol), 'Pupil day table missing day_index-like column.');
    Pup.day_index = Pup.(dcol);
end

pcol = pickVar(Pup, {'pupil_mean','pupil','mean_pupil','pupil_diam','pupil_diameter','pupil_area','mean_area'});
assert(~isempty(pcol), 'Could not find pupil value column in pupil day table.');
Pup.pupil_mean = double(Pup.(pcol));
end

function PupN = normalizePupilByBaselineDay(Pup, baselineDay)
% baselineDay (e.g., 3) per mouse
Gm = findgroups(Pup.mouse_key);
base = splitapply(@(d,v) nanmean(v(d==baselineDay)), double(Pup.day_index), double(Pup.pupil_mean), Gm);
PupN = Pup;
PupN.pupil_norm = PupN.pupil_mean ./ base(Gm);
PupN.pupil_delta = PupN.pupil_mean - base(Gm);
end

function fig = plotPupilDayCompare(PupN)
% requested: day4-5 vs day12-13, normalized by day3 baseline
isD45   = PupN.day_index==4 | PupN.day_index==5;
isD1213 = PupN.day_index==12 | PupN.day_index==13;

% mouse-level mean
G1 = findgroups(PupN.mouse_key, isD45);
mD45 = splitapply(@(x) nanmean(x), PupN.pupil_norm, G1);
G2 = findgroups(PupN.mouse_key, isD1213);
mD1213 = splitapply(@(x) nanmean(x), PupN.pupil_norm, G2);

% align by mouse keys
mice = unique(PupN.mouse_key,'stable');
v45 = nan(numel(mice),1); v1213 = v45;
for i=1:numel(mice)
    rr = PupN.mouse_key==mice(i);
    v45(i)   = nanmean(PupN.pupil_norm(rr & isD45));
    v1213(i) = nanmean(PupN.pupil_norm(rr & isD1213));
end

fig = figure('Color','w','Position',[80 80 700 420]);
hold on; box on;
for i=1:numel(mice)
    if isfinite(v45(i)) && isfinite(v1213(i))
        plot([1 2],[v45(i) v1213(i)],'-','LineWidth',1.2);
    end
end
scatter(1 + 0.06*(rand(sum(isfinite(v45)),1)-0.5), v45(isfinite(v45)), 30, 'filled','MarkerFaceAlpha',0.6);
scatter(2 + 0.06*(rand(sum(isfinite(v1213)),1)-0.5), v1213(isfinite(v1213)), 30, 'filled','MarkerFaceAlpha',0.6);

xlim([0.6 2.4]); set(gca,'XTick',[1 2],'XTickLabel',{'D4-5 (water PR)','D12-13 (morphine PR)'});
ylabel('Pupil normalized ( / baseline D3 )');
title('Pupil day-level: D4-5 vs D12-13 (normalized by D3 baseline)');
grid on;
end

%% ================= EVENT-ALIGNED PUPIL (OPTIONAL) =================

function p = findAnyEventAlignedPupil(runDir)
% try common names first
cands = { ...
    'pupil_event_aligned.csv', ...
    'event_aligned_pupil.csv', ...
    'pupil_reward_aligned.csv', ...
    'pupil_lick_aligned.csv' ...
};
p = '';
for i=1:numel(cands)
    pp = findOneFile(runDir, cands{i});
    if ~isempty(pp), p = pp; return; end
end
dd = dir(fullfile(runDir,'**','*pupil*event*.csv'));
if isempty(dd), return; end
[~,ix]=max([dd.datenum]);
p = fullfile(dd(ix).folder, dd(ix).name);
end

function [BoutTbl, fig1, fig2] = eventAlignedPupilBoutQC(E, cohort, preSec, postSec, boutGap, transitionDays)
% This is a "best-effort" template. It expects that E has enough info to:
%   - identify mouse_key/day
%   - identify event times and pupil time series samples around events
% If E is already in "event-aligned matrix" form (one row per event, columns are t_-2...t_+2),
% we can summarize directly.
%
% Supported formats:
%   A) Wide: columns like t_-2.00 ... t_+2.00 or p_-2...p_+2
%   B) Long: columns [mouse_key, day_index, event_id, t_rel, pupil] and optionally event_type/reward_flag
%
% Here we handle (A) robustly; (B) needs your exact columns -> you can extend.

E.mouse_key = normalizeMouseKey(E.mouse_key);
if ~ismember('day_index', E.Properties.VariableNames)
    dcol = pickVar(E, {'day','Day','dayIdx','day_index'});
    assert(~isempty(dcol), 'Event pupil table missing day_index-like column.');
    E.day_index = E.(dcol);
end
E = attachCohortAndPeriods(E, cohort, transitionDays);

% detect wide pupil columns
cols = string(E.Properties.VariableNames);
pcols = cols(startsWith(cols,'p_') | startsWith(cols,'pupil_') | startsWith(cols,'P_'));
tcols = cols(startsWith(cols,'t_') | startsWith(cols,'T_'));

isWide = ~isempty(pcols) || (~isempty(tcols) && any(startsWith(cols,'mean_')));

if ~isWide
    error('Event-aligned pupil file is not in a recognized "wide" format. Please share column names.');
end

% If columns are like p_-2.00 etc:
if isempty(pcols)
    % fallback: columns starting with 'pupil_' assumed to be pupil
    pcols = cols(contains(cols,'pupil'));
end
% sort by suffix numeric if possible
pcols = sortWideByNumericSuffix(pcols);

Pmat = E{:, pcols};
Pmat = double(Pmat);

% baseline = mean of pre window (assume first half is pre)
n = size(Pmat,2);
preN = max(1, floor(n * (preSec/(preSec+postSec+eps))));
base = mean(Pmat(:,1:preN),2,'omitnan');
Pnorm = Pmat - base;

% reward vs nonreward split if reward flag exists
rcol = pickVar(E, {'reward','reward_flag','is_reward','Reward','rewarded'});
if ~isempty(rcol)
    rew = logical(E.(rcol));
else
    rew = false(height(E),1);
end

% period focus: passive D6-10 vs post D12-13 etc (your focus)
isPassiveD6_10 = (E.ActPass=="Passive") & (E.day_index>=6 & E.day_index<=10);
isPostD12_13   = (E.day_index>=12 & E.day_index<=13);

% summarize: peak delta and AUC
peak = max(Pnorm,[],2,'omitnan');
auc  = mean(Pnorm,2,'omitnan'); % average delta over window

BoutTbl = table(E.mouse_key, E.day_index, E.ActPass, E.PairID, E.period, rew, peak, auc, ...
    'VariableNames', {'mouse_key','day_index','ActPass','PairID','period','is_reward','pupil_peak_delta','pupil_mean_delta'});

% fig1: reward vs nonreward distributions
fig1 = figure('Color','w','Position',[80 80 900 420]);
hold on; box on;
grp = categorical(rew, [false true], {'nonreward','reward'});
for i=1:2
    rr = (grp==categories(grp){i});
    scatter(i + 0.25*(rand(sum(rr),1)-0.5), peak(rr), 20, 'filled', 'MarkerFaceAlpha',0.35);
    plot([i-0.2 i+0.2], [median(peak(rr),'omitnan') median(peak(rr),'omitnan')], 'k-', 'LineWidth',2);
end
set(gca,'XTick',[1 2],'XTickLabel',categories(grp));
ylabel('Pupil peak delta (baseline=pre window)');
title('Event-aligned pupil: reward vs nonreward (peak delta)');
grid on;

% fig2: passive D6-10 vs post D12-13 (peak)
fig2 = figure('Color','w','Position',[80 80 900 420]);
hold on; box on;
v1 = peak(isPassiveD6_10);
v2 = peak(isPostD12_13);
scatter(1 + 0.25*(rand(numel(v1),1)-0.5), v1, 20, 'filled', 'MarkerFaceAlpha',0.35);
scatter(2 + 0.25*(rand(numel(v2),1)-0.5), v2, 20, 'filled', 'MarkerFaceAlpha',0.35);
plot([1-0.2 1+0.2],[median(v1,'omitnan') median(v1,'omitnan')],'k-','LineWidth',2);
plot([2-0.2 2+0.2],[median(v2,'omitnan') median(v2,'omitnan')],'k-','LineWidth',2);
set(gca,'XTick',[1 2],'XTickLabel',{'Passive D6-10 (forced)','D12-13 (post morphine PR)'});
ylabel('Pupil peak delta');
title('Event-aligned pupil: focus comparison');
grid on;
end

function colsSorted = sortWideByNumericSuffix(cols)
% cols like "p_-2.00" or "p_0.10"
cols = string(cols);
num = nan(size(cols));
for i=1:numel(cols)
    tok = regexp(cols(i), '[-+]?\d+\.?\d*', 'match','once');
    if ~isempty(tok), num(i)=str2double(tok); end
end
[~,ord] = sort(num);
colsSorted = cols(ord);
end
