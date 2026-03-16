function analyze_passive_active_dashboardnew()
% Fully revised, fully functional dashboard for your NEW cohort + NEW timeline.
%
% Key updates:
% 1) Hard-coded cohort mapping (cage/color/sex/group + pair IDs), no guessing.
% 2) Epochs: Pre(3–5), During(6–10), Post(11–13), Withdrawal(14–16), Reexposure(17–18).
% 3) Runs analyses twice:
%       - include ALL days
%       - EXCLUDE “transition/suspect” days (customizable)
% 4) Adds lick bout metrics + cumulative lick curves across session time (minute-binned),
%    per mouse (day-by-day), and group aggregate (Passive vs Active).
%
% Outputs:
%   <runDir>\figs\dashboard\INCL_ALL_DAYS\
%   <runDir>\figs\dashboard\EXCL_TRANSITION_DAYS\

%% ---------------- user settings ----------------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% Habituation days (ignored in epoch analyses & cumulative plots)
habituationDays = [1 2];

% “Transition / suspect” days to optionally exclude
% (You mentioned day4, day6, day11, day14 are less reliable; include day3 & day17 too by default)
transitionDays = [3 4 6 11 14 17];

% Lick bout definition: new bout starts when gap between lick onsets > boutGap_s
boutGap_s = 1.0;   % adjust if you use a different bout criterion

% Cumulative lick curves: minute bin width
cumBin_min = 1;    % 1-minute bins

% Max minutes to display in cumulative curves (leave empty to auto)
cumMaxMin = [];    % e.g., 60

%% ---------------- locate latest run + load ----------------
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder,d(ix).name);
csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);

fprintf('Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% Normalize core types
T = ensureString(T,'mouse_key');
T = ensureString(T,'day_name');
T = ensureString(T,'Session_Paradigm');

% TTL numerics handling (keep raw; threshold later)
if ismember('Lick_TTL',T.Properties.VariableNames),     T.Lick_TTL(isnan(T.Lick_TTL)) = 0; end
if ismember('Injector_TTL',T.Properties.VariableNames), T.Injector_TTL(isnan(T.Injector_TTL)) = 0; end

% Ensure day_index exists
assert(ismember('day_index',T.Properties.VariableNames),'CSV must contain day_index.');

%% ---------------- cohort mapping (NO guessing) ----------------
Cohort = buildCohortTable();  % cage/color/sex/group/pair/alive_last_day
[T, mapReport] = attachCohortInfo(T, Cohort);

% Hard check: if anything unmatched, stop (you said: do not auto-detect/guess)
if any(mapReport.unmatchedMouseKeysCount > 0)
    disp(mapReport);
    error('Unmatched mouse_key(s) exist in CSV. Fix parsing or update cohort table.');
end

% Drop any rows beyond alive_last_day (e.g., 6099 white after day13 if present)
dropDead = (T.day_index > T.alive_last_day);
if any(dropDead)
    fprintf('Dropping %d rows beyond alive_last_day.\n', nnz(dropDead));
    T(dropDead,:) = [];
end

%% ---------------- compute per-session + per-day metrics (with caching) ----------------
[S, D] = fast_session_day_metrics_v2(T, runDir, boutGap_s);

% Keep only analysis days (3–18)
keepDays = (D.day_index >= 3 & D.day_index <= 18);
D = D(keepDays,:);

% Group is fixed from cohort map
D.Group = categorical(string(D.group_fixed), {'Active','Passive'});

% Epoch label
D.Epoch = epochOfDay5(double(D.day_index));

%% ---------------- output dirs ----------------
outBase = fullfile(runDir,'figs','dashboard');
if ~exist(outBase,'dir'), mkdir(outBase); end

%% ---------------- run analyses in 2 modes ----------------
modes = { ...
    struct('name','INCL_ALL_DAYS',          'excludeDays', []), ...
    struct('name','EXCL_TRANSITION_DAYS',   'excludeDays', unique([habituationDays transitionDays])) ...
    };

for mi = 1:numel(modes)
    mode = modes{mi};
    outDir = fullfile(outBase, mode.name);
    if ~exist(outDir,'dir'), mkdir(outDir); end

    fprintf('\n===== MODE: %s =====\n', mode.name);

    % Filter D for this mode
    Dm = D;
    if ~isempty(mode.excludeDays)
        Dm = Dm(~ismember(Dm.day_index, mode.excludeDays), :);
    end

    % Recompute epoch (in case all days of some epoch removed)
    Dm.Epoch = epochOfDay5(double(Dm.day_index));

    % ---------------- metrics list ----------------
    baseMetrics = { ...
        'RequirementLast', ...
        'lick_freq_per_min','lick_n', ...
        'lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
        'lick_bout_freq_per_min','lick_bout_n','lick_bout_meanDur_s','lick_bout_meanSize','lick_bout_medianIBI_s', ...
        'rew_freq_per_min','rew_n', ...
        'rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s', ...
        'pupil_mean'};

    baseLabels  = { ...
        'Requirement', ...
        'Lick freq (/min)','Lick count', ...
        'Lick mean dur (s)','Lick total dur (s)','Lick median IEI (s)', ...
        'Lick bout freq (/min)','Lick bout count','Lick bout mean dur (s)','Licks per bout (mean)','Bout median IBI (s)', ...
        'Reward freq (/min)','Reward count', ...
        'Reward mean dur (s)','Reward total dur (s)','Reward median IRI (s)', ...
        'Pupil mean (px)'};

    % Optional extras from your CSV (Immersion/TST/HOT)
    extraM = {}; extraL = {};
    immCol = colLike(Dm,'Immersion_Latency_s');
    if ~isempty(immCol), extraM{end+1}=immCol; extraL{end+1}='Immersion latency (s)'; end

    tstCols = Dm.Properties.VariableNames(contains(lower(Dm.Properties.VariableNames),'tst_') & contains(lower(Dm.Properties.VariableNames),'pct'));
    hotCols = Dm.Properties.VariableNames(contains(lower(Dm.Properties.VariableNames),'hot_') & contains(lower(Dm.Properties.VariableNames),'pct'));
    tstPick = pickPctColumn(tstCols, 'nonmoving'); if ~isempty(tstPick), extraM{end+1}=tstPick; extraL{end+1}='TST nonmoving (%)'; end
    hotPick = pickPctColumn(hotCols, 'nonmoving'); if ~isempty(hotPick), extraM{end+1}=hotPick; extraL{end+1}='HOT nonmoving (%)'; end

    metrics = [baseMetrics, extraM];
    labels  = [baseLabels,  extraL];

    keepIdx = ismember(metrics, Dm.Properties.VariableNames);
    metrics = metrics(keepIdx);
    labels  = labels(keepIdx);

    % ---------------- main plots + stats ----------------
    allStats = table();

    for k = 1:numel(metrics)
        ycol = metrics{k};
        ylab = labels{k};

        W = perMouseEpochTable_v2(Dm, ycol);   % one row per mouse×epoch (median across days)
        if isempty(W), continue; end

        % Dual-panel epoch lines (Passive panel / Active panel)
        fh = figure('Color','w','Position',[80 60 980 880]);
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
        dualPanelEpochLines_v2(W, ylab);
        printpng(fh, fullfile(outDir, sprintf('dual_epoch_%s.png', safeName(ycol)))); close(fh);

        % Change-from-Pre per epoch (strip+box) with p-values (Passive vs Active per epoch)
        fh = figure('Color','w','Position',[100 80 980 520]);
        changeFromPreStrip_v2(W, ylab, 'none');
        printpng(fh, fullfile(outDir, sprintf('delta_from_pre_%s.png', safeName(ycol)))); close(fh);

        % Within-group repeated-measures (Friedman across epochs) + N annotation
        fh = figure('Color','w','Position',[100 90 980 420]);
        rmFriedmanWithinGroupFigure_v2(W, ylab);
        printpng(fh, fullfile(outDir, sprintf('rm_friedman_%s.png', safeName(ycol)))); close(fh);

        % Stats tables
        ST  = runStats_v2(W, ycol);         % within-group signrank (Pre vs each epoch), between-group ranksum per epoch
        STd = runDeltaStats_v2(W, ycol);    % between-group ranksum on deltas per epoch
        STrm = runFriedmanStats_v2(W, ycol);

        ST.metric  = repmat(string(ycol), height(ST),1);
        ST.label   = repmat(string(ylab), height(ST),1);
        STd.metric = repmat(string(ycol), height(STd),1);
        STd.label  = repmat(string(ylab), height(STd),1);
        STrm.metric= repmat(string(ycol), height(STrm),1);
        STrm.label = repmat(string(ylab), height(STrm),1);

        allStats = [allStats; ST; STd; STrm]; %#ok<AGROW>
    end

    % Global multiple-comparison corrections across ALL tests for ALL metrics
    if ~isempty(allStats) && ismember('p', allStats.Properties.VariableNames)
        allStats.p_FDR  = fdrBH(allStats.p);
        allStats.p_Holm = holmBonferroni(allStats.p);
        writetable(allStats, fullfile(outDir,'stats_summary.csv'));
        disp(allStats(:,{'metric','test','effect','level','p','p_FDR','p_Holm','N'}));
    end

    % ---------------- cumulative lick curves ----------------
    % Uses RAW T (not D), but respects this mode's day exclusions
    outCum = fullfile(outDir,'cumulative_lick');
    if ~exist(outCum,'dir'), mkdir(outCum); end

    dayExcl = unique([habituationDays mode.excludeDays]);
    buildCumulativeLickPlots(T, outCum, dayExcl, cumBin_min, cumMaxMin);

    % ---------------- focused comparison: Day11–13 only (Post epoch) ----------------
    % (Pupil, immersion, TST, HOT) are already in Post epoch; also export a focused CSV
    postDays = [11 12 13];
    if ~isempty(mode.excludeDays), postDays = setdiff(postDays, mode.excludeDays); end
    Dpost = D(D.day_index>=11 & D.day_index<=13, :);
    if ~isempty(mode.excludeDays)
        Dpost = Dpost(~ismember(Dpost.day_index, mode.excludeDays), :);
    end
    if ~isempty(Dpost)
        colsFocus = {'mouse_key','day_index','Group','pupil_mean'};
        if ismember('Immersion_Latency_s', Dpost.Properties.VariableNames), colsFocus{end+1}='Immersion_Latency_s'; end
        if ~isempty(tstPick), colsFocus{end+1}=tstPick; end
        if ~isempty(hotPick), colsFocus{end+1}=hotPick; end
        colsFocus = colsFocus(ismember(colsFocus, Dpost.Properties.VariableNames));
        writetable(Dpost(:,colsFocus), fullfile(outDir,'post_day11_13_focus.csv'));
    end

    fprintf('MODE done. Outputs:\n  %s\n', outDir);
end

fprintf('\nALL DONE. Base output:\n  %s\n', outBase);
end

%% ===================== cohort definition =====================
function Cohort = buildCohortTable()
% cage (4-digit), color (subfolder), sex, group_fixed, pair_id, alive_last_day
rows = { ...
% 6100: orange+red passive vs black active
'6100','red',    'F','Passive',1, inf; ...
'6100','orange', 'F','Passive',1, inf; ...
'6100','black',  'F','Active', 1, inf; ...
% 0911 pair 1: red active vs orange passive
'0911','red',    'F','Active', 2, inf; ...
'0911','orange', 'F','Passive',2, inf; ...
% 0911 pair 2: white active vs black passive
'0911','white',  'F','Active', 3, inf; ...
'0911','black',  'F','Passive',3, inf; ...
% 0910: black active vs orange+red passive
'0910','black',  'M','Active', 4, inf; ...
'0910','orange', 'M','Passive',4, inf; ...
'0910','red',    'M','Passive',4, inf; ...
% 6099 pair 1: orange active vs red passive
'6099','orange', 'M','Active', 5, inf; ...
'6099','red',    'M','Passive',5, inf; ...
% 6099 pair 2: black active vs white passive (white died day1–13 only)
'6099','black',  'M','Active', 6, inf; ...
'6099','white',  'M','Passive',6, 13; ...
};
Cohort = cell2table(rows, 'VariableNames', ...
    {'cage','color','sex','group_fixed','pair_id','alive_last_day'});
Cohort.cage  = string(Cohort.cage);
Cohort.color = lower(string(Cohort.color));
Cohort.sex   = upper(string(Cohort.sex));
Cohort.group_fixed = string(Cohort.group_fixed);
end

function [T, report] = attachCohortInfo(T, Cohort)
% Parse mouse_key into cage/color and map to Cohort. No guessing beyond parsing.

mk = lower(strtrim(string(T.mouse_key)));
mk = regexprep(mk,'[_\-]+',' ');

% extract 4 digits cage
cage = regexp(mk,'(?<!\d)(\d{4})(?!\d)','tokens','once');
cage = string(cellfun(@(x) iff(isempty(x), "", x{1}), cage, 'UniformOutput', false));

% extract color token (black/red/orange/white)
color = strings(size(mk));
colors = ["black","red","orange","white"];
for i=1:numel(colors)
    hit = contains(mk, " "+colors(i)+" ") | startsWith(mk, colors(i)+" ") | endsWith(mk, " "+colors(i)) | contains(mk, colors(i));
    color(hit & color=="") = colors(i);
end

T.cage  = cage;
T.color = color;

% Map
keyT = lower(T.cage) + "|" + lower(T.color);
keyC = lower(Cohort.cage) + "|" + lower(Cohort.color);

[tf, loc] = ismember(keyT, keyC);

% Report by mouse_key
uKeys = unique(string(T.mouse_key),'stable');
unmatched = false(size(uKeys));
for i=1:numel(uKeys)
    mask = string(T.mouse_key)==uKeys(i);
    unmatched(i) = any(~tf(mask));
end

report = table();
report.unmatchedMouseKeysCount = sum(unmatched);
if report.unmatchedMouseKeysCount > 0
    report.unmatchedMouseKeys = {uKeys(unmatched)};
end

% Attach cohort columns
T.sex         = strings(height(T),1);
T.group_fixed = strings(height(T),1);
T.pair_id     = nan(height(T),1);
T.alive_last_day = inf(height(T),1);

T.sex(tf)         = Cohort.sex(loc(tf));
T.group_fixed(tf) = Cohort.group_fixed(loc(tf));
T.pair_id(tf)     = Cohort.pair_id(loc(tf));
T.alive_last_day(tf)= Cohort.alive_last_day(loc(tf));

% If anything not mapped, leave blank and force failure in caller
end

function y = iff(cond, a, b)
if cond, y = a; else, y = b; end
end

%% ===================== epoch mapping =====================
function E = epochOfDay5(d)
% Habituation: day1–2 (ignored in analysis)
% Pre: day3–5
% During: day6–10
% Post: day11–13
% Withdrawal: day14–16
% Reexposure: day17–18
e = strings(size(d));
e(d>=3  & d<=5 )  = "Pre";
e(d>=6  & d<=10)  = "During";
e(d>=11 & d<=13)  = "Post";
e(d>=14 & d<=16)  = "Withdrawal";
e(d>=17 & d<=18)  = "Reexposure";
E = categorical(e, ["Pre","During","Post","Withdrawal","Reexposure"], 'Ordinal', true);
end

%% ===================== per-session + per-day metrics (cached) =====================
function [S, D] = fast_session_day_metrics_v2(T, runDir, boutGap_s)
cacheMat = fullfile(runDir, 'S_D_cache_v2.mat');
if exist(cacheMat,'file')
    L = load(cacheMat,'S','D','cache_hash','boutGap_s_cached');
    if isfield(L,'cache_hash') && isequal(L.cache_hash, local_hash(T)) && ...
       isfield(L,'boutGap_s_cached') && isequal(L.boutGap_s_cached, boutGap_s)
        fprintf('Loaded S/D from cache: %s\n', cacheMat);
        S = L.S; D = L.D; return;
    end
end

tic;

need = {'mouse_key','day_index','day_name','session_idx', ...
        'Diameter_px','Lick_TTL','Injector_TTL', ...
        'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps', ...
        'Session_Paradigm','RequirementLast', ...
        'group_fixed','pair_id','sex','alive_last_day','cage','color'};

V = T.Properties.VariableNames;
keepExtras = V( contains(V,'Immersion_Latency','IgnoreCase',true) | ...
                contains(V,'TST_','IgnoreCase',true)            | ...
                contains(V,'HOT_','IgnoreCase',true) );

need = unique([need, keepExtras], 'stable');
need = intersect(need, V, 'stable');
T = T(:, need);

% Types
T.mouse_key = categorical(string(T.mouse_key));
T.day_name  = categorical(string(T.day_name));
if ismember('Session_Paradigm', T.Properties.VariableNames)
    T.Session_Paradigm = categorical(string(T.Session_Paradigm));
end

% TTL threshold
if ismember('Lick_TTL', T.Properties.VariableNames)
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;
end
if ismember('Injector_TTL', T.Properties.VariableNames)
    T.Injector_TTL(isnan(T.Injector_TTL)) = 0; T.Injector_TTL = T.Injector_TTL > 0.5;
end

% Group by session
[g, keys_mouse, keys_day, keys_dayname, keys_sess] = findgroups( ...
    T.mouse_key, T.day_index, T.day_name, T.session_idx);

nG = max(g);
fprintf('Computing per-session metrics for %d sessions...\n', nG);

S = table();
S.mouse_key   = removecats(keys_mouse);
S.day_index   = double(keys_day);
S.day_name    = removecats(keys_dayname);
S.session_idx = double(keys_sess);

% fixed cohort columns (same within session)
S.group_fixed = strings(nG,1);
S.pair_id     = nan(nG,1);
S.sex         = strings(nG,1);
S.cage        = strings(nG,1);
S.color       = strings(nG,1);
S.alive_last_day = inf(nG,1);

if ismember('group_fixed',T.Properties.VariableNames), S.group_fixed = splitapply(@modeString, string(T.group_fixed), g); end
if ismember('pair_id',T.Properties.VariableNames),     S.pair_id     = splitapply(@(x) mode(x(~isnan(x))), T.pair_id, g); end
if ismember('sex',T.Properties.VariableNames),         S.sex         = splitapply(@modeString, string(T.sex), g); end
if ismember('cage',T.Properties.VariableNames),        S.cage        = splitapply(@modeString, string(T.cage), g); end
if ismember('color',T.Properties.VariableNames),       S.color       = splitapply(@modeString, string(T.color), g); end
if ismember('alive_last_day',T.Properties.VariableNames), S.alive_last_day = splitapply(@(x) mode(x(~isnan(x))), T.alive_last_day, g); end

S.RequirementLast = nan(nG,1);
S.SessionMinutes  = nan(nG,1);
S.Session_Paradigm= strings(nG,1);

vars = { ...
'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
'lick_bout_n','lick_bout_freq_per_min','lick_bout_meanDur_s','lick_bout_meanSize','lick_bout_medianIBI_s', ...
'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s', ...
'pupil_mean'};

for v = vars, S.(v{1}) = nan(nG,1); end
for j=1:numel(keepExtras), S.(keepExtras{j}) = nan(nG,1); end

tb_all = pickTimebase_fast(T);

lastTick = tic; step = max(1, floor(nG/50));
for k = 1:nG
    idx = (g == k);
    tb  = tb_all(idx);
    dur = finiteRange_fast(tb);
    S.SessionMinutes(k) = dur/60;

    if ismember('RequirementLast',T.Properties.VariableNames)
        S.RequirementLast(k) = mean(double(T.RequirementLast(idx)),'omitnan');
    end
    if ismember('Session_Paradigm', T.Properties.VariableNames)
        sp = T.Session_Paradigm(idx);
        sp = string(sp); sp = sp(~ismissing(sp));
        if isempty(sp), sp=""; else, sp=sp(1); end
        S.Session_Paradigm(k) = sp;
    end

    if ismember('Diameter_px',T.Properties.VariableNames)
        S.pupil_mean(k) = mean(double(T.Diameter_px(idx)),'omitnan');
    end

    if ismember('Lick_TTL',T.Properties.VariableNames)
        ttl = logical(T.Lick_TTL(idx));
        [n,md,td,iei,onTimes] = eventMetrics_fast_v2(tb, ttl);
        S.lick_n(k)=n; S.lick_meanDur_s(k)=md; S.lick_totalDur_s(k)=td; S.lick_medianIEI_s(k)=iei;
        if S.SessionMinutes(k)>0, S.lick_freq_per_min(k)=n/S.SessionMinutes(k); end

        % Bout metrics from lick onsets
        [bn, bmd, bsz, bibi] = boutMetrics_fromOnsets(onTimes, boutGap_s);
        S.lick_bout_n(k) = bn;
        S.lick_bout_meanDur_s(k) = bmd;
        S.lick_bout_meanSize(k)  = bsz;
        S.lick_bout_medianIBI_s(k)= bibi;
        if S.SessionMinutes(k)>0, S.lick_bout_freq_per_min(k)=bn/S.SessionMinutes(k); end
    end

    if ismember('Injector_TTL',T.Properties.VariableNames)
        [n,md,td,iri] = eventMetrics_fast_v2(tb, logical(T.Injector_TTL(idx)));
        S.rew_n(k)=n; S.rew_meanDur_s(k)=md; S.rew_totalDur_s(k)=td; S.rew_medianIRI_s(k)=iri;
        if S.SessionMinutes(k)>0, S.rew_freq_per_min(k)=n/S.SessionMinutes(k); end
    end

    for j=1:numel(keepExtras)
        col = keepExtras{j};
        if ismember(col, T.Properties.VariableNames)
            S.(col)(k) = mean(double(T.(col)(idx)),'omitnan');
        end
    end

    if mod(k, step)==0 && toc(lastTick)>0.5
        fprintf('  %d/%d (%.0f%%)\n', k, nG, 100*k/nG); lastTick = tic;
    end
end
fprintf('Per-session metrics done in %.1f s.\n', toc);

% Per-day collapse (median across sessions)
fprintf('Collapsing to per-day medians...\n');
[g2, mk2, di2, dn2, gf2] = findgroups(S.mouse_key, S.day_index, S.day_name, string(S.group_fixed));
D = table(removecats(mk2), double(di2), removecats(dn2), categorical(gf2), ...
          'VariableNames',{'mouse_key','day_index','day_name','group_fixed'});

% Carry cohort meta (mode over sessions)
D.pair_id = splitapply(@(x) mode(x(~isnan(x))), S.pair_id, g2);
D.sex     = splitapply(@modeString, string(S.sex), g2);
D.cage    = splitapply(@modeString, string(S.cage), g2);
D.color   = splitapply(@modeString, string(S.color), g2);
D.alive_last_day = splitapply(@(x) mode(x(~isnan(x))), S.alive_last_day, g2);

list = unique([{'RequirementLast','SessionMinutes','Session_Paradigm'}, vars, keepExtras'], 'stable');
list = intersect(list, S.Properties.VariableNames, 'stable');

for v = list
    if strcmp(v{1},'Session_Paradigm')
        D.Session_Paradigm = splitapply(@modeString, string(S.Session_Paradigm), g2);
    else
        D.(v{1}) = splitapply(@(x) median(x,'omitnan'), S.(v{1}), g2);
    end
end

cache_hash = local_hash(T);
boutGap_s_cached = boutGap_s;
try
    save(cacheMat,'S','D','cache_hash','boutGap_s_cached','-v7.3');
    fprintf('Cached S/D to: %s\n', cacheMat);
catch
end
end

function s = modeString(x)
x = string(x); x = x(x~="" & ~ismissing(x));
if isempty(x), s=""; return; end
[u,~,ic] = unique(x);
counts = accumarray(ic,1);
[~,ix] = max(counts);
s = u(ix);
end

%% ===================== cumulative lick plots =====================
function buildCumulativeLickPlots(T, outDir, excludeDays, binMin, maxMin)
% Build:
% 1) Per mouse: day-by-day cumulative lick curve (one line per day)
% 2) Group aggregate (Passive vs Active): pooled across all included days
% 3) Group aggregate by day: mean curve per day per group (optional)

assert(ismember('Lick_TTL',T.Properties.VariableNames),'No Lick_TTL in table.');
assert(ismember('day_index',T.Properties.VariableNames),'No day_index in table.');
assert(ismember('group_fixed',T.Properties.VariableNames),'No group_fixed (cohort mapping missing).');

T = T;
T.mouse_key = categorical(string(T.mouse_key));
T.group_fixed = categorical(string(T.group_fixed), {'Active','Passive'});

% Keep analysis days 3–18 and exclude requested
keep = (T.day_index>=3 & T.day_index<=18) & ~ismember(T.day_index, excludeDays);
T = T(keep,:);

if isempty(T)
    fprintf('No data left for cumulative plots after exclusions.\n');
    return;
end

% Timebase
tb = pickTimebase_fast(T);
T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
T.Lick_TTL = T.Lick_TTL > 0.5;

% Group by session
if ~ismember('session_idx',T.Properties.VariableNames)
    error('CSV must contain session_idx for cumulative curves.');
end

[g, mk, di, gf, si] = findgroups(T.mouse_key, T.day_index, T.group_fixed, T.session_idx);

nG = max(g);
C = table();
C.mouse_key = removecats(mk);
C.day_index = double(di);
C.Group     = gf;
C.session_idx = double(si);
C.t_min = cell(nG,1);
C.cum_lick = cell(nG,1);

% build per session curve (minute binned)
for k = 1:nG
    idx = (g==k);
    t = tb(idx);
    ttl = T.Lick_TTL(idx);

    [~,~,~,~,onTimes] = eventMetrics_fast_v2(t, ttl);
    if isempty(onTimes)
        C.t_min{k} = (0:binMin:binMin)';  % minimal
        C.cum_lick{k} = 0;
        continue;
    end

    t0 = min(onTimes);
    xMin = (onTimes - t0)/60;

    % decide max
    if isempty(maxMin)
        m = ceil(max(xMin));
    else
        m = maxMin;
    end
    edges = 0:binMin:(m+binMin);
    counts = histcounts(xMin, edges);
    cumc = cumsum(counts);

    C.t_min{k} = edges(1:end-1)';
    C.cum_lick{k} = cumc(:);
end

% --------- (1) Per mouse: day-by-day lines ---------
outMouse = fullfile(outDir,'per_mouse');
if ~exist(outMouse,'dir'), mkdir(outMouse); end

mice = unique(C.mouse_key,'stable');
for i=1:numel(mice)
    Cm = C(C.mouse_key==mice(i),:);
    if isempty(Cm), continue; end

    % Aggregate within day: if multiple sessions, sum counts by aligning on t grid
    days = unique(Cm.day_index,'sorted');

    fh = figure('Color','w','Position',[80 80 980 560]); hold on;
    cmap = lines(numel(days));
    for d = 1:numel(days)
        Cd = Cm(Cm.day_index==days(d),:);
        [tGrid, ySum] = sumCurves(Cd.t_min, Cd.cum_lick);
        plot(tGrid, ySum, '-', 'LineWidth',1.8, 'Color',cmap(d,:));
    end
    xlabel('Minutes from session start');
    ylabel('Cumulative licks');
    title(sprintf('Cumulative licks by day — %s', string(mice(i))));
    grid on; box on;
    lg = legend("Day "+string(days),'Location','eastoutside'); %#ok<NASGU>
    printpng(fh, fullfile(outMouse, sprintf('cum_lick_%s.png', safeName(string(mice(i))))));
    close(fh);
end

% --------- (2) Group aggregate pooled across all included days ---------
fh = figure('Color','w','Position',[90 90 980 520]); hold on;

for G = ["Passive","Active"]
    Cg = C(string(C.Group)==G,:);
    if isempty(Cg), continue; end
    [tGrid, yMat] = stackCurves(Cg.t_min, Cg.cum_lick);
    yMean = mean(yMat,2,'omitnan');
    ySem  = std(yMat,0,2,'omitnan') ./ sqrt(sum(isfinite(yMat),2));

    if G=="Passive"
        plot(tGrid, yMean, '-', 'LineWidth',2.4);
    else
        plot(tGrid, yMean, '-', 'LineWidth',2.4);
    end
    fill([tGrid; flipud(tGrid)], [yMean-ySem; flipud(yMean+ySem)], [0 0 0], ...
        'FaceAlpha',0.10,'EdgeColor','none');
end
xlabel('Minutes from session start');
ylabel('Cumulative licks');
title('Cumulative licks (pooled across days, mean ± SEM)');
grid on; box on;
legend({'Passive','Passive SEM','Active','Active SEM'},'Location','eastoutside');
printpng(fh, fullfile(outDir, 'cum_lick_group_pooled.png'));
close(fh);

end

function [tGrid, ySum] = sumCurves(tCells, yCells)
% sum multiple cumulative curves onto a common grid (pad with last value)
[tGrid, yMat] = stackCurves(tCells, yCells);
ySum = sum(yMat,2,'omitnan');
end

function [tGrid, yMat] = stackCurves(tCells, yCells)
% Build common grid and pad each curve with last value to the grid end
n = numel(tCells);
maxT = 0;
for i=1:n
    if isempty(tCells{i}), continue; end
    maxT = max(maxT, max(tCells{i}));
end
if maxT==0
    tGrid = (0:1:1)'; yMat = zeros(numel(tGrid), n); return;
end
% assume all have same bin step as first non-empty
dt = [];
for i=1:n
    if numel(tCells{i})>=2
        dt = median(diff(tCells{i})); break;
    end
end
if isempty(dt) || ~isfinite(dt) || dt<=0, dt=1; end
tGrid = (0:dt:maxT)';

yMat = nan(numel(tGrid), n);
for i=1:n
    t = tCells{i}; y = yCells{i};
    if isempty(t) || isempty(y), continue; end
    % align by interpolation on the grid using step-wise (previous) hold
    yPad = padToGrid_hold(t, y, tGrid);
    yMat(:,i) = yPad;
end
end

function yG = padToGrid_hold(t, y, tGrid)
t = double(t(:)); y = double(y(:));
[tu, ia] = unique(t,'stable');
yu = y(ia);
yG = nan(numel(tGrid),1);
% for each grid point, take last observed cumulative value at or before that time
j = 1;
last = yu(1);
for i=1:numel(tGrid)
    while j <= numel(tu) && tu(j) <= tGrid(i)
        last = yu(j);
        j = j + 1;
    end
    yG(i) = last;
end
end

%% ===================== plotting & stats =====================
function COL = palette()
COL.passive = [0 0.45 0.74];       % blue-ish
COL.active  = [0.85 0.33 0.10];    % orange-ish
COL.grey    = 0.7*[1 1 1];
COL.epochOrder = {'Pre','During','Post','Withdrawal','Reexposure'};
COL.epochMark  = {'o','^','s','d','v'};
end

function W = perMouseEpochTable_v2(D, ycols)
if ischar(ycols) || isstring(ycols), ycols = cellstr(string(ycols)); end

if numel(ycols)==1
    y = ycols{1};
    G = groupsummary(D(:,{'mouse_key','Group','Epoch',y}), {'mouse_key','Group','Epoch'}, 'median', y);
    W = renamevars(G, "median_"+y, "value");
else
    key = unique(D(:,{'mouse_key','Group','Epoch'}),'rows','stable');
    W = key;
    for j=1:numel(ycols)
        y = ycols{j};
        G = groupsummary(D(:,{'mouse_key','Group','Epoch',y}), {'mouse_key','Group','Epoch'}, 'median', y);
        W = outerjoin(W, renamevars(G, "median_"+y, y), 'Keys',{'mouse_key','Group','Epoch'}, 'MergeKeys',true);
    end
end
end

function dualPanelEpochLines_v2(W, ylab)
COL = palette(); epochs = COL.epochOrder;

nexttile; hold on; title(sprintf('Passive — %s', ylab));
drawPanel_v2(W(W.Group=="Passive",:), epochs, COL, ylab);

nexttile; hold on; title(sprintf('Active — %s', ylab));
drawPanel_v2(W(W.Group=="Active",:), epochs, COL, ylab);
end

function drawPanel_v2(Wg, epochs, COL, ylab)
if isempty(Wg)
    text(0.02,0.95,'N=0 mice','Units','normalized');
    set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs);
    xlim([0.7 numel(epochs)+0.3]); grid on; box on; ylabel(ylab); return;
end

mice = unique(Wg.mouse_key,'stable');
for i=1:numel(mice)
    s = Wg(Wg.mouse_key==mice(i),:);
    [x,ord] = sort(categorical(s.Epoch,epochs,'Ordinal',true));
    y = s.value(ord);
    plot(double(x), y, '-o', 'Color',COL.grey, 'MarkerSize',3, 'LineWidth',0.9, ...
        'MarkerFaceColor',COL.grey*0.9, 'MarkerEdgeColor','none');
end

M = groupsummary(Wg, 'Epoch', 'mean', 'value');
E = groupsummary(Wg, 'Epoch', @(x) std(x,'omitnan')./sqrt(sum(isfinite(x))), 'value');
x = double(categorical(M.Epoch,epochs,'Ordinal',true));
y = M.mean_value;
e = E.fun1_value;

fill([x; flipud(x)], [y-e; flipud(y+e)], [0 0 0], 'FaceAlpha',0.12, 'EdgeColor','none');
plot(x, y, 'k-', 'LineWidth',2.0);
errorbar(x, y, e, 'k','LineWidth',1.1,'CapSize',8);

set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs);
xlim([0.7 numel(epochs)+0.3]); grid on; box on; ylabel(ylab);
text(0.02,0.95,sprintf('N=%d mice', numel(mice)),'Units','normalized','FontWeight','bold');
end

function changeFromPreStrip_v2(W, ylab, adjustMethod)
if nargin<3 || isempty(adjustMethod), adjustMethod = 'none'; end
COL = palette();
epochs = COL.epochOrder;

W = W(~isnan(W.value),:);

% baseline per mouse (Pre)
base = groupsummary(W(W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
base.Properties.VariableNames{'median_value'} = 'base';
W = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key', 'MergeKeys',true);
W.delta = W.value - W.base;

hold on
groups = {'Passive','Active'};
off = [-0.14 0.14];

% legend handles
hG(1) = scatter(nan,nan,42,'o','filled','MarkerFaceColor',COL.passive,'MarkerEdgeColor','k','DisplayName','Passive');
hG(2) = scatter(nan,nan,42,'o','filled','MarkerFaceColor',COL.active, 'MarkerEdgeColor','k','DisplayName','Active');

topY = -inf(1,numel(epochs));
for g=1:numel(groups)
    G = W(W.Group==groups{g},:);
    thisCol = COL.passive; if g==2, thisCol = COL.active; end

    for e=1:numel(epochs)
        Ge = G(string(G.Epoch)==epochs{e} & isfinite(G.delta),:);
        if isempty(Ge), continue; end
        x = e + off(g);

        b = boxchart(repmat(x,height(Ge),1), Ge.delta, 'BoxWidth',0.20, ...
            'MarkerStyle','none', 'BoxFaceAlpha',0.18);
        if isprop(b,'BoxFaceColor'), b.BoxFaceColor = thisCol; end
        if isprop(b,'BoxEdgeColor'), b.BoxEdgeColor = thisCol; end
        scatter(repmat(x,height(Ge),1), Ge.delta, 26, 'filled', ...
            'MarkerFaceColor',thisCol, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.75);

        topY(e) = max(topY(e), max(Ge.delta));
    end
end

set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs);
xlim([0.5 numel(epochs)+0.5]); yline(0,'k:'); grid on; box on;
ylabel(sprintf('%s – Pre', ylab));
legend(hG,'Location','best');
title('Change from Pre (per mouse)');

% Passive vs Active ranksum per epoch (on deltas), skip Pre
pRaw  = nan(1,numel(epochs));
valid = false(1,numel(epochs));
for e=1:numel(epochs)
    if strcmp(epochs{e},'Pre'), continue; end
    Pa = W(W.Group=="Passive" & string(W.Epoch)==epochs{e} & isfinite(W.delta),:);
    Ac = W(W.Group=="Active"  & string(W.Epoch)==epochs{e} & isfinite(W.delta),:);
    x = Pa.delta; y = Ac.delta;
    if numel(x)>=2 && numel(y)>=2
        pRaw(e) = ranksum(x,y);
        valid(e)= true;
    end
end

pAdj = pRaw;
switch lower(string(adjustMethod))
    case "holm", if any(valid), pAdj(valid) = holmBonferroni(pRaw(valid)); end
    case "fdr",  if any(valid), pAdj(valid) = fdrBH(pRaw(valid)); end
    otherwise
end

yl = ylim; rngY = diff(yl); pad = 0.06*rngY;
for e=1:numel(epochs)
    if ~valid(e), continue; end
    x1 = e + off(1); x2 = e + off(2);
    yb = max(topY(e), yl(2)-2*pad) + pad;
    plot([x1 x1 x2 x2],[yb-0.5*pad yb yb yb-0.5*pad],'k-','LineWidth',1);
    tag = upper(char(adjustMethod)); if strcmp(tag,'NONE'), tag = 'raw'; end
    txt = sprintf('%s (%s %.3g)', starStr(pAdj(e)), tag, pAdj(e));
    text(mean([x1 x2]), yb+0.2*pad, txt, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',9, 'FontWeight','bold');
    yl(2) = max(yl(2), yb+0.9*pad);
end
ylim(yl);
end

function rmFriedmanWithinGroupFigure_v2(W, ylab)
COL = palette(); epochs = COL.epochOrder;
groups = ["Passive","Active"];

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gi=1:numel(groups)
    G = groups(gi);
    nexttile; hold on;
    title(sprintf('%s — %s', ylab, G));

    Tg = W(W.Group==G,:);
    mice = unique(Tg.mouse_key,'stable');
    if isempty(mice)
        text(0.05,0.9,'N=0','Units','normalized'); axis off; continue;
    end

    % wide matrix: mice x epochs
    X = nan(numel(mice), numel(epochs));
    for i=1:numel(mice)
        s = Tg(Tg.mouse_key==mice(i),:);
        for e=1:numel(epochs)
            r = s(string(s.Epoch)==epochs{e},:);
            if ~isempty(r), X(i,e) = r.value(1); end
        end
    end

    % plot box per epoch
    for e=1:numel(epochs)
        xe = X(:,e);
        xe = xe(isfinite(xe));
        if isempty(xe), continue; end
        boxchart(repmat(e,numel(xe),1), xe, 'BoxWidth',0.22, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
        scatter(repmat(e,numel(xe),1), xe, 22, 'k','filled', 'MarkerFaceAlpha',0.55);
        text(e, max(xe), sprintf('n=%d', numel(xe)), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontWeight','bold','FontSize',9);
    end

    % Friedman p (complete cases only)
    goodRows = all(isfinite(X),2);
    p = NaN;
    if sum(goodRows) >= 3 && size(X,2) >= 3
        p = friedman(X(goodRows,:), 1, 'off');
    end

    set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs);
    xlim([0.5 numel(epochs)+0.5]); grid on; box on;
    ylabel(ylab);
    if isfinite(p)
        subtitle(sprintf('Friedman across epochs (complete mice): p=%.3g, N=%d', p, sum(goodRows)));
    else
        subtitle(sprintf('Friedman n/a (need >=3 epochs and >=3 complete mice).'));
    end
end
end

%% ===================== stats =====================
function ST = runStats_v2(W, ycol)
epochs = {'Pre','During','Post','Withdrawal','Reexposure'};
T = W(:,{'mouse_key','Group','Epoch','value'});
T.mouse_key = categorical(T.mouse_key);
T.Group     = categorical(T.Group);
T.Epoch     = categorical(T.Epoch, epochs, 'Ordinal',true);

rows = {};

% Within-group: paired signrank (Pre vs each other epoch)
for G = ["Passive","Active"]
    A = groupsummary(T(T.Group==G & T.Epoch=="Pre",:), 'mouse_key','median','value');
    A.Properties.VariableNames{'median_value'}='pre';
    for E = ["During","Post","Withdrawal","Reexposure"]
        B = groupsummary(T(T.Group==G & T.Epoch==E,:),'mouse_key','median','value');
        B.Properties.VariableNames{'median_value'}='other';
        M = outerjoin(A(:,{'mouse_key','pre'}), B(:,{'mouse_key','other'}), 'Keys','mouse_key','MergeKeys',true);
        good = isfinite(M.pre) & isfinite(M.other);
        p = NaN; N = nnz(good);
        if N >= 3
            p = signrank(M.pre(good), M.other(good));
        end
        rows(end+1,:) = {'WithinNP', 'Group='+string(G), 'Pre|'+string(E), p, N}; %#ok<AGROW>
    end
end

% Between-group: ranksum at each epoch
for E = ["Pre","During","Post","Withdrawal","Reexposure"]
    Pa = T(T.Group=="Passive" & T.Epoch==E & isfinite(T.value),:);
    Ac = T(T.Group=="Active"  & T.Epoch==E & isfinite(T.value),:);
    x = Pa.value; y = Ac.value;
    p = NaN; N = numel(x)+numel(y);
    if numel(x)>=2 && numel(y)>=2
        p = ranksum(x,y);
    end
    rows(end+1,:) = {'BetweenNP','Epoch='+string(E),'Passive vs Active', p, N}; %#ok<AGROW>
end

ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
ST.metric = repmat(string(ycol), height(ST),1);
end

function ST = runDeltaStats_v2(W, ycol)
epochs = {'Pre','During','Post','Withdrawal','Reexposure'};

% baseline per mouse (Pre)
base = groupsummary(W(W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
base.Properties.VariableNames{'median_value'}='base';
T = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key','MergeKeys',true);
T.delta = T.value - T.base;

rows = {};
for E = ["During","Post","Withdrawal","Reexposure"]
    Te = T(string(T.Epoch)==string(E) & isfinite(T.delta),:);
    Pa = Te(Te.Group=="Passive",:);
    Ac = Te(Te.Group=="Active", :);
    x = Pa.delta; y = Ac.delta;
    p = NaN; N = numel(x)+numel(y);
    if numel(x)>=2 && numel(y)>=2
        p = ranksum(x,y);
    end
    rows(end+1,:) = {'DeltaBetweenNP','Epoch='+string(E),'Passive vs Active (delta)', p, N}; %#ok<AGROW>
end

ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
ST.metric = repmat(string(ycol), height(ST),1);
end

function STrm = runFriedmanStats_v2(W, ycol)
COL = palette(); epochs = COL.epochOrder;
rows = {};

for G = ["Passive","Active"]
    Tg = W(W.Group==G,:);
    mice = unique(Tg.mouse_key,'stable');
    X = nan(numel(mice), numel(epochs));
    for i=1:numel(mice)
        s = Tg(Tg.mouse_key==mice(i),:);
        for e=1:numel(epochs)
            r = s(string(s.Epoch)==epochs{e},:);
            if ~isempty(r), X(i,e) = r.value(1); end
        end
    end

    goodRows = all(isfinite(X),2);
    p = NaN; N = sum(goodRows);
    if N >= 3 && numel(epochs) >= 3
        p = friedman(X(goodRows,:), 1, 'off');
    end
    rows(end+1,:) = {'FriedmanRM','Group='+string(G),'All epochs', p, N}; %#ok<AGROW>
end

STrm = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
STrm.metric = repmat(string(ycol), height(STrm),1);
end

%% ===================== event & bout helpers =====================
function tb = pickTimebase_fast(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
tb = nan(height(T),1);
for i=1:numel(cands)
    if ismember(cands{i}, T.Properties.VariableNames)
        v = double(T.(cands{i}));
        if any(isfinite(v)), tb = v; return; end
    end
end
end

function r = finiteRange_fast(x)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), r = 0; else, r = max(x)-min(x); end
end

function [n, meanDur, totalDur, medianIEI, onTimes] = eventMetrics_fast_v2(t, ttl)
% Returns onset times too (for lick bout metrics and cumulative curves)
t = double(t(:)); ttl = logical(ttl(:));
good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);

onTimes = [];
if numel(t)<2
    n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return;
end

dt = diff(t);
md = median(dt(isfinite(dt)));
if ~isfinite(md), md = 1/30; end

% monotonic enforcement (light)
t(2:end) = max(t(2:end), t(1:end-1)+md*0.5);

d = diff([false; ttl; false]);
on  = find(d==1);
off = find(d==-1)-1;

n = numel(on);
if n==0
    meanDur=NaN; totalDur=0; medianIEI=NaN; return;
end

edges = [t; t(end)+md];
segDur = sum(edges(off+1) - edges(on), 2, 'omitnan');
meanDur  = mean(segDur,'omitnan');
totalDur = sum(segDur,'omitnan');

onTimes = t(on);
if n>=2
    medianIEI = median(diff(onTimes),'omitnan');
else
    medianIEI = NaN;
end
end

function [bn, meanBoutDur, meanBoutSize, medianIBI] = boutMetrics_fromOnsets(onTimes, boutGap_s)
% Bout = consecutive onsets separated by <= boutGap_s
onTimes = double(onTimes(:));
onTimes = onTimes(isfinite(onTimes));
if numel(onTimes) < 1
    bn=0; meanBoutDur=NaN; meanBoutSize=NaN; medianIBI=NaN; return;
end
if numel(onTimes)==1
    bn=1; meanBoutDur=0; meanBoutSize=1; medianIBI=NaN; return;
end

gap = diff(onTimes);
isNewBout = [true; gap > boutGap_s];

boutStartIdx = find(isNewBout);
boutEndIdx = [boutStartIdx(2:end)-1; numel(onTimes)];
bn = numel(boutStartIdx);

boutDur = onTimes(boutEndIdx) - onTimes(boutStartIdx);
boutSize = boutEndIdx - boutStartIdx + 1;

meanBoutDur  = mean(boutDur,'omitnan');
meanBoutSize = mean(boutSize,'omitnan');

if bn >= 2
    IBI = onTimes(boutStartIdx(2:end)) - onTimes(boutEndIdx(1:end-1));
    medianIBI = median(IBI,'omitnan');
else
    medianIBI = NaN;
end
end

%% ===================== misc helpers =====================
function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames)
    T2.(nm) = repmat("",height(T2),1);
    return;
end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end

function pick = pickPctColumn(cols, preferKey)
pick = '';
if isempty(cols), return; end
cols = cellstr(cols);
pctCols = cols(contains(lower(cols),'pct'));
if isempty(pctCols), return; end
pref = pctCols(contains(lower(pctCols), lower(preferKey)));
if ~isempty(pref), pick = pref{1}; else, pick = pctCols{1}; end
end

function nm = colLike(T, key)
nm = '';
if istable(T)
    V = T.Properties.VariableNames;
    hit = find(strcmpi(V,key) | contains(lower(V), lower(key)), 1, 'first');
    if ~isempty(hit), nm = V{hit}; end
end
end

function s = safeName(nm)
s = regexprep(string(nm),'[^a-zA-Z0-9]+','_');
s = char(s);
end

function printpng(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',180);
catch
    print(fh, fn, '-dpng','-r180');
end
end

function s = starStr(p)
if ~isfinite(p), s='n/a'; return; end
if p < 1e-4, s='****';
elseif p < 1e-3, s='***';
elseif p < 1e-2, s='**';
elseif p < 0.05, s='*';
else, s='n.s.';
end
end

function q = fdrBH(p)
p = double(p);
m = numel(p);
[ps,idx] = sort(p);
q = nan(size(p));
ranks = (1:m)';
adj = ps.*m./ranks;
for i=m-1:-1:1
    adj(i)=min(adj(i),adj(i+1));
end
q(idx) = min(adj,1);
end

function p_holm = holmBonferroni(p)
p = double(p(:));
[ps,idx] = sort(p);
m = numel(p);
adj = (m - (1:m)' + 1) .* ps;
for k = 2:m
    adj(k) = max(adj(k), adj(k-1));
end
p_holm = zeros(size(p));
p_holm(idx) = min(adj,1);
end

function h = local_hash(T)
% Hash based on mouse_key/day_index/session_idx + row count.
try
    mk = categorical(string(T.mouse_key));
    di = double(T.day_index);
    if ismember('session_idx',T.Properties.VariableNames)
        si = double(T.session_idx);
    else
        si = zeros(height(T),1);
    end
    raw = [uint32(double(grp2idx(mk))), uint32(di), uint32(si), uint32(height(T))];
    h = uint64(sum(uint64(raw(:)).*1664525 + 1013904223));
catch
    h = uint64(now*1e6);
end
end
