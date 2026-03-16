function run_motivation_analysis_new()
% Motivation analysis from longitudinal CSV (UPDATED for new cohort + new timeline)
%
% NEW in this version:
%  1) Explicit cohort mapping from your list (cage+color -> sex, Active/Passive, PairID).
%  2) Phase labels by day:
%        Pre        = D3–5
%        During     = D6–10
%        Post       = D11–13
%        Withdrawal = D14–16
%        ReExposure = D17–18
%     (Habituation D1–2 is ignored by default in phase plots)
%  3) Two parallel outputs:
%        A) INCLUDE transition days
%        B) EXCLUDE transition days (less reliable): D4, D6, D11, D14, D17
%  4) Group categories changed to: Active / Passive (instead of ActiveOnly / HadPassive)
%
% Outputs under newest run_*:
%   run_*/figs/motivation/includeTransitions/TrialsMotiv.csv
%   run_*/figs/motivation/includeTransitions/SessionsMotiv.csv
%   run_*/figs/motivation/excludeTransitions/TrialsMotiv.csv
%   run_*/figs/motivation/excludeTransitions/SessionsMotiv.csv
%   + PNGs (+ bin_counts.csv)

%% ---------- locate latest run_* & read CSV ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
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

csvPath  = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s', csvPath);
fprintf('Reading: %s\n', csvPath);

T = readtable(csvPath,'VariableNamingRule','preserve');

%% ---------- hygiene + timebase ----------
T = ensureString(T,'mouse_key');

% session index
if ~ismember('session_idx',T.Properties.VariableNames), T.session_idx = ones(height(T),1); end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end

% day_index (robust fallback)
if ~ismember('day_index',T.Properties.VariableNames) || all(~isfinite(T.day_index))
    [G, ~, si] = findgroups(string(T.mouse_key), double(T.session_idx));
    ord = splitapply(@(x) tiedrank(x), si, G);  % 1,2,3,… within mouse
    T.day_index = ord;
else
    if ~isnumeric(T.day_index), T.day_index = double(T.day_index); end
end

% passive flag (optional)
if ismember('isPassive',T.Properties.VariableNames) && ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% licking TTL
assert(ismember('Lick_TTL',T.Properties.VariableNames),'Lick_TTL missing.');
T.Lick_TTL(isnan(T.Lick_TTL)) = 0;

% pick a usable timebase
tb = pickTimebase_local(T);
assert(any(isfinite(tb)),'No usable timebase found (CamTime_rel_s, PlotTime_s_30fps, CamTime_s, or Frame/30).');

%% ---------- NEW: cohort map (cage+color -> sex/group/pair) ----------
cohort = build_cohort_map_local();  % containers.Map("cage_color" -> struct)
pairs  = build_pairs_local();       % containers.Map("cage_color" -> "Pair#")

% Add parsed fields to raw T (best effort, non-fatal)
[ckey, cageN, color] = parse_mouse_key_local(string(T.mouse_key));
T.cohort_key = ckey;
T.cage       = cageN;
T.color      = color;

% Attach Sex / Group / PairID if in cohort
T.sex = repmat("", height(T), 1);
T.group = repmat("", height(T), 1);    % "Active" / "Passive"
T.pair_id = repmat("", height(T), 1);

for i=1:height(T)
    k = string(T.cohort_key(i));
    if strlength(k)==0, continue; end
    if isKey(cohort, k)
        s = cohort(k);
        T.sex(i)     = string(s.sex);
        T.group(i)   = string(s.group);
        if isKey(pairs, k), T.pair_id(i) = string(pairs(k)); end
    end
end

% Fallback: if cohort mapping missing, infer group from isPassive / Session_Paradigm
missingGroup = (strlength(string(T.group))==0);
if any(missingGroup)
    fprintf('Warning: %d rows missing cohort group; using fallback isPassive/Session_Paradigm.\n', nnz(missingGroup));
    g = repmat("Active", height(T),1);
    if ismember('isPassive',T.Properties.VariableNames)
        ip = double(T.isPassive);
        g(ip==1) = "Passive";
    elseif ismember('Session_Paradigm',T.Properties.VariableNames)
        p = T.Session_Paradigm;
        if iscategorical(p), p = string(p); end
        g(contains(lower(string(p)),'passive')) = "Passive";
    end
    T.group(missingGroup) = g(missingGroup);
end

%% ---------- TRIAL table ----------
Trials = build_trial_table_motivation(T, tb);

% Attach cohort-level labels to Trials (mouse-level; stable within mouse)
Trials = attach_cohort_to_trials_local(Trials, T);

% Phase labels (timeline)
Trials.phase = phase_label_local(Trials.day_index);

%% ---------- SESSION table ----------
Sessions = aggregate_session_motivation(Trials, T, tb);
Sessions.phase = phase_label_local(Sessions.day_index);

%% ---------- Run two outputs: include vs exclude transition days ----------
baseOutDir = fullfile(runDir,'figs','motivation');
if ~exist(baseOutDir,'dir'), mkdir(baseOutDir); end

transitionDays = [4 6 11 14 17];   % user-defined "less reliable"
dayMin = 3; dayMax = 18;           % we focus on D3–18 for phase comparisons

run_one_variant(true,  fullfile(baseOutDir,'includeTransitions'));
run_one_variant(false, fullfile(baseOutDir,'excludeTransitions'));

fprintf('Done.\nOutputs under:\n  %s\n', baseOutDir);

%% =================== nested runner ===================
    function run_one_variant(includeTransitions, outDir)
        if ~exist(outDir,'dir'), mkdir(outDir); end
        tag = 'INCLUDE transition days';
        if ~includeTransitions, tag = 'EXCLUDE transition days'; end
        fprintf('\n--- %s ---\nOutDir: %s\n', tag, outDir);

        % Filter to timeline of interest
        keepTrials = (Trials.day_index>=dayMin & Trials.day_index<=dayMax);
        keepSess   = (Sessions.day_index>=dayMin & Sessions.day_index<=dayMax);

        if ~includeTransitions
            keepTrials = keepTrials & ~ismember(Trials.day_index, transitionDays);
            keepSess   = keepSess   & ~ismember(Sessions.day_index, transitionDays);
        end

        Tr = Trials(keepTrials,:);
        Se = Sessions(keepSess,:);

        % Coverage report
        write_bin_coverage_local(Se, outDir);

        % Save tables
        writetable(Tr, fullfile(outDir,'TrialsMotiv.csv'));
        writetable(Se, fullfile(outDir,'SessionsMotiv.csv'));

        % Plots (phase-based)
        plot_requirement_progress(Tr, outDir);
        plot_session_motivation_summaries_BY_PHASE(Se, outDir);

        % Quick dynamics + alt figs (still useful)
        plot_time_to_RLminus1_vs_day(Tr, outDir, true);
        plot_early_late_quick(Tr, outDir, 10);
        alt_plot_delta_vs_day_simple(Tr, outDir, 10, 3);
        make_alt_motivation_figs(Tr, outDir, 10);

        % Survival + scatter (uses session table)
        plot_time_to_RL_survival(Se, outDir);
        plot_RL_vs_TimeToRL(Se, outDir);

        % Correlations (session numeric features)
        corr_session_features(Se, outDir);
    end
end

%% =================== cohort builders ===================
function cohort = build_cohort_map_local()
% Key: "cage_color" (e.g., "6100_red")
% Values: struct('sex',"F"/"M",'group',"Active"/"Passive")
cohort = containers.Map('KeyType','char','ValueType','any');

% From user list:
% 6100red f_s   (ignore trailing "_s" if it appears; treat as female)
% 6100orange f_p
% 6100black f_a
% 0911red f_a
% 0911orange f_p
% 0911black f_p
% 0911white f_a
% 0910 red m_p
% 0910 orange m_p
% 0910 black m_a
% 6099 red_m_p
% 6099 orange_m_a
% 6099 black_m_a
% 6099 white died (not in list; if present in data, you can add it here)

add('6100','red'   ,'F','Passive');  % user wrote f_s; treated as passive? (pair says red passive)
add('6100','orange','F','Passive');
add('6100','black' ,'F','Active');

add('0911','red'   ,'F','Active');
add('0911','orange','F','Passive');
add('0911','black' ,'F','Passive');
add('0911','white' ,'F','Active');

add('0910','red'   ,'M','Passive');
add('0910','orange','M','Passive');
add('0910','black' ,'M','Active');

add('6099','red'   ,'M','Passive');
add('6099','orange','M','Active');
add('6099','black' ,'M','Active');
% If 6099 white exists in your dataset (passive; died D13), uncomment:
% add('6099','white','M','Passive');

    function add(cage,color,sex,group)
        k = sprintf('%s_%s', string(cage), lower(string(color)));
        cohort(char(k)) = struct('sex',string(sex),'group',string(group));
    end
end

function pairs = build_pairs_local()
% Optional: Pair labeling for bookkeeping (not required for analysis)
% Key: "cage_color" -> "Pair#"
pairs = containers.Map('KeyType','char','ValueType','char');

% Pair definitions from user:
% 1) 6100 orange, red (passive) vs 6100 black (active)
setPair('6100',{'orange','red','black'},'Pair1');

% 2) 0911 red (active) vs 0911 orange (passive)
setPair('0911',{'red','orange'},'Pair2');

% 3) 0911 white (active) vs 0911 black (passive)
setPair('0911',{'white','black'},'Pair3');

% 4) 0910 black (active) vs 0910 orange, red (passive)
setPair('0910',{'black','orange','red'},'Pair4');

% 5) 6099 orange (active) vs 6099 red (passive)
setPair('6099',{'orange','red'},'Pair5');

% 6) 6099 black (active) vs 6099 white (passive; died D13)
setPair('6099',{'black','white'},'Pair6');

    function setPair(cage, colors, pid)
        for i=1:numel(colors)
            k = sprintf('%s_%s', string(cage), lower(string(colors{i})));
            pairs(char(k)) = char(pid);
        end
    end
end

function Trials = attach_cohort_to_trials_local(Trials, Traw)
% Attach stable cohort fields to trial table by mouse_key via mode of raw rows
Trials.sex   = repmat("", height(Trials), 1);
Trials.group = categorical(repmat("Active",height(Trials),1), ["Active","Passive"]);
Trials.pair_id = repmat("", height(Trials), 1);
Trials.cage  = nan(height(Trials),1);
Trials.color = repmat("", height(Trials), 1);

for i=1:height(Trials)
    mk = string(Trials.mouse_key(i));
    r  = (string(Traw.mouse_key)==mk);
    if ~any(r), continue; end

    % cage/color
    if ismember('cage',Traw.Properties.VariableNames)
        v = double(Traw.cage(r)); v = v(isfinite(v));
        if ~isempty(v), Trials.cage(i) = mode(v); end
    end
    if ismember('color',Traw.Properties.VariableNames)
        c = string(Traw.color(r)); c = c(strlength(c)>0);
        if ~isempty(c), Trials.color(i) = mode(c); end
    end

    % sex
    if ismember('sex',Traw.Properties.VariableNames)
        s = string(Traw.sex(r)); s = s(strlength(s)>0);
        if ~isempty(s), Trials.sex(i) = mode(s); end
    end

    % group
    if ismember('group',Traw.Properties.VariableNames)
        g = string(Traw.group(r)); g = g(strlength(g)>0);
        if ~isempty(g)
            gg = mode(g);
            if strcmpi(gg,'Passive'), Trials.group(i) = categorical("Passive",["Active","Passive"]);
            else,                    Trials.group(i) = categorical("Active", ["Active","Passive"]);
            end
        end
    end

    % pair
    if ismember('pair_id',Traw.Properties.VariableNames)
        p = string(Traw.pair_id(r)); p = p(strlength(p)>0);
        if ~isempty(p), Trials.pair_id(i) = mode(p); end
    end
end

% Keep legacy-compatible name used in plots
Trials.GroupMouse = Trials.group;
end

%% =================== PHASE LABEL ===================
function c = phase_label_local(di)
% Timeline:
%  Pre        D3-5
%  During     D6-10
%  Post       D11-13
%  Withdrawal D14-16
%  ReExposure D17-18
cats = ["Pre","During","Post","Withdrawal","ReExposure","<other>"];
lab = repmat("<other>", size(di));

lab(di>=3  & di<=5 )  = "Pre";
lab(di>=6  & di<=10)  = "During";
lab(di>=11 & di<=13)  = "Post";
lab(di>=14 & di<=16)  = "Withdrawal";
lab(di>=17 & di<=18)  = "ReExposure";

c = categorical(lab, cats, 'Ordinal', true);
end

function write_bin_coverage_local(Sessions, outDir)
phaseCats = categories(Sessions.phase);
grpCats   = categories(Sessions.GroupMouse);

[~,~,ip] = unique(Sessions.phase);
[~,~,ig] = unique(Sessions.GroupMouse);

counts = accumarray([ip ig], 1, [numel(phaseCats) numel(grpCats)], @sum, 0);
BinCoverage = array2table(counts, 'VariableNames', regexprep(grpCats,'\W','_'), 'RowNames', phaseCats);
writetable(BinCoverage, fullfile(outDir,'bin_counts.csv'),'WriteRowNames',true);
disp('Coverage (sessions per phase × group):'); disp(BinCoverage);
end

%% =================== TRIAL-LEVEL builder ===================
function Trials = build_trial_table_motivation(T, tb)
assert(ismember('Trial',T.Properties.VariableNames),'Trial column missing');
assert(ismember('TrialRequirement',T.Properties.VariableNames),'TrialRequirement missing');

tt = double(tb(:));
trial = double(T.Trial(:));
req   = double(T.TrialRequirement(:));
lick  = logical(T.Lick_TTL(:));

mouse = string(T.mouse_key(:));
day   = double(T.day_index(:));
sess  = double(T.session_idx(:));

% pupil columns optional
hasPupil = ismember('Diameter_px',T.Properties.VariableNames);
if hasPupil, pupil = double(T.Diameter_px(:)); else, pupil = nan(height(T),1); end

good = isfinite(tt) & isfinite(trial) & isfinite(req) & ~isnan(lick);
tt=tt(good); trial=trial(good); req=req(good); lick=lick(good);
mouse=mouse(good); day=day(good); sess=sess(good); pupil=pupil(good);

% sort by mouse/day/session/time
[~,ord] = sortrows(table(mouse,day,sess,tt), {'mouse','day','sess','tt'});
tt=tt(ord); trial=trial(ord); req=req(ord); lick=lick(ord);
mouse=mouse(ord); day=day(ord); sess=sess(ord); pupil=pupil(ord);

G = findgroups(mouse, day, sess, trial);
u = unique(G);

rows = {};
for g = u(:)'
    idx = (G==g);
    mk = mouse(find(idx,1,'first'));
    di = day(find(idx,1,'first'));
    si = sess(find(idx,1,'first'));
    tr = trial(find(idx,1,'first'));

    ri = req(idx);
    r  = nanmode(ri);                       % frame-level mode within trial
    r_int = NaN; if isfinite(r) && r>=1, r_int = round(r); end

    t  = tt(idx);
    y  = lick(idx);
    on = ttl_onsets_local(t,y);             % lick onsets (s)
    on = unique(on,'stable');

    % trial anchors
    t_start = min(t);                       % trial start from clock
    lat1 = NaN; n_before = NaN; t_reward = NaN; t_to_reward = NaN; rate_to_reward = NaN; medIEI = NaN;

    if ~isempty(on)
        lat1 = max(0, on(1)-t_start);
        if isfinite(r_int) && numel(on) >= r_int
            t_reward    = on(r_int);        % reward time = requirement-th lick
            n_before    = r_int;
            t_to_reward = t_reward - on(1);
            if t_to_reward>0, rate_to_reward = n_before / t_to_reward; else, t_to_reward = NaN; end
        else
            n_before = numel(on);
        end
        if numel(on) >= 2
            iei = diff(on); iei = iei(iei>=0.02 & iei<=10);
            if ~isempty(iei), medIEI = median(iei,'omitnan'); end
        end
    end

    overshoot  = NaN; efficiency = NaN; vigor = NaN;
    if isfinite(r_int) && isfinite(n_before)
        overshoot  = max(0, n_before - r_int);
        efficiency = r_int / max(1, n_before);
    end
    if isfinite(medIEI) && medIEI>0, vigor = 1/medIEI; end

    % optional pupil summaries
    preP = NaN; inP = NaN; dP = NaN;
    if hasPupil && ~isempty(on)
        preMask = t >= (on(1)-1) & t < on(1);
        if any(preMask), preP = mean(pupil(preMask),'omitnan'); end
        if isfinite(t_reward), inMask = t >= on(1) & t <= t_reward;
        else,                  inMask = t >= on(1) & t <= (on(1)+2);
        end
        if any(inMask), inP = mean(pupil(inMask),'omitnan'); end
        if isfinite(preP) && isfinite(inP), dP = inP - preP; end
    end

    % NOTE: t_trial_s (duration) will be filled AFTER we know the next trial start
    rows(end+1,:) = {mk, di, si, tr, r, r_int, n_before, overshoot, efficiency, ...
                     NaN, rate_to_reward, lat1, medIEI, vigor, t_reward, t_start, ...
                     preP, inP, dP}; %#ok<AGROW>
end

Trials = cell2table(rows, 'VariableNames', ...
    {'mouse_key','day_index','session_idx','trial','req','req_i','n_licks_before_reward','overshoot','efficiency', ...
     't_trial_s','rate_in_trial_hz','latency_firstlick_s','iei_median_trial_s','vigor_hz','t_reward_s','t_start_s', ...
     'pupil_pre','pupil_in','pupil_delta'});

% Fill true trial DURATIONS: start(n) -> start(n+1)
Trials = sortrows(Trials, {'mouse_key','day_index','session_idx','trial','t_start_s'});

[Gsess, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
SessEnd = splitapply(@(ts, tr) max([ts; tr],[],'omitnan'), Trials.t_start_s, Trials.t_reward_s, Gsess);

for g = 1:max(Gsess)
    idx = find(Gsess==g);
    if numel(idx)>=2
        Trials.t_trial_s(idx(1:end-1)) = Trials.t_start_s(idx(2:end)) - Trials.t_start_s(idx(1:end-1));
    end
    last = idx(end);
    tend = SessEnd(g);
    if isfinite(tend)
        dt = tend - Trials.t_start_s(last);
        if dt>0, Trials.t_trial_s(last) = dt; end
    end
end

% ITI to next trial (reward_n -> next trial start)
Trials.ITI_to_next_s = nan(height(Trials),1);
for i=1:height(Trials)-1
    same = Trials.mouse_key(i)==Trials.mouse_key(i+1) & ...
           Trials.day_index(i)==Trials.day_index(i+1) & ...
           Trials.session_idx(i)==Trials.session_idx(i+1);
    if same && isfinite(Trials.t_reward_s(i)) && isfinite(Trials.t_start_s(i+1))
        Trials.ITI_to_next_s(i) = Trials.t_start_s(i+1) - Trials.t_reward_s(i);
    end
end
end

%% =================== SESSION-LEVEL aggregator ===================
function Sessions = aggregate_session_motivation(Trials, T, tb)
keyS = unique(Trials(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
Sessions = keyS;

% Group from Trials (mode)
Sessions.GroupMouse = categorical(repmat("Active",height(Sessions),1),["Active","Passive"]);
for i=1:height(Sessions)
    r = Trials.mouse_key==Sessions.mouse_key(i) & Trials.day_index==Sessions.day_index(i) & Trials.session_idx==Sessions.session_idx(i);
    if any(r)
        gm = mode(categorical(Trials.GroupMouse(r)));
        if any(strcmpi(string(gm),'Passive'))
            Sessions.GroupMouse(i) = categorical("Passive",["Active","Passive"]);
        else
            Sessions.GroupMouse(i) = categorical("Active", ["Active","Passive"]);
        end
    end
end

% RequirementLast from raw table (if present)
Sessions.RequirementLast = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key), string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
    if any(r) && ismember('RequirementLast', T.Properties.VariableNames)
        v = double(T.RequirementLast(r)); v = v(isfinite(v));
        if ~isempty(v), Sessions.RequirementLast(i) = max(v); end
    end
end

% Session duration (min) from timebase
Sessions.session_min = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key),string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
    times = double(tb(r)); times = times(isfinite(times));
    if ~isempty(times), Sessions.session_min(i) = (max(times)-min(times))/60; end
end

% Per-session fills
Sessions.TrialsCompleted             = groupsummary_fill_local(Trials, Sessions, @(tr) nnz(isfinite(tr.t_reward_s)));
Sessions.MaxRequirementAchieved      = groupsummary_fill_local(Trials, Sessions, @(tr) max(tr.req,[],'omitnan'));
Sessions.TrialVelocity_per_min       = Sessions.TrialsCompleted ./ max(1e-6, Sessions.session_min);

% Licks per min from raw stream
Sessions.LicksPerMin = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key),string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
    if any(r) && Sessions.session_min(i)>0
        on = ttl_onsets_local(double(tb(r)), logical(T.Lick_TTL(r)));
        Sessions.LicksPerMin(i) = numel(on) / Sessions.session_min(i);
    end
end

Sessions.RequirementVelocity_per_min = groupsummary_fill_local(Trials, Sessions, @fit_req_velocity);
Sessions.efficiency_median           = groupsummary_fill_local(Trials, Sessions, @(tr) median(tr.efficiency,'omitnan'));
Sessions.ttrial_median_s             = groupsummary_fill_local(Trials, Sessions, @(tr) median(tr.t_trial_s,'omitnan'));
Sessions.overshoot_median            = groupsummary_fill_local(Trials, Sessions, @(tr) median(tr.overshoot,'omitnan'));
Sessions.median_ITI_s                = groupsummary_fill_local(Trials, Sessions, @(tr) median(tr.ITI_to_next_s,'omitnan'));
Sessions.vigor_median_hz             = groupsummary_fill_local(Trials, Sessions, @(tr) median(tr.vigor_hz,'omitnan'));

[deff, dtt]                          = groupsummary_fill2_local(Trials, Sessions, @early_late_deltas_local);
Sessions.dEff_late_minus_early       = deff;
Sessions.dTtrial_late_minus_early    = dtt;

[effH, tH]                           = groupsummary_fill2_local(Trials, Sessions, @high_req_summaries_local);
Sessions.eff_med_highReq             = effH;
Sessions.ttrial_med_highReq          = tH;

[idle, streakN]                      = groupsummary_fill2_local(Trials, Sessions, @idle_and_streaks_local);
Sessions.IdleFraction                = idle;
Sessions.ITI_streaks_gt30            = streakN;

% Time-to-RL & speed indices
Sessions.TimeToRequirementLast_min      = nan(height(Sessions),1);
Sessions.TimeTo50pctRequirement_min     = nan(height(Sessions),1);
Sessions.SatiationSpeedIndex_RL_per_min = nan(height(Sessions),1);

for i=1:height(Sessions)
    rS = Trials.mouse_key==Sessions.mouse_key(i) & ...
         Trials.day_index==Sessions.day_index(i) & ...
         Trials.session_idx==Sessions.session_idx(i);
    tr = Trials(rS,:);
    RL = Sessions.RequirementLast(i);
    if isfinite(RL) && ~isempty(tr)
        tRL = min(tr.t_reward_s(tr.req==RL),[],'omitnan');
        if isfinite(tRL)
            Sessions.TimeToRequirementLast_min(i) = tRL/60;
            Sessions.SatiationSpeedIndex_RL_per_min(i) = RL / max(1e-6, Sessions.TimeToRequirementLast_min(i));
        end
        halfRL = ceil(0.5*RL);
        tHalf = min(tr.t_reward_s(tr.req>=halfRL),[],'omitnan');
        if isfinite(tHalf)
            Sessions.TimeTo50pctRequirement_min(i) = tHalf/60;
        end
    end
end

% Attach ONLY requested raw features
Sessions = attach_session_means_from_raw_local(Sessions, T);

%% ------- nested helpers -------
    function out = groupsummary_fill_local(TrialsTbl, SessionsTbl, funHandle)
        out = nan(height(SessionsTbl),1);
        for ii=1:height(SessionsTbl)
            r = TrialsTbl.mouse_key==SessionsTbl.mouse_key(ii) & ...
                TrialsTbl.day_index==SessionsTbl.day_index(ii) & ...
                TrialsTbl.session_idx==SessionsTbl.session_idx(ii);
            tr = TrialsTbl(r,:);
            if ~isempty(tr)
                try
                    out(ii) = funHandle(tr);
                catch
                    out(ii) = NaN;
                end
            end
        end
    end

    function [a,b] = groupsummary_fill2_local(TrialsTbl, SessionsTbl, funHandle)
        a = nan(height(SessionsTbl),1); b = nan(height(SessionsTbl),1);
        for ii = 1:height(SessionsTbl)
            r = TrialsTbl.mouse_key==SessionsTbl.mouse_key(ii) & ...
                TrialsTbl.day_index==SessionsTbl.day_index(ii) & ...
                TrialsTbl.session_idx==SessionsTbl.session_idx(ii);
            tr = TrialsTbl(r,:);
            if isempty(tr), continue; end
            try
                [a(ii), b(ii)] = funHandle(tr);
            catch ME
                warning('groupsummary_fill2_local failed for %s d%d s%d: %s', ...
                    string(SessionsTbl.mouse_key(ii)), SessionsTbl.day_index(ii), ...
                    SessionsTbl.session_idx(ii), ME.message);
                a(ii) = NaN; b(ii) = NaN;
            end
        end
    end

    function [dEff, dT] = early_late_deltas_local(tr)
        N = 10;
        tr = sortrows(tr, {'trial','t_start_s'});
        k = min(N, floor(height(tr)/2));
        if k<1, dEff=NaN; dT=NaN; return; end
        goodE = find(isfinite(tr.efficiency));
        goodT = find(isfinite(tr.t_trial_s));
        if numel(goodE)<2*k || numel(goodT)<2*k, dEff=NaN; dT=NaN; return; end
        dEff = mean(tr.efficiency(goodE(end-k+1:end)),'omitnan') - ...
               mean(tr.efficiency(goodE(1:k)),'omitnan');
        dT   = nansum(tr.t_trial_s(goodT(end-k+1:end))) - ...
               nansum(tr.t_trial_s(goodT(1:k)));
    end

    function [effH, tH] = high_req_summaries_local(tr)
        rq = tr.req; rq = rq(isfinite(rq));
        if isempty(rq), effH=NaN; tH=NaN; return; end
        thr = quantile(rq,.67);
        m   = tr.req>=thr & isfinite(tr.t_trial_s);
        effH = median(tr.efficiency(m),'omitnan');
        tH   = median(tr.t_trial_s(m),'omitnan');
    end

    function [idleFrac, streakN] = idle_and_streaks_local(tr)
        idleFrac=NaN; streakN=NaN;
        if ~any(isfinite(tr.t_reward_s)), return; end
        tStart = min(tr.t_start_s,[],'omitnan');
        tEnd   = max(tr.t_reward_s,[],'omitnan');
        trial_time = nansum(tr.t_trial_s);
        total_time = tEnd - tStart;
        if isfinite(total_time) && total_time>0
            idleFrac = max(0, 1 - trial_time/max(1e-6,total_time));
        end
        iti = tr.ITI_to_next_s; iti = iti(isfinite(iti));
        streakN = nnz(iti > 30);
    end

    function SessionsOut = attach_session_means_from_raw_local(SessionsIn, Traw)
        rawList = {'TST_Frames_Non_moving','HOT_Frames_Non_moving','Diameter_px','Immersion_Latency_s'};
        aggFun  = {@median,                 @median,                 @mean,        @median};
        SessionsOut = SessionsIn;
        for k = 1:numel(rawList)
            nm = rawList{k};
            if ~ismember(nm, Traw.Properties.VariableNames), continue; end
            vals = nan(height(SessionsOut),1);
            for ii=1:height(SessionsOut)
                r = strcmp(string(Traw.mouse_key),string(SessionsOut.mouse_key(ii))) & ...
                    Traw.day_index==SessionsOut.day_index(ii) & Traw.session_idx==SessionsOut.session_idx(ii);
                x = double(Traw.(nm)(r)); x = x(isfinite(x));
                if isempty(x), continue; end
                vals(ii) = aggFun{k}(x);
            end
            short = regexprep(nm,'[^a-zA-Z0-9]+','_');
            SessionsOut.(sprintf('raw_%s',short)) = vals;
        end
    end
end

%% =================== PLOTS (phase-based summaries) ===================
function plot_session_motivation_summaries_BY_PHASE(S, outDir)
varsToPlot = {'RequirementLast','TrialsCompleted','TrialVelocity_per_min','RequirementVelocity_per_min', ...
              'LicksPerMin','efficiency_median','ttrial_median_s','overshoot_median','median_ITI_s','vigor_median_hz', ...
              'TimeToRequirementLast_min','TimeTo50pctRequirement_min','SatiationSpeedIndex_RL_per_min'};
labels = containers.Map( ...
varsToPlot, ...
{'RequirementLast','Trials completed','Trials/min','Requirement slope (/min)', ...
 'Licks/min','Efficiency (median)','Trial time (s, median)','Overshoot (median)','ITI (s, median)','Vigor (1/IEI, Hz)', ...
 'Time to RL (min)','Time to 50% RL (min)','Speed index (RL/min)'} );

phases = categories(S.phase);

for k=1:numel(varsToPlot)
    v = varsToPlot{k};
    if ~ismember(v, S.Properties.VariableNames), continue; end
    fig = figure('Color','w','Position',[80 80 1180 460]);
    t = tiledlayout(1,numel(phases),'TileSpacing','compact','Padding','compact');
    for b=1:numel(phases)
        nexttile; hold on
        row = (S.phase==phases{b});
        Xa = S.(v)(row & S.GroupMouse=="Active");  Xa = Xa(isfinite(Xa));
        Xp = S.(v)(row & S.GroupMouse=="Passive"); Xp = Xp(isfinite(Xp));

        if isempty(Xa) && isempty(Xp)
            text(0.5,0.5,'no data in this phase','Units','normalized','HorizontalAlignment','center');
            grid on; box on; set(gca,'XTick',[]); title(phases{b}); ylabel(labels(v));
            continue
        end

        try
            swarmchart(ones(size(Xa))*1, Xa, 16, 'filled','MarkerFaceAlpha',.65);
            swarmchart(ones(size(Xp))*2, Xp, 16, 'filled','MarkerFaceAlpha',.65);
        catch
            jitter = @(n) (rand(n,1)-.5)*.25;
            scatter(1+jitter(numel(Xa)),Xa,16,'filled');
            scatter(2+jitter(numel(Xp)),Xp,16,'filled');
        end

        if ~isempty(Xa)
            m1 = mean(Xa,'omitnan'); s1 = std(Xa,'omitnan')/sqrt(max(1,numel(Xa)));
            errorbar(1, m1, s1, 'k_','LineWidth',1.2);
        end
        if ~isempty(Xp)
            m2 = mean(Xp,'omitnan'); s2 = std(Xp,'omitnan')/sqrt(max(1,numel(Xp)));
            errorbar(2, m2, s2, 'k_','LineWidth',1.2);
        end
        set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'}); title(phases{b});
        ylabel(labels(v)); grid on; box on
    end
    title(t, sprintf('%s by group × phase', labels(v)));
    savepng_local(fig, fullfile(outDir, sprintf('byGroupPhase_%s.png', v))); close(fig);
end
end

%% =================== Existing plots/helpers (unchanged logic) ===================
function plot_requirement_progress(Trials, outDir)
fig = figure('Color','w','Position',[80 80 980 620]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% A: example mice (≤6)
uM = unique(Trials.mouse_key,'stable');
rng(1); pick = uM(randperm(numel(uM), min(6,numel(uM))));
nexttile; hold on
for i=1:numel(pick)
    r = Trials.mouse_key==pick(i) & isfinite(Trials.t_reward_s);
    tr = Trials(r, {'t_reward_s','req'});
    if ~isempty(tr), plot(tr.t_reward_s/60, tr.req, '-', 'LineWidth',1.2); end
end
xlabel('Time (min)'); ylabel('Requirement'); title('Requirement vs time (examples)');
grid on; box on

% B: all sessions light gray
nexttile; hold on
[G, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
Treq = splitapply(@(t,r) {sortrows([t(:) r(:)],1)}, Trials.t_reward_s, Trials.req, G);
for j=1:numel(Treq)
    M = Treq{j}; if isempty(M), continue; end
    plot(M(:,1)/60, M(:,2), '-', 'Color',[.85 .85 .85]);
end
xlabel('Time (min)'); ylabel('Requirement'); title('All sessions');
grid on; box on

% C: max requirement per session by group
nexttile; hold on
[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
Rmax = splitapply(@(r) max(r,[],'omitnan'), Trials.req, Gs);
Gm   = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);
cats = categories(Gm);
for c=1:numel(cats)
    x = Rmax(Gm==cats{c}); x=x(isfinite(x));
    histogram(x,'DisplayStyle','stairs','LineWidth',1.6);
end
legend(cats,'Location','best'); xlabel('Max requirement'); ylabel('Sessions');
title('Max requirement per session'); grid on; box on

% D: time to last reward by group
nexttile; hold on
t_end = splitapply(@(t) max(t,[],'omitnan'), Trials.t_reward_s, Gs);
for c=1:numel(cats)
    y = t_end(Gm==cats{c})/60; y=y(isfinite(y));
    swarmchart(c*ones(size(y)), y, 20, 'filled');
    if ~isempty(y)
        errorbar(c, mean(y,'omitnan'), std(y,'omitnan')/sqrt(numel(y)), 'k_', 'LineWidth',1.6);
    end
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('Time to last reward (min)');
title('Session completion time'); grid on; box on

savepng_local(fig, fullfile(outDir,'requirement_progress.png')); close(fig);
end

function plot_time_to_RL_survival(S, outDir)
fig = figure('Color','w','Position',[80 80 720 520]); hold on
cats = categories(S.GroupMouse);
for c = 1:numel(cats)
    row = S.GroupMouse==cats{c};
    t  = S.TimeToRequirementLast_min(row);
    Tm = S.session_min(row);
    tt = t; miss = isnan(t) & isfinite(Tm);
    tt(miss) = Tm(miss);        % end-of-session censoring
    cens = miss;
    [f,x] = ecdf(tt, 'Censoring', cens);
    plot(x, f, 'LineWidth', 2);
end
legend(cats,'Location','southeast');
xlabel('Time to RequirementLast (min)'); ylabel('P(reached RL)');
title('Time-to-RL (survival with censoring)');
grid on; box on
savepng_local(fig, fullfile(outDir,'survival_time_to_RL.png')); close(fig);
end

function plot_RL_vs_TimeToRL(S, outDir)
fig = figure('Color','w','Position',[80 80 660 520]); hold on
cats = categories(S.GroupMouse);
mk = {'o','^','s','d'};
for c = 1:numel(cats)
    row = S.GroupMouse==cats{c} & isfinite(S.RequirementLast) & isfinite(S.TimeToRequirementLast_min);
    scatter(S.TimeToRequirementLast_min(row), S.RequirementLast(row), 36, 'filled', mk{min(c,numel(mk))});
end
xlabel('Time to RequirementLast (min)'); ylabel('RequirementLast');
title('Height vs time to finish'); legend(cats,'Location','best');
grid on; box on
savepng_local(fig, fullfile(outDir,'scatter_RL_vs_TimeToRL.png')); close(fig);
end

function corr_session_features(Sessions, outDir)
Snum = Sessions(:, vartype('numeric'));
keep = varfun(@(x) nnz(isfinite(x))>10 & nanstd(x)>0, Snum, 'OutputFormat','uniform');
Smat = Snum(:, find(keep));

X = table2array(Smat);
names = Smat.Properties.VariableNames;
if size(X,2) < 2
    warning('Not enough numeric columns for correlation.');
    return;
end

R = corr(X, 'Rows','pairwise', 'Type','Pearson');

fig = figure('Color','w','Position',[60 60 980 800]);

useAbsForClustering = true;
A = R;
if useAbsForClustering, A = abs(A); end
A(isnan(A)) = 0;

D = 1 - A; D = (D + D')/2; D(1:size(D,1)+1:end) = 0;
dvec = squareform(D, 'tovector');

Z   = linkage(dvec, 'average');
ord = optimalleaforder(Z, dvec);

R     = R(ord, ord);
names = names(ord);
imagesc(R,[-1 1]); axis square; colorbar; colormap(parula);
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',45, ...
         'YTick',1:numel(names),'YTickLabel',names);
title('Pearson correlation (session-level features)');
for i=1:numel(names)
    for j=1:numel(names)
        text(j,i,sprintf('%.2f',R(i,j)),'HorizontalAlignment','center','Color','k');
    end
end
savepng_local(fig, fullfile(outDir,'corr_session_features.png')); close(fig);
end

%% =================== ALT figures and utilities (kept) ===================
% (These are copied from your prior version with NO logic change, except they now
%  rely on Trials.GroupMouse categories Active/Passive, which is already true.)

function make_alt_motivation_figs(Trials, outDir, N)
if nargin<3, N = 10; end
if ~exist(outDir,'dir'), mkdir(outDir); end
plot_early_late_quick_alt(Trials, outDir, N);
[dEff, dT, day, grp] = simple_deltas_from_trials_alt(Trials, N);
plot_delta_histos_alt(dEff, dT, grp, outDir);
plot_delta_vs_day_simple_from_vectors(dEff, dT, day, grp, outDir);
plot_highreq_simple_alt(Trials, outDir);
plot_idle_and_streaks_simple_alt(Trials, outDir);
end

function plot_early_late_quick(Trials, outDir, N)
if nargin<3, N=10; end
[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
S = table();
S.mouse_key  = splitapply(@(x) x(1), Trials.mouse_key, Gs);
S.day_index  = splitapply(@(x) x(1), Trials.day_index, Gs);
S.session_idx= splitapply(@(x) x(1), Trials.session_idx, Gs);
S.group      = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

Ee  = nan(height(S),1); El  = Ee;
TTe = Ee; TTl = Ee;
ITIe= Ee; ITIl= Ee;
VGe = Ee; VGl= Ee;

for i=1:height(S)
    r = Gs==i; tr = Trials(r,:);
    tr = sortrows(tr, {'trial','t_start_s'});
    k  = min(N, floor(height(tr)/2));
    if k<1, continue; end

    finiteE = find(isfinite(tr.efficiency));
    finiteT = find(isfinite(tr.t_trial_s));
    finiteI = find(isfinite(tr.ITI_to_next_s));
    finiteV = find(isfinite(tr.vigor_hz));

    if numel(finiteE)>=2*k
        Ee(i) = mean(tr.efficiency(finiteE(1:k)),'omitnan');
        El(i) = mean(tr.efficiency(finiteE(end-k+1:end)),'omitnan');
    end
    if numel(finiteT)>=2*k
        TTe(i) = nansum(tr.t_trial_s(finiteT(1:k)));
        TTl(i) = nansum(tr.t_trial_s(finiteT(end-k+1:end)));
    end
    if numel(finiteI)>=2*k
        ITIe(i)= median(tr.ITI_to_next_s(finiteI(1:k)),'omitnan');
        ITIl(i)= median(tr.ITI_to_next_s(finiteI(end-k+1:end)),'omitnan');
    end
    if numel(finiteV)>=2*k
        VGe(i) = median(tr.vigor_hz(finiteV(1:k)),'omitnan');
        VGl(i) = median(tr.vigor_hz(finiteV(end-k+1:end)),'omitnan');
    end
end

S.EarlyEff=Ee;  S.LateEff=El;
S.EarlyT =TTe;  S.LateT =TTl;
S.EarlyITI=ITIe;S.LateITI=ITIl;
S.EarlyV=VGe;   S.LateV=VGl;

fig = figure('Color','w','Position',[80 80 1200 420]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

col = containers.Map({'Active','Passive'},{[0 0 0],[0 .4 1]});

p = {@(s) [s.EarlyEff s.LateEff]  , 'Efficiency';
     @(s) [s.EarlyT   s.LateT]    , 'Time for N trials (s)';
     @(s) [s.EarlyITI s.LateITI]  , 'ITI (median, s)';
     @(s) [s.EarlyV   s.LateV]    , 'Vigor (Hz)'};

for k=1:4
    nexttile; hold on; title(p{k,2});
    for g = categories(S.group)'
        row = S.group==g{1};
        X = p{k,1}(S(row,:));
        X = X(all(isfinite(X),2),:);
        if isempty(X), continue; end
        plot([1 2], [X(:,1) X(:,2)]', '-', 'Color', col(g{1}), 'LineWidth', 0.8);
        m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')/sqrt(size(X,1));
        errorbar([1 2], m, s, 'k_', 'LineWidth',1.6);
    end
    xlim([0.8 2.2]); set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'}); grid on; box on
end
sgtitle(sprintf('Early vs Late (N=%d trials) per session', N));
if ~exist(outDir,'dir'), mkdir(outDir); end
savepng_local(fig, fullfile(outDir, sprintf('early_late_dynamics_N%d.png', N))); close(fig);
writetable(S, fullfile(outDir, sprintf('EarlyLateBlocks_N%d.csv',N)));
end

% --- The remaining ALT helpers are unchanged (copied from your prior script) ---
function plot_early_late_quick_alt(Trials, outDir, N)
[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
S = table();
S.group = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

K = height(S);
Ee = nan(K,1); El = nan(K,1);
TTe= nan(K,1); TTl = nan(K,1);
ITIe=nan(K,1); ITIl= nan(K,1);
VGe= nan(K,1); VGl = nan(K,1);

for i=1:K
    tr = sortrows(Trials(Gs==i,:), {'trial','t_start_s'});
    k  = min(N, floor(height(tr)/2));
    if k<1, continue; end

    iE = find(isfinite(tr.efficiency));
    iT = find(isfinite(tr.t_trial_s));
    iI = find(isfinite(tr.ITI_to_next_s));
    iV = find(isfinite(tr.vigor_hz));

    if numel(iE)>=2*k
        Ee(i) = mean(tr.efficiency(iE(1:k)),'omitnan');
        El(i) = mean(tr.efficiency(iE(end-k+1:end)),'omitnan');
    end
    if numel(iT)>=2*k
        TTe(i) = nansum(tr.t_trial_s(iT(1:k)));
        TTl(i) = nansum(tr.t_trial_s(iT(end-k+1:end)));
    end
    if numel(iI)>=2*k
        ITIe(i) = median(tr.ITI_to_next_s(iI(1:k)),'omitnan');
        ITIl(i) = median(tr.ITI_to_next_s(iI(end-k+1:end)),'omitnan');
    end
    if numel(iV)>=2*k
        VGe(i) = median(tr.vigor_hz(iV(1:k)),'omitnan');
        VGl(i) = median(tr.vigor_hz(iV(end-k+1:end)),'omitnan');
    end
end

S.EarlyEff=Ee; S.LateEff=El;
S.EarlyT=TTe;  S.LateT=TTl;
S.EarlyITI=ITIe; S.LateITI=ITIl;
S.EarlyV=VGe;   S.LateV=VGl;

fig = figure('Color','w','Position',[80 80 1200 420]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
labs = {'Efficiency','Time for N trials (s)','ITI (median, s)','Vigor (Hz)'};
getXY = {@(r)[r.EarlyEff r.LateEff], @(r)[r.EarlyT r.LateT], @(r)[r.EarlyITI r.LateITI], @(r)[r.EarlyV r.LateV]};
cats = categories(S.group);
for p=1:4
    nexttile; hold on; title(labs{p})
    for c=1:numel(cats)
        R = S(S.group==cats{c},:);
        X = getXY{p}(R); X = X(all(isfinite(X),2),:);
        if isempty(X), continue; end
        plot([1 2],[X(:,1) X(:,2)]','-','Color',[.75 .75 .75])
        m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')/sqrt(size(X,1));
        errorbar([1 2],m,s,'k_','LineWidth',1.6)
    end
    set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'}); xlim([.8 2.2]); grid on; box on
end
sgtitle(sprintf('ALT: Early vs Late (N=%d) per session',N));
savepng_local(fig, fullfile(outDir, sprintf('ALT_early_late_dynamics_N%d.png', N))); close(fig);
end

function [dEff,dT,day,grp] = simple_deltas_from_trials_alt(Trials, N)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
dEff = nan(K,1); dT = nan(K,1);
day  = nan(K,1);
grp  = categorical(repmat("Active",K,1),["Active","Passive"]);
for i=1:K
    tr = sortrows(Trials(Gs==i,:),{'trial','t_start_s'});
    day(i) = tr.day_index(1);
    grp(i) = mode(tr.GroupMouse);
    k = min(N, floor(height(tr)/2)); if k<1, continue; end
    iE = find(isfinite(tr.efficiency));
    iT = find(isfinite(tr.t_trial_s));
    if numel(iE)>=2*k
        dEff(i) = mean(tr.efficiency(iE(end-k+1:end)),'omitnan') - ...
                  mean(tr.efficiency(iE(1:k)),'omitnan');
    end
    if numel(iT)>=2*k
        dT(i) = nansum(tr.t_trial_s(iT(end-k+1:end))) - ...
                nansum(tr.t_trial_s(iT(1:k)));
    end
end
end

function plot_delta_histos_alt(dEff, dT, grp, outDir)
cats = categories(grp);
fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; title('\DeltaEfficiency (late-early)')
for c=1:numel(cats)
    v = dEff(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    histogram(v,'DisplayStyle','stairs','LineWidth',1.6);
end
xline(0,'k:'); legend(cats,'Location','best'); grid on; box on

nexttile; hold on; title('\DeltaTime_N (late-early)  (s)')
for c=1:numel(cats)
    v = dT(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    histogram(v,'DisplayStyle','stairs','LineWidth',1.6);
end
xline(0,'k:'); legend(cats,'Location','best'); grid on; box on

savepng_local(fig, fullfile(outDir,'ALT_delta_histos.png')); close(fig);
end

function plot_delta_vs_day_simple_from_vectors(dEff, dT, day, grp, outDir)
groups = categories(grp);
col.Active  = [0 114 189]/255;
col.Passive = [217 83 25]/255;

fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    function panel(xDay, yVal, ttl, ylab, ax)
        axes(ax); cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for g=1:numel(groups)
            gname = groups{g};
            cg = col.(gname);
            mask = (string(grp)==gname) & isfinite(xDay) & isfinite(yVal);
            scatter(xDay(mask), yVal(mask), 18, 'filled', ...
                'MarkerFaceColor', cg, 'MarkerEdgeColor', 'none');
            if nnz(mask)>=1
                d = xDay(mask); v = yVal(mask);
                [Gb, days] = findgroups(round(d));
                mu = splitapply(@(z) mean(z,'omitnan'), v, Gb);
                [days,ix] = sort(days); mu = mu(ix);
                plot(days, mu, '-', 'Color', cg, 'LineWidth', 2);
                if numel(mu)>2
                    mu_s = movmean(mu, 3, 'omitnan');
                    plot(days, mu_s, '-', 'Color', cg, 'LineWidth', 3, 'HandleVisibility','off');
                end
            end
        end
        yline(0,'k:'); title(ttl); xlabel('day\_index'); ylabel(ylab);
        legend(groups,'Location','best');
    end

ax1=subplot(1,2,1);
panel(day, dEff, '\DeltaEfficiency vs day', '\DeltaEff', ax1);
ax2=subplot(1,2,2);
panel(day, dT,   '\DeltaTime_N vs day (s)', '\DeltaTime (s)', ax2);

savepng_local(fig, fullfile(outDir,'ALT_delta_vs_day_SIMPLE.png')); close(fig);
end

function plot_highreq_simple_alt(Trials, outDir)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
effH = nan(K,1); tH = nan(K,1); grp = categorical(repmat("Active",K,1),["Active","Passive"]);
for i=1:K
    tr = Trials(Gs==i,:);
    grp(i) = mode(tr.GroupMouse);
    rq = tr.req; rq = rq(isfinite(rq));
    if ~isempty(rq)
        thr = quantile(rq,.67); m = tr.req>=thr;
        if any(m)
            effH(i) = median(tr.efficiency(m),'omitnan');
            tH(i)   = median(tr.t_trial_s(m),'omitnan');
        end
    end
    if ~isfinite(effH(i)), effH(i) = median(tr.efficiency,'omitnan'); end
    if ~isfinite(tH(i)),   tH(i)   = median(tr.t_trial_s,'omitnan'); end
end

cats = categories(grp);
fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; title('Efficiency (high-req or median) by group')
for c=1:numel(cats)
    v = effH(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('Efficiency'); grid on; box on

nexttile; hold on; title('Trial time (s, high-req or median) by group')
for c=1:numel(cats)
    v = tH(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('t\_trial (s)'); grid on; box on

savepng_local(fig, fullfile(outDir,'ALT_highreq_simple.png')); close(fig);
end

function plot_idle_and_streaks_simple_alt(Trials, outDir)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
idle = nan(K,1); streaks = nan(K,1); grp = categorical(repmat("Active",K,1),["Active","Passive"]);

for i=1:K
    tr = sortrows(Trials(Gs==i,:),{'trial','t_start_s'});
    grp(i) = mode(tr.GroupMouse);

    if any(isfinite(tr.t_reward_s))
        tStart = min(tr.t_start_s,[],'omitnan');
        tEnd   = max(tr.t_reward_s,[],'omitnan');
        trial_time = nansum(tr.t_trial_s);
        total_time = tEnd - tStart;
        if isfinite(total_time) && total_time>0
            idle(i) = max(0, 1 - trial_time/max(1e-6,total_time));
        end
    end
    iti = tr.ITI_to_next_s; iti = iti(isfinite(iti));
    if ~isempty(iti), streaks(i) = nnz(iti > 30); end
end

cats = categories(grp);
fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; title('Idle fraction by group (per session)')
for c=1:numel(cats)
    v = idle(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('Idle fraction'); grid on; box on

nexttile; hold on; title('Streaks ITI>30s by group (per session)')
for c=1:numel(cats)
    v = streaks(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('# streaks'); grid on; box on

savepng_local(fig, fullfile(outDir,'ALT_idle_and_streaks.png')); close(fig);
end

function y = fit_req_velocity(tr)
y = NaN;
M = tr(~isnan(tr.t_reward_s) & isfinite(tr.req), {'t_reward_s','req'});
if ~isempty(M)
    x = M.t_reward_s/60; yv = M.req;
    if numel(x)>=3
        p = polyfit(x,yv,1);
        y = p(1);
    end
end
end

function alt_plot_delta_vs_day_simple(Trials, outDir, N, win)
if nargin<3, N=10; end
if nargin<4, win=3; end

[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
S = table();
S.mouse_key   = splitapply(@(x) x(1), Trials.mouse_key, Gs);
S.day_index   = splitapply(@(x) x(1), Trials.day_index, Gs);
S.session_idx = splitapply(@(x) x(1), Trials.session_idx, Gs);
S.group       = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

Ee = nan(height(S),1); El = Ee;
TTe= Ee; TTl= Ee;

for i=1:height(S)
    r = (Gs==i);
    tr = Trials(r,:);
    tr = sortrows(tr,{'trial','t_start_s'});
    k = min(N, floor(height(tr)/2));
    if k<1, continue; end

    idxE = find(isfinite(tr.efficiency));
    idxT = find(isfinite(tr.t_trial_s));
    if numel(idxE)>=2*k
        Ee(i) = mean(tr.efficiency(idxE(1:k)),'omitnan');
        El(i) = mean(tr.efficiency(idxE(end-k+1:end)),'omitnan');
    end
    if numel(idxT)>=2*k
        TTe(i) = nansum(tr.t_trial_s(idxT(1:k)));
        TTl(i) = nansum(tr.t_trial_s(idxT(end-k+1:end)));
    end
end

S.dEff  = El - Ee;
S.dTime = TTl - TTe;

col.Active  = [0 114 189]/255;
col.Passive = [217 83 25]/255;
groups = {'Active','Passive'};

    function plot_metric(ax, yvals, ttl, ylab)
        axes(ax); cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for g=1:numel(groups)
            gname = groups{g};
            cg = col.(gname);
            row = (string(S.group)==gname) & isfinite(yvals);
            scatter(S.day_index(row), yvals(row), 22, 'filled', ...
                'MarkerFaceColor', cg, 'MarkerEdgeColor', 'none');

            [Gd, days] = findgroups(S.day_index(row));
            mu = splitapply(@(v) mean(v,'omitnan'), yvals(row), Gd);
            [days, ix] = sort(days); mu = mu(ix);
            plot(days, mu, '-', 'Color', cg, 'LineWidth', 2);

            if numel(mu)>=2
                mu_s = movmean(mu, win, 'omitnan');
                plot(days, mu_s, '-', 'Color', cg, 'LineWidth', 3, 'HandleVisibility','off');
            end
        end
        yline(0,'k:');
        title(ttl); xlabel('day\_index'); ylabel(ylab);
        legend(groups,'Location','best');
    end

fig = figure('Color','w','Position',[80 80 1200 480]);
ax1 = subplot(1,2,1);
plot_metric(ax1, S.dEff,  '\DeltaEfficiency vs day', '\DeltaEff (late - early)');
ax2 = subplot(1,2,2);
plot_metric(ax2, S.dTime, '\DeltaTime_N vs day (s)', '\DeltaTime_N (s, late - early)');

if ~exist(outDir,'dir'), mkdir(outDir); end
savepng_local(fig, fullfile(outDir,'ALT_delta_vs_day_SIMPLE.png'));
close(fig);
end

function plot_time_to_RLminus1_vs_day(Trials, outDir, useLastSessionPerDay)
if nargin < 3, useLastSessionPerDay = true; end
if ~exist(outDir,'dir'), mkdir(outDir); end

[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
S = table();
S.mouse_key   = splitapply(@(x) x(1), Trials.mouse_key, Gs);
S.day_index   = splitapply(@(x) x(1), Trials.day_index, Gs);
S.session_idx = splitapply(@(x) x(1), Trials.session_idx, Gs);
S.GroupMouse  = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

K = height(S);
S.RL             = nan(K,1);
S.TimeToRLm1_min = nan(K,1);

for i = 1:K
    r  = (Gs == i);
    tr = sortrows(Trials(r,:),{'trial','t_start_s'});

    good = isfinite(tr.t_reward_s) & isfinite(tr.req);
    if ~any(good), continue; end
    rq  = double(tr.req(good));
    tR  = double(tr.t_reward_s(good));
    t0  = min(tr.t_start_s,[],'omitnan');

    RL = floor(max(rq));
    S.RL(i) = RL;

    if isfinite(RL) && RL >= 2
        idx = (rq == (RL-1));
        if ~any(idx), idx = (rq >= (RL-1)); end
        if any(idx)
            t_first = min(tR(idx),[],'omitnan');
            if isfinite(t_first) && isfinite(t0)
                S.TimeToRLm1_min(i) = max(0, (t_first - t0)/60);
            end
        end
    end
end

if useLastSessionPerDay
    [Gmd, ~] = findgroups(S.mouse_key, S.day_index);
    keep = false(height(S),1);
    for g = 1:max(Gmd)
        ii = find(Gmd==g);
        if isempty(ii), continue; end
        [~,j] = max(S.session_idx(ii));
        keep(ii(j)) = true;
    end
    S = S(keep,:);
end

cols = struct('Active',[0 114 189]/255, 'Passive',[217 83 25]/255);
cats = categories(S.GroupMouse);

fig = figure('Color','w','Position',[80 80 900 620]); hold on; grid on; box on

for c = 1:numel(cats)
    cg = cats{c};
    row = (S.GroupMouse==cg) & isfinite(S.TimeToRLm1_min) & isfinite(S.day_index);
    if ~any(row), continue; end
    scatter(S.day_index(row), S.TimeToRLm1_min(row), 28, ...
        'filled','MarkerFaceColor',cols.(cg),'MarkerEdgeColor','none', ...
        'DisplayName',cg);
end

for c = 1:numel(cats)
    cg = cats{c};
    row = (S.GroupMouse==cg) & isfinite(S.TimeToRLm1_min) & isfinite(S.day_index);
    if ~any(row), continue; end
    d  = double(S.day_index(row));
    y  = double(S.TimeToRLm1_min(row));

    [Gd, days] = findgroups(d);
    mu  = splitapply(@(v) mean(v,'omitnan'), y, Gd);
    sem = splitapply(@(v) std(v,0,'omitnan')/sqrt(nnz(isfinite(v))), y, Gd);

    [days,ix] = sort(days); mu = mu(ix); sem = sem(ix);
    errorbar(days, mu, sem, 'o-', 'Color', cols.(cg), 'LineWidth', 2, ...
        'MarkerFaceColor','w', 'MarkerSize',5, 'CapSize',0, ...
        'HandleVisibility','off');
end

xlabel('day\_index');
ylabel('Time to (RL - 1) (min)');
title('Time to reach RequirementLast - 1 (per mouse/day)');
legend('Location','best');

xlim([min(S.day_index,[],'omitnan')-0.5, max(S.day_index,[],'omitnan')+0.5]);
ylim([0, max(S.TimeToRLm1_min,[],'omitnan')*1.05]);

savepng_local(fig, fullfile(outDir,'ALT_time_to_RLminus1_vs_day.png'));
close(fig);
end

%% =================== parsing helpers ===================
function [cohort_key, cageN, color] = parse_mouse_key_local(mouse_key)
% Robust parsing of mouse_key into:
%  - cageN: numeric 4-digit cage number if found
%  - color: lower-case color token if found (letters)
%  - cohort_key: "cage_color" if both found else ""
mk = string(mouse_key);
cohort_key = repmat("", size(mk));
cageN = nan(size(mk));
color = repmat("", size(mk));

for i=1:numel(mk)
    s = lower(strtrim(mk(i)));
    if strlength(s)==0, continue; end

    % Find first 4-digit block
    mC = regexp(s, '(\d{4})', 'tokens','once');
    if isempty(mC), continue; end
    cageN(i) = str2double(mC{1});

    % Remove separators and digits, then find a color-like token
    % Accept letters after the 4-digit (e.g., "6100red", "6100_red", "6100-red")
    mK = regexp(s, '\d{4}\s*[_\-\s]*([a-z]+)', 'tokens','once');
    if isempty(mK)
        % fallback: any letter token
        mK = regexp(s, '([a-z]+)', 'tokens','once');
    end
    if ~isempty(mK)
        color(i) = string(mK{1});
    end

    if isfinite(cageN(i)) && strlength(color(i))>0
        cohort_key(i) = sprintf('%04d_%s', cageN(i), color(i));
    end
end
end

%% =================== basic helpers ===================
function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end

function on = ttl_onsets_local(t, ttl, min_gap)
if nargin < 3, min_gap = 0; end
t   = double(t(:));
ttl = ttl(:);
if ~islogical(ttl), ttl = ttl > 0.5; end
good = isfinite(t) & ~isnan(ttl);
t = t(good); ttl = logical(ttl(good));
if numel(t)<2, on = []; return; end
d = diff([false; ttl; false]);
on = t(d(1:end-1)==1);
if min_gap>0 && numel(on)>1
    keep = [true; diff(on) > min_gap];
    on = on(keep);
end
end

function m = nanmode(v)
v = v(isfinite(v));
if isempty(v), m = NaN; else, m = mode(v); end
end

function tb = pickTimebase_local(T)
tb = nan(height(T),1);
if ismember('CamTime_rel_s',T.Properties.VariableNames)
    v = double(T.CamTime_rel_s); if any(isfinite(v)), tb=v; return; end
end
if ismember('PlotTime_s_30fps',T.Properties.VariableNames)
    v = double(T.PlotTime_s_30fps); if any(isfinite(v)), tb=v; return; end
end
if ismember('CamTime_s',T.Properties.VariableNames)
    v = double(T.CamTime_s); if any(isfinite(v)), tb=v; return; end
end
if ismember('Frame',T.Properties.VariableNames)
    v = double(T.Frame); if any(isfinite(v)), tb=v/30; return; end
end
end

function savepng_local(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end
