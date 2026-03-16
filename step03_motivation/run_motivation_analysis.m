function run_motivation_analysis()
% Stand-alone motivation analysis from longitudinal CSV
% Outputs under newest run_*:
%   run_*/figs/motivation/TrialsMotiv.csv
%   run_*/figs/motivation/SessionsMotiv.csv
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

% ---------- hygiene + timebase ----------
T = ensureString(T,'mouse_key');

% session index
if ~ismember('session_idx',T.Properties.VariableNames), T.session_idx = ones(height(T),1); end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end

% day_index (robust fallback so plots never go empty just because of missing days)
if ~ismember('day_index',T.Properties.VariableNames) || all(~isfinite(T.day_index))
    % fallback: per-mouse session order as day index
    [G, mk, si] = findgroups(string(T.mouse_key), double(T.session_idx));
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

%% ---------- mouse group (EVER passive => HadPassive) ----------
HadPassive = classifyHadPassive_mouseLevel(T);
mice = unique(string(T.mouse_key),'stable');
hpmap = containers.Map(mice, repmat("ActiveOnly",numel(mice),1));
for i=1:numel(mice)
    if isKey(HadPassive, mice(i)) && HadPassive(mice(i))
        hpmap(mice(i)) = "HadPassive";
    end
end

%% ---------- TRIAL table ----------
Trials = build_trial_table_motivation(T, tb);

% tag group + time bin
Trials.GroupMouse = repmat(categorical("ActiveOnly",["ActiveOnly","HadPassive"]), height(Trials),1);
for i=1:height(Trials)
    mk = string(Trials.mouse_key(i));
    if isKey(hpmap, mk) && hpmap(mk)=="HadPassive"
        Trials.GroupMouse(i) = categorical("HadPassive",["ActiveOnly","HadPassive"]);
    end
end
Trials.bin5 = bin5_label(Trials.day_index);


%% ---------- SESSION table ----------
Sessions = aggregate_session_motivation(Trials, T, tb);

% coverage report (helps explain any empty panels)
outDir = fullfile(runDir,'figs','motivation');

if ~exist(outDir,'dir'), mkdir(outDir); end

plot_time_to_RLminus1_vs_day(Trials, outDir, true);



binCats = categories(Sessions.bin5);
grpCats = categories(Sessions.GroupMouse);
[~,~,ib] = unique(Sessions.bin5);
[~,~,ig] = unique(Sessions.GroupMouse);
counts = accumarray([ib ig], 1, [numel(binCats) numel(grpCats)]);
BinCoverage = array2table(counts,'VariableNames',regexprep(grpCats,'\W','_'),'RowNames',binCats);
writetable(BinCoverage, fullfile(outDir,'bin_counts.csv'),'WriteRowNames',true);
disp('Bin coverage (sessions per bin × group):'); disp(BinCoverage);

% quick early/late paired lines
plot_early_late_quick(Trials, outDir, 10);

% simple Δ vs day (robust)
alt_plot_delta_vs_day_simple(Trials, outDir, 10, 3);  % N=10, 3-day window

% save outputs
writetable(Trials,   fullfile(outDir,'TrialsMotiv.csv'));
writetable(Sessions, fullfile(outDir,'SessionsMotiv.csv'));

%% ---------- plots ----------
plot_requirement_progress(Trials, outDir);

% Show only the stable-by-bin summaries; the 6 empty ones are replaced by ALT figs below
plot_session_motivation_summaries_STABLE(Sessions, outDir);

% Non-binned, simple alternatives (replace the formerly empty 6)
make_alt_motivation_figs(Trials, outDir, 10);

plot_time_to_RL_survival(Sessions, outDir);
plot_RL_vs_TimeToRL(Sessions, outDir);

%% ---------- Pearson correlations ----------
corr_session_features(Sessions, outDir);

fprintf('Done.\nOutputs:\n  %s\n  %s\n  %s\n', ...
    fullfile(outDir,'TrialsMotiv.csv'), fullfile(outDir,'SessionsMotiv.csv'), fullfile(outDir,'bin_counts.csv'));
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

% ===== Fill true trial DURATIONS from the clock: start(n) -> start(n+1) =====
Trials = sortrows(Trials, {'mouse_key','day_index','session_idx','trial','t_start_s'});

% per-session end (best-effort): use max of start/reward times we saw
[Gsess, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
SessEnd = splitapply(@(ts, tr) max([ts; tr],[],'omitnan'), Trials.t_start_s, Trials.t_reward_s, Gsess);

% assign ends using the *next* trial's start; last trial -> session end
for g = 1:max(Gsess)
    idx = find(Gsess==g);
    if numel(idx)>=2
        Trials.t_trial_s(idx(1:end-1)) = Trials.t_start_s(idx(2:end)) - Trials.t_start_s(idx(1:end-1));
    end
    % last trial duration to session end (if positive), else leave NaN
    last = idx(end);
    tend = SessEnd(g);
    if isfinite(tend)
        dt = tend - Trials.t_start_s(last);
        if dt>0, Trials.t_trial_s(last) = dt; end
    end
end

% ITI to next trial (reward_n -> first lick_{n+1})
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
Sessions.bin5 = bin5_label(Sessions.day_index);

Sessions.GroupMouse = categorical(repmat("ActiveOnly",height(Sessions),1),["ActiveOnly","HadPassive"]);
% fill group from Trials (mode)
for i=1:height(Sessions)
    r = Trials.mouse_key==Sessions.mouse_key(i) & Trials.day_index==Sessions.day_index(i) & Trials.session_idx==Sessions.session_idx(i);
    if any(r), Sessions.GroupMouse(i) = mode(Trials.GroupMouse(r)); end
end

% RequirementLast from raw table
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

% ------- per-session fills using nested helpers -------
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

[deff, dtt]                           = groupsummary_fill2_local(Trials, Sessions, @early_late_deltas_local);
Sessions.dEff_late_minus_early        = deff;
Sessions.dTtrial_late_minus_early     = dtt;

[effH, tH]                             = groupsummary_fill2_local(Trials, Sessions, @high_req_summaries_local);
Sessions.eff_med_highReq               = effH;
Sessions.ttrial_med_highReq            = tH;

[idle, streakN]                        = groupsummary_fill2_local(Trials, Sessions, @idle_and_streaks_local);
Sessions.IdleFraction                  = idle;
Sessions.ITI_streaks_gt30              = streakN;

% ---- Time-to-RL & speed indices ----
Sessions.TimeToRequirementLast_min   = nan(height(Sessions),1);
Sessions.TimeTo50pctRequirement_min  = nan(height(Sessions),1);
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

% ---- Attach ONLY the requested raw features (session-level) ----
Sessions = attach_session_means_from_raw_local(Sessions, T);

% ===== nested helpers (scoped to aggregate_session_motivation) =====
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
    % first/last N trials within the session (use fewer if needed)
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
        % Only these raw features are added to the session table:
        %   TST_Frames_Non_moving (median)
        %   HOT_Frames_Non_moving (median)
        %   Diameter_px           (mean)
        %   Immersion_Latency_s   (median)
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

%% ============== PLOTS & CORR ===================
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

% ---------- NEW: stable subset of by-bin summaries (skip the empty 6) ----------
function plot_session_motivation_summaries_STABLE(S, outDir)
varsToPlot = {'RequirementLast','TrialsCompleted','TrialVelocity_per_min','RequirementVelocity_per_min', ...
              'LicksPerMin','efficiency_median','ttrial_median_s','overshoot_median','median_ITI_s','vigor_median_hz', ...
              'TimeToRequirementLast_min','TimeTo50pctRequirement_min','SatiationSpeedIndex_RL_per_min'};
labels = containers.Map( ...
varsToPlot, ...
{'RequirementLast','Trials completed','Trials/min','Requirement slope (/min)', ...
 'Licks/min','Efficiency (median)','Trial time (s, median)','Overshoot (median)','ITI (s, median)','Vigor (1/IEI, Hz)', ...
 'Time to RL (min)','Time to 50% RL (min)','Speed index (RL/min)'} );

bins = {'D3-5','D6-8','D9-11','D12-14','D15-16','<undef>'};

for k=1:numel(varsToPlot)
    v = varsToPlot{k};
    if ~ismember(v, S.Properties.VariableNames), continue; end
    fig = figure('Color','w','Position',[80 80 1180 460]);
    t = tiledlayout(1,numel(bins),'TileSpacing','compact','Padding','compact');
    for b=1:numel(bins)
        nexttile; hold on
        row = string(S.bin5)==bins{b};
        Xa = S.(v)(row & S.GroupMouse=="ActiveOnly"); Xa = Xa(isfinite(Xa));
        Xp = S.(v)(row & S.GroupMouse=="HadPassive");  Xp = Xp(isfinite(Xp));

        if isempty(Xa) && isempty(Xp)
            text(0.5,0.5,'no data in this bin','Units','normalized','HorizontalAlignment','center');
            grid on; box on; set(gca,'XTick',[]); title(bins{b}); ylabel(labels(v));
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
        set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'}); title(bins{b});
        ylabel(labels(v)); grid on; box on
    end
    title(t, sprintf('%s by group × time bin', labels(v)));
    savepng_local(fig, fullfile(outDir, sprintf('byGroupBin_%s.png', v))); close(fig);
end
end

function c = bin5_label(di)
cats = ["D3-5","D6-8","D9-11","D12-14","D15-16","<undef>"];  % include <undef>
lab = repmat("<undef>", size(di));
lab(di>=3  & di<=5 )  = "D3-5";
lab(di>=6  & di<=8 )  = "D6-8";
lab(di>=9  & di<=11)  = "D9-11";
lab(di>=12 & di<=14)  = "D12-14";
lab(di>=15 & di<=25)  = "D15-16";
c = categorical(lab, cats, 'Ordinal', true);
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

%% =================== helpers ===================
function plot_early_late_quick(Trials, outDir, N)
if nargin<3, N=10; end
[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);

% compute session-level early/late blocks
S = table();
S.mouse_key  = splitapply(@(x) x(1), Trials.mouse_key, Gs);
S.day_index  = splitapply(@(x) x(1), Trials.day_index, Gs);
S.session_idx= splitapply(@(x) x(1), Trials.session_idx, Gs);
S.group      = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

Ee  = nan(height(S),1); El  = Ee;   % efficiency
TTe = Ee; TTl = Ee;                 % trial time (s), SUM over N trials
ITIe= Ee; ITIl= Ee;                 % ITI to next (s), MEDIAN
VGe = Ee; VGl= Ee;                  % vigor (Hz), MEDIAN

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

% ---- plot: per-session paired lines, color by group (black=ActiveOnly, blue=HadPassive) ----
fig = figure('Color','w','Position',[80 80 1200 420]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

col = containers.Map({'ActiveOnly','HadPassive'},{[0 0 0],[0 .4 1]});  % black vs blue

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
        % group mean ± SEM (black bar)
        m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')/sqrt(size(X,1));
        errorbar([1 2], m, s, 'k_', 'LineWidth',1.6);
    end
    xlim([0.8 2.2]); set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'}); grid on; box on
end
sgtitle(sprintf('Early vs Late (N=%d trials) per session', N));
if ~exist(outDir,'dir'), mkdir(outDir); end
savepng_local(fig, fullfile(outDir, sprintf('early_late_dynamics_N%d.png', N))); close(fig);

% optional: write the table for quick QA
writetable(S, fullfile(outDir, sprintf('EarlyLateBlocks_N%d.csv',N)));
end

function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end

function on = ttl_onsets_local(t, ttl, min_gap)
% Detect rising edges (onsets) in a TTL-like boolean vector.
% t      : time vector (seconds)
% ttl    : logical or numeric (thresholded >0.5 if numeric)
% min_gap: optional min separation (s) to collapse very-close onsets
if nargin < 3, min_gap = 0; end
t   = double(t(:));
ttl = ttl(:);
if ~islogical(ttl), ttl = ttl > 0.5; end
good = isfinite(t) & ~isnan(ttl);
t = t(good); ttl = logical(ttl(good));
if numel(t)<2, on = []; return; end
d = diff([false; ttl; false]);   % +1 at rising, -1 at falling
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

function hpmap = classifyHadPassive_mouseLevel(T)
% Return a containers.Map: mouse_id -> true/false (ever passive)
mk = string(T.mouse_key);
mice = unique(mk,'stable');
hpmap = containers.Map(mice, repmat(false,numel(mice),1));
haveIP = ismember('isPassive',T.Properties.VariableNames);
havePar= ismember('Session_Paradigm',T.Properties.VariableNames);
for i = 1:numel(mice)
    r = mk==mice(i);
    flag = false;
    if haveIP
        ip = double(T.isPassive(r)); flag = flag | any(ip==1, 'all');
    end
    if ~flag && havePar
        p = T.Session_Paradigm(r);
        if iscategorical(p), p = string(p); end
        flag = flag | any(contains(lower(string(p)),'passive'), 'all');
    end
    hpmap(mice(i)) = flag;
end
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

%% =================== ALT: SIMPLE FIGURES (no time bins) ===================
function make_alt_motivation_figs(Trials, outDir, N)
if nargin<3, N = 10; end
if ~exist(outDir,'dir'), mkdir(outDir); end

% 1) Early vs Late paired lines (per session) – lighter/simple version
plot_early_late_quick_alt(Trials, outDir, N);

% 2) Build simple per-session deltas from Trials
[dEff, dT, day, grp] = simple_deltas_from_trials_alt(Trials, N);

% 3) Δ histograms by group
plot_delta_histos_alt(dEff, dT, grp, outDir);

% 4) Δ vs day scatter with moving mean (toolbox-free)
plot_delta_vs_day_simple_from_vectors(dEff, dT, day, grp, outDir);

% 5) High-req medians (efficiency & trial time) by group (or median fallback)
plot_highreq_simple_alt(Trials, outDir);

% 6) Idle fraction & # long ITIs by group
plot_idle_and_streaks_simple_alt(Trials, outDir);
end

% ---------- 1) Early vs Late paired lines (ALT) ----------
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
    k  = min(N, floor(height(tr)/2));   % auto-shrink if needed
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

% ---------- 2) Per-session Δs built from Trials ----------
function [dEff,dT,day,grp] = simple_deltas_from_trials_alt(Trials, N)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
dEff = nan(K,1); dT = nan(K,1);
day  = nan(K,1);
grp  = categorical(repmat("ActiveOnly",K,1),["ActiveOnly","HadPassive"]);
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

% ---------- 3) Δ histograms by group ----------
function plot_delta_histos_alt(dEff, dT, grp, outDir)
cats = categories(grp);
fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ΔEff
nexttile; hold on; title('\DeltaEfficiency (late-early)')
for c=1:numel(cats)
    v = dEff(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    histogram(v,'DisplayStyle','stairs','LineWidth',1.6);
end
xline(0,'k:'); legend(cats,'Location','best'); grid on; box on

% ΔTime(N)
nexttile; hold on; title('\DeltaTime_N (late-early)  (s)')
for c=1:numel(cats)
    v = dT(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    histogram(v,'DisplayStyle','stairs','LineWidth',1.6);
end
xline(0,'k:'); legend(cats,'Location','best'); grid on; box on

savepng_local(fig, fullfile(outDir,'ALT_delta_histos.png')); close(fig);
end

% ---------- 4) Δ vs day (simple, toolbox-free) ----------
function plot_delta_vs_day_simple_from_vectors(dEff, dT, day, grp, outDir)
groups = categories(grp);
col.ActiveOnly = [0 114 189]/255;   % MATLAB blue
col.HadPassive = [217 83 25]/255;   % MATLAB orange

fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% helper
    function panel(xDay, yVal, ttl, ylab, ax)
        axes(ax); cla(ax); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for g=1:numel(groups)
            gname = groups{g};
            cg = col.(gname);
            mask = (string(grp)==gname) & isfinite(xDay) & isfinite(yVal);
            scatter(xDay(mask), yVal(mask), 18, 'filled', ...
                'MarkerFaceColor', cg, 'MarkerEdgeColor', 'none');
            if nnz(mask)>=1
                % bin by integer day
                d = xDay(mask); v = yVal(mask);
                [Gb, days] = findgroups(round(d));
                mu = splitapply(@(z) mean(z,'omitnan'), v, Gb);
                [days,ix] = sort(days); mu = mu(ix);
                plot(days, mu, '-', 'Color', cg, 'LineWidth', 2);
                % moving mean for a touch of smoothing
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

% ---------- 5) High-req simple scatters (fallback if none) ----------
function plot_highreq_simple_alt(Trials, outDir)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
effH = nan(K,1); tH = nan(K,1); grp = categorical(repmat("ActiveOnly",K,1),["ActiveOnly","HadPassive"]);
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
    % fallback so the plot isn’t empty for this session
    if ~isfinite(effH(i)), effH(i) = median(tr.efficiency,'omitnan'); end
    if ~isfinite(tH(i)),   tH(i)   = median(tr.t_trial_s,'omitnan'); end
end

cats = categories(grp);
fig=figure('Color','w','Position',[80 80 1000 380]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Efficiency
nexttile; hold on; title('Efficiency (high-req or median) by group')
for c=1:numel(cats)
    v = effH(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('Efficiency'); grid on; box on

% Trial time
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

% ---------- 6) Idle fraction & # long ITIs by group ----------
function plot_idle_and_streaks_simple_alt(Trials, outDir)
[Gs,~]=findgroups(Trials.mouse_key,Trials.day_index,Trials.session_idx);
K = max(Gs);
idle = nan(K,1); streaks = nan(K,1); grp = categorical(repmat("ActiveOnly",K,1),["ActiveOnly","HadPassive"]);

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

% Idle fraction
nexttile; hold on; title('Idle fraction by group (per session)')
for c=1:numel(cats)
    v = idle(grp==cats{c}); v=v(isfinite(v));
    if isempty(v), continue; end
    swarmchart(c*ones(size(v)), v, 18, 'filled');
    errorbar(c, mean(v,'omitnan'), std(v,'omitnan')/sqrt(numel(v)), 'k_', 'LineWidth',1.6);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel('Idle fraction'); grid on; box on

% Streaks ITI>30s
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

% ---------- small utility fits ----------
function y = fit_req_velocity(tr)
y = NaN;
M = tr(~isnan(tr.t_reward_s) & isfinite(tr.req), {'t_reward_s','req'});
if ~isempty(M)
    x = M.t_reward_s/60; yv = M.req;
    if numel(x)>=3
        p = polyfit(x,yv,1);
        y = p(1); % slope per minute
    end
end
end

% (legacy helpers kept for reference)
function [dEff, dT] = early_late_deltas(tr)
N = 10;
dEff = NaN; dT = NaN;
if height(tr) < 2*N, return; end
tr = sortrows(tr, {'trial','t_start_s'});
goodE = isfinite(tr.efficiency);
goodT = isfinite(tr.t_trial_s);
if nnz(goodE) < 2*N || nnz(goodT) < 2*N, return; end
idxE = find(goodE);
idxT = find(goodT);
firstE = idxE(1:N); lastE = idxE(end-N+1:end);
firstT = idxT(1:N); lastT = idxT(end-N+1:end);
dEff = mean(tr.efficiency(lastE), 'omitnan') - mean(tr.efficiency(firstE), 'omitnan');
dT   = nansum(tr.t_trial_s(lastT)) - nansum(tr.t_trial_s(firstT));
end

function [effH, tH] = high_req_summaries(tr)
effH = NaN; tH = NaN;
rq = tr.req; rq = rq(isfinite(rq));
if isempty(rq), return; end
thr = quantile(rq, .67);
maskH = tr.req >= thr & isfinite(tr.t_trial_s);
effH = median(tr.efficiency(maskH),'omitnan');
tH   = median(tr.t_trial_s(maskH),'omitnan');
end

% simple Δ vs day plots (called earlier)
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

S.dEff  = El - Ee;         % Late − Early
S.dTime = TTl - TTe;       % Late − Early (s)

col.ActiveOnly = [0 114 189]/255;
col.HadPassive = [217 83 25]/255;
groups = {'ActiveOnly','HadPassive'};

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
% Time to reach (RequirementLast - 1) vs day_index.
% • Each dot = one mouse/day (optionally last session of that day).
% • Overlays per-day mean ± SEM line for ActiveOnly and HadPassive.

if nargin < 3, useLastSessionPerDay = true; end
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------- Per-session RL and time to RL-1
[Gs, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
S = table();
S.mouse_key   = splitapply(@(x) x(1), Trials.mouse_key, Gs);
S.day_index   = splitapply(@(x) x(1), Trials.day_index, Gs);
S.session_idx = splitapply(@(x) x(1), Trials.session_idx, Gs);
S.GroupMouse  = splitapply(@(g) mode(categorical(g)), Trials.GroupMouse, Gs);

K = height(S);
S.RL             = nan(K,1);
S.TimeToRLm1_min = nan(K,1);   % minutes from session start to first reward at req = RL-1

for i = 1:K
    r  = (Gs == i);
    tr = sortrows(Trials(r,:),{'trial','t_start_s'});

    good = isfinite(tr.t_reward_s) & isfinite(tr.req);
    if ~any(good), continue; end
    rq  = double(tr.req(good));
    tR  = double(tr.t_reward_s(good));
    t0  = min(tr.t_start_s,[],'omitnan');

    RL = floor(max(rq));              % RequirementLast achieved in this session
    S.RL(i) = RL;

    if isfinite(RL) && RL >= 2
        idx = (rq == (RL-1));         % exact RL-1 first
        if ~any(idx), idx = (rq >= (RL-1)); end   % robust fallback
        if any(idx)
            t_first = min(tR(idx),[],'omitnan');
            if isfinite(t_first) && isfinite(t0)
                S.TimeToRLm1_min(i) = max(0, (t_first - t0)/60);
            end
        end
    end
end

%% ------- Keep last session per mouse/day (optional)
if useLastSessionPerDay
    [Gmd, ~] = findgroups(S.mouse_key, S.day_index);
    keep = false(height(S),1);
    for g = 1:max(Gmd)
        ii = find(Gmd==g);
        if isempty(ii), continue; end
        [~,j] = max(S.session_idx(ii));  % last session that day
        keep(ii(j)) = true;
    end
    S = S(keep,:);
end

%% ------- Scatter + per-day mean±SEM lines
cols = struct('ActiveOnly',[0 114 189]/255, 'HadPassive',[217 83 25]/255);
cats = categories(S.GroupMouse);

fig = figure('Color','w','Position',[80 80 900 620]); hold on; grid on; box on

% Scatter points
for c = 1:numel(cats)
    cg = cats{c};
    row = (S.GroupMouse==cg) & isfinite(S.TimeToRLm1_min) & isfinite(S.day_index);
    if ~any(row), continue; end
    scatter(S.day_index(row), S.TimeToRLm1_min(row), 28, ...
        'filled','MarkerFaceColor',cols.(cg),'MarkerEdgeColor','none', ...
        'DisplayName',cg);
end

% Mean ± SEM per day lines
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

    % errorbar with line overlay
    hE = errorbar(days, mu, sem, 'o-', 'Color', cols.(cg), 'LineWidth', 2, ...
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


function [idleFrac, streakN] = idle_and_streaks(tr)
idleFrac = NaN; streakN = NaN;
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
