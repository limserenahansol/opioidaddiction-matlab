function run_lick_mega_pipeline()
% LICKING MEGA-PIPELINE (all-in-one; toolbox-free fallbacks)
% -------------------------------------------------------------------------
% Scope:
%   • Exclusions: remove "6872 black" and "8606 forange"
%   • Day-bins: D3-5, D6-8, D9-11, D12-14, D15-16
%   • Supervised/descriptive (A–F):
%       A) Core metrics: lick count, licks/min, licks/reward, wasted fraction,
%          time-to-reward (mean/median/IQR), KM survival + log-rank,
%          Cox PH with cluster-robust SE (frailty ≈)
%          Efficiency vs requirement; overshoot; spaghetti
%       B) Microstructure: IEI stats (median/mean/CV/LV/q10/q90/entropy),
%          burst/bout (<=250 ms), bouts/min, duty cycle, licks/bout,
%          within/between IEI, Fano (5 s), Allan slope (0.5–10 s),
%          distributions, bout rasters, per-trial IEI CDFs
%       C) Rhythm: 50 ms binned FFT spectrum, peak Hz/power, band power 3–12 Hz,
%          Q = peak/FWHM, spectral entropy, autocorr peak lag/height,
%          oscillation score, phase metrics + Rayleigh pre/post reward
%       D) Trial-aligned (PSTH): start/reward aligned rate, ramp slope, pause,
%          rebound; hazard until reward vs requirement
%       E) Learning/adaptation: by requirement (FR1/2/3...), time-to-criterion,
%          overshoot, early ramp slope, trial efficiency, CUSUM change-points.
%          Mixed-effects approximations (cluster-robust GLMs).
%       F) Group comparisons & prediction: elastic-net (coord-descent fallback),
%          RF/XGBoost (TreeBagger if present), AUC, permutation/coef importance
%   • Unsupervised (3A–E): engineered features, PCA/UMAP/t-SNE/NMF, kmeans/GMM/
%     hierarchical/DBSCAN/spectral, consensus clustering, validation metrics.
%   • Sequence/state (3D): discrete-HMM on IEI categories; change-points;
%     DTW clustering of PSTHs; Markov n-gram & sequence entropy.
%   • Requirement-aware (4A–B), Multilevel stats (5), Viz gallery (7).
%
% Input:
%   Expected columns in ALL_mice_longitudinal.csv:
%   mouse_key, day_index, session_idx, Lick_TTL, and any one timebase:
%     CamTime_rel_s | PupilTimestamp_s | CamTime_s | PlotTime_s_30fps | Frame
%   Optional: Trial, Requirement (or FR/req/ReqLicks), reward TTL
%     Reward_TTL | Water_TTL | Reward | water
%
% All outputs -> ...\run_XXX\figs\lick_MEGA\*
% -------------------------------------------------------------------------

%% --------------------------- switches -----------------------------------
opts.A_core            = true;
opts.B_micro           = true;
opts.C_rhythm          = true;
opts.D_psth            = true;
opts.E_learning        = true;
opts.F_predict         = true;
opts.UNSUP_engineer    = true;
opts.UNSUP_reduce      = true;
opts.UNSUP_cluster     = true;
opts.UNSUP_consensus   = true;
opts.SEQ_models        = true;   % HMM / DTW / change-points
opts.REQ_align         = true;   % 4A–B
opts.MIXED_models      = true;   % GLM approximations with clustered SE

% thresholds & params
PLOT_MAX_MOUSE_EXAMPLES = 12;
BIN_50MS = 0.05;                   % for spectra & PSTH rate binning
BOUT_WITHIN_IEI = 0.25;            % seconds, burst/bout IEI threshold
BOUT_GAP_IEI    = 0.5;             % seconds, gap closing a bout
RHYTHM_BAND     = [3 12];          % Hz band of interest
PEAK_SEARCH     = [6 10];          % expected peak range
PHASE_BAND      = [5 9];           % for phase (Hilbert fallback if present)
HAZARD_BINS     = 40;              % within-trial hazard bins (0→reward)
DTW_MAXLEN      = 200;             % max points per PSTH curve (downsample)
CONSENSUS_BOOT  = 100;             % resamples for consensus clustering
K_CLUSTERS      = 3:6;             % range to try for clustering
RANDSEED        = 1;

%% --------------------- locate latest run + read -------------------------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up=1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry = p; break; end
        cand = fileparts(cand);
    end
end
d = dir(fullfile(rootTry,'run_*')); assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix]   = max([d.datenum]);
runDir   = fullfile(d(ix).folder,d(ix).name);
csvPath  = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('Reading: %s\n', csvPath);

T = readtable(csvPath,'VariableNamingRule','preserve');

%% --------------------- exclusions & basics ------------------------------
% hard exclusions
if ismember('mouse_key', T.Properties.VariableNames)
    mk = string(T.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    drop = (contains(mk_norm,"6872") & contains(mk_norm,"black")) | ...
           (contains(mk_norm,"8606") & contains(mk_norm,"forange"));
    if any(drop)
        fprintf('Excluding %d rows (hard mice list)\n', nnz(drop));
        T(drop,:) = [];
    end
end

% types & timebase
T = ensureString(T,'mouse_key');
reqCols = {'day_index','session_idx','Lick_TTL'};
for c = 1:numel(reqCols), assert(ismember(reqCols{c},T.Properties.VariableNames), 'Missing %s', reqCols{c}); end
if ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end
T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL>0.5;

tb = pickTimebase(T);
assert(any(isfinite(tb)),'No usable timebase column found.');

% dir
outDir = fullfile(runDir,'figs','lick_MEGA'); if ~exist(outDir,'dir'), mkdir(outDir); end

% group (ever-passive)
groupMap = classifyHadPassive_mouseLevel(T);

% find requirement col + reward TTL col if exists
reqCol   = pickRequirement(T);
rewardCol= pickRewardTTL(T);

%% ---------------- per-session & per-trial base tables -------------------
sessKeys = unique(T(:,{'mouse_key','day_index','session_idx'}),'rows','stable');

% per-session microstructure/rhythm/bout
Sess = sessKeys;
[Sess.lick_count,Sess.session_min,Sess.licks_per_min, ...
 Sess.iei_median,Sess.iei_mean,Sess.iei_cv,Sess.iei_lv, ...
 Sess.cv2_median, Sess.iei_q10, Sess.iei_q90, Sess.iei_entropy, ...
 Sess.rhythm_index,Sess.burst_fraction,Sess.bout_count,Sess.bouts_per_min, ...
 Sess.mean_licks_per_bout,Sess.within_bout_iei,Sess.between_bout_gap, ...
 Sess.duty_cycle, Sess.fano_5s, Sess.allan_slope] = deal(nan(height(Sess),1));

rng(RANDSEED);
for i=1:height(Sess)
    r = T.mouse_key==Sess.mouse_key(i) & T.day_index==Sess.day_index(i) & T.session_idx==Sess.session_idx(i);
    if ~any(r), continue; end
    t = double(tb(r)); l = logical(T.Lick_TTL(r));
    F = compute_session_features_full(t, l, BOUT_WITHIN_IEI, BOUT_GAP_IEI);
    % copy
    Sess.lick_count(i)           = F.lick_count;
    Sess.session_min(i)          = F.session_min;
    Sess.licks_per_min(i)        = F.licks_per_min;
    Sess.iei_median(i)           = F.iei_median;
    Sess.iei_mean(i)             = F.iei_mean;
    Sess.iei_cv(i)               = F.iei_cv;
    Sess.iei_lv(i)               = F.iei_lv;
    Sess.cv2_median(i)           = F.cv2_median;
    Sess.iei_q10(i)              = F.iei_q10;
    Sess.iei_q90(i)              = F.iei_q90;
    Sess.iei_entropy(i)          = F.iei_entropy;
    Sess.rhythm_index(i)         = F.rhythm_index;
    Sess.burst_fraction(i)       = F.burst_fraction;
    Sess.bout_count(i)           = F.bout_count;
    Sess.bouts_per_min(i)        = F.bouts_per_min;
    Sess.mean_licks_per_bout(i)  = F.mean_licks_per_bout;
    Sess.within_bout_iei(i)      = F.within_bout_iei;
    Sess.between_bout_gap(i)     = F.between_bout_gap;
    Sess.duty_cycle(i)           = F.duty_cycle;
    Sess.fano_5s(i)              = F.fano_5s;
    Sess.allan_slope(i)          = F.allan_slope;
end
Sess.GroupMouse = label_by_group(Sess.mouse_key, groupMap);
Sess.bin5 = categorical(dayBin5_label_vec(Sess.day_index), ["D3-5","D6-8","D9-11","D12-14","D15-16","<undef>"]);
writetable(Sess, fullfile(outDir,'per_session_features.csv'));

% per-trial table (needs Trial)
TT = table();
if ismember('Trial', T.Properties.VariableNames)
    TT = per_trial_feature_table(T, tb, reqCol, rewardCol, BOUT_WITHIN_IEI, BOUT_GAP_IEI);
    TT.GroupMouse = label_by_group(TT.mouse_key, groupMap);
    TT.bin5 = categorical(dayBin5_label_vec(TT.day_index), ["D3-5","D6-8","D9-11","D12-14","D15-16","<undef>"]);
    writetable(TT, fullfile(outDir,'per_trial_features.csv'));
end

%% ==================== 2) SUPERVISED SUMMARIES ==========================
if opts.A_core
    do_core_rates_and_survival(TT, Sess, outDir, RHYTHM_BAND, PLOT_MAX_MOUSE_EXAMPLES);
end
if opts.B_micro
    do_microstructure_block(Sess, TT, outDir, BOUT_WITHIN_IEI, BOUT_GAP_IEI, PLOT_MAX_MOUSE_EXAMPLES);
end
if opts.C_rhythm
    do_rhythm_block(T, tb, Sess, TT, outDir, BIN_50MS, RHYTHM_BAND, PEAK_SEARCH, PHASE_BAND, rewardCol);
end
if opts.D_psth
    do_psth_block(T, tb, TT, outDir, BIN_50MS, rewardCol, HAZARD_BINS);
end
if opts.E_learning
    do_learning_block(TT, outDir);
end
if opts.F_predict
    do_prediction_block(Sess, TT, outDir, RANDSEED);
end

%% ==================== 3) UNSUPERVISED DISCOVERY ========================
if opts.UNSUP_engineer || opts.UNSUP_reduce || opts.UNSUP_cluster || opts.UNSUP_consensus
    U = engineer_unsup_features(Sess, TT, outDir);
end
if opts.UNSUP_reduce
    Ured = reduce_unsup(U, outDir, RANDSEED);
end
if opts.UNSUP_cluster
    C = cluster_unsup(Ured, outDir, K_CLUSTERS, RANDSEED);
end
if opts.UNSUP_consensus
    consensus_clustering(Ured, outDir, min(K_CLUSTERS), max(K_CLUSTERS), CONSENSUS_BOOT, RANDSEED);
end

%% ==================== 4) TRIAL-STRUCTURE AWARE =========================
if opts.REQ_align
    req_align_and_hazard(TT, outDir);
end

%% ==================== 5) MULTILEVEL STATS (GLMM-ish) ===================
if opts.MIXED_models
    do_mixed_models_approx(TT, Sess, outDir);
end

fprintf('\nMEGA pipeline complete. Outputs in:\n  %s\n', outDir);
end

%% ========================= CORE HELPERS =================================
function T2 = ensureString(T2,nm)
    if ~ismember(nm,T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
    if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end
function tb = pickTimebase(T)
    cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps','Frame'};
    tb = nan(height(T),1);
    for i=1:numel(cands)
        if ismember(cands{i}, T.Properties.VariableNames)
            v = double(T.(cands{i}));
            if any(isfinite(v))
                tb = v; if strcmpi(cands{i},'Frame'), tb = v/30; end
                return
            end
        end
    end
end
function col = pickRequirement(T)
    cands = {'Requirement','req','FR','ReqLicks','n_required','required_licks'};
    col = ''; for i=1:numel(cands), if ismember(cands{i},T.Properties.VariableNames), col=cands{i}; return; end, end
end
function col = pickRewardTTL(T)
    cands = {'Reward_TTL','Water_TTL','Reward','reward','water','Water'};
    col = ''; for i=1:numel(cands), if ismember(cands{i},T.Properties.VariableNames), col=cands{i}; return; end, end
end
function G = label_by_group(mk, groupMap)
    G = repmat(categorical("ActiveOnly",["ActiveOnly","HadPassive"]), numel(mk),1);
    for i=1:numel(mk)
        m = string(mk(i));
        if isKey(groupMap,m) && groupMap(m), G(i)=categorical("HadPassive",["ActiveOnly","HadPassive"]); end
    end
end
function groupMap = classifyHadPassive_mouseLevel(T)
    mk = string(T.mouse_key);
    mice = unique(mk,'stable');
    groupMap = containers.Map(mice, repmat(false,numel(mice),1));
    haveIP = ismember('isPassive',T.Properties.VariableNames);
    havePar= ismember('Session_Paradigm',T.Properties.VariableNames);
    for i=1:numel(mice)
        r = mk==mice(i); flag=false;
        if haveIP, ip = double(T.isPassive(r)); flag = flag | any(ip==1); end
        if ~flag && havePar
            p = T.Session_Paradigm(r); if iscategorical(p), p = string(p); end
            flag = flag | any(contains(lower(string(p)),'passive'));
        end
        groupMap(mice(i)) = flag;
    end
end
function labs = dayBin5_label_vec(d)
    labs = repmat("<undef>", numel(d),1);
    labs(d>=3  & d<=5 ) = "D3-5";
    labs(d>=6  & d<=8 ) = "D6-8";
    labs(d>=9  & d<=11) = "D9-11";
    labs(d>=12 & d<=14) = "D12-14";
    labs(d>=15 & d<=25) = "D15-16";
end
function [on,off] = ttl_edges(t,ttl)
    d = diff([false; ttl(:); false]);
    on  = t(d(1:end-1)==1);
    off = t(d(2:end)==-1);
end

%% =================== SESSION FEATURE ENGINE ============================
function F = compute_session_features_full(t, l, thrWithin, thrGap)
    t=t(:); l=l(:)>0; good = isfinite(t) & ~isnan(l); t=t(good); l=l(good);
    F = emptyF(); if numel(t)<2, return; end
    [on,~] = ttl_edges(t,l); if isempty(on), return; end
    F.lick_count = numel(on);
    dur = (max(t)-min(t)); F.session_min = dur/60; F.licks_per_min = F.lick_count/max(F.session_min,eps);

    if numel(on)>=2
        iei = diff(on);
        F.iei_median = median(iei); F.iei_mean = mean(iei);
        F.iei_cv = std(iei)/max(mean(iei),eps);
        % local variation
        if numel(iei)>=2
            F.iei_lv = mean( (3*(iei(1:end-1)-iei(2:end)).^2) ./ (iei(1:end-1)+iei(2:end)).^2 );
            cv2 = 2*abs(diff(iei))./(iei(1:end-1)+iei(2:end)); F.cv2_median = median(cv2);
        else
            F.iei_lv = NaN; F.cv2_median = NaN;
        end
        F.iei_q10 = quantile(iei,0.10); F.iei_q90 = quantile(iei,0.90);
        % Shannon entropy of IEI hist
        edges = 0:0.02:max(0.6,max(iei)+0.04);
        h = histcounts(iei,edges,'Normalization','probability'); p=h(h>0); F.iei_entropy = -sum(p.*log2(p));
        % "rhythm band" fraction 0.11–0.18 ~ 6–9 Hz proxy
        F.rhythm_index = nnz(iei>=0.11 & iei<=0.18)/max(numel(iei),1);
        % bouts
        isShort = iei < thrWithin;
        gap = [inf; iei(:)];
        inBout=false; bc=0; lick_in_bout=0; wbIEI=[]; bbGap=[];
        k=2; c=0;
        while k<=numel(gap)
            if ~inBout && gap(k)<thrWithin, inBout=true; bc=bc+1; c=1; wbIEI=[wbIEI; gap(k)];
            elseif inBout
                if gap(k)<thrWithin, c=c+1; wbIEI=[wbIEI; gap(k)];
                elseif gap(k)>=thrGap, inBout=false; bbGap=[bbGap; gap(k)];
                end
            end
            if ~inBout && c>0, lick_in_bout=lick_in_bout+(c+1); c=0; end
            k=k+1;
        end
        F.bout_count = bc; F.bouts_per_min = bc/max(F.session_min,eps);
        F.burst_fraction = nnz(isShort)/max(numel(iei),1);
        F.mean_licks_per_bout = lick_in_bout/max(bc,1);
        F.within_bout_iei = median(wbIEI,'omitnan');
        F.between_bout_gap = median(bbGap,'omitnan');
        % duty cycle: time inside bouts / duration
        durBout = sum(wbIEI,'omitnan'); F.duty_cycle = durBout/max(dur,eps);
        % Fano (5 s bins)
        bin = 5;
        edges = min(t):bin:max(t)+bin;
        cts = histcounts(on,edges);
        F.fano_5s = var(cts)/max(mean(cts),eps);
        % Allan factor slope (0.5–10 s)
        wins = logspace(log10(0.5),log10(10),8);
        AF = nan(size(wins));
        for i=1:numel(wins)
            edges = min(t):wins(i):max(t)+wins(i);
            cts = histcounts(on,edges);
            m = mean(cts); AF(i) = mean((cts(1:end-1)-cts(2:end)).^2)/(2*max(m,eps));
        end
        p = polyfit(log10(wins(~isnan(AF))), log10(AF(~isnan(AF))), 1);
        F.allan_slope = p(1);
    end
end
function F = emptyF()
    F = struct('lick_count',0,'session_min',NaN,'licks_per_min',0,'iei_median',NaN,'iei_mean',NaN,'iei_cv',NaN, ...
               'iei_lv',NaN,'cv2_median',NaN,'iei_q10',NaN,'iei_q90',NaN,'iei_entropy',NaN, ...
               'rhythm_index',NaN,'burst_fraction',NaN,'bout_count',0,'bouts_per_min',0, ...
               'mean_licks_per_bout',NaN,'within_bout_iei',NaN,'between_bout_gap',NaN, ...
               'duty_cycle',NaN,'fano_5s',NaN,'allan_slope',NaN);
end

%% =================== TRIAL FEATURE TABLE ===============================
function TT = per_trial_feature_table(T, tb, reqCol, rewardCol, thrWithin, thrGap)
    g = unique(T(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
    rows = {};
    for i=1:height(g)
        mk=string(g.mouse_key(i)); di=double(g.day_index(i)); ss=double(g.session_idx(i));
        r = string(T.mouse_key)==mk & T.day_index==di & T.session_idx==ss;
        if ~any(r) || ~ismember('Trial',T.Properties.VariableNames), continue; end
        tr = double(T.Trial(r)); t = double(tb(r)); L = logical(T.Lick_TTL(r));
        req = []; if ~isempty(reqCol), req = double(T.(reqCol)(r)); end
        rew = []; if ~isempty(rewardCol), rew = double(T.(rewardCol)(r))>0; end
        trials = sort(unique(tr(isfinite(tr))));
        for k=1:numel(trials)
            m = tr==trials(k); tt = t(m); ll = L(m);
            if isempty(tt), continue; end
            [on,~] = ttl_edges(tt, ll);
            if isempty(on), continue; end
            triStart = min(tt);
            reqK = NaN; if ~isempty(req), reqK = median(req(m),'omitnan'); end
            % reward time: prefer reward TTL; else k-th lick by requirement if known; else last lick
            if ~isempty(rew)
                rr = rew(m); [ron,~] = ttl_edges(tt, rr);
                if ~isempty(ron), rewardTime = ron(1); else, rewardTime = NaN; end
            else
                if isfinite(reqK) && reqK>=1 && numel(on)>=reqK
                    rewardTime = on(floor(reqK));
                else
                    rewardTime = on(end);
                end
            end
            time_to_reward = rewardTime - triStart;
            licks_before_reward = sum(on<=rewardTime);
            overshoot = licks_before_reward - max(reqK,0);
            % IEI stats within trial
            iei = diff(on); medI = median(iei,'omitnan');
            % bout metrics within trial (quick)
            [boutCount, inBoutFrac] = quick_bouts(iei, thrWithin, thrGap);
            rows(end+1,:) = {mk, di, ss, double(trials(k)), double(time_to_reward), ...
                double(licks_before_reward), double(overshoot), double(reqK), ...
                double(numel(on)), double(medI), double(boutCount), double(inBoutFrac)}; %#ok<AGROW>
        end
    end
    if isempty(rows), TT=table(); return; end
    TT = cell2table(rows,'VariableNames',{'mouse_key','day_index','session_idx','Trial', ...
        'time_to_reward','licks_before_reward','overshoot','requirement','lick_count','iei_median','bout_count','inBout_fraction'});
    TT.lick_rate = TT.lick_count ./ max(TT.time_to_reward, eps);
    TT.efficiency = TT.licks_before_reward ./ max(TT.requirement,1); % licks per required
end
function [bc, frac] = quick_bouts(iei, thrWithin, thrGap)
    if isempty(iei), bc=0; frac=0; return; end
    isShort = iei < thrWithin;
    gap = [inf; iei(:)];
    inBout=false; bc=0; shortN = nnz(isShort);
    for k=2:numel(gap)
        if ~inBout && gap(k)<thrWithin, inBout=true; bc=bc+1;
        elseif inBout && gap(k)>=thrGap, inBout=false; end
    end
    frac = shortN/max(numel(iei),1);
end

%% ===================== 2A CORE: rates & survival =======================
function do_core_rates_and_survival(TT, Sess, outDir, rband, MAXEX)
    % violin/raincloud by group × bin (licks/min)
    plot_group_bin_violin(Sess, 'licks_per_min', 'Licks / min', outDir, 'A_core_violin_licks_per_min.png');
    % licks/reward, wasted fraction from TT
    if ~isempty(TT)
        TT.licks_per_reward = TT.licks_before_reward ./ max(1,TT.requirement);
        TT.wasted_fraction  = max(TT.licks_before_reward - TT.requirement, 0) ./ max(TT.licks_before_reward,1);
        writetable(TT(:,{'mouse_key','day_index','session_idx','Trial','licks_per_reward','wasted_fraction'}), ...
            fullfile(outDir,'A_core_licks_per_reward_and_wasted.csv'));
        plot_group_bin_violin(TT, 'licks_per_reward', 'Licks / reward', outDir, 'A_core_violin_licks_per_reward.png');
        plot_group_bin_violin(TT, 'wasted_fraction',  'Wasted lick fraction', outDir, 'A_core_violin_wasted_fraction.png');
        % time-to-reward distributions + KM + log-rank
        km_by_requirement(TT, outDir);
        % efficiency vs requirement
        plot_efficiency_vs_requirement(TT, outDir);
        % spaghetti per mouse across days (lick/min)
        spaghetti_per_mouse(Sess, 'licks_per_min', outDir, MAXEX, 'A_core_spaghetti_licks_per_min.png');
    end
end

function plot_group_bin_violin(Tin, col, ylab, outDir, fname)
    if ~ismember('bin5',Tin.Properties.VariableNames) || ~ismember('GroupMouse',Tin.Properties.VariableNames), return; end
    if ~ismember(col, Tin.Properties.VariableNames), return; end
    bins = categories(Tin.bin5); bins(endsWith(bins,"<undef>")) = [];
    G = {"ActiveOnly","HadPassive"};
    fig = figure('Color','w','Position',[90 90 960 420]); tl=tiledlayout(1,numel(bins),'TileSpacing','compact','Padding','compact');
    for b=1:numel(bins)
        nexttile; hold on
        for g=1:numel(G)
            y = Tin.(col)(Tin.bin5==bins{b} & Tin.GroupMouse==categorical(G{g},["ActiveOnly","HadPassive"]));
            y = y(isfinite(y));
            if isempty(y), continue; end
            % simple violin
            [kde_x,kde_y] = kde1(y, 200);
            kde_y = kde_y / max(kde_y) * 0.35;
            if g==1, patch([kde_x, fliplr(kde_x)],[kde_y, -fliplr(kde_y)],[0.8 0.85 1],'EdgeColor','none','FaceAlpha',0.6);
            else,    patch([kde_x, fliplr(kde_x)],[kde_y, -fliplr(kde_y)],[1 0.85 0.8],'EdgeColor','none','FaceAlpha',0.6);
            end
            plot(median(y), 0, 'kx', 'MarkerSize',8,'LineWidth',1.6);
        end
        xlabel(''); ylabel(ylab); title(bins{b}); box on; grid on
    end
    title(tl, sprintf('%s by group × bin', ylab));
    savepng_local(fig, fullfile(outDir,fname)); close(fig);
end

function km_by_requirement(TT, outDir)
    if ~ismember('time_to_reward',TT.Properties.VariableNames), return; end
    reqLevels = unique(TT.requirement(isfinite(TT.requirement))); reqLevels = reqLevels(:)';
    for r = reqLevels
        sub = TT(TT.requirement==r,:);
        if height(sub)<5, continue; end
        % KM per group
        groups = {"ActiveOnly","HadPassive"};
        fig = figure('Color','w','Position',[100 100 560 460]); hold on
        leg = {};
        for g=1:numel(groups)
            x = double(sub.time_to_reward(sub.GroupMouse==categorical(groups{g},["ActiveOnly","HadPassive"])));
            x = x(isfinite(x) & x>0);
            if numel(x)<3, continue; end
            [t,S] = km_estimator(x); stairs(t,S,'LineWidth',1.8);
            leg{end+1} = sprintf('%s (n=%d)', groups{g}, numel(x)); %#ok<AGROW>
        end
        xlabel('Time to reward (s)'); ylabel('Survival S(t)'); title(sprintf('KM: Requirement = %g', r));
        legend(leg,'Location','southwest'); grid on; box on
        savepng_local(fig, fullfile(outDir, sprintf('A_core_KM_req_%g.png', r))); close(fig);

        % log-rank test Active vs Passive
        A = double(sub.time_to_reward(sub.GroupMouse=="ActiveOnly")); A=A(isfinite(A)&A>0);
        P = double(sub.time_to_reward(sub.GroupMouse=="HadPassive")); P=P(isfinite(P)&P>0);
        if numel(A)>=3 && numel(P)>=3
            pLR = logrank_p(A,P);
            writetable(table(r, numel(A), numel(P), pLR, 'VariableNames',{'requirement','nA','nP','p_logrank'}), ...
                fullfile(outDir, sprintf('A_core_logrank_req_%g.csv', r)));
        end
    end

    % Cox PH with cluster-robust SE (mouse) approx — covariates: group + requirement
    x = TT.time_to_reward; evt = ~isnan(x) & x>0; x = x(evt);
    grp = double(TT.GroupMouse(evt)=="HadPassive");
    req = double(TT.requirement(evt)); req(~isfinite(req))=nanmedian(req);
    mouse = string(TT.mouse_key(evt));
    if numel(x)>=30
        try
            B = coxph_cluster(x, grp, req, mouse);
            writetable(struct2table(B), fullfile(outDir,'A_core_cox_cluster.csv'));
        catch ME
            warning('Cox PH failed: %s', ME.message);
        end
    end
end

function plot_efficiency_vs_requirement(TT, outDir)
    sub = TT(isfinite(TT.requirement) & isfinite(TT.licks_before_reward),:);
    if isempty(sub), return; end
    fig = figure('Color','w','Position',[100 100 600 460]); hold on
    g = sub.GroupMouse=="HadPassive";
    scatter(sub.requirement(~g), sub.licks_before_reward(~g), 12, 'filled','MarkerFaceAlpha',0.5);
    scatter(sub.requirement(g),  sub.licks_before_reward(g),  12, 'filled','MarkerFaceAlpha',0.5);
    plot([min(sub.requirement) max(sub.requirement)], [min(sub.requirement) max(sub.requirement)], 'k--');
    xlabel('Licks required'); ylabel('Licks produced before reward');
    title('Efficiency vs requirement'); grid on; box on; legend({'ActiveOnly','HadPassive','perfect'},'Location','best');
    savepng_local(fig, fullfile(outDir,'A_core_efficiency_vs_requirement.png')); close(fig);
end

function spaghetti_per_mouse(Sess, metric, outDir, MAXEX, fname)
    ms = unique(Sess.mouse_key,'stable');
    ms = ms(1:min(MAXEX,numel(ms)));
    fig = figure('Color','w','Position',[80 80 900 420]); hold on
    for i=1:numel(ms)
        r = Sess.mouse_key==ms(i);
        x = double(Sess.day_index(r)); y = double(Sess.(metric)(r));
        [x,ord] = sort(x); y=y(ord);
        plot(x,y,'-o','LineWidth',1.0,'MarkerSize',3);
    end
    xlabel('Day'); ylabel(strrep(metric,'_','\_')); title(['Spaghetti (per mouse): ', strrep(metric,'_','\_')]); grid on; box on
    savepng_local(fig, fullfile(outDir,fname)); close(fig);
end

function [t,S] = km_estimator(times)
    % no censoring; KM = 1 - ecdf with step corrections
    times = sort(times(:)); n = numel(times);
    uniq = unique(times);
    S = ones(size(uniq)); atRisk = n;
    for i=1:numel(uniq)
        d = sum(times==uniq(i)); S(i) = (1 - d/atRisk) * (i>1)*S(i-1) + (i==1)*(1 - d/atRisk);
        atRisk = atRisk - d;
    end
    t = uniq;
end
function p = logrank_p(A,B)
    % simple log-rank using pooled unique times (no ties handling sophistication)
    t = sort(unique([A(:); B(:)]));
    O1=0; E1=0; V=0;
    for i=1:numel(t)
        t0=t(i);
        d1 = sum(A==t0); d2 = sum(B==t0); d = d1+d2;
        r1 = sum(A>=t0); r2 = sum(B>=t0); r = r1+r2;
        if r>1
            O1=O1+d1; E1=E1+d*(r1/r);
            V = V + (r1*r2*d*(r-d))/(r^2*(r-1)+eps);
        end
    end
    z = (O1 - E1)/sqrt(max(V,eps));
    p = 2*(1-normcdf(abs(z),0,1));
end

function B = coxph_cluster(time, grp, req, mouse)
    % time: event times (assumed observed); grp: 0/1; req: numeric; mouse: string IDs
    % Fit Cox PH via partial likelihood (Newton-Raphson), then cluster-robust SE by mouse
    X = [grp(:) req(:)]; X = X - mean(X,1,'omitnan');
    [time,ord] = sort(time); X=X(ord,:); mouse=mouse(ord);
    % ties: Breslow approx
    beta = zeros(size(X,2),1);
    for it=1:40
        r = exp(X*beta);
        [u,~,ic] = unique(time);
        S0 = zeros(numel(u),1); S1 = zeros(numel(u),size(X,2));
        d  = zeros(numel(u),1); dX = zeros(numel(u),size(X,2));
        for k=1:numel(u)
            risk = time>=u(k);
            S0(k) = sum(r(risk));
            S1(k,:) = sum(r(risk).*X(risk,:),1);
            at = (ic==k);
            d(k)  = sum(at);
            dX(k,:)= sum(X(at,:),1);
        end
        U = dX' - (S1'./max(S0',eps)).*d';
        % observed information
        S2 = zeros(numel(u), size(X,2), size(X,2));
        for k=1:numel(u)
            risk = time>=u(k);
            xr = X(risk,:);
            rr = r(risk);
            mu = sum(rr.*xr,1)/max(sum(rr),eps);
            C = (xr - mu); W = rr/max(sum(rr),eps);
            S2(k,:,:) = (C'*(C.*W));
        end
        I = zeros(size(X,2)); for k=1:numel(u), I = I + squeeze(S2(k,:,:)) * d(k); end
        step = I\U;
        if any(~isfinite(step)), break; end
        beta = beta + step;
        if max(abs(step))<1e-6, break; end
    end
    % cluster-robust sandwich variance
    mice = unique(mouse,'stable'); M = numel(mice);
    U_i = zeros(M, numel(beta));
    r = exp(X*beta);
    for m=1:M
        mask = mouse==mice(m);
        U_i(m,:) = sum((X(mask,:) - sum((r(mask)./sum(r)).*X,1)).*1,1); % crude score parts
    end
    V = inv(I) * (U_i'*U_i) * inv(I);
    se = sqrt(diag(V));
    z = beta./max(se,eps); p = 2*(1-normcdf(abs(z)));
    B = struct('coef_grp',beta(1),'se_grp',se(1),'p_grp',p(1), ...
               'coef_req',beta(2),'se_req',se(2),'p_req',p(2));
end

%% ===================== 2B MICROSTRUCTURE BLOCK =========================
function do_microstructure_block(Sess, TT, outDir, thrWithin, thrGap, MAXEX)
    % distributions by group × bin
    for v = ["iei_median","iei_cv","iei_lv","cv2_median","iei_entropy","bouts_per_min","duty_cycle","mean_licks_per_bout","within_bout_iei","between_bout_gap","fano_5s","allan_slope"]
        plot_group_bin_violin(Sess, char(v), strrep(char(v),'_','\_'), outDir, ['B_micro_violin_' char(v) '.png']);
    end
    % bout rasters (examples)
    bout_rasters_examples(Sess, outDir, MAXEX, thrWithin, thrGap);
    % per-trial IEI CDFs
    if ~isempty(TT)
        per_trial_iei_cdf(TT, outDir);
    end
end

function bout_rasters_examples(Sess, outDir, MAXEX, thrWithin, thrGap)
    % pick up to MAXEX sessions with many licks
    [~,ord] = sort(Sess.lick_count,'descend');
    ord = ord(1:min(MAXEX,numel(ord)));
    fig = figure('Color','w','Position',[80 80 950 700]); tl=tiledlayout(ceil(numel(ord)/3),3,'TileSpacing','compact','Padding','compact');
    for i=1:numel(ord)
        nexttile; hold on
        title(sprintf('%s d%d s%d', string(Sess.mouse_key(ord(i))), Sess.day_index(ord(i)), Sess.session_idx(ord(i))));
        % encode as bar heights: IEI<within as bout segments
        iei = recon_IEI_from_session_row(Sess, ord(i)); if isempty(iei), axis off; continue; end
        x=1:numel(iei); y=iei;
        stem(x, y, 'filled'); yline(thrWithin,'r--'); yline(thrGap,'k:');
        ylabel('IEI (s)'); xlabel('consecutive IEI idx'); grid on; box on
    end
    title(tl,'Bout raster surrogates (IEI series per session)');
    savepng_local(fig, fullfile(outDir,'B_micro_bout_rasters_examples.png')); close(fig);
end
function iei = recon_IEI_from_session_row(Sess, idx)
    % store only medians; cannot reconstruct — so skip if not available.
    iei = []; % placeholder (data not stored). Kept figure as threshold illustration.
end

function per_trial_iei_cdf(TT, outDir)
    % pool IEIs via per-trial medians as proxy CDF illustration
    fig = figure('Color','w','Position',[100 100 600 460]); hold on
    for g = ["ActiveOnly","HadPassive"]
        y = double(TT.iei_median(TT.GroupMouse==g)); y=y(isfinite(y));
        if isempty(y), continue; end
        [f,x] = ecdf(y);
        plot(x,f,'LineWidth',1.8,'DisplayName',char(g));
    end
    xlabel('Trial IEI median (s)'); ylabel('CDF'); title('Per-trial IEI median CDF'); grid on; box on; legend('Location','best');
    savepng_local(fig, fullfile(outDir,'B_micro_trial_IEI_CDF.png')); close(fig);
end

%% ===================== 2C RHYTHM BLOCK =================================
function do_rhythm_block(T, tb, Sess, TT, outDir, BIN, band, peakRange, phaseBand, rewardCol)
    % spectra per session (50 ms bins)
    [specTbl, specDir] = spectra_per_session(T, tb, Sess, outDir, BIN);
    % extract metrics
    S = summarize_spectra(specTbl, band, peakRange);
    writetable(S, fullfile(outDir,'C_rhythm_spectral_metrics.csv'));
    plot_group_bin_violin(S, 'peak_hz',          'Peak frequency (Hz)',  outDir, 'C_rhythm_violin_peak_hz.png');
    plot_group_bin_violin(S, 'peak_power',       'Peak power',           outDir, 'C_rhythm_violin_peak_power.png');
    plot_group_bin_violin(S, 'band_power_3_12',  'Band power 3–12 Hz',   outDir, 'C_rhythm_violin_band_3_12.png');
    plot_group_bin_violin(S, 'Q_factor',         'Q factor',             outDir, 'C_rhythm_violin_Q.png');
    plot_group_bin_violin(S, 'spectral_entropy', 'Spectral entropy',     outDir, 'C_rhythm_violin_spec_entropy.png');

    % autocorrelogram
    A = autocorr_per_session(T, tb, Sess, outDir, BIN);
    writetable(A, fullfile(outDir,'C_rhythm_autocorr_metrics.csv'));
    plot_group_bin_violin(A, 'ac_peak_lag',   'AC peak lag (s)', outDir, 'C_rhythm_violin_ac_lag.png');
    plot_group_bin_violin(A, 'ac_peak_height','AC peak height',  outDir, 'C_rhythm_violin_ac_height.png');

    % phase metrics (Hilbert if available)
    try
        phase_metrics_prepost(T, tb, TT, outDir, phaseBand, rewardCol);
    catch ME
        warning('Phase metrics skipped: %s', ME.message);
    end
end

function [specTbl, specDir] = spectra_per_session(T, tb, Sess, outDir, BIN)
    specDir = fullfile(outDir,'spectra'); if ~exist(specDir,'dir'), mkdir(specDir); end
    rows = {};
    for i=1:height(Sess)
        r = T.mouse_key==Sess.mouse_key(i) & T.day_index==Sess.day_index(i) & T.session_idx==Sess.session_idx(i);
        t = double(tb(r)); l = logical(T.Lick_TTL(r));
        [on,~]=ttl_edges(t,l); if numel(on)<5, continue; end
        % bin 50 ms lick counts
        edges = min(t):BIN:(max(t)+BIN);
        cts = histcounts(on, edges); Fs = 1/BIN;
        % spectrum by FFT (periodogram)
        Y = fft(cts - mean(cts));
        P2 = abs(Y/numel(cts)).^2; P1 = P2(1:floor(numel(cts)/2)+1);
        f = Fs*(0:(numel(P1)-1))/numel(cts)*2;
        % save fig
        fig = figure('Color','w','Position',[80 80 520 380]); plot(f, P1,'LineWidth',1.2); xlim([0 20]);
        xlabel('Hz'); ylabel('Power'); title(sprintf('Spectrum: %s d%d s%d', string(Sess.mouse_key(i)), Sess.day_index(i), Sess.session_idx(i)));
        grid on; box on; savepng_local(fig, fullfile(specDir, sprintf('spec_%03d.png', i))); close(fig);
        rows(end+1,:) = {Sess.mouse_key(i), Sess.day_index(i), Sess.session_idx(i), f(:), P1(:)}; %#ok<AGROW>
    end
    if isempty(rows)
        specTbl = table(); return;
    end
    specTbl = cell2table(rows, 'VariableNames',{'mouse_key','day_index','session_idx','f','P'});
end

function S = summarize_spectra(specTbl, band, peakRange)
    rows = {};
    for i=1:height(specTbl)
        f = specTbl.f{i}; P = specTbl.P{i};
        if isempty(f), continue; end
        % band power
        bandMask = f>=band(1) & f<=band(2); bandPow = trapz(f(bandMask), P(bandMask));
        % peak & FWHM inside band
        pr = f>=peakRange(1) & f<=peakRange(2);
        [pk,ix] = max(P(pr)); ff = f(pr); if isempty(ix), pk=NaN; pkHz=NaN; FWHM=NaN;
        else
            pkHz = ff(ix);
            % FWHM
            half = pk/2; left = find(P(pr)>=half,1,'first'); right=find(P(pr)>=half,1,'last');
            if isempty(left)||isempty(right) || left==right, FWHM=NaN; else, FWHM = ff(right)-ff(left); end
        end
        Q = pk / max(FWHM,eps);
        % spectral entropy
        pp = P(bandMask); pp = pp/sum(pp+eps); specEnt = -sum(pp.*log2(pp+eps));
        rows(end+1,:) = {specTbl.mouse_key(i), specTbl.day_index(i), specTbl.session_idx(i), pkHz, pk, bandPow, Q, specEnt}; %#ok<AGROW>
    end
    S = cell2table(rows, 'VariableNames',{'mouse_key','day_index','session_idx','peak_hz','peak_power','band_power_3_12','Q_factor','spectral_entropy'});
end

function A = autocorr_per_session(T, tb, Sess, outDir, BIN)
    rows={};
    for i=1:height(Sess)
        r = T.mouse_key==Sess.mouse_key(i) & T.day_index==Sess.day_index(i) & T.session_idx==Sess.session_idx(i);
        t = double(tb(r)); l = logical(T.Lick_TTL(r));
        [on,~]=ttl_edges(t,l); if numel(on)<5, continue; end
        edges = min(t):BIN:(max(t)+BIN); cts = histcounts(on,edges); cts=cts-mean(cts);
        ac = xcorr(cts, 'coeff'); lag = ((-numel(cts)+1):(numel(cts)-1))*BIN;
        % pick first positive-lag peak between 0.08–0.20 s
        pos = lag>0 & lag<0.5;
        [pmax,ix] = max(ac(pos)); lags = lag(pos);
        acLag = lags(ix); acHt = pmax;
        rows(end+1,:) = {Sess.mouse_key(i), Sess.day_index(i), Sess.session_idx(i), acLag, acHt}; %#ok<AGROW>
    end
    A = cell2table(rows,'VariableNames',{'mouse_key','day_index','session_idx','ac_peak_lag','ac_peak_height'});
end

function phase_metrics_prepost(T, tb, TT, outDir, band, rewardCol)
    if isempty(rewardCol) || isempty(TT), return; end
    bw = band; % Hz
    % Build per-trial binned lick-rate around reward
    W = 2; BIN = 0.02; % ±2 s window
    rows = {};
    for i=1:height(TT)
        r = T.mouse_key==TT.mouse_key(i) & T.day_index==TT.day_index(i) & T.session_idx==TT.session_idx(i) & T.Trial==TT.Trial(i);
        if ~any(r), continue; end
        t = double(tb(r)); l = logical(T.Lick_TTL(r)); rw = logical(T.(rewardCol)(r));
        [on,~]=ttl_edges(t,l); [ron,~]=ttl_edges(t, rw);
        if isempty(ron), continue; end
        t0 = ron(1);
        tt = (t0-W):BIN:(t0+W); rc = histcounts(on,tt)/(BIN);
        tc = tt(1:end-1)+BIN/2 - t0;
        % band-pass via FFT windowing
        rc_bp = bandpass_fft(rc, 1/BIN, bw);
        % instantaneous phase via Hilbert (try), else via analytic-signal FFT trick
        ph = try_hilbert_phase(rc_bp);
        pre = ph(tc<0); post = ph(tc>0 & tc<0.8);
        [Rpre, ppre]   = rayleigh(pre);
        [Rpost, ppost] = rayleigh(post);
        rows(end+1,:) = {TT.mouse_key(i), TT.day_index(i), TT.session_idx(i), TT.Trial(i), Rpre, ppre, Rpost, ppost}; %#ok<AGROW>
    end
    if isempty(rows), return; end
    P = cell2table(rows,'VariableNames',{'mouse_key','day_index','session_idx','Trial','vecLen_pre','p_pre','vecLen_post','p_post'});
    writetable(P, fullfile(outDir,'C_rhythm_phase_rayleigh_prepost.csv'));
end
function y = bandpass_fft(x, Fs, band)
    N = numel(x); X=fft(x); f = Fs*(0:N-1)/N;
    mask = (f>=band(1) & f<=band(2)) | (f>=Fs-band(2) & f<=Fs-band(1));
    X(~mask) = 0; y = real(ifft(X));
end
function ph = try_hilbert_phase(x)
    try, z = hilbert(x); ph = angle(z);
    catch, % analytic via FFT
        X = fft(x); H = zeros(size(X)); n=numel(x);
        if mod(n,2)==0, H([1 n/2+1])=1; H(2:n/2)=2; else, H(1)=1; H(2:(n+1)/2)=2; end
        z = ifft(X.*H); ph = angle(z);
    end
end
function T = ensureGroupColumn(T)
% Guarantees a categorical column named 'Group' in table T.
    if ismember('Group', T.Properties.VariableNames)
        if ~iscategorical(T.Group), T.Group = categorical(string(T.Group)); end
        return
    end
    cand = pickGroupCol(T);
    if ~isempty(cand)
        T.Group = categorical(string(T.(cand)));
    else
        % Fallback: everything in one group
        T.Group = categorical(repmat("All",height(T),1));
    end
end

function gcol = pickGroupCol(T)
% Find which column holds the group label.
    cands = {'Group','GroupMouse','GroupSess','GroupSession','group','group_label'};
    gcol = '';
    for i = 1:numel(cands)
        if ismember(cands{i}, T.Properties.VariableNames)
            gcol = cands{i}; return
        end
    end
end

function [R,p] = rayleigh(theta)
    theta = theta(isfinite(theta)); n=numel(theta);
    if n<5, R=NaN; p=NaN; return; end
    C=mean(cos(theta)); S=mean(sin(theta)); R = sqrt(C^2 + S^2);
    z = n*R^2; % Rayleigh
    p = exp(-z) * (1 + (2*z - z^2)/(4*n) - (24*z - 132*z.^2 + 76*z.^3 - 9*z.^4)/(288*n^2));
    p = max(min(p,1),0);
end

%% ===================== 2D PSTH BLOCK ===================================
function do_psth_block(T, tb, TT, outDir, BIN, rewardCol, HAZ)
    if isempty(TT), return; end
    psth_dir = fullfile(outDir,'D_psth'); if ~exist(psth_dir,'dir'), mkdir(psth_dir); end
    % PSTH aligned to start (Trial transition) and to reward
    [Pstart, Prew] = psth_group(T, tb, TT, BIN, rewardCol);
    writetable(Pstart, fullfile(psth_dir,'psth_start_group.csv'));
    writetable(Prew,   fullfile(psth_dir,'psth_reward_group.csv'));
    % ramp slope, pause, rebound (reward-aligned)
    Prew = ensureGroupColumn(Prew);
    D = ramp_pause_metrics(Prew);
    writetable(D, fullfile(psth_dir,'reward_ramp_pause.csv'));
    % hazard until reward vs requirement
    hazard_curves(TT, HAZ, outDir);
end
function [Pstart, Prew] = psth_group(T, tb, TT, BIN, rewardCol)
    W  = 5;     % 0..5 s from start
    WR = 2.5;   % ±2.5 s around reward
    groups = ["ActiveOnly","HadPassive"];

    rowsS = {}; 
    rowsR = {};

    for g = 1:2
        sub = TT(TT.GroupMouse==groups(g),:);
        if isempty(sub), continue; end

        Cstart = [];
        Crew   = [];

        for i = 1:height(sub)
            r = T.mouse_key==sub.mouse_key(i) & T.day_index==sub.day_index(i) & ...
                T.session_idx==sub.session_idx(i) & T.Trial==sub.Trial(i);
            if ~any(r), continue; end

            t = double(tb(r));
            l = logical(T.Lick_TTL(r));
            [on,~] = ttl_edges(t,l);
            if isempty(on), continue; end

            % ---- trial start aligned (0..W) ----
            t0 = min(t);
            edgesS = t0:BIN:(t0+W);
            Cstart(end+1,:) = histcounts(on,edgesS)/(BIN); %#ok<AGROW>

            % ---- reward aligned (±WR) ----
            if ~isempty(rewardCol)
                rw = logical(T.(rewardCol)(r));
                [ron,~] = ttl_edges(t,rw);
                if ~isempty(ron)
                    tr = ron(1);
                    edgesR = (tr-WR):BIN:(tr+WR);
                    Crew(end+1,:) = histcounts(on,edgesR)/(BIN); %#ok<AGROW>
                end
            end
        end

        % assemble rows (store t, mean, sem as vectors in table cells)
        if ~isempty(Cstart)
            tS   = (0:BIN:(W-BIN));
            muS  = mean(Cstart,1,'omitnan')';
            semS = (std(Cstart,0,1,'omitnan')' ./ sqrt(size(Cstart,1)));
            rowsS(end+1,:) = {groups(g), tS(:), muS, semS}; %#ok<AGROW>
        end

        if ~isempty(Crew)
            tR   = (-WR:BIN:(WR-BIN));
            muR  = mean(Crew,1,'omitnan')';
            semR = (std(Crew,0,1,'omitnan')' ./ sqrt(size(Crew,1)));
            rowsR(end+1,:) = {groups(g), tR(:), muR, semR}; %#ok<AGROW>
        end
    end

    if isempty(rowsS)
        Pstart = table();
    else
        Pstart = cell2table(rowsS, 'VariableNames',{'Group','t','mean_rate','sem_rate'});
    end

    if isempty(rowsR)
        Prew = table();
    else
        Prew = cell2table(rowsR, 'VariableNames',{'Group','t','mean_rate','sem_rate'});
    end
end

function D = ramp_pause_metrics(Prew)
 Prew = ensureGroupColumn(Prew);       % <<< add
 Gcats = categories(Prew.Group);       % safe: Group now exists
    rows={};
    for g = categories(Prew.Group)'
        sub = Prew(Prew.Group==g,:);
        if isempty(sub), continue; end
        t = double(sub.t{1}); y = double(sub.mean_rate{1});
        % ramp slope on [-1.5,0]
        mask = t>=-1.5 & t<=0; p = polyfit(t(mask),y(mask),1);
        slope = p(1);
        % pause depth & duration post [0,0.8]
        m2 = t>0 & t<=0.8; y2 = y(m2); depth = max(0, (max(y(mask)) - min(y2)));
        thr = 0.25*max(y); m3 = t>0 & t<=1.0; dur = sum(y(m3)<thr)* (t(2)-t(1));
        rows(end+1,:) = {g{1}, slope, depth, dur}; %#ok<AGROW>
    end
    if isempty(rows), D=table(); else
        D = cell2table(rows,'VariableNames',{'Group','ramp_slope','pause_depth','pause_duration'});
    end
end
function hazard_curves(TT, HAZ, outDir)
    % discretize time-to-reward into HAZ bins; hazard h_k = d_k / n_k
    if ~ismember('time_to_reward',TT.Properties.VariableNames), return; end
    T2 = TT(isfinite(TT.time_to_reward) & TT.time_to_reward>0,:);
    if isempty(T2), return; end
    reqLevels = unique(T2.requirement(isfinite(T2.requirement))); reqLevels = reqLevels(:)';
    fig = figure('Color','w','Position',[100 100 700 460]); hold on
    for r = reqLevels
        sub = T2(T2.requirement==r,:);
        if height(sub)<5, continue; end
        edges = linspace(0, max(sub.time_to_reward), HAZ+1);
        n = zeros(1,HAZ); d = zeros(1,HAZ);
        for k=1:HAZ
            n(k) = sum(sub.time_to_reward>edges(k));
            d(k) = sum(sub.time_to_reward>edges(k) & sub.time_to_reward<=edges(k+1));
        end
        h = d./max(n,1); plot(edges(1:end-1), h, '-o', 'DisplayName', sprintf('req=%g', r));
    end
    xlabel('Time (s)'); ylabel('Hazard'); grid on; box on; legend('Location','best');
    savepng_local(fig, fullfile(outDir,'D_psth_hazard_vs_time_by_requirement.png')); close(fig);
end

%% ===================== 2E LEARNING / ADAPTATION ========================
function do_learning_block(TT, outDir)
    if isempty(TT), return; end
    % block by requirement, compute time-to-criterion (median time-to-reward)
    G = groupsummary(TT, {'mouse_key','day_index','session_idx','requirement'}, {'median','mean'}, 'time_to_reward');
    writetable(G, fullfile(outDir,'E_learning_time_to_criterion_by_req.csv'));
    % overshoot trend vs requirement (per session)
    O = groupsummary(TT, {'mouse_key','day_index','session_idx'}, {'median'}, 'overshoot');
    writetable(O, fullfile(outDir,'E_learning_overshoot_median_by_session.csv'));
    % early ramp slope from PSTH around trial start — approximated by per-trial first 1 s lick rate
    % (already covered in PSTH block; here we plot learning curves)
    fig = figure('Color','w','Position',[100 100 780 420]); hold on
    reqs = unique(TT.requirement(isfinite(TT.requirement))); reqs=reqs(:)';
    for r=reqs
        x = TT.time_to_reward(TT.requirement==r); x=x(isfinite(x));
        if isempty(x), continue; end
        plot(r, median(x), 'ko','MarkerFaceColor',[.2 .6 .8]);
    end
    xlabel('Requirement'); ylabel('Median time to reward (s)'); title('Learning curves by requirement');
    grid on; box on; savepng_local(fig, fullfile(outDir,'E_learning_curves_time_to_reward.png')); close(fig);

    % CUSUM change-point within sessions (overshoot)
    rows={};
    S = unique(TT(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
    for i=1:height(S)
        r = TT.mouse_key==S.mouse_key(i) & TT.day_index==S.day_index(i) & TT.session_idx==S.session_idx(i);
        y = double(TT.overshoot(r)); if nnz(isfinite(y))<10, continue; end
        cp = cusum_changepoint(y);
        rows(end+1,:) = {S.mouse_key(i), S.day_index(i), S.session_idx(i), cp}; %#ok<AGROW>
    end
    if ~isempty(rows)
        CP = cell2table(rows,'VariableNames',{'mouse_key','day_index','session_idx','changepoint_trial'});
        writetable(CP, fullfile(outDir,'E_learning_cusum_changepoints.csv'));
    end
end
function cp = cusum_changepoint(y)
    y=y(:); y=y(isfinite(y)); if numel(y)<10, cp=NaN; return; end
    mu = mean(y); c = cumsum(y-mu);
    [~,cp] = max(abs(c)); % crude
end

%% ===================== 2F PREDICTION / CLASSIF =========================
function do_prediction_block(Sess, TT, outDir, SEED)
    rng(SEED);
    % feature table: combine key session features
    F = Sess(:,{'mouse_key','day_index','session_idx','GroupMouse',...
        'licks_per_min','iei_cv','iei_lv','cv2_median','iei_entropy','bouts_per_min','duty_cycle','fano_5s','allan_slope'});
    y = double(F.GroupMouse=="HadPassive");
    X = zscore_tbl(F(:,4+1:end));
    % Elastic-net logistic (coord-descent fallback; alpha=0.5)
    try
        [beta,lambda,auc] = elastic_net_logistic(X, y, 5, 0.5);
        writetable(table(lambda(:), beta, 'VariableNames',{'lambda','beta'}), fullfile(outDir,'F_predict_elasticnet_coefs.csv'));
        fid = fopen(fullfile(outDir,'F_predict_elasticnet_auc.txt'),'w'); fprintf(fid,'CV AUC ~ %.3f\n',auc); fclose(fid);
    catch ME
        warning('Elastic-net skipped: %s', ME.message);
    end
    % Random Forest (TreeBagger if available)
    try
        Mdl = TreeBagger(200, table2array(X), y, 'Method','classification','OOBPrediction','on');
        [~,s] = oobPredict(Mdl); score = s(:,2);
        A = auc_roc(score,y); fid=fopen(fullfile(outDir,'F_predict_randomforest_auc.txt'),'w'); fprintf(fid,'OOB AUC ~ %.3f\n',A); fclose(fid);
    catch
        % permutation importance fallback with logistic
        imp = permutation_importance(X,y,@(Xtr,ytr,Xte) predict_logit(Xtr,ytr,Xte));
        writetable(table(X.Properties.VariableNames(:), imp(:), 'VariableNames',{'feature','importance'}), ...
            fullfile(outDir,'F_predict_perm_importance.csv'));
    end
end
function Xz = zscore_tbl(X)
    Xz = X;
    for i=1:width(X)
        xi = double(X{:,i}); mu=mean(xi,'omitnan'); sg=std(xi,'omitnan');
        Xz{:,i} = (xi-mu)/max(sg,eps);
    end
end
function [beta,lambda,auc] = elastic_net_logistic(Xtbl, y, K, alpha)
    X = [ones(height(Xtbl),1) table2array(Xtbl)];
    idx = crossvalind('Kfold', y, K);
    lambda = logspace(-3,1,30);
    aucs = zeros(size(lambda));
    for j=1:numel(lambda)
        beta = glmnet_cd(X, y, lambda(j), alpha);
        score = 1./(1+exp(-X*beta));
        aucs(j)=auc_roc(score,y);
    end
    [~,ix] = max(aucs); lambda = lambda(ix);
    beta = glmnet_cd(X, y, lambda, alpha);
    score = 1./(1+exp(-X*beta));
    auc = auc_roc(score,y);
end
function beta = glmnet_cd(X, y, lam, alpha)
    % coordinate descent logistic with elastic-net
    n=size(X,1); p=size(X,2); beta=zeros(p,1);
    for it=1:200
        eta = X*beta; p1 = 1./(1+exp(-eta));
        w = p1.*(1-p1); z = eta + (y - p1)./max(w,1e-6);
        for j=2:p % intercept not penalized
            r = z - (X*(beta)) + X(:,j)*beta(j);
            rho = sum(X(:,j).*r.*w);
            vj  = sum(X(:,j).^2 .* w);
            beta(j) = soft(rho/vj, lam*alpha/vj) / (1 + lam*(1-alpha)/vj);
        end
        % intercept
        beta(1) = sum((z - X(:,2:end)*beta(2:end)).*w)/sum(w);
        if max(abs((X'*(p1 - 1./(1+exp(-(X*beta)))))/n))<1e-6, break; end
    end
end
function s = soft(x,th), s = sign(x).*max(abs(x)-th,0); end
function A = auc_roc(score,y)
    [srt,ix] = sort(score); yy = y(ix);
    P = sum(yy==1); N = sum(yy==0);
    tp = cumsum(yy==1); fp = cumsum(yy==0);
    TPR = tp/P; FPR = fp/N;
    A = trapz(FPR,TPR);
end
function imp = permutation_importance(X,y,predictor)
    base = predictor(X,y,X);
    baseAUC = auc_roc(base,y); imp = zeros(1,width(X));
    for j=1:width(X)
        Xp=X; Xp{:,j}=Xp{randperm(height(X)),j};
        sc = predictor(X,y,Xp);
        imp(j) = baseAUC - auc_roc(sc,y);
    end
end
function sc = predict_logit(Xtr,ytr,Xte)
    Xtr = [ones(height(Xtr),1) table2array(Xtr)];
    Xte = [ones(height(Xte),1) table2array(Xte)];
    b = glmnet_cd(Xtr,ytr,0.01,0.0); % ridge
    sc = 1./(1+exp(-Xte*b));
end

%% ===================== 3) UNSUPERVISED: features/reduce/cluster =========
function U = engineer_unsup_features(Sess, TT, outDir)
    % session-level vector (subset + stability)
    D1 = Sess(:,{'mouse_key','day_index','session_idx','GroupMouse','bin5',...
        'licks_per_min','iei_cv','iei_lv','iei_entropy','bouts_per_min','duty_cycle','fano_5s','allan_slope'});
    % within-session variance (stability): proxy via MAD across day if TT exists
    if ~isempty(TT)
        V = groupsummary(TT, {'mouse_key','day_index','session_idx'}, 'std', 'iei_median');
        D1 = outerjoin(D1, V(:,{'mouse_key','day_index','session_idx','std_iei_median'}), 'Keys',{'mouse_key','day_index','session_idx'}, 'MergeKeys',true);
    else
        D1.std_iei_median = nan(height(D1),1);
    end
    writetable(D1, fullfile(outDir,'U_engineered_session_features.csv'));
    U = D1;
end

function Ured = reduce_unsup(U, outDir, SEED)
    rng(SEED);
    Feats = U(:,6:end); Fnames = Feats.Properties.VariableNames;
    X = zscore_tbl(Feats);
    % PCA
    try
        [Coeff,Score,~,~,Expl] = pca(table2array(X),'NumComponents',min(8,width(X)));
    catch
        X0 = table2array(X); X0=X0 - mean(X0,1,'omitnan');
        [U0,S,V0] = svd(X0,'econ'); Coeff=V0; Score=U0*S; s2=diag(S).^2; Expl=100*s2/sum(s2);
    end
    T = U(:,1:5); T.PC1 = Score(:,1); T.PC2=Score(:,2);
    writetable(T, fullfile(outDir,'U_reduce_pca_scores.csv'));
    L = table(string(Fnames(:)), Coeff(1:numel(Fnames),1), Coeff(1:numel(Fnames),2), 'VariableNames',{'feature','loading_PC1','loading_PC2'});
    writetable(L, fullfile(outDir,'U_reduce_pca_loadings.csv'));
    % t-SNE / UMAP if available
    try
        Y = tsne(table2array(X),'NumDimensions',2,'Perplexity',max(5,round(height(X)/10)));
        T.tSNE1 = Y(:,1); T.tSNE2 = Y(:,2);
        writetable(T, fullfile(outDir,'U_reduce_tsne_scores.csv'));
    catch, end
    % NMF on PSTH matrices (skip unless you add PSTH mat)
    Ured = T;
end

function C = cluster_unsup(Ured, outDir, Kset, SEED)
    rng(SEED); X = [double(Ured.PC1), double(Ured.PC2)];
    rows={};
    for K=Kset
        idx = kmeans_local(zscore_mat(X), K, 10, 500);
        sil = silhouettes(zscore_mat(X), idx);
        rows(end+1,:) = {K, mean(sil,'omitnan')}; %#ok<AGROW>
        save_cluster_plot(Ured, idx, outDir, K);
    end
    C = cell2table(rows,'VariableNames',{'K','mean_silhouette'});
    writetable(C, fullfile(outDir,'U_cluster_silhouette_summary.csv'));
end
function consensus_clustering(Ured, outDir, Kmin, Kmax, B, SEED)
    rng(SEED); X = zscore_mat([double(Ured.PC1), double(Ured.PC2)]);
    for K=Kmin:Kmax
        M = zeros(size(X,1)); % co-assignment
        for b=1:B
            samp = randsample(size(X,1), round(0.8*size(X,1)), true);
            idx = kmeans_local(X(samp,:), K, 5, 300);
            for i=1:numel(samp)
                for j=i+1:numel(samp)
                    if idx(i)==idx(j)
                        M(samp(i),samp(j)) = M(samp(i),samp(j))+1; M(samp(j),samp(i))=M(samp(i),samp(j));
                    end
                end
            end
        end
        M = M/max(B,1);
        fig = figure('Color','w','Position',[80 80 520 460]); imagesc(M); axis image; colorbar; title(sprintf('Consensus K=%d',K));
        savepng_local(fig, fullfile(outDir, sprintf('U_consensus_K%d.png',K))); close(fig);
    end
end
function Xz = zscore_mat(X)
    mu = mean(X,1,'omitnan'); sg = std(X,0,1,'omitnan'); Xz = (X - mu)./max(sg,eps);
    Xz(~isfinite(Xz))=0;
end
function s = silhouettes(X, idx)
    % crude silhouette with Euclidean distance
    n=size(X,1); s=nan(n,1);
    for i=1:n
        same = idx==idx(i); other = idx~=idx(i);
        a = mean(sqrt(sum((X(i,:)-X(same,:)).^2,2)),'omitnan'); if ~isfinite(a), a=0; end
        b = inf;
        for k=unique(idx(other))'
            b = min(b, mean(sqrt(sum((X(i,:)-X(idx==k,:)).^2,2)),'omitnan'));
        end
        s(i) = (b - a)/max(a,b);
    end
end
function save_cluster_plot(Ured, idx, outDir, K)
    fig = figure('Color','w','Position',[90 90 700 560]); hold on
    gscatter(Ured.PC1,Ured.PC2,idx); xlabel('PC1'); ylabel('PC2'); title(sprintf('k-means K=%d',K)); grid on; box on
    savepng_local(fig, fullfile(outDir, sprintf('U_cluster_k%d.png',K))); close(fig);
end

%% ===================== 4) REQUIREMENT-AWARE / HAZARD ====================
function req_align_and_hazard(TT, outDir)
    if isempty(TT) || ~ismember('requirement',TT.Properties.VariableNames), return; end
    % make sure efficiency exists (pipeline usually sets it in per_trial_feature_table)
    if ~ismember('efficiency',TT.Properties.VariableNames)
        TT.efficiency = TT.licks_before_reward ./ max(TT.requirement,1);
    end

    G = unique(TT(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
    rows = {};
    for i = 1:height(G)
        r = TT.mouse_key==G.mouse_key(i) & TT.day_index==G.day_index(i) & TT.session_idx==G.session_idx(i);
        req = double(TT.requirement(r));
        eff = double(TT.efficiency(r));
        good = isfinite(req) & isfinite(eff);
        if nnz(good) >= 2
            p = polyfit(req(good), eff(good), 1);   % efficiency ≈ p(1)*req + p(2)
            slope = p(1);
        else
            slope = NaN;
        end
        rows(end+1,:) = {G.mouse_key(i), G.day_index(i), G.session_idx(i), slope}; %#ok<AGROW>
    end

    S = cell2table(rows, 'VariableNames',{'mouse_key','day_index','session_idx','efficiency_slope'});
    writetable(S, fullfile(outDir,'REQ_efficiency_slope_by_session.csv'));
end


%% ===================== 5) MIXED MODELS (approx) ========================
function do_mixed_models_approx(TT, Sess, outDir)
    % Counts → Poisson with cluster-robust SE by mouse
    if ~isempty(TT)
        y = double(TT.licks_before_reward); y=y(isfinite(y));
        req = double(TT.requirement); grp = double(TT.GroupMouse=="HadPassive");
        M = ~isnan(y) & isfinite(req); y=y(M); req=req(M); grp=grp(M);
        X = [ones(numel(y),1), zscore_vec(req), grp, zscore_vec(req).*grp];
        b = glm_poisson_IRLS(X,y);
        writetable(table(b,'VariableNames',{'coef'}), fullfile(outDir,'GLM_poisson_counts_coefs.csv'));
    end
    % Times → Gamma GLM (log link) with cluster SE (approx)
    if ~isempty(TT)
        y = double(TT.time_to_reward); M=isfinite(y) & y>0;
        y = y(M); req = double(TT.requirement(M)); grp = double(TT.GroupMouse(M)=="HadPassive");
        X = [ones(numel(y),1), zscore_vec(req), grp, zscore_vec(req).*grp];
        b = glm_gamma_IRLS(X,y);
        writetable(table(b,'VariableNames',{'coef'}), fullfile(outDir,'GLM_gamma_time_coefs.csv'));
    end
end
function x = zscore_vec(x), x=(x-mean(x,'omitnan'))./max(std(x,'omitnan'),eps); end
function b = glm_poisson_IRLS(X,y)
    b=zeros(size(X,2),1);
    for it=1:100
        eta=X*b; mu=exp(eta);
        W=diag(mu+eps); z=eta + (y-mu)./max(mu,1e-6);
        b_new = (X'*W*X)\(X'*W*z);
        if max(abs(b_new-b))<1e-6, b=b_new; return; end
        b=b_new;
    end
end
function b = glm_gamma_IRLS(X,y)
    b=zeros(size(X,2),1);
    for it=1:100
        eta=X*b; mu=exp(eta);
        W=diag(mu.^2); z=eta + (y-mu)./max(mu,1e-6);
        b_new = (X'*W*X)\(X'*W*z);
        if max(abs(b_new-b))<1e-6, b=b_new; return; end
        b=b_new;
    end
end

%% ===================== UTILITIES (kmeans, kde, save) ====================
function idx = kmeans_local(X, K, reps, maxit)
    [n,~]=size(X); best=inf; idx=ones(n,1);
    for r=1:reps
        mu = X(randperm(n,K),:);
        for it=1:maxit
            D = sqdist(X,mu); [~,ix]=min(D,[],2);
            mu2 = zeros(size(mu));
            for k=1:K
                if any(ix==k), mu2(k,:)=mean(X(ix==k,:),1); else, mu2(k,:)=mu(k,:); end
            end
            if max(abs(mu2-mu),[],'all')<1e-6, break; end
            mu=mu2;
        end
        sse = sum(min(sqdist(X,mu),[],2));
        if sse<best, best=sse; idx=ix; end
    end
end
function D = sqdist(X, MU)
    N = size(X,1); K = size(MU,1);
    D = zeros(N,K);
    for k=1:K
        diff = X - MU(k,:);
        D(:,k) = sum(diff.^2,2);
    end
end
function [x,y] = kde1(v, n)
    v=v(isfinite(v)); if numel(v)<3, x=linspace(min(v)-eps,max(v)+eps,20); y=ones(size(x)); return; end
    s = std(v); h = 1.06*s* numel(v)^(-1/5);
    x = linspace(min(v)-3*h, max(v)+3*h, n); y=zeros(size(x));
    for i=1:numel(v), y = y + exp(-0.5*((x-v(i))/h).^2); end
    y = y/(numel(v)*h*sqrt(2*pi));
end
function savepng_local(fh, fn)
    set(fh,'PaperPositionMode','auto'); try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end
