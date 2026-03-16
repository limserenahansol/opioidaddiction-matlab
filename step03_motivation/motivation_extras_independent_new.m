function motivation_extras_independent_new()
% motivation_extras_independent (REVISED FOR NEW COHORT + PERIODS)
%
% One-shot extras:
% (1) Peri-event pupil clusters (reward- or lick-locked), Active vs Passive
% (2) Within-mouse regressions: pupil_delta ~ behavioral predictors (+VIF/cond)
% (3) Lick "burstiness" with learned thresholds from ACTIVE sessions (GMM on log IEI)
% (4) Tail-immersion latency vs TrialsCompleted (per mouse)
% (5) NEW: Period summaries (pre/during/post/withdrawal/reexposure),
%          computed both including vs excluding "switch days".
%
% Outputs land in:   run_*/figs/motivation/

%% ===================== USER CONFIG =====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
% Use cohort map (recommended). If true, Active/Passive comes from map,
% overriding T.isPassive if present.
USE_COHORT_MAP = true;

% Only analyze these days (ignore habituation D1-2)
VALID_DAYS = 3:18;

% "Less reliable" days you may want to exclude for period comparisons
% (you mentioned 4,6,11,14; reexposure switch is day 17)
UNRELIABLE_DAYS = [4 6 11 14 17];

% ---- exact-name exclusions (case-insensitive, trim spaces) ----
% NOTE: these should match canonical mouse keys like "6872 black"
EXCL_PUPIL = keynorm(["6872 black","7597 black","8606 forange"]);  % pupil-only
EXCL_ALL   = keynorm(["6872 black","8606 forange"]);               % all other analyses
%% =======================================================

%% ---- locate latest run_* & read CSV ----
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

% ---- hygiene ----
T = ensureString(T,'mouse_key');
if ~ismember('session_idx',T.Properties.VariableNames), T.session_idx = ones(height(T),1); end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end

% robust day_index fallback
if ~ismember('day_index',T.Properties.VariableNames) || all(~isfinite(T.day_index))
    [G, ~, si] = findgroups(string(T.mouse_key), double(T.session_idx));
    ord = splitapply(@(x) tiedrank(x), si, G);
    T.day_index = ord;
elseif ~isnumeric(T.day_index)
    T.day_index = double(T.day_index);
end

% canonicalize mouse_key (handles underscores / missing space / leading zeros)
T.mouse_key = canonical_mouse_key(T.mouse_key);

% sanity echo
allMice = unique(keynorm(T.mouse_key));
fprintf('Excl (pupil): %s\n', strjoin(intersect(allMice,EXCL_PUPIL), ', '));
fprintf('Excl (other): %s\n', strjoin(intersect(allMice,EXCL_ALL),   ', '));

% passive flag tidy (if present)
if ismember('isPassive',T.Properties.VariableNames) && ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% lick TTL
assert(ismember('Lick_TTL',T.Properties.VariableNames),'Lick_TTL missing.');
T.Lick_TTL(isnan(T.Lick_TTL)) = 0;

% timebase
tb = pickTimebase_local(T);
assert(any(isfinite(tb)),'No usable timebase found (CamTime_rel_s, PlotTime_s_30fps, CamTime_s, or Frame/30).');

% output dir
outDir = fullfile(runDir,'figs','motivation');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- cohort map (your new cohort) ----
cohort = build_new_cohort_map();  % canonical mouse_key inside
% attach cohort columns to T (Group, Sex, PairID) if matched
T = attach_cohort_to_rows(T, cohort);

% Decide Active/Passive row mask source
if USE_COHORT_MAP
    if ~ismember('Group', T.Properties.VariableNames)
        error('USE_COHORT_MAP=true but cohort Group column missing after attach.');
    end
    fprintf('Grouping source: COHORT MAP (override)\n');
else
    fprintf('Grouping source: CSV (isPassive / Session_Paradigm)\n');
end

%% ---- filter to valid days ----
keepDay = isfinite(T.day_index) & ismember(T.day_index, VALID_DAYS);
T = T(keepDay,:);

%% ---- build minimal Trials & Sessions needed by the analyses ----
Trials   = build_trial_table_min(T, tb);           % reward time, vigor, pupil_delta, ITI, etc.
Sessions = build_session_table_min(T, tb, Trials, cohort, USE_COHORT_MAP); % session_min, TrialsCompleted, Group

% apply "other analyses" exclusions
Trials   = Trials(~ismember(keynorm(Trials.mouse_key),   EXCL_ALL), :);
Sessions = Sessions(~ismember(keynorm(Sessions.mouse_key), EXCL_ALL), :);

% attach session-level median Immersion_Latency_s if present
Sessions = attach_session_median(Sessions, T, 'Immersion_Latency_s', 'raw_Immersion_Latency_s');

% persist for reference (optional)
try
    writetable(Trials,   fullfile(outDir,'Trials_for_extras.csv'));
    writetable(Sessions, fullfile(outDir,'Sessions_for_extras.csv'));
catch
end

%% ---- (1) Peri-event pupil clusters (Active vs Passive) ----
if ismember('Diameter_px',T.Properties.VariableNames)
    analyze_peri_event_clusters(T, Trials, tb, outDir, 'reward', 'Passive', EXCL_PUPIL);
    analyze_peri_event_clusters(T, Trials, tb, outDir, 'lick',   'Passive', EXCL_PUPIL);
    analyze_peri_event_clusters(T, Trials, tb, outDir, 'reward', 'Active',  EXCL_PUPIL);
    analyze_peri_event_clusters(T, Trials, tb, outDir, 'lick',   'Active',  EXCL_PUPIL);
else
    warning('No Diameter_px in T; skipping peri-event pupil analyses.');
end

%% ---- (2) Within-mouse regression: pupil vs behavior (+ collinearity) ----
regress_pupil_diameter_per_mouse(Trials, Sessions, outDir);

%% ---- (3) Burstiness (GMM thresholds learned from ACTIVE) ----
[burstTbl, ~] = analyze_lick_bursts_GMM(T, tb, Sessions, outDir, EXCL_ALL, cohort, USE_COHORT_MAP);

% sanity checks of learned thresholds
burst_sanity_check_freqs(T, tb, Sessions, outDir, EXCL_ALL);

% NEW: period-based plots + summaries (incl vs excl switch days)
period_summaries_all(burstTbl, Sessions, outDir, UNRELIABLE_DAYS);

% Also keep legacy-style across-days plot (still useful)
plot_bursts_by_day(burstTbl, Sessions, outDir);

%% ---- (4) Tail immersion latency vs TrialsCompleted ----
plot_immersion_vs_trials(Sessions, outDir);

fprintf('\nDone. Extras written to: %s\n', outDir);
end

%% =================== COHORT MAP (NEW) ===================
function cohort = build_new_cohort_map()
% Returns a table with canonical mouse_key, Sex, Group, PairID

rows = {
    % cage 6100 (female): active = black, passive = orange/red
    "6100 black",  "F", "Active",  1
    "6100 orange", "F", "Passive", 1
    "6100 red",    "F", "Passive", 1

    % cage 0911 (female): pair A: red active vs orange passive
    "0911 red",    "F", "Active",  2
    "0911 orange", "F", "Passive", 2

    % cage 0911 (female): pair B: white active vs black passive
    "0911 white",  "F", "Active",  3
    "0911 black",  "F", "Passive", 3

    % cage 0910 (male): active black vs passive orange/red
    "0910 black",  "M", "Active",  4
    "0910 orange", "M", "Passive", 4
    "0910 red",    "M", "Passive", 4

    % cage 6099 (male): pair A: orange active vs red passive
    "6099 orange", "M", "Active",  5
    "6099 red",    "M", "Passive", 5

    % cage 6099 (male): pair B: black active vs white passive (white died day13)
    "6099 black",  "M", "Active",  6
    "6099 white",  "M", "Passive", 6
};

cohort = cell2table(rows, 'VariableNames', {'mouse_key','Sex','Group','PairID'});
cohort.mouse_key = canonical_mouse_key(cohort.mouse_key);
cohort.Sex   = string(cohort.Sex);
cohort.Group = categorical(string(cohort.Group), ["Active","Passive"]);
cohort.PairID = double(cohort.PairID);
end

function T = attach_cohort_to_rows(T, cohort)
% Adds Group, Sex, PairID columns to T if mouse_key matches cohort
mk = string(T.mouse_key);
T.Group  = categorical(repmat("Unknown",height(T),1), ["Active","Passive","Unknown"]);
T.Sex    = string(repmat("",height(T),1));
T.PairID = nan(height(T),1);

[tf, loc] = ismember(mk, string(cohort.mouse_key));
T.Group(tf)  = categorical(string(cohort.Group(loc(tf))), ["Active","Passive","Unknown"]);
T.Sex(tf)    = string(cohort.Sex(loc(tf)));
T.PairID(tf) = cohort.PairID(loc(tf));
end

%% =================== (1) Peri-event pupil clusters ===================
function analyze_peri_event_clusters(T, Trials, tb, outDir, eventKind, groupName, excludeMice)
% eventKind: 'reward' or 'lick'
% groupName: 'Active' or 'Passive'
if nargin<5, eventKind='reward'; end
if nargin<6, groupName='Passive'; end
if ~ismember('Diameter_px',T.Properties.VariableNames), return; end

win = [-2 6];  baseW = [-0.5 0];  fs = 20;

miceAll = unique(string(T.mouse_key),'stable');
if nargin>=7 && ~isempty(excludeMice)
    miceAll = miceAll(~ismember(keynorm(miceAll), keynorm(excludeMice)));
end

% filter mice by cohort Group if present
if ismember('Group', T.Properties.VariableNames)
    wanted = categorical(groupName, ["Active","Passive","Unknown"]);
    mice = string(unique(T.mouse_key(T.Group==wanted),'stable'));
else
    % fallback: infer from isPassive / Session_Paradigm if no cohort columns
    mice = miceAll;
end

% ---- per-mouse peri-event trace (mean across events) ----
S = table();  % mouse_key, trace
for i = 1:numel(mice)
    mk = mice(i);
    rowT   = string(T.mouse_key)==mk;
    times  = double(tb(rowT));
    pupil  = double(T.Diameter_px(rowT));
    if ~any(isfinite(pupil)), continue; end

    switch lower(eventKind)
        case 'reward'
            rowTr = string(Trials.mouse_key)==mk & isfinite(Trials.t_reward_s);
            ev = Trials.t_reward_s(rowTr);
        otherwise % 'lick'
            L  = logical(T.Lick_TTL(rowT));
            ev = ttl_onsets_local(times, L);
    end
    ev = ev(isfinite(ev)); if isempty(ev), continue; end

    tgrid = win(1):1/fs:win(2);
    M = nan(numel(ev), numel(tgrid));
    for k=1:numel(ev)
        t0 = ev(k);
        mask = times >= t0+win(1) & times <= t0+win(2);
        if nnz(mask) < 5, continue; end
        tseg = times(mask) - t0;   yseg = pupil(mask);
        b = tseg>=baseW(1) & tseg<=baseW(2);
        if any(b), yseg = yseg - mean(yseg(b),'omitnan'); end
        [tuniq, ia] = unique(tseg);
        if numel(tuniq)>=5
            M(k,:) = interp1(tuniq, yseg(ia), tgrid, 'linear', NaN);
        end
    end
    tr = mean(M,1,'omitnan');
    if any(isfinite(tr)), S(end+1,:) = table(mk, {tr}, 'VariableNames',{'mouse_key','trace'}); end %#ok<AGROW>
end

if height(S) < 3
    warning('Not enough mice for %s clustering (%s).', groupName, eventKind);
    return
end

% ---- cluster z-scored waveforms ----
X = cat(1, S.trace{:});
X = fillmissing(X,'linear',2,'EndValues','nearest');
X = zscore(X,0,2);
goodRows = all(isfinite(X),2);  S = S(goodRows,:);  X = X(goodRows,:);
bestK=2; bestSil=-Inf;
for K=2:4
    idxCand = kmeans(X,K,'Replicates',20,'MaxIter',500,'Distance','sqeuclidean','Display','off');
    s = mean(silhouette(X,idxCand));
    if s>bestSil, bestSil=s; bestK=K; end
end
idx = kmeans(X,bestK,'Replicates',50,'MaxIter',1000,'Distance','sqeuclidean','Display','off');

% ---- plot cluster means ± SEM ----
tgrid = linspace(win(1),win(2),size(X,2)).';
fig = figure('Color','w','Position',[60 60 1100 380]);
tiledlayout(1,bestK,'TileSpacing','compact');
for k=1:bestK
    nexttile; hold on
    Y  = X(idx==k,:);     nK = size(Y,1);
    mu = mean(Y,1,'omitnan')';
    se = (nK<2) * zeros(size(mu)) + (nK>=2) * (std(Y,0,1,'omitnan')'/sqrt(nK));
    fill([tgrid; flipud(tgrid)],[mu-se; flipud(mu+se)],[0.85 0.9 1],'EdgeColor','none');
    plot(tgrid,mu,'b','LineWidth',2);
    yline(0,'k-'); xline(0,'k:');
    title(sprintf('Cluster %d (n=%d)',k,nK));
    xlabel('Time from event (s)'); ylabel('\Delta pupil (z)'); box on; grid on
end
sgtitle(sprintf('Peri-event pupil clusters — %s (%s)', eventKind, groupName));
fn = fullfile(outDir, sprintf('pupil_clusters_%s_%s.png', lower(groupName), eventKind));
savepng_local(fig, fn); close(fig);

% membership CSV
Tmem = table(string(S.mouse_key), idx(:), 'VariableNames',{'mouse_key','cluster'});
writetable(Tmem, fullfile(outDir, sprintf('pupil_clusters_%s_%s_membership.csv',lower(groupName),eventKind)));
end

%% =================== (2) Within-mouse pupil regressions ===================
function regress_pupil_diameter_per_mouse(Trials, Sessions, outDir)
Yname = 'pupil_delta';
Xnames = {'req','vigor_hz','t_trial_s','overshoot','ITI_to_next_s'};

keep = isfinite(Trials.(Yname));
for k=1:numel(Xnames), keep = keep & isfinite(Trials.(Xnames{k})); end
D = Trials(keep, [{'mouse_key','day_index','session_idx',Yname} Xnames]);

mice = unique(D.mouse_key,'stable');
B = nan(numel(mice), numel(Xnames));
VIFmat = nan(numel(mice), numel(Xnames));
condX  = nan(numel(mice),1);

for i=1:numel(mice)
    rows = D.mouse_key==mice(i);
    if nnz(rows) < 30, continue; end
    y = double(D.(Yname)(rows));
    X = zeros(nnz(rows), numel(Xnames));
    for k=1:numel(Xnames), X(:,k) = double(D.(Xnames{k})(rows)); end
    % z-score within mouse
    X = zscore(X); y = zscore(y);
    X = [ones(size(X,1),1) X];
    try
        beta = X\y; B(i,:) = beta(2:end)';  % drop intercept
        % collinearity checks on predictors
        R = corr(X(:,2:end), 'Rows','pairwise');
        if any(~isfinite(R(:)))
            R = corrcoef(fillmissing(X(:,2:end),'linear'));
        end
        VIFmat(i,:) = diag(inv(R));
        condX(i)    = cond(R);
    catch
    end
end

BetaTbl = table(string(mice), B, condX, 'VariableNames',{'mouse_key','betas','cond_corr'});
for k=1:numel(Xnames)
    BetaTbl.(sprintf('beta_%s',Xnames{k})) = B(:,k);
    BetaTbl.(sprintf('VIF_%s',Xnames{k}))  = VIFmat(:,k);
end
writetable(BetaTbl, fullfile(outDir,'reg_pupil_betas_per_mouse.csv'));

% ---- clustered heatmap (robust) ----
rows = all(isfinite(B),2);
nR = sum(rows);
if nR>=2
    XH = B(rows,:);
    m = mean(XH,2); s = std(XH,0,2); s(s==0) = eps;
    XH = (XH - m) ./ s;

    try
        Z   = linkage(XH,'average','correlation');
        ord = optimalleaforder(Z, pdist(XH,'correlation'));
    catch
        ord = 1:size(XH,1);
    end

    fig = figure('Color','w','Position',[60 60 880 520]);
    imagesc(XH(ord,:), [-2 2]); axis tight; colorbar; colormap(parula);
    yt   = 1:nR;
    ylab = string(mice(rows));
    ylab = ylab(ord);
    set(gca,'YTick',yt,'YTickLabel',cellstr(ylab));
    set(gca,'XTick',1:numel(Xnames),'XTickLabel',Xnames,'XTickLabelRotation',45);
    title('Per-mouse betas: pupil\_delta ~ predictors (z-scored)');
    savepng_local(fig, fullfile(outDir,'reg_pupil_betas_heatmap.png')); close(fig);
elseif nR==1
    XH = B(rows,:);
    fig = figure('Color','w','Position',[400 200 420 240]);
    imagesc(XH, [-2 2]); axis tight; colorbar; colormap(parula);
    set(gca,'YTick',1,'YTickLabel',cellstr(string(mice(rows))));
    set(gca,'XTick',1:numel(Xnames),'XTickLabel',Xnames,'XTickLabelRotation',45);
    title('Per-mouse betas (only 1 mouse with valid betas)');
    savepng_local(fig, fullfile(outDir,'reg_pupil_betas_heatmap.png')); close(fig);
end

bad = any(VIFmat>10,2) | condX>30;
if any(bad)
    writetable(BetaTbl(bad,:), fullfile(outDir,'reg_pupil_high_collinearity.csv'));
end
end

%% =================== (3) Burstiness with GMM-learned thresholds ===================
function [burstTbl, params] = analyze_lick_bursts_GMM(T, tb, Sessions, outDir, excludeMice, cohort, useCohort)
% Learn IEI thresholds from ACTIVE sessions; then compute burst metrics per session.

% exclusions
allow = true(height(T),1);
if nargin>=5 && ~isempty(excludeMice)
    allow = ~ismember(keynorm(T.mouse_key), keynorm(excludeMice));
end

% ---- Active-only rows for learning ----
isActiveRow = allow;
if useCohort && ismember('Group', T.Properties.VariableNames)
    isActiveRow = isActiveRow & (T.Group==categorical("Active",["Active","Passive","Unknown"]));
else
    if ismember('isPassive',T.Properties.VariableNames)
        isActiveRow = isActiveRow & (T.isPassive==0 | isnan(T.isPassive));
    elseif ismember('Session_Paradigm',T.Properties.VariableNames)
        p = T.Session_Paradigm; if iscategorical(p), p = string(p); end
        isActiveRow = isActiveRow & ~contains(lower(string(p)),'passive');
    end
end

[Gsess, ~] = findgroups(string(T.mouse_key(isActiveRow)), ...
                        double(T.day_index(isActiveRow)), ...
                        double(T.session_idx(isActiveRow)));
IEI = [];
tbA = double(tb(isActiveRow));
Lck = logical(T.Lick_TTL(isActiveRow));
for g=1:max(Gsess)
    r = find(Gsess==g);
    on = ttl_onsets_local(tbA(r), Lck(r));
    if numel(on)>=3
        d = diff(on);
        d  = d(d>=0.02 & d<=10);
        IEI = [IEI; d(:)]; %#ok<AGROW>
    end
end
assert(numel(IEI)>50,'Too few Active IEIs to learn thresholds.');

% ---- GMM on log(IEI), pick K=2..3 by BIC ----
x  = log(IEI);
opts = statset('MaxIter',1000);
GM = cell(1,3); BIC = inf(1,3);
for K=2:3
    try
        GM{K} = fitgmdist(x, K, 'Replicates',12, 'Options',opts, 'RegularizationValue',1e-6);
        BIC(K) = GM{K}.BIC;
    catch
    end
end
[~,Kbest] = min(BIC(2:3)); Kbest = Kbest+1;
g = GM{Kbest};

muS  = g.mu(:)';                      % 1×K
sigS = sqrt(squeeze(g.Sigma))';       % 1×K
wS   = g.ComponentProportion(:)';     % 1×K
[muS, ord] = sort(muS); sigS = sigS(ord); wS = wS(ord);

thr_intra_log = fzero(@(z) wS(1)*normpdf(z,muS(1),sigS(1)) - ...
                           wS(2)*normpdf(z,muS(2),sigS(2)), mean(muS(1:2)));
if Kbest==3 && wS(3)>0.03 && (muS(3)-muS(2))>0.2
    thr_gap_log = fzero(@(z) wS(2)*normpdf(z,muS(2),sigS(2)) - ...
                             wS(3)*normpdf(z,muS(3),sigS(3)), mean(muS(2:3)));
else
    thr_gap_log = muS(end) + 1.645*sigS(end);
end
thr_intra = exp(thr_intra_log);
thr_gap   = exp(thr_gap_log);
if thr_gap < thr_intra, thr_gap = max(thr_intra*1.5, thr_intra+0.1); end

% ---- diagnostics plot ----
fig = figure('Color','w','Position',[80 80 840 420]); hold on
edges = [0:0.02:1.5, 2:0.1:10];
histogram(IEI, edges, 'Normalization','pdf', 'FaceColor',[.85 .85 .9], 'EdgeColor','none');
zcol   = linspace(min(x), max(x), 400)';                     % 400×1
Zcomp  = normpdf(zcol, muS, sigS);                           % 400×K
pdfLog = Zcomp * wS(:);                                      % 400×1
pdfIEI = pdfLog ./ exp(zcol);                                % change-of-variables
plot(exp(zcol), pdfIEI, 'k-', 'LineWidth', 1.8);
xline(thr_intra, 'b--', 'LineWidth', 1.5);
xline(thr_gap,   'r--', 'LineWidth', 1.5);
set(gca,'XScale','log'); xlabel('IEI (s, log scale)'); ylabel('PDF');
title(sprintf('Learned thresholds (ACTIVE): intra=%.3fs, gap=%.3fs, K=%d',thr_intra,thr_gap,Kbest));
legend({'IEI histogram','GMM pdf','thr_{intra}','thr_{gap}'}, 'Location','northwest');
savepng_local(fig, fullfile(outDir,'IEI_GMM_active.png')); close(fig);

% ---- per-session burst metrics using learned thresholds ----
burstTbl = Sessions(:,{'mouse_key','day_index','session_idx','Group','Sex','PairID'});
[burstTbl.bursts_per_min, burstTbl.licks_per_burst, ...
 burstTbl.rate_within_burst_hz, burstTbl.frac_licks_in_bursts, ...
 burstTbl.n_bursts] = deal(nan(height(Sessions),1));

for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key),string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i) & ...
        allow;
    t  = double(tb(r)); L = logical(T.Lick_TTL(r));
    on = ttl_onsets_local(t,L);
    if numel(on)<3 || ~isfinite(Sessions.session_min(i)), continue; end

    [bID, keep] = segment_bursts_with_thresholds(on, thr_intra, thr_gap);
    keep = keep & bID > 0;

    K = max(bID(keep));
    burstTbl.n_bursts(i) = K;
    if ~isfinite(K) || K < 1, continue; end

    licksPer = accumarray(bID(keep), 1, [K 1]);

    rates = nan(K,1);
    for k = 1:K
        li = find(bID==k & keep);
        if numel(li) >= 2
            dur = on(li(end)) - on(li(1));
            if dur > 0
                rates(k) = (numel(li)-1)/dur;
            end
        end
    end

    burstTbl.bursts_per_min(i)       = K / Sessions.session_min(i);
    burstTbl.licks_per_burst(i)      = mean(licksPer,'omitnan');
    burstTbl.rate_within_burst_hz(i) = mean(rates,'omitnan');
    burstTbl.frac_licks_in_bursts(i) = nansum(licksPer) / numel(on);
end

writetable(burstTbl, fullfile(outDir,'burst_metrics_per_session.csv'));
params = struct('thr_intra_s',thr_intra,'thr_gap_s',thr_gap,'Kbest',Kbest, ...
                'mu_log',muS,'sigma_log',sigS,'w',wS, ...
                'plot',fullfile(outDir,'IEI_GMM_active.png'));
save(fullfile(outDir,'burst_thresholds_active.mat'),'-struct','params');
end

function [bID, keep] = segment_bursts_with_thresholds(on, thr_intra, thr_gap)
iei = diff(on); N = numel(on);
bID = zeros(N,1); keep = false(N,1);
if N<2, return; end
k = 0; i = 1; keep(1)=true;
while i < N
    if iei(i) > thr_gap
        k = k + 1;
        bID(i+1) = k; keep(i+1)=true; i=i+1; continue
    end
    if k==0
        k=1; bID(i)=k; bID(i+1)=k; keep([i i+1])=true; i=i+1; continue
    end
    if iei(i) <= thr_intra
        bID(i+1) = k; keep(i+1)=true; i=i+1;
    else
        k = k + 1;
        bID(i+1) = k; keep(i+1)=true; i=i+1;
    end
end
for kk=1:max(bID)
    if nnz(bID==kk) < 2, keep(bID==kk)=false; end
end
end

%% =================== (4) Immersion latency vs TrialsCompleted ===================
function plot_immersion_vs_trials(S, outDir)
if ~ismember('raw_Immersion_Latency_s', S.Properties.VariableNames)
    warning('No raw_Immersion_Latency_s attached to Sessions.'); return;
end
mice = unique(S.mouse_key,'stable');
R = table(string(mice), nan(numel(mice),1), 'VariableNames',{'mouse_key','rho_within_mouse'});
for i=1:numel(mice)
    r = S.mouse_key==mice(i) & isfinite(S.day_index);
    if nnz(r)<3, continue; end
    d = S.day_index(r);
    y1= S.raw_Immersion_Latency_s(r);
    y2= S.TrialsCompleted(r);

    ok = isfinite(y1) & isfinite(y2);
    if nnz(ok)>=3, R.rho_within_mouse(i) = corr(y1(ok),y2(ok),'Type','Spearman'); end

    fig = figure('Color','w','Position',[80 80 700 420]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile; plot(d,y1,'-o','LineWidth',1.4); ylabel('Immersion latency (s)'); grid on; title(sprintf('Mouse %s',string(mice(i))));
    nexttile; plot(d,y2,'-o','LineWidth',1.4); ylabel('Trials completed'); xlabel('Day'); grid on;
    savepng_local(fig, fullfile(outDir, sprintf('immersion_vs_trials_%s.png', string(mice(i))))); close(fig);
end
writetable(R, fullfile(outDir,'immersion_trials_correlations.csv'));
end

%% =================== NEW: PERIOD SUMMARIES ===================
function period_summaries_all(burstTbl, Sessions, outDir, unreliableDays)
% Runs summaries twice:
%   (A) including switch/unreliable days
%   (B) excluding them
subA = fullfile(outDir,'period_summary_includingSwitchDays');
subB = fullfile(outDir,'period_summary_excludingSwitchDays');
if ~exist(subA,'dir'), mkdir(subA); end
if ~exist(subB,'dir'), mkdir(subB); end

period_summaries_one(burstTbl, Sessions, subA, false, unreliableDays);
period_summaries_one(burstTbl, Sessions, subB, true,  unreliableDays);
end

function period_summaries_one(burstTbl, Sessions, outDir, excludeUnreliable, unreliableDays)
% Attach period label
S = Sessions;
S.period = period_label(S.day_index, excludeUnreliable, unreliableDays);

% attach burst metrics by join
B = burstTbl;
B.period = period_label(B.day_index, excludeUnreliable, unreliableDays);

% drop excluded
S = S(S.period~="<exclude>",:);
B = B(B.period~="<exclude>",:);

% --------- Session metrics summary by period × group ---------
sessMetrics = {};
if ismember('TrialsCompleted', S.Properties.VariableNames), sessMetrics{end+1} = 'TrialsCompleted'; end %#ok<AGROW>
if ismember('RequirementLast',  S.Properties.VariableNames), sessMetrics{end+1} = 'RequirementLast';  end %#ok<AGROW>
if ismember('session_min',      S.Properties.VariableNames), sessMetrics{end+1} = 'session_min';      end %#ok<AGROW>
if ismember('raw_Immersion_Latency_s', S.Properties.VariableNames), sessMetrics{end+1} = 'raw_Immersion_Latency_s'; end %#ok<AGROW>

sessSummary = summarize_table_by_period_group(S, sessMetrics, 'Sessions');
writetable(sessSummary, fullfile(outDir,'session_metrics_by_period_group.csv'));

plot_metric_panels_by_period_group(S, sessMetrics, 'Session', outDir);

% --------- Burst metrics summary by period × group ---------
burstMetrics = {'bursts_per_min','frac_licks_in_bursts','licks_per_burst','rate_within_burst_hz','n_bursts'};
burstSummary = summarize_table_by_period_group(B, burstMetrics, 'Burst');
writetable(burstSummary, fullfile(outDir,'burst_metrics_by_period_group.csv'));

plot_metric_panels_by_period_group(B, burstMetrics, 'Burst', outDir);
end

function Tsum = summarize_table_by_period_group(T, metrics, label)
% T must contain: period (string/categorical), Group (categorical)
Tsum = table();
if ~ismember('Group', T.Properties.VariableNames)
    return;
end
periods = unique(string(T.period),'stable');
groups  = categories(T.Group);
for mi = 1:numel(metrics)
    m = metrics{mi};
    if ~ismember(m, T.Properties.VariableNames), continue; end
    for pi = 1:numel(periods)
        for gi = 1:numel(groups)
            rows = string(T.period)==periods(pi) & T.Group==groups{gi};
            x = double(T.(m)(rows));
            x = x(isfinite(x));
            if isempty(x)
                mu=NaN; se=NaN; n=0;
            else
                if strcmp(m,'n_bursts')
                    mu = nansum(x); se = NaN; n = numel(x);
                else
                    mu = mean(x,'omitnan'); se = std(x,'omitnan')/sqrt(max(1,numel(x))); n = numel(x);
                end
            end
            Tsum = [Tsum; table(string(label), string(m), string(periods(pi)), string(groups{gi}), mu, se, n, ...
                'VariableNames',{'table','metric','period','group','value','sem','N'})]; %#ok<AGROW>
        end
    end
end
end

function plot_metric_panels_by_period_group(T, metrics, label, outDir)
if ~ismember('Group', T.Properties.VariableNames), return; end
periods = unique(string(T.period),'stable');
groups  = categories(T.Group);

for mi = 1:numel(metrics)
    m = metrics{mi};
    if ~ismember(m, T.Properties.VariableNames), continue; end

    fig = figure('Color','w','Position',[80 80 1200 360]);
    tl = tiledlayout(1, numel(periods), 'TileSpacing','compact','Padding','compact');
    for pi = 1:numel(periods)
        nexttile; hold on
        for gi = 1:numel(groups)
            rows = string(T.period)==periods(pi) & T.Group==groups{gi};
            x = double(T.(m)(rows));
            x = x(isfinite(x));
            if isempty(x), continue; end
            try
                swarmchart(gi*ones(size(x)), x, 18, 'filled','MarkerFaceAlpha',.65);
            catch
                jitter = @(n) (rand(n,1)-.5)*.25;
                scatter(gi+jitter(numel(x)), x, 18, 'filled');
            end
            if strcmp(m,'n_bursts')
                % totals are not a mean metric; still show mean marker for visual reference
                errorbar(gi, mean(x,'omitnan'), std(x,'omitnan')/sqrt(max(1,numel(x))), 'k_', 'LineWidth',1.1);
            else
                errorbar(gi, mean(x,'omitnan'), std(x,'omitnan')/sqrt(max(1,numel(x))), 'k_', 'LineWidth',1.1);
            end
        end
        set(gca,'XTick',1:numel(groups),'XTickLabel',groups);
        title(periods(pi));
        ylabel(m); grid on; box on
    end
    title(tl, sprintf('%s metric by period × group — %s', label, m));
    savepng_local(fig, fullfile(outDir, sprintf('%s_metric_byPeriod_%s.png', lower(label), m))); close(fig);
end
end

function p = period_label(day_index, excludeUnreliable, unreliableDays)
% Returns string labels for your timeline:
% pre: 3-5, during: 6-10, post: 11-13, withdrawal: 14-16, reexposure: 17-18
d = double(day_index);
p = repmat("<exclude>", size(d));

% by default exclude days outside 3-18
inRange = isfinite(d) & d>=3 & d<=18;
p(~inRange) = "<exclude>";

% unreliable day exclusion (optional)
if excludeUnreliable && ~isempty(unreliableDays)
    bad = ismember(d, unreliableDays);
else
    bad = false(size(d));
end

% assign labels where not bad
ok = inRange & ~bad;
p(ok & d>=3  & d<=5 )  = "pre_D3-5";
p(ok & d>=6  & d<=10)  = "during_D6-10";
p(ok & d>=11 & d<=13)  = "post_D11-13";
p(ok & d>=14 & d<=16)  = "withdrawal_D14-16";
p(ok & d>=17 & d<=18)  = "reexposure_D17-18";
end

%% =================== Minimal builders & helpers ===================
function Trials = build_trial_table_min(T, tb)
assert(ismember('Trial',T.Properties.VariableNames),'Trial missing');
assert(ismember('TrialRequirement',T.Properties.VariableNames),'TrialRequirement missing');

tt = double(tb(:));
trial = double(T.Trial(:));
req   = double(T.TrialRequirement(:));
lick  = logical(T.Lick_TTL(:));
mouse = string(T.mouse_key(:));
day   = double(T.day_index(:));
sess  = double(T.session_idx(:));
hasPupil = ismember('Diameter_px',T.Properties.VariableNames);
if hasPupil, pupil = double(T.Diameter_px(:)); else, pupil = nan(height(T),1); end

good = isfinite(tt) & isfinite(trial) & isfinite(req) & ~isnan(lick);
tt=tt(good); trial=trial(good); req=req(good); lick=lick(good);
mouse=mouse(good); day=day(good); sess=sess(good); pupil=pupil(good);

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
    r  = nanmode(ri);
    r_int = NaN; if isfinite(r) && r>=1, r_int = round(r); end

    t  = tt(idx);
    y  = lick(idx);
    on = ttl_onsets_local(t,y); on = unique(on,'stable');

    t_start = min(t);
    lat1 = NaN; n_before = NaN; t_reward = NaN; t_to_reward = NaN; rate_to_reward = NaN; medIEI = NaN;

    if ~isempty(on)
        lat1 = max(0, on(1)-t_start);
        if isfinite(r_int) && numel(on) >= r_int
            t_reward    = on(r_int);
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

    rows(end+1,:) = {mk, di, si, tr, r, n_before, overshoot, efficiency, ...
                     NaN, rate_to_reward, lat1, medIEI, vigor, t_reward, t_start, preP, inP, dP}; %#ok<AGROW>
end

Trials = cell2table(rows, 'VariableNames', ...
    {'mouse_key','day_index','session_idx','trial','req','n_licks_before_reward','overshoot','efficiency', ...
     't_trial_s','rate_in_trial_hz','latency_firstlick_s','iei_median_trial_s','vigor_hz','t_reward_s','t_start_s', ...
     'pupil_pre','pupil_in','pupil_delta'});

% Fill true trial durations from next-start / session end
Trials = sortrows(Trials, {'mouse_key','day_index','session_idx','trial','t_start_s'});
[Gsess, ~] = findgroups(Trials.mouse_key, Trials.day_index, Trials.session_idx);
SessEnd = splitapply(@(ts,tr) max([ts; tr],[],'omitnan'), Trials.t_start_s, Trials.t_reward_s, Gsess);
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

% ITI to next trial
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

function Sessions = build_session_table_min(T, tb, Trials, cohort, useCohort)
keyS = unique(Trials(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
Sessions = keyS;

% Attach cohort labels at mouse-level
Sessions.Group  = categorical(repmat("Unknown",height(Sessions),1), ["Active","Passive","Unknown"]);
Sessions.Sex    = string(repmat("",height(Sessions),1));
Sessions.PairID = nan(height(Sessions),1);

if useCohort
    [tf, loc] = ismember(string(Sessions.mouse_key), string(cohort.mouse_key));
    Sessions.Group(tf)  = categorical(string(cohort.Group(loc(tf))), ["Active","Passive","Unknown"]);
    Sessions.Sex(tf)    = string(cohort.Sex(loc(tf)));
    Sessions.PairID(tf) = cohort.PairID(loc(tf));
else
    % fallback if you really want to use CSV labeling
    if ismember('isPassive',T.Properties.VariableNames)
        % mouse-level: if mouse ever isPassive==1, call it Passive else Active
        mice = unique(string(T.mouse_key),'stable');
        mp = containers.Map(mice, repmat("Active",numel(mice),1));
        for i=1:numel(mice)
            r = string(T.mouse_key)==mice(i);
            ip = double(T.isPassive(r));
            if any(ip==1,'all'), mp(mice(i))="Passive"; end
        end
        for i=1:height(Sessions)
            mk = string(Sessions.mouse_key(i));
            if isKey(mp,mk)
                Sessions.Group(i) = categorical(string(mp(mk)), ["Active","Passive","Unknown"]);
            end
        end
    end
end

% RequirementLast if present
Sessions.RequirementLast = nan(height(Sessions),1);
if ismember('RequirementLast', T.Properties.VariableNames)
    for i=1:height(Sessions)
        r = strcmp(string(T.mouse_key), string(Sessions.mouse_key(i))) & ...
            T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
        v = double(T.RequirementLast(r)); v = v(isfinite(v));
        if ~isempty(v), Sessions.RequirementLast(i) = max(v); end
    end
end

% session duration
Sessions.session_min = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key),string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
    times = double(tb(r)); times = times(isfinite(times));
    if ~isempty(times), Sessions.session_min(i) = (max(times)-min(times))/60; end
end

% trials completed
Sessions.TrialsCompleted = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = Trials.mouse_key==Sessions.mouse_key(i) & ...
        Trials.day_index==Sessions.day_index(i) & ...
        Trials.session_idx==Sessions.session_idx(i);
    tr = Trials(r,:);
    Sessions.TrialsCompleted(i) = nnz(isfinite(tr.t_reward_s));
end
end

function Sessions = attach_session_median(Sessions, T, colName, outName)
if nargin<4, outName = ['raw_' regexprep(colName,'\W','_')]; end
if ~ismember(colName, T.Properties.VariableNames), return; end
vals = nan(height(Sessions),1);
for i=1:height(Sessions)
    r = strcmp(string(T.mouse_key),string(Sessions.mouse_key(i))) & ...
        T.day_index==Sessions.day_index(i) & T.session_idx==Sessions.session_idx(i);
    x = double(T.(colName)(r)); x = x(isfinite(x));
    if ~isempty(x), vals(i) = median(x); end
end
Sessions.(outName) = vals;
end

%% ---------- plotting helper: bursts across days ----------
function plot_bursts_by_day(burstTbl, Sessions, outDir)
vars   = {'bursts_per_min','frac_licks_in_bursts','licks_per_burst','rate_within_burst_hz'};
labels = {'Bursts/min','Fraction of licks in bursts','Licks/burst','Within-burst rate (Hz)'};

% Ensure group present
if ~ismember('Group', burstTbl.Properties.VariableNames)
    burstTbl = outerjoin(burstTbl, Sessions(:,{'mouse_key','day_index','session_idx','Group'}), ...
                         'Keys',{'mouse_key','day_index','session_idx'}, 'MergeKeys',true);
end
G = categories(burstTbl.Group);

for v = 1:numel(vars)
    fig = figure('Color','w','Position',[100 100 1000 360]); hold on
    for gi = 1:numel(G)
        rows = burstTbl.Group==G{gi} & isfinite(burstTbl.(vars{v})) & isfinite(burstTbl.day_index);
        d    = burstTbl.day_index(rows);
        x    = burstTbl.(vars{v})(rows);
        [ud,~,id] = unique(d);
        mu = accumarray(id, x, [], @mean, NaN);
        plot(ud, mu, '-o','LineWidth',1.6);
        scatter(d + (gi-1)*0.04, x, 14, 'filled', 'MarkerFaceAlpha', 0.3);
    end
    xlabel('Day'); ylabel(labels{v});
    legend(G,'Location','northwest'); grid on; box on
    title(sprintf('Burst metric across days — %s', labels{v}));
    savepng_local(fig, fullfile(outDir, sprintf('burst_byDay_%s.png', vars{v}))); close(fig);
end
end

%% ---------- misc utilities ----------
function savepng_local(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end

function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
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

function on = ttl_onsets_local(t, ttl)
ttl = logical(ttl(:)); t = double(t(:));
good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);
if numel(t) < 2, on = []; return; end
d = diff([false; ttl; false]); on = t(d(1:end-1)==1);
end

function m = nanmode(v)
v = v(isfinite(v));
if isempty(v), m = NaN; else, m = mode(v); end
end

function s = keynorm(x)
s = lower(strtrim(string(x)));
end

function mk = canonical_mouse_key(x)
% Robustly canonicalize mouse keys to: "#### color"
% Handles:
%   "0910 red", "910 red", "0910red", "6099 red_m_p", "6100orange f_p"
% Strategy:
%   - find first number block (3-5 digits) => cage, pad to 4
%   - find first color token after that (letters only)
%   - output "#### color" (lower)
sx = lower(strtrim(string(x)));
mk = strings(size(sx));
for i=1:numel(sx)
    s = sx(i);
    if strlength(s)==0
        mk(i) = "";
        continue
    end
    % remove common separators
    s2 = replace(s, ["_", "-", ",", ";"], " ");
    s2 = regexprep(s2, '\s+', ' ');

    % cage digits
    tok = regexp(s2, '(\d{3,5})', 'tokens', 'once');
    if isempty(tok)
        mk(i) = strtrim(s2);
        continue
    end
    cage = str2double(tok{1});
    if ~isfinite(cage)
        mk(i) = strtrim(s2);
        continue
    end
    cageStr = sprintf('%04d', cage);

    % remove everything up to end of cage digits occurrence
    idx = regexp(s2, tok{1}, 'once');
    rest = strtrim(extractAfter(s2, idx + strlength(tok{1}) - 1));

    % color token = first letters-only word
    colTok = regexp(rest, '([a-z]+)', 'tokens', 'once');
    if isempty(colTok)
        % maybe it was glued like "0910red"
        glued = regexp(s2, '^\d{3,5}([a-z]+)', 'tokens', 'once');
        if ~isempty(glued), colTok = glued; end
    end
    if isempty(colTok)
        mk(i) = strtrim(s2);
        continue
    end
    color = string(colTok{1});
    mk(i) = strtrim(cageStr + " " + color);
end
end

%% =================== Burst sanity checks (unchanged logic) ===================
function burst_sanity_check_freqs(T, tb, Sessions, outDir, excludeMice)
mat = fullfile(outDir,'burst_thresholds_active.mat');
assert(exist(mat,'file')>0, 'Missing %s (run analyze_lick_bursts_GMM first).', mat);
S = load(mat);

allow = true(height(T),1);
if nargin>=5 && ~isempty(excludeMice)
    allow = ~ismember(keynorm(T.mouse_key), keynorm(excludeMice));
end

[Gsess, ~] = findgroups(string(T.mouse_key(allow)), double(T.day_index(allow)), double(T.session_idx(allow)));
IEI = [];
tbA = double(tb(allow)); L = logical(T.Lick_TTL(allow));
for g=1:max(Gsess)
    r = find(Gsess==g);
    on = ttl_onsets_local(tbA(r), L(r));
    if numel(on)>=3
        d = diff(on);
        d  = d(d>=0.02 & d<=10);
        IEI = [IEI; d(:)]; %#ok<AGROW>
    end
end
IEI = IEI(isfinite(IEI));
if numel(IEI) < 50
    warning('Too few IEIs for sanity plots.'); return;
end

thr_intra = S.thr_intra_s;
thr_gap   = S.thr_gap_s;
inside    = IEI <= thr_intra;
outside   = IEI  > thr_intra;
outside_nogap = IEI > thr_intra & IEI < thr_gap;

Hz_in  = 1./IEI(inside);
Hz_out = 1./IEI(outside);

fig = figure('Color','w','Position',[80 80 820 420]); hold on
edgesHz = 0:0.5:20;
histogram(Hz_in,  edgesHz, 'Normalization','probability','FaceAlpha',0.7);
histogram(Hz_out, edgesHz, 'Normalization','probability','FaceAlpha',0.7);
xlabel('Instantaneous lick frequency (Hz)'); ylabel('Probability');
title(sprintf('Inside vs outside bursts (thr_{intra}=%.3fs, thr_{gap}=%.3fs)', thr_intra, thr_gap));
legend({'Inside burst (IEI \le thr_{intra})','Outside burst (IEI > thr_{intra})'},'Location','northeast');
grid on; box on
savepng_local(fig, fullfile(outDir,'sanity_freq_inside_vs_outside.png')); close(fig);

fig = figure('Color','w','Position',[80 80 820 420]); hold on
edgesIEI = [0:0.02:1.5, 2:0.1:10];
histogram(IEI(inside),  edgesIEI, 'Normalization','pdf','FaceAlpha',0.7);
histogram(IEI(outside), edgesIEI, 'Normalization','pdf','FaceAlpha',0.7);
set(gca,'XScale','log'); xlabel('IEI (s, log)'); ylabel('PDF');
title('IEI distributions: inside vs outside bursts'); grid on; box on
legend({'Inside','Outside'},'Location','northeast');
xline(thr_intra,'b--','thr_{intra}');
xline(thr_gap,'r--','thr_{gap}');
savepng_local(fig, fullfile(outDir,'sanity_iei_inside_vs_outside.png')); close(fig);

sumtbl = table();
sumtbl.metric     = ["Hz_inside"; "Hz_outside"];
sumtbl.N          = [numel(Hz_in); numel(Hz_out)];
sumtbl.mean       = [mean(Hz_in,'omitnan');  mean(Hz_out,'omitnan')];
sumtbl.median     = [median(Hz_in,'omitnan'); median(Hz_out,'omitnan')];
sumtbl.std        = [std(Hz_in,'omitnan');    std(Hz_out,'omitnan')];
x = Hz_in; y = Hz_out;
pooled = sqrt(((numel(x)-1)*var(x,'omitnan') + (numel(y)-1)*var(y,'omitnan'))/max(1,(numel(x)+numel(y)-2)));
cohen_d = (mean(x,'omitnan')-mean(y,'omitnan')) / max(pooled,eps);
try, p_rs = ranksum(x,y); catch, p_rs = NaN; end
try, [~,p_ks] = kstest2(x,y); catch, p_ks = NaN; end
sumtbl.cohen_d   = [cohen_d; cohen_d];
sumtbl.p_ranksum = [p_rs; p_rs];
sumtbl.p_ks      = [p_ks; p_ks];

Textra = table(nnz(inside), nnz(outside), nnz(outside_nogap), thr_intra, thr_gap, 'VariableNames', ...
    {'n_IEI_inside','n_IEI_outside','n_IEI_outside_noGap','thr_intra_s','thr_gap_s'});
writetable(sumtbl, fullfile(outDir,'sanity_freq_inside_vs_outside_summary.csv'));
writetable(Textra, fullfile(outDir,'sanity_freq_inside_vs_outside_counts.csv'));
end
