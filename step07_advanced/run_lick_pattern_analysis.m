function run_lick_master_analysis()
% MASTER licking analysis (standalone, no toolboxes required).
% - Builds per-trial lick rate
% - PCA + kmeans (unsupervised)
% - Cluster behavior across day and across trial (mean±SEM)
% - Trend classifications: trials-only, days-only, per-day trial trend, and combined
% - Only raw p-values are reported; no FDR anywhere.
%
% Expected columns in ALL_mice_longitudinal.csv:
% mouse_key, day_index, session_idx, Lick_TTL, and any one timebase:
%   CamTime_rel_s | PupilTimestamp_s | CamTime_s | PlotTime_s_30fps | Frame
% Optional: Trial

%% ---------- locate latest run + read ----------
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

% ---- optional hard exclusions; edit or remove if you want everything
if ismember('mouse_key', T.Properties.VariableNames)
    mk      = string(T.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    drop = (contains(mk_norm,"6872") & contains(mk_norm,"black")) | ...
           (contains(mk_norm,"8606") & contains(mk_norm,"forange"));
    if any(drop)
        fprintf('Excluding %d rows (hard list)\n', nnz(drop)); T(drop,:) = [];
    end
end

% ---- types & basics
T = ensureString(T,'mouse_key');
reqCols = {'day_index','session_idx','Lick_TTL'};
for c = 1:numel(reqCols), assert(ismember(reqCols{c},T.Properties.VariableNames), 'Missing %s', reqCols{c}); end
if ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end
T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;

% ---- timebase (seconds)
tb = pickTimebase(T);
assert(any(isfinite(tb)),'No usable timebase column found.');

outDir = fullfile(runDir,'figs','lick_patterns_MASTER');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---------- per-session + per-trial metrics ----------
sessKeys = unique(T(:,{'mouse_key','day_index','session_idx'}),'rows','stable');

% per-session: lick_per_min + IEI, CV2, rhythm, bursts/bouts
Sess = sessKeys;
[Sess.lick_count,Sess.session_min,Sess.lick_per_min, ...
 Sess.iei_median,Sess.iei_mean,Sess.iei_cv,Sess.cv2_median, ...
 Sess.rhythm_index,Sess.burst_fraction,Sess.bout_count,Sess.bout_rate_per_min] ...
   = deal(nan(height(Sess),1));

for i=1:height(Sess)
    r = T.mouse_key==Sess.mouse_key(i) & T.day_index==Sess.day_index(i) & T.session_idx==Sess.session_idx(i);
    if ~any(r), continue; end
    t  = double(tb(r));
    l  = logical(T.Lick_TTL(r));
    F = compute_session_features(t,l);
    Sess.lick_count(i)        = F.lick_count;
    Sess.session_min(i)       = F.session_min;
    Sess.lick_per_min(i)      = F.lick_per_min;
    Sess.iei_median(i)        = F.iei_median;
    Sess.iei_mean(i)          = F.iei_mean;
    Sess.iei_cv(i)            = F.iei_cv;
    Sess.cv2_median(i)        = F.cv2_median;
    Sess.rhythm_index(i)      = F.rhythm_index;
    Sess.burst_fraction(i)    = F.burst_fraction;
    Sess.bout_count(i)        = F.bout_count;
    Sess.bout_rate_per_min(i) = F.bout_rate_per_min;
end
writetable(Sess, fullfile(outDir,'per_session_features.csv'));

% per-trial lick rate (Hz)  = #licks_in_trial / trial_duration
TT = table(); TTc = table();
if ismember('Trial', T.Properties.VariableNames)
    TT = per_trial_rates(T, tb);
    writetable(TT, fullfile(outDir,'per_trial_rates.csv'));
end

%% ---------- PCA + k-means (unsupervised) ----------
% Use simple, interpretable features: lick_per_min, cv2_median, rhythm_index, bout_rate_per_min
feat = Sess(:,{'mouse_key','day_index','session_idx','lick_per_min','cv2_median','rhythm_index','bout_rate_per_min'});

% If trials exist, add mean per-trial lick rate as another feature
if ~isempty(TT)
    mTR = groupsummary(TT, {'mouse_key','day_index','session_idx'}, 'mean','lick_rate');
    feat = outerjoin(feat, mTR(:,{'mouse_key','day_index','session_idx','mean_lick_rate'}), ...
                     'Keys',{'mouse_key','day_index','session_idx'}, 'MergeKeys',true);
else
    feat.mean_lick_rate = feat.lick_per_min;
end

X = [double(feat.lick_per_min), double(feat.cv2_median), ...
     double(feat.rhythm_index), double(feat.bout_rate_per_min), ...
     double(feat.mean_lick_rate)];
good = all(isfinite(X),2);
feat = feat(good,:); X = X(good,:);
Xz = zscore_nan(X);

% PCA
try
    [~,S,~,~,expl] = pca(Xz);
catch
    Xc = Xz - mean(Xz,1,'omitnan'); [U,~,~] = svd(Xc,'econ'); S = U(:,1:2); expl=[NaN NaN];
end

% kmeans (fallback to local impl if toolbox missing)
K = 3;
try
    idx = kmeans(Xz, K, 'Replicates',10, 'MaxIter',500, 'Display','off');
catch
    idx = kmeans_local(Xz, K, 10, 500);
end
feat.cluster = idx;

% PCA scatter by cluster
fig = figure('Color','w','Position',[80 80 720 560]); hold on
mk = {'o','s','^','d','p','h'};
for k=1:K
    mask = feat.cluster==k;
    scatter(S(mask,1), S(mask,2), 28, 'filled', 'Marker', mk{mod(k-1,numel(mk))+1}, ...
        'DisplayName', sprintf('Cluster %d',k), 'MarkerFaceAlpha',0.8);
end
xlabel(sprintf('PC1 (%.1f%%)', expl(1))); ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
title('Sessions: PCA colored by k-means cluster'); grid on; box on; legend('Location','bestoutside');
savepng(fig, fullfile(outDir,'pca_kmeans_clusters.png'));

% Clusters across DAY (mean±SEM lick_per_min vs day)
fig = figure('Color','w','Position',[90 90 860 420]); hold on
days = unique(feat.day_index); days = days(:)';
for k=1:K
    y = nan(size(days)); e = nan(size(days));
    for j=1:numel(days)
        rr = feat.cluster==k & feat.day_index==days(j);
        v  = feat.lick_per_min(rr);
        y(j)= mean(v,'omitnan'); e(j)= sem(v);
    end
    shaded_line(days, y, e, sprintf('Cluster %d',k));
end
xlabel('Day'); ylabel('Licks / min'); title('Clusters: mean ± SEM across day');
grid on; box on; legend('Location','best'); applyRobustYLim(gca);
savepng(fig, fullfile(outDir,'clusters_across_day_meanSEM.png'));

% Clusters across TRIAL (mean±SEM lick_rate vs trial)
if ~isempty(TT)
    key = TT(:,{'mouse_key','day_index','session_idx'});
    key = innerjoin(key, feat(:,{'mouse_key','day_index','session_idx','cluster'}), 'Keys',{'mouse_key','day_index','session_idx'});
    TTc = innerjoin(TT, key, 'Keys',{'mouse_key','day_index','session_idx'});
    trials = unique(TTc.Trial); trials = trials(:)';
    fig = figure('Color','w','Position',[90 90 860 420]); hold on
    for k=1:K
        y = nan(size(trials)); e = y;
        for j=1:numel(trials)
            v = TTc.lick_rate(TTc.cluster==k & TTc.Trial==trials(j));
            y(j)= mean(v,'omitnan'); e(j)= sem(v);
        end
        shaded_line(trials, y, e, sprintf('Cluster %d',k));
    end
    xlabel('Trial'); ylabel('Lick rate (Hz)'); title('Clusters: mean ± SEM across trial');
    grid on; box on; legend('Location','best'); applyRobustYLim(gca);
    savepng(fig, fullfile(outDir,'clusters_across_trial_meanSEM.png'));
end

%% ---------- Trend classifications ----------
% 1) Across trials (pooled over days)
TrialClass = table();
if ~isempty(TTc)
    um = unique(TTc.mouse_key,'stable'); rows={};
    for i=1:numel(um)
        m = um(i);
        r = TTc.mouse_key==m & isfinite(TTc.Trial) & isfinite(TTc.lick_rate);
        [rho,p] = spearman(double(TTc.Trial(r)), double(TTc.lick_rate(r)));
        rows(end+1,:) = {m, rho, p, trend_label(rho,p)}; %#ok<AGROW>
    end
    TrialClass = cell2table(rows,'VariableNames',{'mouse_key','rho_trials','p_trials','class_trials'});
    writetable(TrialClass, fullfile(outDir,'class_trials_overall.csv'));
end

% 2) Across days (trial-agnostic): per-day median lick_per_min
um = unique(Sess.mouse_key,'stable'); rows={};
for i=1:numel(um)
    m = um(i);
    r = Sess.mouse_key==m & isfinite(Sess.day_index) & isfinite(Sess.lick_per_min);
    if ~any(r), continue; end
    D = groupsummary(Sess(r,{'day_index','lick_per_min'}), 'day_index','median','lick_per_min');
    [rho,p] = spearman(double(D.day_index), double(D.median_lick_per_min));
    rows(end+1,:) = {m, rho, p, trend_label(rho,p)}; %#ok<AGROW>
end
DayClass = cell2table(rows,'VariableNames',{'mouse_key','rho_days','p_days','class_days'});
writetable(DayClass, fullfile(outDir,'class_days_overall.csv'));

% 3) Per-day trial trend (one label per mouse×day)
TrialPerDay = table();
if ~isempty(TTc)
    gd = unique(TTc(:,{'mouse_key','day_index'}),'rows','stable'); rows={};
    for i=1:height(gd)
        m = gd.mouse_key(i); d = gd.day_index(i);
        r = TTc.mouse_key==m & TTc.day_index==d & isfinite(TTc.Trial) & isfinite(TTc.lick_rate);
        if nnz(r)<3, continue; end
        [rho,p] = spearman(double(TTc.Trial(r)), double(TTc.lick_rate(r)));
        rows(end+1,:) = {m, d, rho, p, trend_label(rho,p)}; %#ok<AGROW>
    end
    TrialPerDay = cell2table(rows,'VariableNames',{'mouse_key','day_index','rho_trial_per_day','p_trial_per_day','class_trial_per_day'});
    writetable(TrialPerDay, fullfile(outDir,'class_trials_per_day.csv'));
end

% 4) Combined 2-D view and merged label
Comb = outerjoin(TrialClass, DayClass, 'Keys','mouse_key','MergeKeys',true);
Comb.combined_label = strings(height(Comb),1);
for i=1:height(Comb)
    a = safeStr(Comb.class_trials(i)); b = safeStr(Comb.class_days(i));
    if a=="Faster" && b=="Faster", Comb.combined_label(i)="Faster";
    elseif a=="Slower" && b=="Slower", Comb.combined_label(i)="Slower";
    else, Comb.combined_label(i)="NoChange";
    end
end
writetable(Comb, fullfile(outDir,'class_combined_day_and_trial.csv'));

% Scatter map of rho (day) vs rho (trial)
fig = figure('Color','w','Position',[90 90 560 500]); hold on
scatter(Comb.rho_days, Comb.rho_trials, 40, 'filled');
xlabel('Spearman \rho (day vs lick/min)'); ylabel('Spearman \rho (trial vs lick rate)');
title('Mouse trend map (day-trend vs trial-trend)'); grid on; box on; applyRobustYLim(gca);
savepng(fig, fullfile(outDir,'trend_map_day_vs_trial.png'));

fprintf('Done. Outputs in:\n  %s\n', outDir);
end

%% ======================= helpers (all inline) =======================
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
            if any(isfinite(v)), tb = strcmpi(cands{i},'Frame') * (v/30) + ~strcmpi(cands{i},'Frame') * v; return; end
        end
    end
end
function [nL, durMin] = session_counts(t,l)
    t = double(t(:)); l = logical(l(:));
    good = isfinite(t) & ~isnan(l); t=t(good); l=l(good);
    d  = diff([false; l; false]); on = t(d(1:end-1)==1);
    nL = numel(on); durMin = max(t)-min(t); durMin = durMin/60;
end
function F = compute_session_features(t, ttl)
    t = double(t(:)); ttl = logical(ttl(:));
    good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);
    F = struct('lick_count',0,'session_min',NaN,'lick_per_min',0,'iei_median',NaN,'iei_mean',NaN,'iei_cv',NaN,'cv2_median',NaN,'rhythm_index',NaN,'burst_fraction',NaN,'bout_count',0,'bout_rate_per_min',0);
    if numel(t)<2, return; end
    d  = diff([false; ttl; false]); on = t(d(1:end-1)==1); if isempty(on), return; end
    dur = (max(t)-min(t)); F.session_min = dur/60; nL = numel(on); F.lick_count = nL; F.lick_per_min = nL/max(F.session_min,eps);
    if nL>=2
        iei = diff(on);
        F.iei_median = median(iei,'omitnan');
        F.iei_mean   = mean(iei,'omitnan');
        F.iei_cv     = std(iei,'omitnan')/max(mean(iei,'omitnan'),eps);
        cv2 = 2*abs(diff(iei))./(iei(1:end-1)+iei(2:end)); F.cv2_median = median(cv2,'omitnan');
        band = iei>=0.11 & iei<=0.18; F.rhythm_index = nnz(band)/max(numel(iei),1);
        thr_within = 0.20; thr_gap = 0.50;
        isShort = iei < thr_within; F.burst_fraction = nnz(isShort)/max(numel(iei),1);
        gap = [inf; iei(:)]; inBout=false; bc=0;
        for k=2:numel(gap)
            if ~inBout && gap(k)<thr_within, inBout=true; bc=bc+1;
            elseif inBout && gap(k)>=thr_gap, inBout=false; end
        end
        F.bout_count = bc; F.bout_rate_per_min = bc/max(F.session_min,eps);
    end
end
function TT = per_trial_rates(T, tb)
    g = unique(T(:,{'mouse_key','day_index','session_idx'}),'rows','stable');
    out = {};
    for i=1:height(g)
        mk=g.mouse_key(i); di=g.day_index(i); ss=g.session_idx(i);
        rows = T.mouse_key==mk & T.day_index==di & T.session_idx==ss;
        if ~any(rows) || ~ismember('Trial',T.Properties.VariableNames), continue; end
        Trial = double(T.Trial(rows)); good = isfinite(Trial) & Trial>0;
        if ~any(good), continue; end
        ti = double(tb(rows)); ti = ti(good); tr = Trial(good); lk = logical(T.Lick_TTL(rows)); lk = lk(good);
        trials = sort(unique(tr));
        for k=1:numel(trials)
            m = tr==trials(k); tK=ti(m); lK=lk(m);
            on = lick_onsets(tK,lK); dur = max(tK)-min(tK);
            n   = numel(on); rate = n/max(dur,eps);
            out(end+1,:) = {mk, di, ss, double(trials(k)), double(dur), double(n), double(rate)}; %#ok<AGROW>
        end
    end
    if isempty(out), TT=table(); return; end
    TT = cell2table(out,'VariableNames',{'mouse_key','day_index','session_idx','Trial','trial_duration_s','lick_count','lick_rate'});
end
function on = lick_onsets(t,ttl)
    d = diff([false; ttl(:); false]); on = double(t(d(1:end-1)==1));
end
function x = zscore_nan(x)
    mu = mean(x,1,'omitnan'); sg = std(x,0,1,'omitnan'); x = (x - mu)./max(sg,eps); x(~isfinite(x))=0;
end
function semv = sem(v), v=v(isfinite(v)); if isempty(v), semv=NaN; else, semv=std(v)/sqrt(numel(v)); end, end
function shaded_line(x,y,e,lab)
    x = x(:)'; y=y(:)'; e=e(:)'; if numel(x)<2, return; end
    fill([x fliplr(x)], [y-e fliplr(y+e)], [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.55); hold on
    plot(x,y,'-o','LineWidth',1.6,'MarkerSize',4,'DisplayName',lab);
end
function savepng(fh, fn)
    set(fh,'PaperPositionMode','auto'); try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
    close(fh);
end
function applyRobustYLim(ax)
    ln = findobj(ax,'Type','line'); vals=[]; for i=1:numel(ln), vals=[vals; get(ln(i),'YData')]; end
    vals=vals(isfinite(vals)); if numel(vals)<4, return; end
    lo=prctile(vals,2); hi=prctile(vals,98); if isfinite(lo)&&isfinite(hi)&&hi>lo, pad=0.08*(hi-lo); ylim(ax,[lo-pad, hi+pad]); end
end
function [rho,p] = spearman(x,y)
    x=x(:); y=y(:); good=isfinite(x)&isfinite(y); x=x(good); y=y(good);
    if numel(x)<3, rho=NaN; p=NaN; return; end
    try, [rho,p]=corr(x,y,'Type','Spearman','Rows','complete');
    catch, xr=tiedrank(x); yr=tiedrank(y); [rho,p]=corr(xr,yr,'Rows','complete'); end
end
function lab = trend_label(rho,p)
    if ~isfinite(rho)||~isfinite(p)||p>=0.05, lab="NoChange"; elseif rho>0, lab="Faster"; else, lab="Slower"; end
end
function s = safeStr(x), if ismissing(x) || strlength(string(x))==0, s=""; else, s=string(x); end, end
function idx = kmeans_local(X, K, reps, maxit)
    % very small fallback k-means (Lloyd)
    [n,~]=size(X); best=inf; idx=ones(n,1);
    for r=1:reps
        mu = X(randperm(n,K),:);
        for it=1:maxit
            D = pdist2(X,mu); [~,ix]=min(D,[],2);
            mu2 = zeros(size(mu));
            for k=1:K
                if any(ix==k), mu2(k,:)=mean(X(ix==k,:),1); else, mu2(k,:)=mu(k,:); end
            end
            if max(abs(mu2-mu),[],'all')<1e-6, break; end
            mu=mu2;
        end
        sse = sum(min(pdist2(X,mu).^2,[],2));
        if sse<best, best=sse; idx=ix; end
    end
end

