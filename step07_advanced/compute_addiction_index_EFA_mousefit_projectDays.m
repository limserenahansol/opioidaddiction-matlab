function compute_addiction_index_EFA_mousefit_projectDays_v2()
% compute_addiction_index_EFA_mousefit_projectDays_v2
%
% Paper-style:
%   1) Fit EFA using mouse-level aggregated features (one row per mouse).
%   2) Choose nf by increasing until >= 85% common variance explained.
%   3) Project mouse×day data onto learned loadings -> factor scores per mouse×day.
%   4) Addiction Index (AI) per mouse×day + plots (spaghetti + phase comparisons).
%
% Key fixes vs prior version:
%   - Excludes day1,2 entirely (uses day>=3 only).
%   - Uses MOUSE-LEVEL grouping (MouseGroup) for Active vs Passive.
%   - Avoids size mismatch by carrying mu/sd only for final kept features.
%   - Robust feature reduction for small nMice.
%
% Requires Statistics and Machine Learning Toolbox (factoran).

%% ----------------- USER SETTINGS -----------------
runDir = 'K:\addiction_concate_Dec_2025\longitudinal_outputs\run_009'; % adjust
outDir = fullfile(runDir,'figs','addiction_index_EFA');
if ~exist(outDir,'dir'), mkdir(outDir); end

cacheMat = fullfile(runDir,'S_D_cache.mat');
assert(exist(cacheMat,'file')>0, 'Missing %s', cacheMat);
S = load(cacheMat);
assert(isfield(S,'D'), 'S_D_cache.mat must contain table D');
D = S.D;

% Exclude day 1-2 (HARD REQUIREMENT)
assert(ismember('day_index', D.Properties.VariableNames), 'D must contain day_index');
D = D(double(D.day_index) >= 3, :);

% Transition days (optional)
transitionDays = [3 6 11 14];
includeTransitionDays = true; % set false if you want to exclude them
if ~includeTransitionDays
    D = D(~ismember(double(D.day_index), transitionDays), :);
end

% Period definition (edit if needed)
D.Period = periodOfDay(double(D.day_index));
D = D(~isundefined(D.Period), :);

% ----------------- MOUSE GROUP (Active vs Passive) -----------------
% Recommended: Fill these with your cohort mouse IDs (mouse_key values).
% If left empty, the script will try to infer group from existing columns.
activeKeys  = string([]);   % e.g. ["6100","0911","..."]
passiveKeys = string([]);   % e.g. ["6099","0910","..."]

D = attachMouseGroup(D, activeKeys, passiveKeys);
% After this, D.MouseGroup exists and is categorical {'Active','Passive'} where known.

% ----------------- Feature list -----------------
candidateFeatures = { ...
    'RequirementLast','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
    'bout_n','bout_meanDur_s','rew_freq_per_min','rew_meanDur_s','pupil_mean', ...
    'lick_totalDur_s','bout_totalDur_s','rew_totalDur_s', ...
    'Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min' ...
    };

features = candidateFeatures(ismember(candidateFeatures, D.Properties.VariableNames));
assert(numel(features) >= 4, 'Not enough features found in D to run EFA.');

% Missingness thresholds
maxMissingPerFeature = 0.40; % drop feature if >40% missing at mouse-level
maxMissingPerMouse   = 0.40; % drop mouse if >40% missing across retained features

% Factor selection
targetExplainedCommon = 0.85;  % add factors until >=85% common variance explained
maxFactorsHard = 12;

rotationMethod = 'varimax';

% Optional sign alignment reference (mouse-level)
signAlignRef = 'RequirementLast';
doSignAlign = ismember(signAlignRef, features);

%% ----------------- 1) Mouse-level aggregation -----------------
[mice, mouseTbl, Xmouse_raw, features_kept] = buildMouseLevelMatrix(D, features, ...
    maxMissingPerFeature, maxMissingPerMouse);

features = features_kept; %#ok<NASGU>

% Standardize mouse-level (store mu/sd for FINAL features)
[XmouseZ, muFeat, sdFeat] = zscore_with_store(Xmouse_raw); % muFeat/sdFeat correspond to final features

% Robust pre-EFA feature reduction due to small nMice
nMice = size(XmouseZ,1);
pFeat = size(XmouseZ,2);
fprintf('\nMouse-level matrix: nMice=%d, pFeat=%d\n', nMice, pFeat);

% Rule of thumb cap: p <= max(6, nMice-3)
pCap = max(6, nMice - 3);
if pFeat > pCap
    fprintf('Too many features for EFA stability (p=%d > %d). Reducing features...\n', pFeat, pCap);
    [XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features_final] = reduce_features_by_variance( ...
        XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features, pCap);
else
    features_final = features;
end

% Drop near-constant (in case)
[XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features_final] = drop_near_constant_features( ...
    XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features_final, 1e-6);

% Drop highly correlated duplicates
[XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features_final] = drop_high_corr_features( ...
    XmouseZ, Xmouse_raw, muFeat, sdFeat, mouseTbl, features_final, 0.98);

% Update after filtering
nMice = size(XmouseZ,1);
pFeat = size(XmouseZ,2);
maxFactors = min(maxFactorsHard, max(1, pFeat-1));

fprintf('After filtering: nMice=%d, pFeat=%d, maxFactors=%d\n', nMice, pFeat, maxFactors);

%% ----------------- 2) Choose number of factors by >=85% common variance explained -----------------
fprintf('\nSelecting # factors to reach >= %.0f%% common variance explained...\n', 100*targetExplainedCommon);

best = struct('nf',[],'Lambda',[],'Psi',[],'T',[],'stats',[],'commonExpl',[],'failMsg','');
for nf = 1:maxFactors
    [ok, Lambda, Psi, Tmat, stats, msg] = try_factoran_robust(XmouseZ, nf, rotationMethod);
    if ~ok
        fprintf('  nf=%d FAILED: %s\n', nf, msg);
        continue
    end

    commonExpl = sum(1 - Psi) / numel(Psi);
    fprintf('  nf=%d OK -> common variance explained ~ %.1f%%\n', nf, 100*commonExpl);

    best.nf = nf; best.Lambda = Lambda; best.Psi = Psi; best.T = Tmat; best.stats = stats;
    best.commonExpl = commonExpl; best.failMsg = '';

    if commonExpl >= targetExplainedCommon
        break
    end
end

assert(~isempty(best.nf), 'Could not fit EFA for any nf.');

nf     = best.nf;
Lambda = best.Lambda;
Psi    = best.Psi;
fprintf('Selected nf=%d (common variance explained ~ %.1f%%)\n', nf, 100*best.commonExpl);

%% ----------------- 3) Mouse-level factor scores + AI (one per mouse) -----------------
Fmouse = score_factors_regression(XmouseZ, Lambda, Psi); % nMice x nf

% Optional sign alignment (make factors positively correlate with ref variable)
if doSignAlign && ismember(signAlignRef, features_final)
    ref = double(mouseTbl.(signAlignRef));
    refZ = zscore_safe_vector(ref);
    for k = 1:nf
        r = corr(Fmouse(:,k), refZ, 'rows','complete');
        if isfinite(r) && r < 0
            Fmouse(:,k) = -Fmouse(:,k);
            Lambda(:,k) = -Lambda(:,k);
        end
    end
end

% Paper-like AI per mouse: mean of z-scored factor scores
[FmouseZ, muF, sdF] = zscore_with_store(Fmouse);
AI_mouse_raw = mean(FmouseZ, 2, 'omitnan');

% Scale to 0-1 for plotting
AI_mouse = rescale01(AI_mouse_raw);
mouseTbl.AI_mouse = AI_mouse;

%% ----------------- 4) Project mouse×day -> AI per mouse×day -----------------
Dsub = D(:, intersect([{'mouse_key','day_index','Period','MouseGroup'}, features_final], ...
    D.Properties.VariableNames, 'stable'));

% Ensure all final features exist
for j=1:numel(features_final)
    f = features_final{j};
    if ~ismember(f, Dsub.Properties.VariableNames)
        Dsub.(f) = NaN(height(Dsub),1);
    end
end

Xday_raw = Dsub{:, features_final};

% Impute missing with mouse-level feature medians (muFeat is median from zscore_with_store)
for j=1:size(Xday_raw,2)
    col = Xday_raw(:,j);
    col(~isfinite(col)) = muFeat(j);
    Xday_raw(:,j) = col;
end

% Standardize using mouse-level mu/sd (fixed basis)
XdayZ = (Xday_raw - muFeat) ./ sdFeat;
XdayZ(~isfinite(XdayZ)) = 0;

Fday = score_factors_regression(XdayZ, Lambda, Psi); % nRows x nf

% Put day-level factors on same scale as mouse-level factor z-scoring
FdayZ = (Fday - muF) ./ sdF;
FdayZ(~isfinite(FdayZ)) = 0;

AI_day_raw = mean(FdayZ, 2, 'omitnan');
AI_day = rescale01(AI_day_raw);

Dsub.AI = AI_day;

%% ----------------- 5) Plots -----------------
plot_spaghetti_AI(Dsub, outDir);              % all mice + by MouseGroup
plot_phase_comparison_AI(Dsub, outDir);       % phase comparison, by MouseGroup

%% ----------------- 6) Save outputs -----------------
save(fullfile(outDir,'AI_results.mat'), ...
    'mouseTbl','Dsub','features_final','nf','Lambda','Psi','muFeat','sdFeat','muF','sdF', ...
    'best','includeTransitionDays','transitionDays','targetExplainedCommon');

writetable(mouseTbl, fullfile(outDir,'AI_mouse_summary.csv'));
writetable(Dsub, fullfile(outDir,'AI_mouse_day.csv'));

fprintf('Saved AI outputs to:\n  %s\n', outDir);

end

%% ===================== Helpers =====================

function D = attachMouseGroup(D, activeKeys, passiveKeys)
% Creates D.MouseGroup as mouse-level cohort membership, not row-level period/trial label.

assert(ismember('mouse_key', D.Properties.VariableNames), 'D must contain mouse_key');

% If an existing mouse-level group column exists, use it
candCols = {'MouseGroup','mouseGroup','CohortGroup','cohortGroup','GroupType','groupType'};
for c = 1:numel(candCols)
    if ismember(candCols{c}, D.Properties.VariableNames)
        g = string(D.(candCols{c}));
        g = lower(strtrim(g));
        mg = strings(height(D),1);
        mg(contains(g,'passive')) = "Passive";
        mg(contains(g,'active'))  = "Active";
        mg(mg=="") = missing;
        D.MouseGroup = categorical(mg, ["Active","Passive"]);
        return
    end
end

% Otherwise infer from row-level columns if present (ever-passive rule)
rowG = strings(height(D),1);
if ismember('Group', D.Properties.VariableNames)
    rowG = string(D.Group);
elseif ismember('group', D.Properties.VariableNames)
    rowG = string(D.group);
elseif ismember('Condition', D.Properties.VariableNames)
    rowG = string(D.Condition);
end
rowG = lower(strtrim(rowG));

mice = unique(D.mouse_key,'stable');
mouseGroup = strings(numel(mice),1);

for i=1:numel(mice)
    rows = (D.mouse_key==mice(i));
    g = rowG(rows);
    g = g(g~="" & ~ismissing(g));

    % ever rule
    if any(contains(g,'passive'))
        mouseGroup(i) = "Passive";
    elseif any(contains(g,'active'))
        mouseGroup(i) = "Active";
    else
        mouseGroup(i) = "";
    end
end

% Deterministic override using provided keys (RECOMMENDED)
if ~isempty(activeKeys) || ~isempty(passiveKeys)
    for i=1:numel(mice)
        if any(mice(i)==passiveKeys), mouseGroup(i) = "Passive"; end
        if any(mice(i)==activeKeys),  mouseGroup(i) = "Active";  end
    end
end

% Attach to rows
map = containers.Map(cellstr(string(mice)), cellstr(string(mouseGroup)));
mg = strings(height(D),1);
for r=1:height(D)
    k = char(string(D.mouse_key(r)));
    if isKey(map,k)
        mg(r) = string(map(k));
    else
        mg(r) = "";
    end
end
mg(mg=="") = missing;
D.MouseGroup = categorical(mg, ["Active","Passive"]);
end

function P = periodOfDay(d)
% Adjust if your periods differ
p = strings(size(d));
p(d>=3  & d<=5 )  = "Pre";
p(d>=6  & d<=10)  = "During";
p(d>=11 & d<=13)  = "Post";
p(d>=14 & d<=16)  = "Withdrawal";
p(d>=17 & d<=18)  = "Re-exposure";
P = categorical(p, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

function [mice, mouseTbl, Xmouse_raw, features_kept] = buildMouseLevelMatrix(D, features, maxMissingPerFeature, maxMissingPerMouse)

mice = unique(D.mouse_key,'stable');
mouseTbl = table(mice, 'VariableNames', {'mouse_key'});

% Carry MouseGroup to mouseTbl (mode; should be constant)
if ismember('MouseGroup', D.Properties.VariableNames)
    mg = strings(numel(mice),1);
    for i=1:numel(mice)
        rows = D.mouse_key==mice(i);
        v = string(D.MouseGroup(rows));
        v = v(~ismissing(v));
        if isempty(v), mg(i) = ""; else, mg(i) = string(mode(categorical(v))); end
    end
    mg(mg=="") = missing;
    mouseTbl.MouseGroup = categorical(mg, ["Active","Passive"]);
end

% Mouse-level median per feature across all available days (>=3 already filtered)
for j=1:numel(features)
    f = features{j};
    x = nan(numel(mice),1);
    for i=1:numel(mice)
        rows = D.mouse_key==mice(i);
        v = double(D.(f)(rows));
        x(i) = median(v,'omitnan');
    end
    mouseTbl.(f) = x;
end

% Drop features with too much missing (mouse-level)
featKeep = true(1,numel(features));
for j=1:numel(features)
    f = features{j};
    missFrac = mean(~isfinite(double(mouseTbl.(f))));
    if missFrac > maxMissingPerFeature
        featKeep(j) = false;
        fprintf('Dropping feature %s (mouse-level missing=%.1f%%)\n', f, 100*missFrac);
    end
end
features_kept = features(featKeep);
assert(numel(features_kept) >= 4, 'After missingness filtering, not enough features for EFA.');

Xmouse_raw = mouseTbl{:, features_kept};

% Drop mice with too much missing across kept features
missPerMouse = mean(~isfinite(Xmouse_raw),2);
keepMouse = missPerMouse <= maxMissingPerMouse;
if any(~keepMouse)
    fprintf('Dropping %d mice due to missingness > %.0f%%\n', nnz(~keepMouse), 100*maxMissingPerMouse);
end

mouseTbl = mouseTbl(keepMouse,:);
Xmouse_raw = mouseTbl{:, features_kept};

% Impute remaining missing in mouse-level with feature medians
for j=1:size(Xmouse_raw,2)
    col = Xmouse_raw(:,j);
    med = median(col,'omitnan');
    col(~isfinite(col)) = med;
    Xmouse_raw(:,j) = col;
end

end

function [Z, mu, sd] = zscore_with_store(X)
% Robust standardization:
% - mu = median (omitnan)
% - sd = std (omitnan), guarded against 0
X = double(X);
mu = median(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(sd==0 | ~isfinite(sd)) = 1;
Z = (X - mu) ./ sd;
Z(~isfinite(Z)) = 0;
end

function z = zscore_safe_vector(x)
x = double(x(:));
mu = median(x,'omitnan');
sd = std(x,0,'omitnan');
if ~isfinite(sd) || sd==0, sd = 1; end
z = (x - mu) ./ sd;
z(~isfinite(z)) = 0;
end

function y = rescale01(x)
x = double(x);
good = isfinite(x);
if ~any(good)
    y = 0.5*ones(size(x));
    return
end
lo = min(x(good));
hi = max(x(good));
if hi==lo
    y = 0.5*ones(size(x));
else
    y = (x - lo) / (hi - lo);
end
y(~isfinite(y)) = 0.5;
end

function F = score_factors_regression(Xz, Lambda, Psi)
% Regression factor scores with fixed loadings
Xz = double(Xz);
Psi = double(Psi(:));
k = size(Lambda,2);
invPsi = diag(1 ./ max(Psi, 1e-6));
A = (Lambda' * invPsi * Lambda + eye(k));
W = (A \ (Lambda' * invPsi)); % k x p
F = (W * Xz')'; % n x k
F(~isfinite(F)) = 0;
end

function [ok, Lambda, Psi, Tmat, stats, msg] = try_factoran_robust(Xz, nf, rotationMethod)
ok = false; Lambda=[]; Psi=[]; Tmat=[]; stats=[]; msg='';

msg1 = ''; msg2 = ''; msg3 = '';

% Attempt 1: raw
try
    [Lambda, Psi, Tmat, stats] = factoran(Xz, nf, 'rotate', rotationMethod, 'scores','regression', 'maxit', 3000);
    ok = true; return
catch ME
    msg1 = ME.message;
end

% Attempt 2: jitter
try
    Xj = Xz + 1e-6*randn(size(Xz));
    [Lambda, Psi, Tmat, stats] = factoran(Xj, nf, 'rotate', rotationMethod, 'scores','regression', 'maxit', 3000);
    ok = true; return
catch ME
    msg2 = ME.message;
end

% Attempt 3: correlation-matrix with ridge
try
    R = corrcoef(Xz);
    R = (R + R')/2;
    ridge = 1e-3;
    R = R + ridge*eye(size(R,1));
    [Lambda, Psi, Tmat, stats] = factoran(R, nf, 'nobs', size(Xz,1), 'rotate', rotationMethod, 'scores','regression', 'maxit', 3000);
    ok = true; return
catch ME
    msg3 = ME.message;
end

msg = sprintf('raw: %s | jitter: %s | R+ridge: %s', msg1, msg2, msg3);
end

function [Xz, Xraw, mu, sd, mouseTbl, feat] = reduce_features_by_variance(Xz, Xraw, mu, sd, mouseTbl, feat, pCap)
v = var(Xz,0,1);
[~,ord] = sort(v,'descend','MissingPlacement','last');
ord = ord(1:min(pCap, numel(ord)));

Xz   = Xz(:,ord);
Xraw = Xraw(:,ord);
mu   = mu(ord);
sd   = sd(ord);
feat = feat(ord);

% also keep those columns in mouseTbl
keepNames = feat;
keepTbl = mouseTbl(:, [{'mouse_key'} intersect(mouseTbl.Properties.VariableNames, {'MouseGroup'}, 'stable') keepNames]);
mouseTbl = keepTbl;
end

function [Xz, Xraw, mu, sd, mouseTbl, feat] = drop_near_constant_features(Xz, Xraw, mu, sd, mouseTbl, feat, tolVar)
v = var(Xz,0,1);
keep = (v > tolVar) & isfinite(v);
if any(~keep)
    fprintf('Dropping %d near-constant features:\n', nnz(~keep));
    disp(feat(~keep));
end
Xz   = Xz(:,keep);
Xraw = Xraw(:,keep);
mu   = mu(keep);
sd   = sd(keep);
feat = feat(keep);

keepNames = feat;
keepTbl = mouseTbl(:, [{'mouse_key'} intersect(mouseTbl.Properties.VariableNames, {'MouseGroup'}, 'stable') keepNames]);
mouseTbl = keepTbl;
end

function [Xz, Xraw, mu, sd, mouseTbl, feat] = drop_high_corr_features(Xz, Xraw, mu, sd, mouseTbl, feat, thr)
if size(Xz,2) < 2, return; end
R = corrcoef(Xz);
R(1:size(R,1)+1:end) = 0;
toDrop = false(1,size(Xz,2));

for j=1:size(R,1)
    if toDrop(j), continue; end
    dup = find(abs(R(j,:)) > thr);
    dup = dup(~toDrop(dup));
    if ~isempty(dup)
        % drop all dup indices (they are "others"); keep j
        toDrop(dup) = true;
    end
end
toDrop = toDrop & (1:numel(toDrop))~=find(~toDrop,1,'first'); % safety (keep at least one)

if any(toDrop)
    fprintf('Dropping %d highly-correlated features (|r|>%.2f):\n', nnz(toDrop), thr);
    disp(feat(toDrop));
end

keep = ~toDrop;
Xz   = Xz(:,keep);
Xraw = Xraw(:,keep);
mu   = mu(keep);
sd   = sd(keep);
feat = feat(keep);

keepNames = feat;
keepTbl = mouseTbl(:, [{'mouse_key'} intersect(mouseTbl.Properties.VariableNames, {'MouseGroup'}, 'stable') keepNames]);
mouseTbl = keepTbl;
end

function plot_spaghetti_AI(Dsub, outDir)
% All mice
fh = figure('Color','w','Position',[80 80 920 520]); hold on;

mice = unique(Dsub.mouse_key,'stable');
for i=1:numel(mice)
    r = Dsub.mouse_key==mice(i);
    dd = double(Dsub.day_index(r));
    [dd,ord] = sort(dd);
    ai = double(Dsub.AI(r)); ai = ai(ord);
    plot(dd, ai, '-', 'LineWidth',0.9);
end

days = unique(double(Dsub.day_index));
meanAI = nan(size(days));
for j=1:numel(days)
    meanAI(j) = mean(double(Dsub.AI(double(Dsub.day_index)==days(j))), 'omitnan');
end
plot(days, meanAI, 'k-', 'LineWidth',2.5);

xlabel('Day'); ylabel('AI (0–1)');
title('Addiction Index trajectories (all mice)');
ylim([0 1]);
grid off; box on;

exportgraphics(fh, fullfile(outDir,'AI_spaghetti_all.png'), 'Resolution',180);
close(fh);

% By MouseGroup (mouse cohort), same y-axis across panels
if ismember('MouseGroup', Dsub.Properties.VariableNames)
    groups = categories(Dsub.MouseGroup);
    groups = groups(~cellfun(@isempty,groups));

    fh = figure('Color','w','Position',[80 80 980 520]);
    tiledlayout(1, numel(groups), 'TileSpacing','compact','Padding','compact');

    for gi=1:numel(groups)
        ax = nexttile; hold(ax,'on');
        G = groups{gi};
        Dg = Dsub(Dsub.MouseGroup==G,:);
        miceg = unique(Dg.mouse_key,'stable');

        for i=1:numel(miceg)
            r = Dg.mouse_key==miceg(i);
            dd = double(Dg.day_index(r));
            [dd,ord] = sort(dd);
            ai = double(Dg.AI(r)); ai = ai(ord);
            plot(ax, dd, ai, '-', 'LineWidth',0.9);
        end

        days = unique(double(Dg.day_index));
        meanAI = nan(size(days));
        for j=1:numel(days)
            meanAI(j) = mean(double(Dg.AI(double(Dg.day_index)==days(j))), 'omitnan');
        end
        plot(ax, days, meanAI, 'k-', 'LineWidth',2.5);

        title(ax, sprintf('%s', G));
        xlabel(ax,'Day'); ylabel(ax,'AI (0–1)');
        ylim(ax,[0 1]);
        grid(ax,'off'); box(ax,'on');
    end

    exportgraphics(fh, fullfile(outDir,'AI_spaghetti_byMouseGroup.png'), 'Resolution',180);
    close(fh);
end
end

function plot_phase_comparison_AI(Dsub, outDir)
% Aggregate to mouse×Period median to avoid day-count bias
assert(ismember('Period',Dsub.Properties.VariableNames), 'Dsub must contain Period');
assert(ismember('AI',Dsub.Properties.VariableNames), 'Dsub must contain AI');

useGroup = ismember('MouseGroup', Dsub.Properties.VariableNames) && ~all(ismissing(string(Dsub.MouseGroup)));

if useGroup
    G = groupsummary(Dsub(:,{'mouse_key','MouseGroup','Period','AI'}), ...
        {'mouse_key','MouseGroup','Period'}, 'median', 'AI');
    G.Properties.VariableNames{'median_AI'} = 'AI_med';
else
    % no group -> still do period comparison per mouse
    G = groupsummary(Dsub(:,{'mouse_key','Period','AI'}), {'mouse_key','Period'}, 'median','AI');
    G.Properties.VariableNames{'median_AI'} = 'AI_med';
end

periods = categories(Dsub.Period);

fh = figure('Color','w','Position',[90 90 980 520]); hold on;

if useGroup
    groups = categories(G.MouseGroup);
    off = linspace(-0.15,0.15,max(1,numel(groups)));

    for gi=1:numel(groups)
        Gi = G(G.MouseGroup==groups{gi},:);

        for pi=1:numel(periods)
            P = Gi(Gi.Period==periods{pi} & isfinite(Gi.AI_med),:);
            if isempty(P), continue; end

            x = pi + off(gi);
            boxchart(repmat(x,height(P),1), P.AI_med, ...
                'BoxWidth',0.22, 'MarkerStyle','none', 'BoxFaceAlpha',0.20);
            scatter(repmat(x,height(P),1), P.AI_med, 24, 'k', 'filled', 'MarkerFaceAlpha',0.65);
        end
    end
else
    for pi=1:numel(periods)
        P = G(G.Period==periods{pi} & isfinite(G.AI_med),:);
        if isempty(P), continue; end
        boxchart(repmat(pi,height(P),1), P.AI_med, 'BoxWidth',0.35, 'MarkerStyle','none', 'BoxFaceAlpha',0.20);
        scatter(repmat(pi,height(P),1), P.AI_med, 24, 'k', 'filled', 'MarkerFaceAlpha',0.65);
    end
end

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlim([0.5 numel(periods)+0.5]);
ylabel('AI (0–1)');
title('AI phase comparison (mouse-level median per phase)');
ylim([0 1]);
grid off; box on;

exportgraphics(fh, fullfile(outDir,'AI_phase_comparison.png'), 'Resolution',180);
close(fh);

% Optional: Passive vs Active stats within each period
if useGroup && numel(categories(G.MouseGroup))==2
    fh = figure('Color','w','Position',[90 90 1200 520]);
    tiledlayout(1,numel(periods),'TileSpacing','compact','Padding','compact');

    for pi=1:numel(periods)
        ax = nexttile; hold(ax,'on');
        Pa = G(G.MouseGroup=="Passive" & G.Period==periods{pi},:);
        Ac = G(G.MouseGroup=="Active"  & G.Period==periods{pi},:);

        x = Pa.AI_med; y = Ac.AI_med;
        p = NaN;
        if numel(x)>=2 && numel(y)>=2
            p = ranksum(x,y);
        end

        scatter(ax, ones(size(x)), x, 30, 'filled');
        scatter(ax, 2*ones(size(y)), y, 30, 'filled');
        set(ax,'XTick',[1 2],'XTickLabel',{'Passive','Active'});
        ylim(ax,[0 1]);
        title(ax, sprintf('%s  p=%.3g', string(periods{pi}), p));
        grid(ax,'off'); box(ax,'on');
    end

    exportgraphics(fh, fullfile(outDir,'AI_phase_stats_panels.png'), 'Resolution',180);
    close(fh);
end
end
