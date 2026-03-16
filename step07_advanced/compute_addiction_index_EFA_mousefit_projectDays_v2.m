function compute_addiction_index_EFA_mousefit_projectDays_v3()
% compute_addiction_index_EFA_mousefit_projectDays_v3
%
% FIXED vs v2:
%   - Forces D.mouse_key and mouseInfo.mouse_key to the SAME type (string) before join
%   - Avoids outerjoin categorical concat bug by joining into a clean key table
%   - Normalizes mouse_key formatting (\, /, - -> _ ; lowercase colors)
%   - Works whether D.mouse_key is categorical / cellstr / string
%
% Also still:
%   - Excludes day 1-2 everywhere
%   - MouseGroup assigned ONLY from your explicit cage/color mapping
%   - Phase-dependent stats within group + group comparison within phase with stars

%% ----------------- USER SETTINGS -----------------
runDir = 'K:\addiction_concate_Dec_2025\longitudinal_outputs\run_009';
outDir = fullfile(runDir,'figs','addiction_index_EFA');
if ~exist(outDir,'dir'), mkdir(outDir); end

cacheMat = fullfile(runDir,'S_D_cache.mat');
assert(exist(cacheMat,'file')>0,'Missing %s. Point to correct runDir or save D.', cacheMat);
S = load(cacheMat);
assert(isfield(S,'D'),'S_D_cache.mat must contain variable D');
D = S.D;

% Exclude habituation days
excludeDays = [1 2];

% Transition day handling (optional)
transitionDays = [3 6 11 14];
includeTransitionDays = true;  % set false to exclude transition days
if ~includeTransitionDays
    D = D(~ismember(double(D.day_index), transitionDays), :);
end

% Factor selection rule
targetExplainedCommon = 0.85;
maxFactorsHardCap = 12;

rotationMethod = 'varimax';
corrDupThr = 0.98;
nearConstVar = 1e-6;
maxMissingPerFeature = 0.40;
maxMissingPerMouse   = 0.40;

signAlignRef = 'RequirementLast';

%% ----------------- 0) Basic checks + normalize keys -----------------
reqCols = {'mouse_key','day_index'};
for i=1:numel(reqCols)
    assert(ismember(reqCols{i}, D.Properties.VariableNames), 'D missing required column: %s', reqCols{i});
end

% Force mouse_key to string and normalize format
D.mouse_key = normalize_mouse_key(D.mouse_key);

D.day_index = double(D.day_index);
D = D(~ismember(D.day_index, excludeDays), :); % HARD EXCLUDE day1-2

% Period assignment
D.Period = periodOfDay(D.day_index);
D = D(~isundefined(D.Period), :);

%% ----------------- 1) Assign MouseGroup + Sex + PairID from explicit mapping -----------------
mouseInfo = build_mouse_map();                 % mouse_key is string
mouseInfo.mouse_key = normalize_mouse_key(mouseInfo.mouse_key);

% Join safely without type-mismatch/categorical issues:
% Build a small mapping table with string keys and categorical values
mapTbl = mouseInfo(:, {'mouse_key','MouseGroup','Sex','PairID'});

% Ensure map columns are categorical (not string) to keep consistent
mapTbl.MouseGroup = categorical(string(mapTbl.MouseGroup), {'Active','Passive'});
mapTbl.Sex        = categorical(string(mapTbl.Sex), {'F','M'});

% Use ismember mapping (more robust than outerjoin for type quirks)
[tf, loc] = ismember(D.mouse_key, mapTbl.mouse_key);

if any(~tf)
    warnMice = unique(D.mouse_key(~tf));
    fprintf('WARNING: %d mice in D have no mapping. They will be DROPPED:\n', numel(warnMice));
    disp(warnMice);
end

% Attach columns
D.MouseGroup = categorical(repmat(missing, height(D), 1), {'Active','Passive'});
D.Sex        = categorical(repmat(missing, height(D), 1), {'F','M'});
D.PairID     = nan(height(D),1);

D.MouseGroup(tf) = mapTbl.MouseGroup(loc(tf));
D.Sex(tf)        = mapTbl.Sex(loc(tf));
D.PairID(tf)     = mapTbl.PairID(loc(tf));

% Drop unmapped rows
D = D(tf,:);

%% ----------------- 2) Features available in D -----------------
candidateFeatures = { ...
    'RequirementLast','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
    'bout_n','bout_meanDur_s','rew_freq_per_min','rew_meanDur_s','pupil_mean', ...
    'lick_totalDur_s','bout_totalDur_s','rew_totalDur_s', ...
    'Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min' ...
    };

features = candidateFeatures(ismember(candidateFeatures, D.Properties.VariableNames));
assert(numel(features) >= 4, 'Not enough feature columns found in D. Found %d.', numel(features));

%% ----------------- 3) Mouse-level aggregation (median across all days; day1-2 excluded) -----------------
mice = unique(D.mouse_key,'stable');
mouseTbl = table(mice, 'VariableNames', {'mouse_key'});

% Add group/sex/pair per mouse
mouseTbl.MouseGroup = categorical(strings(numel(mice),1), {'Active','Passive'});
mouseTbl.Sex        = categorical(strings(numel(mice),1), {'F','M'});
mouseTbl.PairID     = nan(numel(mice),1);

for i=1:numel(mice)
    rows = D.mouse_key==mice(i);
    mouseTbl.MouseGroup(i) = mode(D.MouseGroup(rows));
    mouseTbl.Sex(i)        = mode(D.Sex(rows));
    pr = D.PairID(rows); pr = pr(~isnan(pr));
    if ~isempty(pr), mouseTbl.PairID(i) = pr(1); end
end

% Feature medians per mouse
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
features = features(featKeep);
assert(numel(features) >= 4, 'After missingness filtering, not enough features for EFA.');

Xmouse = mouseTbl{:, features};
missPerMouse = mean(~isfinite(Xmouse),2);
keepMouse = missPerMouse <= maxMissingPerMouse;

if any(~keepMouse)
    fprintf('Dropping %d mice due to missingness > %.0f%%\n', nnz(~keepMouse), 100*maxMissingPerMouse);
end
mouseTbl = mouseTbl(keepMouse,:);
Xmouse   = mouseTbl{:, features};

% Impute remaining missing with feature medians
for j=1:size(Xmouse,2)
    col = Xmouse(:,j);
    med = median(col,'omitnan');
    col(~isfinite(col)) = med;
    Xmouse(:,j) = col;
end

% Standardize (store mu/sd for projecting day rows)
[XmouseZ, muFeat, sdFeat] = zscore_safe(Xmouse);

%% ----------------- 4) Robust pre-processing for EFA stability -----------------
nM = size(XmouseZ,1);
p0 = size(XmouseZ,2);
fprintf('\nMouse-level matrix: nMice=%d, pFeat=%d\n', nM, p0);

pCap = max(6, nM-3);
if p0 > pCap
    fprintf('Reducing features for stability: p=%d > pCap=%d\n', p0, pCap);
    [XmouseZ, features] = reduce_features_for_EFA(XmouseZ, features, pCap);
    % also shrink muFeat/sdFeat accordingly by refitting from Xmouse (safer)
    Xmouse = mouseTbl{:, features};
    [XmouseZ, muFeat, sdFeat] = zscore_safe(Xmouse);
end

[XmouseZ, features] = drop_near_constant(XmouseZ, features, nearConstVar);
[XmouseZ, features] = drop_high_corr(XmouseZ, features, corrDupThr);

% Recompute muFeat/sdFeat in the final feature set
Xmouse = mouseTbl{:, features};
[XmouseZ, muFeat, sdFeat] = zscore_safe(Xmouse);

p = size(XmouseZ,2);
maxFactors = min([maxFactorsHardCap, max(1,p-1)]);
fprintf('After filtering: nMice=%d, pFeat=%d, maxFactors=%d\n', size(XmouseZ,1), p, maxFactors);

%% ----------------- 5) Choose number of factors by >=85% common variance explained -----------------
best = struct('nf',[],'Loadings',[],'Psi',[],'stats',[],'T',[],'commonExpl',[]);
fprintf('\nSelecting # factors to reach >= %.0f%% common variance explained...\n', 100*targetExplainedCommon);

for nf = 1:maxFactors
    [ok, Lambda, Psi, Tmat, stats, failMsg] = try_factoran_robust(XmouseZ, nf, rotationMethod);
    if ~ok
        fprintf('  nf=%d FAILED: %s\n', nf, failMsg);
        continue
    end

    commonExpl = sum(1 - Psi) / numel(Psi);
    fprintf('  nf=%d OK -> common variance explained ~ %.1f%%\n', nf, 100*commonExpl);

    best.nf = nf;
    best.Loadings = Lambda;
    best.Psi = Psi;
    best.T = Tmat;
    best.stats = stats;
    best.commonExpl = commonExpl;

    if commonExpl >= targetExplainedCommon
        break
    end
end

assert(~isempty(best.nf), 'Could not fit EFA for any nf.');
nf     = best.nf;
Lambda = best.Loadings;
Psi    = best.Psi;

fprintf('Selected nf=%d (common variance explained ~ %.1f%%)\n', nf, 100*best.commonExpl);

%% ----------------- 6) Mouse-level factor scores + AI -----------------
Fmouse = score_factors_regression(XmouseZ, Lambda, Psi);

doSignAlign = ismember(signAlignRef, mouseTbl.Properties.VariableNames) && ~isempty(signAlignRef);
if doSignAlign
    ref = double(mouseTbl.(signAlignRef));
    refZ = zscore_safe(ref);
    for k=1:nf
        r = corr(Fmouse(:,k), refZ, 'rows','complete');
        if isfinite(r) && r < 0
            Fmouse(:,k) = -Fmouse(:,k);
            Lambda(:,k) = -Lambda(:,k);
        end
    end
end

FmouseZ = zscore_safe(Fmouse);
AI_mouse = mean(FmouseZ, 2, 'omitnan');
AI_mouse = rescale01(AI_mouse);
mouseTbl.AI_mouse = AI_mouse;

%% ----------------- 7) Project mouse×day -> AI per mouse×day -----------------
Dsub = D(:, intersect([{'mouse_key','day_index','MouseGroup','Period','Sex','PairID'}, features], ...
                      D.Properties.VariableNames, 'stable'));

Xday = double(Dsub{:, features});

muFeat = muFeat(:)'; sdFeat = sdFeat(:)';
for j=1:size(Xday,2)
    col = Xday(:,j);
    col(~isfinite(col)) = muFeat(j);
    Xday(:,j) = col;
end

XdayZ = (Xday - muFeat) ./ sdFeat;
XdayZ(~isfinite(XdayZ)) = 0;

Fday = score_factors_regression(XdayZ, Lambda, Psi);
FdayZ = zscore_safe(Fday);
AI_day = mean(FdayZ, 2, 'omitnan');
AI_day = rescale01(AI_day);
Dsub.AI = AI_day;

%% ----------------- 8) Plots -----------------
plot_spaghetti_byMouseGroup(Dsub, outDir);
plot_phase_withinGroup_stats(Dsub, outDir);     % stars on phase-comparison plots
plot_group_withinPhase_stats(Dsub, outDir);     % stars on group-within-phase plots

%% ----------------- 9) Save outputs -----------------
save(fullfile(outDir,'AI_results.mat'), ...
    'mouseTbl','Dsub','features','nf','Lambda','Psi','muFeat','sdFeat', ...
    'best','includeTransitionDays','transitionDays','targetExplainedCommon');

writetable(mouseTbl, fullfile(outDir,'AI_mouse_summary.csv'));
writetable(Dsub,     fullfile(outDir,'AI_mouse_day.csv'));

fprintf('Saved AI outputs to:\n  %s\n', outDir);

end

%% ===================== Key normalization =====================
function mk = normalize_mouse_key(x)
% Convert categorical/cellstr/string to string, normalize separators, lowercase color part.
if iscategorical(x)
    x = string(x);
elseif iscell(x)
    x = string(x);
else
    x = string(x);
end

x = replace(x,"\","_");
x = replace(x,"/","_");
x = replace(x,"-","_");
x = strtrim(x);

% If it looks like "6100_black", force color to lowercase
mk = x;
for i=1:numel(mk)
    parts = split(mk(i), "_");
    if numel(parts) >= 2
        parts(2) = lower(parts(2));
        mk(i) = strjoin(parts(1:2), "_");
    end
end
end

%% ===================== Cohort mapping =====================
function mouseInfo = build_mouse_map()
keys  = ["6100_red","6100_orange","6100_black", ...
         "0911_red","0911_orange","0911_black","0911_white", ...
         "0910_red","0910_orange","0910_black", ...
         "6099_red","6099_orange","6099_black","6099_white"];

Sex = ["F","F","F",  "F","F","F","F",  "M","M","M",  "M","M","M","M"];
Group = ["Passive","Passive","Active", ...
         "Active","Passive","Passive","Active", ...
         "Passive","Passive","Active", ...
         "Passive","Active","Active","Passive"];
PairID = [1,1,1, 2,2,3,3, 4,4,4, 5,5,6,6]';

mouseInfo = table(string(keys(:)), categorical(Sex(:),{'F','M'}), ...
                  categorical(Group(:),{'Active','Passive'}), PairID, ...
                  'VariableNames',{'mouse_key','Sex','MouseGroup','PairID'});
end

%% ===================== Period mapping =====================
function P = periodOfDay(d)
p = strings(size(d));
p(d>=3  & d<=5 )  = "Pre";
p(d>=6  & d<=10)  = "During";
p(d>=11 & d<=13)  = "Post";
p(d>=14 & d<=16)  = "Withdrawal";
p(d>=17 & d<=18)  = "Re-exposure";
P = categorical(p, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

%% ===================== Math helpers =====================
function [Z, mu, sd] = zscore_safe(X)
X = double(X);
if isvector(X), X = X(:); end
mu = median(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(sd==0 | ~isfinite(sd)) = 1;
Z = (X - mu) ./ sd;
Z(~isfinite(Z)) = 0;
end

function y = rescale01(x)
x = double(x);
finite = isfinite(x);
if ~any(finite)
    y = 0.5*ones(size(x));
    return
end
lo = min(x(finite));
hi = max(x(finite));
if hi==lo
    y = 0.5*ones(size(x));
else
    y = (x - lo) / (hi - lo);
end
end

function F = score_factors_regression(Xz, Lambda, Psi)
Xz = double(Xz);
Psi = double(Psi(:));
k = size(Lambda,2);
invPsi = diag(1 ./ max(Psi, 1e-6));
A = (Lambda' * invPsi * Lambda + eye(k));
W = (A \ (Lambda' * invPsi));
F = (W * Xz')';
F(~isfinite(F)) = 0;
end

%% ===================== Robust factoran wrapper =====================
function [ok, Lambda, Psi, Tmat, stats, msg] = try_factoran_robust(Xz, nf, rotationMethod)
ok = false; Lambda=[]; Psi=[]; Tmat=[]; stats=[]; msg = '';

try
    [Lambda,Psi,Tmat,stats] = factoran(Xz, nf, 'rotate', rotationMethod, 'scores','regression', 'maxit',3000);
    ok = true; return
catch ME
    msg1 = ME.message;
end

try
    Xj = Xz + 1e-6*randn(size(Xz));
    [Lambda,Psi,Tmat,stats] = factoran(Xj, nf, 'rotate', rotationMethod, 'scores','regression', 'maxit',3000);
    ok = true; return
catch ME
    msg2 = ME.message;
end

try
    R = corrcoef(Xz);
    R = (R+R')/2;
    ridge = 1e-3;
    R = R + ridge*eye(size(R,1));
    [Lambda,Psi,Tmat,stats] = factoran(R, nf, 'nobs', size(Xz,1), 'rotate', rotationMethod, 'scores','regression', 'maxit',3000);
    ok = true; return
catch ME
    msg3 = ME.message;
end

msg = sprintf('raw: %s | jitter: %s | R+ridge: %s', msg1, msg2, msg3);
end

function [X, features] = drop_near_constant(X, features, tolVar)
v = var(X,0,1);
keep = (v > tolVar) & isfinite(v);
if any(~keep)
    fprintf('Dropping %d near-constant features:\n', nnz(~keep));
    disp(features(~keep));
end
X = X(:,keep);
features = features(keep);
end

function [X, features] = drop_high_corr(X, features, thr)
R = corrcoef(X);
R(1:size(R,1)+1:end) = 0;
toDrop = false(1,size(X,2));
for j=1:size(R,1)
    if toDrop(j), continue; end
    dup = find(abs(R(j,:)) > thr);
    dup = dup(~toDrop(dup));
    if numel(dup) >= 1
        toDrop(dup) = true;
        toDrop(j) = false;
        dup(dup==j) = [];
        toDrop(dup) = true;
    end
end
if any(toDrop)
    fprintf('Dropping %d highly-correlated features (|r|>%.2f):\n', nnz(toDrop), thr);
    disp(features(toDrop));
end
X = X(:,~toDrop);
features = features(~toDrop);
end

function [X, features] = reduce_features_for_EFA(X, features, pCap)
v = var(X,0,1);
[~,ord] = sort(v,'descend','MissingPlacement','last');
ord = ord(1:pCap);
X = X(:,ord);
features = features(ord);
end

%% ===================== Plot helpers (with stars) =====================
function plot_spaghetti_byMouseGroup(Dsub, outDir)
fh = figure('Color','w','Position',[80 80 1020 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

groups = categories(Dsub.MouseGroup);
for gi=1:numel(groups)
    ax = nexttile; hold(ax,'on');
    G = groups{gi};
    Dg = Dsub(Dsub.MouseGroup==G,:);
    if isempty(Dg)
        title(ax, sprintf('%s',G));
        text(ax, 0.5, 0.5, 'No data for this MouseGroup','HorizontalAlignment','center');
        axis(ax,'off');
        continue
    end

    mice = unique(Dg.mouse_key,'stable');
    for i=1:numel(mice)
        r = Dg.mouse_key==mice(i);
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

    title(ax, sprintf('%s',G));
    xlabel(ax,'Day'); ylabel(ax,'AI (0–1)');
    ylim(ax,[0 1]);
    grid(ax,'off'); box(ax,'on');
end

printpng(fh, fullfile(outDir,'AI_spaghetti_byMouseGroup.png'));
close(fh);
end

function plot_phase_withinGroup_stats(Dsub, outDir)
groups = categories(Dsub.MouseGroup);
periods = categories(Dsub.Period);

fh = figure('Color','w','Position',[90 90 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi=1:numel(groups)
    ax = nexttile; hold(ax,'on');
    G = groups{gi};
    Dg = Dsub(Dsub.MouseGroup==G,:);
    if isempty(Dg)
        title(ax, sprintf('%s (mouse-level median per phase)',G));
        text(ax, 0.5, 0.5, 'No data','HorizontalAlignment','center');
        axis(ax,'off');
        continue
    end

    T = groupsummary(Dg(:,{'mouse_key','Period','AI'}), {'mouse_key','Period'}, 'median','AI');
    T.Properties.VariableNames{'median_AI'} = 'AI_med';

    ymax = 0;
    for pi=1:numel(periods)
        P = T(T.Period==periods{pi},:);
        if isempty(P), continue; end
        x = pi;
        boxchart(ax, repmat(x,height(P),1), P.AI_med, 'BoxWidth',0.28, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
        scatter(ax, repmat(x,height(P),1), P.AI_med, 32, 'k', 'filled', 'MarkerFaceAlpha',0.60);
        ymax = max(ymax, max(P.AI_med));
    end

    set(ax,'XTick',1:numel(periods),'XTickLabel',periods);
    xlim(ax,[0.5 numel(periods)+0.5]);
    ylim(ax,[0 1]);
    ylabel(ax,'AI (0–1)');
    title(ax, sprintf('%s (mouse-level median per phase)',G));
    grid(ax,'off'); box(ax,'on');

    % KW + posthoc stars
    y = []; g = [];
    for pi=1:numel(periods)
        P = T(T.Period==periods{pi},:);
        if isempty(P), continue; end
        y = [y; P.AI_med];
        g = [g; repmat(pi, height(P),1)];
    end

    pKW = NaN;
    try, pKW = kruskalwallis(y, g, 'off'); catch, end
    subtitle(ax, sprintf('Kruskal-Wallis p = %s', fmt_p(pKW)));

    if isfinite(pKW)
        try
            [~,~,st] = kruskalwallis(y, g, 'off');
            C = multcompare(st, 'CType','dunn-sidak', 'Display','off'); % [i j ... p]
            sigPairs = C(C(:,6) < 0.05, [1 2 6]);
            if ~isempty(sigPairs)
                draw_sig_brackets(ax, sigPairs, ymax);
            end
        catch
        end
    end
end

printpng(fh, fullfile(outDir,'AI_phase_comparison_byMouseGroup_withStars.png'));
close(fh);
end

function plot_group_withinPhase_stats(Dsub, outDir)
periods = categories(Dsub.Period);

T = groupsummary(Dsub(:,{'mouse_key','MouseGroup','Period','AI'}), ...
                 {'mouse_key','MouseGroup','Period'}, 'median','AI');
T.Properties.VariableNames{'median_AI'} = 'AI_med';

fh = figure('Color','w','Position',[90 90 1280 520]);
tiledlayout(1,numel(periods),'TileSpacing','compact','Padding','compact');

for pi=1:numel(periods)
    ax = nexttile; hold(ax,'on');
    P = T(T.Period==periods{pi},:);

    A = P(P.MouseGroup=="Active",:);
    B = P(P.MouseGroup=="Passive",:);

    if isempty(A) || isempty(B)
        title(ax, sprintf('%s\n(not enough data)', periods{pi}));
        axis(ax,'off');
        continue
    end

    x1 = 1; x2 = 2;
    boxchart(ax, repmat(x1,height(A),1), A.AI_med, 'BoxWidth',0.35, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
    boxchart(ax, repmat(x2,height(B),1), B.AI_med, 'BoxWidth',0.35, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
    scatter(ax, repmat(x1,height(A),1), A.AI_med, 34, 'k','filled','MarkerFaceAlpha',0.60);
    scatter(ax, repmat(x2,height(B),1), B.AI_med, 34, 'k','filled','MarkerFaceAlpha',0.60);

    p = NaN;
    try
        if height(A)>=2 && height(B)>=2
            p = ranksum(A.AI_med, B.AI_med);
        end
    catch
    end

    set(ax,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
    ylim(ax,[0 1]);
    ylabel(ax,'AI (0–1)');
    title(ax, sprintf('%s', periods{pi}));
    grid(ax,'off'); box(ax,'on');

    yMax = max([A.AI_med; B.AI_med]);
    yb = min(0.96, yMax + 0.06);
    plot(ax, [1 1 2 2], [yb-0.01 yb yb yb-0.01], 'k-', 'LineWidth',1.2);
    text(ax, 1.5, yb+0.02, stars_from_p(p), 'HorizontalAlignment','center', 'FontWeight','bold');
    text(ax, 1.5, yb-0.05, sprintf('p=%s', fmt_p(p)), 'HorizontalAlignment','center');
end

printpng(fh, fullfile(outDir,'AI_group_comparison_withinPhase_withStars.png'));
close(fh);
end

function draw_sig_brackets(ax, sigPairs, ymax)
sigPairs = sortrows(sigPairs, [1 2]);
h = 0.03;
base = min(0.90, ymax + 0.06);
level = 0;
for r=1:size(sigPairs,1)
    i = sigPairs(r,1);
    j = sigPairs(r,2);
    p = sigPairs(r,3);

    level = level + 1;
    y = min(0.96, base + (level-1)*h);

    plot(ax, [i i j j], [y-h/3 y y y-h/3], 'k-', 'LineWidth',1.0);
    text(ax, mean([i j]), y+0.01, stars_from_p(p), 'HorizontalAlignment','center', 'FontWeight','bold');
end
end

function s = stars_from_p(p)
if ~isfinite(p)
    s = 'n.s.';
elseif p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end

function t = fmt_p(p)
if ~isfinite(p)
    t = 'NaN';
elseif p < 1e-4
    t = sprintf('%.2e', p);
else
    t = sprintf('%.4g', p);
end
end

function printpng(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',180);
catch
    print(fh, fn, '-dpng','-r180');
end
end
