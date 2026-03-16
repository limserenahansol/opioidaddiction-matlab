function analyze_addiction_score_efa_dec3()
% Build addiction score (EFA) per mouse-day and plot trajectories + phase comparisons

%% ---------- locate latest run + load ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*')); assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix]   = max([d.datenum]);
runDir   = fullfile(d(ix).folder,d(ix).name);
csvPath  = fullfile(runDir,'ALL_mice_longitudinal.csv');
fprintf('Reading: %s\n', csvPath);

% Prefer cached per-day metrics from dashboard pipeline
cacheMat = fullfile(runDir, 'S_D_cache.mat');
if ~exist(cacheMat,'file')
    error('Missing cache %s. Run analyze_passive_active_dashboard_dec2 first.', cacheMat);
end
L = load(cacheMat,'D');
D = L.D;
if ~ismember('Group', D.Properties.VariableNames)
    D.Group = repmat(categorical("Unknown"), height(D), 1);
end

% Ensure Group is filled per mouse (avoid missing isPassive on some days)
if ismember('isPassive', D.Properties.VariableNames)
    mk = unique(D.mouse_key, 'stable');
    gAll = strings(height(D),1);
    for i=1:numel(mk)
        r = (D.mouse_key==mk(i));
        ip = D.isPassive(r);
        ip = ip(isfinite(ip));
        if any(ip==1)
            gAll(r) = "Passive";
        elseif any(ip==0)
            gAll(r) = "Active";
        else
            gAll(r) = "Unknown";
        end
    end
    D.Group = categorical(gAll);
end

%% ---------- output dir ----------
outDir = fullfile(runDir,'figs','addiction_score_efa');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---------- EFA settings ----------
maxFactors = 8;       % cap for factor count
nPerm      = 100;     % parallel analysis permutations
doShiftNonneg = true; % shift features to non-negative
doLog2p1     = true;  % log2(x+1) transform

%% ---------- assemble feature matrix ----------
% Use all numeric metrics except identifiers (exclude days 1-2)
exclude = ["day_index","session_idx","SessionMinutes","isPassive"];
exclude = [exclude, "mouse_key","day_name","Group","Period","Session_Paradigm"];
vn = string(D.Properties.VariableNames);
numMask = varfun(@isnumeric, D, 'OutputFormat','uniform');
useVars = vn(numMask);
useVars = setdiff(useVars, exclude, 'stable');

% Restrict to day >= 3 for addiction score modeling
D = D(D.day_index >= 3, :);

% Drop columns with too many missing values
X = D{:, useVars};
missRate = mean(~isfinite(X), 1);
keep = missRate <= 0.5;
useVars = useVars(keep);
X = X(:, keep);

if size(X,2) < 3
    error('Not enough numeric metrics for EFA. Found %d columns after filtering.', size(X,2));
end

% Fill missing with column means
for j=1:size(X,2)
    col = X(:,j);
    mu = mean(col(isfinite(col)));
    col(~isfinite(col)) = mu;
    X(:,j) = col;
end

% Shift to non-negative and log2(x+1)
for j=1:size(X,2)
    col = X(:,j);
    if doShiftNonneg
        minv = min(col);
        if isfinite(minv) && minv < 0
            col = col - minv;
        end
    end
    if doLog2p1
        col = log2(col + 1);
    end
    X(:,j) = col;
end

% Z-score globally
Xz = zscore(X);

% Guard against non‑positive‑definite covariance
[Xz, useVars] = dropBadEfaColumns(Xz, useVars);

% Choose number of factors by parallel analysis
nFactors = selectNFactorsParallel(Xz, maxFactors, nPerm);
if nFactors < 1
    error('Not enough independent columns for EFA after filtering.');
end

try
    [Lambda, Psi, T, stats, F] = factoran(Xz, nFactors, 'rotate','varimax', 'scores','regression');
catch
    nFactors = max(1, nFactors-1);
    [Lambda, Psi, T, stats, F] = factoran(Xz, nFactors);
end

varExp = sum(Lambda.^2, 1);
w = varExp ./ max(eps, sum(varExp));
addictionScoreZ = zscore(F * w(:));
addictionScore01 = rescale(addictionScoreZ, 0, 1);

factorScores = nan(height(D), maxFactors);
factorScores(:,1:nFactors) = F(:,1:nFactors);

factorSummary = table(nFactors, numel(useVars), sum(varExp), ...
    'VariableNames', {'n_factors','n_vars_used','varExp_sum'});

%% ---------- build output table ----------
out = D(:, {'mouse_key','day_index','day_name'});
if ismember('Group', D.Properties.VariableNames)
    out.Group = D.Group;
else
    out.Group = repmat("Unknown", height(D),1);
end
out.AddictionScore = addictionScore01;
out.AddictionScore_z = addictionScoreZ;
for k=1:maxFactors
    out.("Factor"+k) = factorScores(:,k);
end
if all(~isfinite(out.AddictionScore))
    warning('AddictionScore is all NaN. Check input metrics and EFA filtering.');
end
writetable(out, fullfile(outDir,'addiction_score_by_mouse_day.csv'));
writetable(factorSummary, fullfile(outDir,'addiction_score_factor_summary.csv'));
writeEfaWorkflowNotes(outDir, useVars, nFactors);

%% ---------- plots ----------
% Spaghetti across days (separate panels)
plotAddictionSpaghetti(out, outDir);
plotAddictionSpaghettiZ(out, outDir);

% Phase comparison (separate panels)
out.Period = periodOfDay(double(out.day_index));
plotAddictionByPhase(out, outDir);
plotAddictionByPhaseZ(out, outDir);
plotAddictionDeltaFromPre(out, outDir);
plotAddictionRawGroupCompare(out, outDir);
writeAddictionStats(out, outDir);

% PR correlation (if RequirementLast exists)
if ismember('RequirementLast', D.Properties.VariableNames)
    plotAddictionVsPR(out, D, outDir);
end

% Additional regressions (restricted to requested days)
if ismember('lick_totalDur_s', D.Properties.VariableNames)
    plotAddictionVsMetric(out, D, 'lick_totalDur_s', 'Lick total duration (s)', outDir);
end
if ismember('lick_freq_per_min', D.Properties.VariableNames)
    plotAddictionVsMetric(out, D, 'lick_freq_per_min', 'Lick frequency (/min)', outDir);
end
plotAddictionVsMetricBatch(out, D, outDir);

% Lick raster plots for selected days (8 passive + 6 active by default)
plotLickRasterPanels(csvPath, outDir, D);

fprintf('Addiction score outputs saved to:\n  %s\n', outDir);
end

%% ========================================================================
function plotLickRasterPanels(csvPath, outDir, D)
days = [4 7 9 12 15 18];
nPassive = 8; nActive = 6;

% Map mouse_key -> Group from D (per-day table)
mk = unique(D.mouse_key, 'stable');
grp = strings(numel(mk),1);
if ismember('Group', D.Properties.VariableNames)
    for i=1:numel(mk)
        r = (D.mouse_key==mk(i));
        if any(r)
            g = string(D.Group(find(r,1,'first')));
            if strlength(g)==0, g="Unknown"; end
            grp(i) = g;
        else
            grp(i) = "Unknown";
        end
    end
else
    % fallback: build Group from raw CSV isPassive if Group not present
    grp(:) = "Unknown";
end

passiveKeys = mk(grp=="Passive");
activeKeys  = mk(grp=="Active");
if isempty(passiveKeys) || isempty(activeKeys), return; end

passiveKeys = passiveKeys(1:min(nPassive, numel(passiveKeys)));
activeKeys  = activeKeys(1:min(nActive,  numel(activeKeys)));
orderedKeys = [passiveKeys; activeKeys];

% Load minimal columns from raw CSV
T = readtable(csvPath, 'VariableNamingRule','preserve');
if ~ismember('Lick_TTL', T.Properties.VariableNames) || ~ismember('day_index', T.Properties.VariableNames)
    return;
end

timeVar = pickTimeVar(T);
if isempty(timeVar), return; end

% Normalize types
if ~isstring(T.mouse_key), T.mouse_key = string(T.mouse_key); end

for di=1:numel(days)
    d = days(di);
    fig = figure('Color','w','Position',[80 80 1100 520]); hold on
    y = 0;
    for i=1:numel(orderedKeys)
        k = orderedKeys(i);
    r = (string(T.mouse_key)==string(k)) & (double(T.day_index)==d);
        t = double(T.(timeVar)(r));
        lk = double(T.Lick_TTL(r));
        if isempty(t), y = y + 1; continue; end
        [t,ord] = sort(t); lk = lk(ord);
        lickTimes = detectRisingEdgesSimple(t, lk);
        y = y + 1;
        if ~isempty(lickTimes)
            plot(lickTimes/60, y*ones(size(lickTimes)), 'k.', 'MarkerSize', 6);
        end
    end
    xlabel('Session time (min)'); ylabel('Mouse (row)');
    title(sprintf('Lick raster (Day %d): top Passive (%d) / bottom Active (%d)', ...
        d, numel(passiveKeys), numel(activeKeys)));
    grid off; box off
    exportgraphics(fig, fullfile(outDir, sprintf('LICK_raster_day%d.png', d)), 'Resolution', 220);
    close(fig);
end
end

function timeVar = pickTimeVar(T)
opts = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
timeVar = '';
for i=1:numel(opts)
    if ismember(opts{i}, T.Properties.VariableNames)
        timeVar = opts{i};
        return;
    end
end
end

function lickTimes = detectRisingEdgesSimple(t, x)
x = x(:);
x(~isfinite(x)) = 0;
x = x > 0.5;
dx = diff([false; x]);
idx = find(dx==1);
lickTimes = t(idx);
lickTimes = lickTimes(isfinite(lickTimes));
end

%% ========================================================================
function plotAddictionSpaghetti(out, outDir)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
cols = struct('Active',[0.9 0.3 0.2], 'Passive',[0.2 0.4 0.9], 'Unknown',[0.4 0.4 0.4]);
groups = ["Passive","Active"];

yl = [min(out.AddictionScore,[],'omitnan'), max(out.AddictionScore,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    if ~isfield(cols, g), cols.(g) = [0.4 0.4 0.4]; end
    mk = unique(out.mouse_key(string(out.Group)==g),'stable');
    for i=1:numel(mk)
        r = (out.mouse_key==mk(i)) & (string(out.Group)==g);
        [d,ord] = sort(out.day_index(r));
        y = double(out.AddictionScore(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',0.8,'MarkerSize',3,'Color',cols.(g), 'HandleVisibility','off');
    end
    dAll = unique(out.day_index(string(out.Group)==g)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(out.AddictionScore(string(out.Group)==g & out.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, '-', 'LineWidth', 2.5, 'Color', cols.(g));
xlabel('Day'); ylabel('Addiction score (0-1)');
    title(sprintf('Addiction score over days (%s)', g));
    if ~isempty(yl), ylim(yl); end
    grid off; box off
end

exportgraphics(fig, fullfile(outDir,'addiction_score_spaghetti.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionSpaghettiZ(out, outDir)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
cols = struct('Active',[0.9 0.3 0.2], 'Passive',[0.2 0.4 0.9], 'Unknown',[0.4 0.4 0.4]);
groups = ["Passive","Active"];

yl = [min(out.AddictionScore_z,[],'omitnan'), max(out.AddictionScore_z,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    if ~isfield(cols, g), cols.(g) = [0.4 0.4 0.4]; end
    mk = unique(out.mouse_key(string(out.Group)==g),'stable');
    for i=1:numel(mk)
        r = (out.mouse_key==mk(i)) & (string(out.Group)==g);
        [d,ord] = sort(out.day_index(r));
        y = double(out.AddictionScore_z(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',0.8,'MarkerSize',3,'Color',cols.(g), 'HandleVisibility','off');
    end
    dAll = unique(out.day_index(string(out.Group)==g)); dAll=sort(dAll);
    mu = nan(size(dAll));
    for j=1:numel(dAll)
        mu(j) = mean(double(out.AddictionScore_z(string(out.Group)==g & out.day_index==dAll(j))), 'omitnan');
    end
    plot(dAll, mu, '-', 'LineWidth', 2.5, 'Color', cols.(g));
    xlabel('Day'); ylabel('Addiction score (z)');
    title(sprintf('Addiction score over days (%s) [z]', g));
    if ~isempty(yl), ylim(yl); end
    grid off; box off
end

exportgraphics(fig, fullfile(outDir,'addiction_score_spaghetti_z.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionByPhase(out, outDir)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
yl = [min(out.AddictionScore,[],'omitnan'), max(out.AddictionScore,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    rG = string(out.Group)==g & ~isundefined(categorical(out.Period));
    if ~any(rG), continue; end
    x = double(categorical(out.Period(rG), periods, 'Ordinal',true));
    y = double(out.AddictionScore(rG));
    b = boxchart(x, y, 'BoxWidth',0.25, 'MarkerStyle','none','BoxFaceAlpha',0.25);
    if g=="Passive", b.BoxFaceColor=[0.2 0.4 0.9]; else, b.BoxFaceColor=[0.9 0.3 0.2]; end
    set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Period'); ylabel('Addiction score (0-1)');
    title(sprintf('Addiction score by phase (%s)', g));
    if ~isempty(yl), ylim(yl); end
    addAddictionSignificanceLines(gca, out, g, periods);
    grid off; box off
end

exportgraphics(fig, fullfile(outDir,'addiction_score_by_phase.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionByPhaseZ(out, outDir)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
yl = [min(out.AddictionScore_z,[],'omitnan'), max(out.AddictionScore_z,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    rG = string(out.Group)==g & ~isundefined(categorical(out.Period));
    if ~any(rG), continue; end
    x = double(categorical(out.Period(rG), periods, 'Ordinal',true));
    y = double(out.AddictionScore_z(rG));
    b = boxchart(x, y, 'BoxWidth',0.25, 'MarkerStyle','none','BoxFaceAlpha',0.25);
    if g=="Passive", b.BoxFaceColor=[0.2 0.4 0.9]; else, b.BoxFaceColor=[0.9 0.3 0.2]; end
    set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
    xlabel('Period'); ylabel('Addiction score (z)');
    title(sprintf('Addiction score by phase (%s) [z]', g));
    if ~isempty(yl), ylim(yl); end
    grid off; box off
end

exportgraphics(fig, fullfile(outDir,'addiction_score_by_phase_z.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionDeltaFromPre(out, outDir)
% Change-from-Pre plot with between-group stats per period
fig = figure('Color','w','Position',[80 80 1000 520]); hold on
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
colors = struct('Passive',[0.2 0.4 0.9],'Active',[0.9 0.3 0.2]);

% baseline per mouse (Pre)
base = groupsummary(out(out.Period=="Pre",:), 'mouse_key', 'median', 'AddictionScore');
base.Properties.VariableNames{'median_AddictionScore'} = 'base';
T = outerjoin(out, base(:,{'mouse_key','base'}), 'Keys','mouse_key','MergeKeys',true);
T.delta = T.AddictionScore - T.base;

off = [-0.18 0.18];
for gi=1:numel(groups)
    g = groups(gi);
    for pi=1:numel(periods)
        r = string(T.Group)==g & T.Period==periods{pi} & isfinite(T.delta);
        if ~any(r), continue; end
        x = pi + off(gi);
        b = boxchart(repmat(x, sum(r),1), T.delta(r), 'BoxWidth',0.25, ...
            'MarkerStyle','none','BoxFaceAlpha',0.25);
        b.BoxFaceColor = colors.(g);
        scatter(repmat(x, sum(r),1), T.delta(r), 22, 'k', 'filled', 'MarkerFaceAlpha',0.7);
    end
end

% Between-group stats per period
yl = ylim; pad = range(yl)*0.08;
for pi=1:numel(periods)
    a = T.delta(string(T.Group)=="Active" & T.Period==periods{pi});
    p = T.delta(string(T.Group)=="Passive" & T.Period==periods{pi});
    a = a(isfinite(a)); p = p(isfinite(p));
    if isempty(a) || isempty(p), continue; end
    pval = ranksum(a,p);
    y = yl(2) - pad;
    plot([pi-0.15 pi+0.15],[y y],'k-');
    text(pi, y+0.02*range(yl), sprintf('%s (p=%.3g)', pvalStars(pval), pval), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Period'); ylabel('Addiction score – Pre');
title('Change from Pre (per mouse)');
legend(groups, 'Location','northwest');
grid off; box off
exportgraphics(fig, fullfile(outDir,'addiction_score_delta_from_pre.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionRawGroupCompare(out, outDir)
% Raw AddictionScore by period: Active vs Passive with p-values
fig = figure('Color','w','Position',[80 80 1000 520]); hold on
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
colors = struct('Passive',[0.2 0.4 0.9],'Active',[0.9 0.3 0.2]);

off = [-0.18 0.18];
for gi=1:numel(groups)
    g = groups(gi);
    for pi=1:numel(periods)
        r = string(out.Group)==g & out.Period==periods{pi} & isfinite(out.AddictionScore);
        if ~any(r), continue; end
        x = pi + off(gi);
        b = boxchart(repmat(x, sum(r),1), out.AddictionScore(r), 'BoxWidth',0.25, ...
            'MarkerStyle','none','BoxFaceAlpha',0.25);
        b.BoxFaceColor = colors.(g);
        scatter(repmat(x, sum(r),1), out.AddictionScore(r), 22, 'k', 'filled', 'MarkerFaceAlpha',0.7);
    end
end

yl = ylim; pad = range(yl)*0.08;
for pi=1:numel(periods)
    a = out.AddictionScore(string(out.Group)=="Active" & out.Period==periods{pi});
    p = out.AddictionScore(string(out.Group)=="Passive" & out.Period==periods{pi});
    a = a(isfinite(a)); p = p(isfinite(p));
    if isempty(a) || isempty(p), continue; end
    pval = ranksum(a,p);
    y = yl(2) - pad;
    plot([pi-0.15 pi+0.15],[y y],'k-');
    text(pi, y+0.02*range(yl), sprintf('%s (p=%.3g)', pvalStars(pval), pval), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Period'); ylabel('Addiction score (0-1)');
title('Addiction score (raw) by period: Active vs Passive');
legend(groups, 'Location','northwest');
grid off; box off
exportgraphics(fig, fullfile(outDir,'addiction_score_raw_by_period.png'), 'Resolution', 220);
close(fig);
end

function addAddictionSignificanceLines(ax, out, groupName, periods)
if nargin < 4, return; end
rG = string(out.Group)==groupName & ~isundefined(categorical(out.Period));
if ~any(rG), return; end

yl = ylim(ax);
pad = range(yl) * 0.04;
y = yl(2);

% Kruskal-Wallis across periods (title only)
try
    pKw = kruskalwallis(out.AddictionScore(rG), out.Period(rG), 'off');
    title(ax, sprintf('%s (KW p=%.3g %s)', ax.Title.String, pKw, pvalStars(pKw)));
catch
end

% Pairwise ranksum with lines
for i=1:numel(periods)
    for j=i+1:numel(periods)
        a = out.AddictionScore(rG & out.Period==periods{i});
        b = out.AddictionScore(rG & out.Period==periods{j});
        a = a(isfinite(a)); b = b(isfinite(b));
        if isempty(a) || isempty(b), continue; end
        p = ranksum(a, b);
        y = y + pad;
        plot(ax, [i i j j], [y-pad/2 y y y-pad/2], 'k-','LineWidth',1);
        text(ax, mean([i j]), y+pad*0.1, sprintf('%s (p=%.3g)', pvalStars(p), p), ...
            'HorizontalAlignment','center','FontWeight','bold');
    end
end
ylim(ax, [yl(1) y+pad]);
end

function writeAddictionStats(out, outDir)
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
fid = fopen(fullfile(outDir,'addiction_score_stats.txt'),'w');
fprintf(fid,'Addiction score stats\n');

% Within-group across periods
for gi=1:numel(groups)
    g = groups(gi);
    rG = string(out.Group)==g & ~isundefined(categorical(out.Period));
    if ~any(rG), continue; end
    try
        pKw = kruskalwallis(out.AddictionScore(rG), out.Period(rG), 'off');
    catch
        pKw = NaN;
    end
    fprintf(fid,'Within %s: KW p=%.4g\n', g, pKw);
    for i=1:numel(periods)
        for j=i+1:numel(periods)
            a = out.AddictionScore(rG & out.Period==periods{i});
            b = out.AddictionScore(rG & out.Period==periods{j});
            a = a(isfinite(a)); b = b(isfinite(b));
            if isempty(a) || isempty(b), continue; end
            p = ranksum(a, b);
            fprintf(fid,'  %s vs %s: p=%.4g\n', periods{i}, periods{j}, p);
        end
    end
end

% Between-group within each period
for pi=1:numel(periods)
    pName = periods{pi};
    a = out.AddictionScore(string(out.Group)=="Active" & out.Period==pName);
    b = out.AddictionScore(string(out.Group)=="Passive" & out.Period==pName);
    a = a(isfinite(a)); b = b(isfinite(b));
    if isempty(a) || isempty(b), continue; end
    p = ranksum(a, b);
    fprintf(fid,'Between groups (%s): p=%.4g\n', pName, p);
end
fclose(fid);
end

function s = pvalStars(p)
if ~isfinite(p)
    s = '';
elseif p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = '';
end
end

function plotAddictionVsPR(out, D, outDir)
% Correlate PR with addiction score (separate panels)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
groups = ["Passive","Active"];
yl = [min(out.AddictionScore,[],'omitnan'), max(out.AddictionScore,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end
allowedDays = [4 5 7 8 9 10 12 13 15 16 17 18];

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    rG = string(out.Group)==g & ismember(double(out.day_index), allowedDays);
    if g == "Passive"
        rG = rG & ~ismember(double(out.day_index), 6:10);
    end
    pr = D.RequirementLast(rG);
    as = out.AddictionScore(rG);
    mask = isfinite(pr) & isfinite(as);
    if any(mask)
        scatter(pr(mask), as(mask), 28, 'filled');
        lsline;
        [r,p] = corr(pr(mask), as(mask), 'Rows','complete');
        title(sprintf('PR vs Addiction score (%s): r=%.2f, p=%.3g', g, r, p));
    else
        title(sprintf('PR vs Addiction score (%s): n/a', g));
    end
    xlabel('PR score (RequirementLast)'); ylabel('Addiction score (0-1)');
    if ~isempty(yl), ylim(yl); end
    grid off; box off
end

exportgraphics(fig, fullfile(outDir,'addiction_score_vs_PR.png'), 'Resolution', 220);
close(fig);
end

function plotAddictionVsMetric(out, D, ycol, ylab, outDir)
fig = figure('Color','w','Position',[80 80 1100 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
groups = ["Passive","Active"];
yl = [min(out.AddictionScore,[],'omitnan'), max(out.AddictionScore,[],'omitnan')];
if any(~isfinite(yl)), yl = []; end
allowedDays = [4 5 7 8 9 10 12 13 15 16 17 18];

for gi=1:numel(groups)
    g = groups(gi);
    nexttile; hold on
    rG = string(out.Group)==g & ismember(double(out.day_index), allowedDays);
    x = D.(ycol)(rG);
    y = out.AddictionScore(rG);
    mask = isfinite(x) & isfinite(y);
    if any(mask)
        scatter(x(mask), y(mask), 28, 'filled');
        lsline;
        [r,p] = corr(x(mask), y(mask), 'Rows','complete');
        title(sprintf('%s vs Addiction score (%s): r=%.2f, p=%.3g', ylab, g, r, p));
    else
        title(sprintf('%s vs Addiction score (%s): n/a', ylab, g));
    end
    xlabel(ylab); ylabel('Addiction score (0-1)');
    if ~isempty(yl), ylim(yl); end
    grid off; box off
end

fname = sprintf('addiction_score_vs_%s.png', regexprep(ycol,'[^a-zA-Z0-9]+','_'));
exportgraphics(fig, fullfile(outDir, fname), 'Resolution', 220);
close(fig);
end

function plotAddictionVsMetricBatch(out, D, outDir)
list = {
  'RequirementLast','RequirementLast'
  'lick_n','lick_n'
  'lick_freq_per_min','lick_freq_per_min'
  'lick_freq_per_10s','lick_freq_per_10s'
  'lick_meanDur_s','lick_meanDur_s'
  'lick_totalDur_s','lick_totalDur_s'
  'lick_medianIEI_s','lick_medianIEI_s'
  'bout_n','bout_n'
  'bout_freq_per_10s','bout_freq_per_10s'
  'bout_freq_per_min','bout_freq_per_min'
  'bout_meanDur_s','bout_meanDur_s'
  'bout_totalDur_s','bout_totalDur_s'
  'rew_n','rew_n'
  'rew_freq_per_min','rew_freq_per_min'
  'rew_meanDur_s','rew_meanDur_s'
  'rew_totalDur_s','rew_totalDur_s'
  'rew_medianIRI_s','rew_medianIRI_s'
  'pupil_mean','pupil_mean'
  'Immersion_Latency_s','Immersion_Latency_s'
  'Requirement_cum','Requirement_cum'
  'Requirement_speed_per_day','Requirement_speed_per_day'
  'Requirement_speed_per_min','Requirement_speed_per_min'
};

for i=1:size(list,1)
    ycol = list{i,1};
    ylab = list{i,2};
    if ismember(ycol, D.Properties.VariableNames)
        plotAddictionVsMetric(out, D, ycol, ylab, outDir);
    end
end
end

function P = periodOfDay(d)
p = strings(size(d));
p(d>=3  & d<=5 )  = "Pre";
p(d>=6  & d<=10)  = "During";
p(d>=11 & d<=13)  = "Post";
p(d>=14 & d<=16)  = "Withdrawal";
p(d>=17 & d<=18)  = "Re-exposure";
P = categorical(p, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

%% ========================================================================
function writeEfaWorkflowNotes(outDir, useVars, nFactors)
fid = fopen(fullfile(outDir,'addiction_score_EFA_workflow.txt'),'w');
fprintf(fid,'Addiction Score (EFA) workflow\n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'Inputs: per-mouse per-day metrics from S_D_cache.mat\n');
fprintf(fid,'Parameters used (raw count): %d\n', numel(useVars));
fprintf(fid,'Variables:\n');
for i=1:numel(useVars)
    fprintf(fid,'  - %s\n', useVars(i));
end
fprintf(fid,'\nPer-day computation:\n');
fprintf(fid,'\nGlobal computation (all mouse-days together):\n');
fprintf(fid,'  1) Fill missing values with column means.\n');
fprintf(fid,'  2) Shift features to non-negative and apply log2(x+1).\n');
fprintf(fid,'  3) Z-score features across all mouse-days.\n');
fprintf(fid,'  4) Drop near-constant / collinear features.\n');
fprintf(fid,'  5) Choose nFactors by parallel analysis (cap=%d).\n', nFactors);
fprintf(fid,'  6) Run EFA (varimax).\n');
fprintf(fid,'  7) Compute factor scores per mouse-day.\n');
fprintf(fid,'  8) Weight factor scores by factor variance and sum to AddictionScore.\n');
fprintf(fid,'  9) Z-score AddictionScore across all mouse-days.\n');
fprintf(fid,' 10) Rescale AddictionScore to 0–1 for plotting/output.\n');
fprintf(fid,'\nOutputs:\n');
fprintf(fid,'  - addiction_score_by_mouse_day.csv\n');
fprintf(fid,'  - addiction_score_factor_summary.csv\n');
fprintf(fid,'  - addiction_score_spaghetti.png\n');
fprintf(fid,'  - addiction_score_by_phase.png\n');
fprintf(fid,'  - addiction_score_vs_PR.png (if RequirementLast exists)\n');
fclose(fid);
end

%% ========================================================================
function nFactors = selectNFactorsParallel(Xz, maxFactors, nPerm)
% Parallel analysis on correlation matrix
if nargin < 3, nPerm = 100; end
if nargin < 2, maxFactors = 8; end

R = corr(Xz, 'rows','pairwise');
eigReal = sort(eig(R), 'descend');
eigReal = eigReal(:);

eigRand = zeros(numel(eigReal), nPerm);
for i=1:nPerm
    Xp = Xz;
    for j=1:size(Xz,2)
        Xp(:,j) = Xp(randperm(size(Xz,1)), j);
    end
    Rp = corr(Xp, 'rows','pairwise');
    eigRand(:,i) = sort(eig(Rp), 'descend');
end

eigP95 = prctile(eigRand', 95)';
keep = eigReal > eigP95;
nFactors = sum(keep);
nFactors = min([nFactors, maxFactors, size(Xz,2)-1, size(Xz,1)-1]);
if nFactors < 1, nFactors = 1; end
end

%% ========================================================================
function [Xz, useVars] = dropBadEfaColumns(Xz, useVars)
% Drop constant / near-constant columns
v = var(Xz, 0, 1, 'omitnan');
keep = v > 1e-8;
Xz = Xz(:, keep);
useVars = useVars(keep);

% Drop perfectly collinear columns using QR
if size(Xz,2) > 1
    [~, R, E] = qr(Xz, 0);
    tol = 1e-10;
    r = sum(abs(diag(R)) > tol);
    keepIdx = sort(E(1:r));
    Xz = Xz(:, keepIdx);
    useVars = useVars(keepIdx);
end

% Ensure covariance is positive definite by diagonal jitter if needed
try
    [~, p] = chol(cov(Xz, 'omitrows'));
catch
    p = 1;
end
if p ~= 0
    Xz = Xz + 1e-6*randn(size(Xz));
end
end
