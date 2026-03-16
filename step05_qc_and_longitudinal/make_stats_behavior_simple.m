function make_stats_behavior_simple()
% Simple behavior-only stats & plots for:
%   - Tail-immersion latency across *all* days available
%   - Hotplate latency (often 4 days)
%   - TST latency (often 4 days)
%
% Mice are grouped at the mouse-level:
%   Passive group  = mice that EVER had any passive session (isPassive==1)
%   Active group   = mice that NEVER had a passive session
%
% For each metric:
%   - per-day (median per mouse) -> group mean ± SEM trajectory
%   - per-day Active vs Passive Wilcoxon rank-sum (p shown even if n.s.)
%   - within-group Kruskal–Wallis across days (overall p in title)
%   - saves PNG + CSV to run_014/figs/stats_behavior_simple

%% ---- paths & read CSV ----
csvPath = 'rebuild_lick_patterns_with_ID\run_004\';
assert(exist(csvPath,'file')>0, 'CSV not found: %s', csvPath);
runDir = fileparts(fileparts(csvPath)); % ...\run_014
outDir = fullfile(runDir, 'figs', 'stats_behavior_simple');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('[behavior-simple] Reading: %s\n', csvPath);
T = readtable(csvPath, 'VariableNamingRule','preserve');

% ensure types
T = ensureString(T,'mouse_key');
if hasVar(T,'day_index') && ~isnumeric(T.day_index),   T.day_index   = double(T.day_index); end
if hasVar(T,'isPassive') && ~isnumeric(T.isPassive),    T.isPassive   = double(T.isPassive); end

%% ---- detect metric columns (first that exists) ----
vTail = firstHit(T, {'TailImmersion_latency_s','TailImmersion_s','Immersion_Latency_s'});
vHot  = firstHit(T, {'Hotplate_latency_s','HotPlate_latency_s'});
vTST  = firstHit(T, {'TST_latency_s','TST_s','TailSuspension_s'});

metrics = {};
labels  = {};
if ~isempty(vTail), metrics{end+1}=vTail; labels{end+1}='Tail immersion latency (s)'; end
if ~isempty(vHot),  metrics{end+1}=vHot;  labels{end+1}='Hotplate latency (s)';       end
if ~isempty(vTST),  metrics{end+1}=vTST;  labels{end+1}='TST latency (s)';            end
assert(~isempty(metrics), 'No Tail-immersion / Hotplate / TST columns found.');

%% ---- mouse-level group (EVER passive?) ----
if ~hasVar(T,'isPassive')
    error('isPassive column not found in CSV; cannot form groups.');
end
Gmouse = groupsummary(T(:,{'mouse_key','isPassive'}), 'mouse_key', @(x) any(x==1), 'isPassive');
Gmouse.Properties.VariableNames(end) = {'hasPassive'}; % logical

%% ---- per-day (median per mouse per day) for each metric ----
for i=1:numel(metrics)
    ycol = metrics{i}; ylbl = labels{i};
    fprintf('[behavior-simple] Metric: %s\n', ycol);

    % Reduce to needed columns and drop NaN rows for this metric
    keep = hasVar(T,ycol) & hasVar(T,'mouse_key') & hasVar(T,'day_index');
    assert(keep, 'Missing required columns.');
    sub  = T(:,{'mouse_key','day_index', ycol, 'isPassive'});
    sub.(ycol) = double(sub.(ycol));
    sub = sub(isfinite(sub.(ycol)) & isfinite(sub.day_index),:);
    if isempty(sub), warning('No data for %s; skipping.', ycol); continue; end

    % per mouse x day -> median
    D = groupsummary(sub, {'mouse_key','day_index'}, 'median', ycol);
    D.Properties.VariableNames(end) = {'y_med'};

    % attach mouse group (ever passive)
    D = outerJoinInto(D, Gmouse, {'mouse_key'}, 'hasPassive', 'mouse_hasPassive');
    D.mouse_hasPassive = double(D.mouse_hasPassive); % 1=Passive-group, 0=Active-group

    % split groups
    A = D(D.mouse_hasPassive==0,:); % active-only mice (all their days)
    P = D(D.mouse_hasPassive==1,:); % had-passive mice (all their days)

    %% ---- Figure: trajectories per group (mean ± SEM) ----
    f = figure('Color','w','Position',[100 100 1100 460]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % ACTIVE subplot
    nexttile; hold on
    plotMeanSEM_byDay(A, 'y_med', 'day_index');
    title(sprintf('ACTIVE group — %s', ylbl));
    xlabel('Day'); ylabel(ylbl); grid on; box on
    p_kw_A = KW_over_days(A, 'y_med', 'day_index');
    if isfinite(p_kw_A), subtitle(sprintf('Within-group (KW): p=%.3g', p_kw_A)); end

    % PASSIVE subplot
    nexttile; hold on
    plotMeanSEM_byDay(P, 'y_med', 'day_index');
    title(sprintf('PASSIVE group — %s', ylbl));
    xlabel('Day'); ylabel(ylbl); grid on; box on
    p_kw_P = KW_over_days(P, 'y_med', 'day_index');
    if isfinite(p_kw_P), subtitle(sprintf('Within-group (KW): p=%.3g', p_kw_P)); end

    % save fig
    savepng(f, fullfile(outDir, sprintf('TRAJ_%s.png', ycol)));
    close(f);

    %% ---- Between-group per day stats + overlay figure ----
    days = union(unique(A.day_index), unique(P.day_index));
    days = days(isfinite(days));
    if isempty(days), continue; end
    days = sort(days);

    % compute group means & SEM per day + p-values (ranksum) per day
    rows = {};
    f2 = figure('Color','w','Position',[100 100 780 460]); hold on
    for d = days'
        Xa = A.y_med(A.day_index==d); Xa = Xa(isfinite(Xa));
        Xp = P.y_med(P.day_index==d); Xp = Xp(isfinite(Xp));

        ma = mean(Xa,'omitnan'); ea = sem(Xa);
        mp = mean(Xp,'omitnan'); ep = sem(Xp);

        % plot with small horizontal offset
        plot([d-0.15 d-0.15],[ma-ea ma+ea],'k-','LineWidth',1.2);
        plot(d-0.15, ma, 'ko','MarkerFaceColor','k','MarkerSize',5);
        plot([d+0.15 d+0.15],[mp-ep mp+ep],'k-','LineWidth',1.2);
        plot(d+0.15, mp, 'ks','MarkerFaceColor','k','MarkerSize',5);

        % between-group p (always displayed)
        p = NaN;
        if ~isempty(Xa) && ~isempty(Xp)
            try, p = ranksum(Xa,Xp); catch, [~,p] = ttest2(Xa,Xp,'Vartype','unequal'); end
        end
        yl = ylim; yTop = yl(2) - 0.02*range(yl);
        text(d, yTop, sprintf('p=%.3g%s', p, sigStar(p,' prefix',' ')), ...
             'HorizontalAlignment','center','FontWeight','bold');

        rows(end+1,:) = {d, numel(Xa), ma, ea, numel(Xp), mp, ep, p}; %#ok<AGROW>
    end
    xlabel('Day'); ylabel(ylbl); grid on; box on
    title(sprintf('Active (circles) vs Passive (squares) — %s', ylbl));
    savepng(f2, fullfile(outDir, sprintf('BETWEEN_byDay_%s.png', ycol))); close(f2);

    % write CSV summary
    Tcsv = cell2table(rows, 'VariableNames', ...
        {'day','n_active','mean_active','sem_active','n_passive','mean_passive','sem_passive','p_ranksum'});
    try
        writetable(Tcsv, fullfile(outDir, sprintf('BETWEEN_byDay_%s.csv', ycol)));
    catch
    end
end

fprintf('[behavior-simple] Done. Outputs: %s\n', outDir);
end

%% -------------------- helpers --------------------
function tf = hasVar(T,var), tf = ismember(var, T.Properties.VariableNames); end
function T = ensureString(T,v), if hasVar(T,v) && ~isstring(T.(v)), T.(v)=string(T.(v)); end, end
function s = sem(x), x=x(isfinite(x)); if isempty(x), s=NaN; else, s=std(x,0,'omitnan')/sqrt(numel(x)); end, end

function nm = firstHit(T, list)
nm = '';
for i=1:numel(list)
    if ismember(list{i}, T.Properties.VariableNames)
        nm = list{i}; return;
    end
end
end

function D = outerJoinInto(D, G, keys, srcName, dstName)
if isempty(G), return; end
if ~ismember(dstName, D.Properties.VariableNames), D.(dstName)=nan(height(D),1); end
kD = buildKey(D, keys); kG = buildKey(G, keys);
[tf,loc] = ismember(kD, kG);
vals = G.(srcName); D.(dstName)(tf) = vals(loc(tf));
end

function k = buildKey(T, keys)
k = repmat("", height(T), 1);
for i=1:numel(keys)
    v = T.(keys{i}); if ~isstring(v), v = string(v); end
    if i==1, k=v; else, k = k + "|" + v; end
end
end

function p = KW_over_days(D, yname, dayname)
p = NaN;
if isempty(D), return; end
y = double(D.(yname)); g = double(D.(dayname));
good = isfinite(y) & isfinite(g);
y = y(good); g = g(good);
if numel(unique(g)) < 2 || numel(y) < 3, return; end
try
    p = kruskalwallis(y, g, 'off');
catch
    p = NaN;
end
end

function plotMeanSEM_byDay(D, yname, dayname)
if isempty(D), axis off; text(0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center'); return; end
ds = groupsummary(D, dayname, 'mean', yname);
es = groupsummary(D, dayname, @(x) std(x,'omitnan')/sqrt(sum(isfinite(x))), yname);
x = double(ds.(dayname)); 
m = double(ds.("mean_"+yname));
e = double(es.("fun1_"+yname)); e(~isfinite(e))=0;
if isempty(x), axis off; return; end
errorbar(x, m, e, 'k-o', 'LineWidth',1.4, 'MarkerFaceColor','k', 'MarkerSize',4);
end

function s = sigStar(p, varargin)
% sigStar(p, 'prefix',' ') -> returns ' ***' / ' **' / ' *' / '' depending on p
prefix = '';
if ~isempty(varargin)
    for k=1:2:numel(varargin), if strcmpi(varargin{k},'prefix'), prefix = varargin{k+1}; end, end
end
if ~isfinite(p), s=''; return; end
if p < 1e-3, s = [prefix '***'];
elseif p < 1e-2, s = [prefix ' **'];
elseif p < 0.05, s = [prefix '  *'];
else, s = '';
end
end

function savepng(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',150); catch, print(fh, fn, '-dpng','-r150'); end
end
