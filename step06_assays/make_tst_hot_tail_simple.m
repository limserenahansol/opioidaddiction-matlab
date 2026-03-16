function make_tst_hot_tail_simple()
% Simple analysis for Tail Immersion, TST, Hotplate.
% - Splits mice: Active-only vs Had-passive (ever isPassive==1)
% - Day-by-day Active vs Passive (ranksum), p printed even if n.s.
% - Within-group across-day one-way ANOVA + Tukey (fallback KW + BH-FDR)
% - Saves figs + CSVs to latest run folder under TST_HOT_TAIL

%% locate latest run + read CSV
tryPath = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
assert(exist(tryPath,'dir')>0,'Cannot find longitudinal_outputs');
d = dir(fullfile(tryPath,'run_*')); assert(~isempty(d),'No run_* folders found.');
[~,idx] = max([d.datenum]); runDir = fullfile(d(idx).folder, d(idx).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('[TST/HOT/TAIL] Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% basic types
if ~isstring(T.mouse_key), T.mouse_key = string(T.mouse_key); end
if ~isstring(T.day_name),  T.day_name  = string(T.day_name);  end
if ismember('day_index',T.Properties.VariableNames) && ~isnumeric(T.day_index)
    T.day_index = double(T.day_index);
end
if ismember('isPassive',T.Properties.VariableNames) && ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% output dir
outDir = fullfile(runDir,'TST_HOT_TAIL'); if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('[TST/HOT/TAIL] Saving to: %s\n', outDir);

% exact column names you specified
colTail = pickVar(T, {'Immersion_Latency_s'});
colTST  = pickVar(T, {'TST_Frames_Non_moving'});
colHOT  = pickVar(T, {'HOT_Frames_Non_moving'});
if isempty(colTail) && isempty(colTST) && isempty(colHOT)
    error('None of Immersion_Latency_s, TST_Frames_Non_moving, HOT_Frames_Non_moving found.');
end

%% mouse grouping: EVER passive?
if ismember('isPassive',T.Properties.VariableNames)
    G = groupsummary(T(:,{'mouse_key','isPassive'}),'mouse_key',@(x)any(x==1),'isPassive');
    G.Properties.VariableNames(end) = {'hadPassive'};  % 1 if mouse had passive at least once
else
    warning('isPassive not in CSV; treating all mice as Active-only.');
    G = unique(T(:,{'mouse_key'})); G.hadPassive = zeros(height(G),1);
end

%% build per-mouse-per-day medians for only needed columns
need = [{'mouse_key','day_index','day_name'}, colTail, colTST, colHOT];
need = need(~cellfun(@isempty,need));
need = need(ismember(need, T.Properties.VariableNames));
D = T(:,need);

keep = false(height(D),1);
if ~isempty(colTail), keep = keep | isfinite(double(D.(colTail))); end
if ~isempty(colTST),  keep = keep | isfinite(double(D.(colTST)));  end
if ~isempty(colHOT),  keep = keep | isfinite(double(D.(colHOT)));  end
D = D(keep,:);

K = unique(D(:,{'mouse_key','day_index','day_name'}),'rows','stable');
varsToFill = setdiff(D.Properties.VariableNames, {'mouse_key','day_index','day_name'});
for v = string(varsToFill)
    K.(v) = nan(height(K),1);
    for i=1:height(K)
        r = D.mouse_key==K.mouse_key(i) & D.day_index==K.day_index(i);
        K.(v)(i) = median(double(D.(v)(r)),'omitnan');
    end
end
K = outerJoinInto(K, G, {'mouse_key'}, 'hadPassive', 'hadPassive');
K.hadPassive = double(K.hadPassive); % 1=Had-passive, 0=Active-only

%% analyze each metric
if ~isempty(colTail)
    analyzeTrajectory(K, colTail, 'Tail immersion latency (s)', outDir, 'TAIL');
end
if ~isempty(colTST)
    analyzeTrajectory(K, colTST, 'TST non-moving (frames)', outDir, 'TST');
end
if ~isempty(colHOT)
    analyzeTrajectory(K, colHOT, 'Hotplate non-moving (frames)', outDir, 'HOT');
end

fprintf('[TST/HOT/TAIL] Done.\n');
end

%% ===== Helpers =====
function varname = pickVar(T, candidates)
varname=''; for i=1:numel(candidates)
    if ismember(candidates{i},T.Properties.VariableNames), varname=candidates{i}; return; end
end, end

function analyzeTrajectory(K, yvar, ylab, outDir, tag)
% Day-by-day Active vs Passive + within-group across-days ANOVA/Tukey
M = K(:,{'mouse_key','day_index','hadPassive'}); M.y = double(K.(yvar));
M = M(isfinite(M.y),:);
if isempty(M), warning('%s: no data', yvar); return; end
days = unique(M.day_index); days = days(:)';

% Trajectory plot (per-group medians + IQR)
fig = figure('Color','w','Position',[120 120 900 520]); hold on
colors = [0 0 0; 0.35 0.35 0.35]; names = {'Active-only','Had-passive'};
for g=0:1
    Sg = M(M.hadPassive==g,:);
    if isempty(Sg), continue; end
    mice = unique(Sg.mouse_key);
    for i=1:numel(mice)
        r = Sg.mouse_key==mice(i);
        [dd,ord] = sort(double(Sg.day_index(r))); yy = double(Sg.y(r)); yy=yy(ord);
        plot(dd, yy, 'o-','Color',[colors(g+1,:) 0.25],'MarkerFaceColor','w','LineWidth',0.8,'MarkerSize',4);
    end
    med = nan(size(days)); q1 = med; q3 = med;
    for k=1:numel(days)
        v = Sg.y(Sg.day_index==days(k));
        med(k)=median(v,'omitnan'); q1(k)=quantile(v,0.25); q3(k)=quantile(v,0.75);
    end
    fill([days fliplr(days)], [q1 fliplr(q3)], colors(g+1,:), 'FaceAlpha',0.12,'EdgeColor','none');
    plot(days, med, '-','Color',colors(g+1,:),'LineWidth',3);
end
grid on; box on; xlabel('Day'); ylabel(ylab);
title(sprintf('%s — daily trajectory (by group)', ylab));
legend(names,'Location','best'); legend boxoff
savepng(fig, fullfile(outDir, sprintf('%s_daily_%s.png', tag, yvar))); close(fig);

% Day-by-day boxplots & p-values (ranksum)
fig = figure('Color','w','Position',[100 100 1200 520]);
tlo = tiledlayout(1, numel(days), 'TileSpacing','compact','Padding','compact');
rows = {};
for k=1:numel(days)
    nexttile; hold on
    Xa = M.y(M.day_index==days(k) & M.hadPassive==0);
    Xp = M.y(M.day_index==days(k) & M.hadPassive==1);
    plotTwoBox(Xa, Xp, {'Active','Passive'}, ylab);
    title(sprintf('Day %d', days(k)));
    p = NaN; if ~isempty(Xa) && ~isempty(Xp)
        try, p = ranksum(Xa,Xp); catch, [~,p]=ttest2(Xa,Xp,'Vartype','unequal'); end
    end
    yl = ylim; addPBar(gca,1,2,yl(2),p,'Active vs Passive');   % always show p
    rows(end+1,:) = {days(k), numel(Xa), numel(Xp), median(Xa,'omitnan'), median(Xp,'omitnan'), p}; %#ok<AGROW>
end
title(tlo, sprintf('%s — Active vs Passive (per day)', ylab));
savepng(fig, fullfile(outDir, sprintf('%s_box_byDay_%s.png', tag, yvar))); close(fig);
Tcsv = cell2table(rows,'VariableNames',{'day','n_active','n_passive','median_active','median_passive','p_ranksum'});
try, writetable(Tcsv, fullfile(outDir, sprintf('%s_groupCompare_%s.csv', tag, yvar))); catch, end

% Within-group across-days ANOVA (+ Tukey) per group
for g=0:1
    Sg = M(M.hadPassive==g,:); if height(Sg)<2, continue; end
    X = double(Sg.y); G = double(Sg.day_index);
    usedKW = false; p_overall = NaN; pairs=[]; padj=[]; gnames=[];
    try
        [p_overall,~,st] = anova1(X, categorical(G), 'off');
        [mc,~,~,gn] = multcompare(st, 'CType','hsd','Display','off');
        pairs = mc(:,1:2); padj = mc(:,6); gnames = gn;
    catch
        usedKW = true;
        try
            p_overall = kruskalwallis(X, categorical(G), 'off');
            [pairs, padj, gnames] = fallbackPairsFDR(X, categorical(G));
        catch, p_overall = NaN; end
    end

    % Boxplot across days + significant pairs
    fig = figure('Color','w','Position',[110 110 980 520]); hold on
    boxplot(X, categorical(G), 'Symbol','k+','Colors','k'); grid on; box on
    ylabel(ylab); xlabel('Day');
    title(sprintf('%s — one-way %s across days (%s)  p=%.3g', ...
          ylab, iff(usedKW,'KW','ANOVA'), iff(g==0,'ACTIVE','PASSIVE'), p_overall));
    yl = ylim; yBase=yl(2); dy=0.05*range(yl); level=0;
    sigIdx = find(padj < 0.05);
    for r = reshape(sigIdx,1,[])
        a=pairs(r,1); b=pairs(r,2); yTop = yBase + (level+1)*dy;
        addSigBar(gca, a, b, yTop, padj(r), sprintf('D%s vs D%s', string(gnames{a}), string(gnames{b})));
        level = level + 1;
    end
    if ~isempty(sigIdx), ylim([yl(1), yBase + (level+3)*dy]); end
    savepng(fig, fullfile(outDir, sprintf('%s_acrossDays_%s_%s.png', tag, iff(g==0,'ACTIVE','PASSIVE'), yvar))); close(fig);

    % write CSV with all pairwise p
    if ~isempty(padj)
        pairStr = strings(numel(padj),1);
        for r=1:numel(padj)
            pairStr(r) = sprintf('D%s vs D%s', string(gnames{pairs(r,1)}), string(gnames{pairs(r,2)}));
        end
        Tmc = table(repmat(string(iff(g==0,'ACTIVE','PASSIVE')),numel(padj),1), pairStr, padj, ...
                    'VariableNames',{'group','pair','p_adj'});
        To = table(string(iff(g==0,'ACTIVE','PASSIVE')), usedKW, p_overall, ...
                   'VariableNames',{'group','used_KW','p_overall'});
        try
            writetable(To, fullfile(outDir, sprintf('%s_acrossDays_%s_%s_overall.csv', tag, iff(g==0,'ACTIVE','PASSIVE'), yvar)));
            writetable(Tmc, fullfile(outDir, sprintf('%s_acrossDays_%s_%s_pairs.csv',   tag, iff(g==0,'ACTIVE','PASSIVE'), yvar)));
        catch, end
    end
end
end

function [pairs, p_adj, names] = fallbackPairsFDR(X, G)
names = cellstr(unique(G,'stable')); k=numel(names);
pairs = nchoosek(1:k,2); p = nan(size(pairs,1),1);
for i=1:size(pairs,1)
    a=pairs(i,1); b=pairs(i,2);
    xa = X(strcmp(cellstr(G), names{a}));
    xb = X(strcmp(cellstr(G), names{b}));
    try, p(i)=ranksum(xa,xb); catch, [~,p(i)] = ttest2(xa,xb,'Vartype','unequal'); end
end
p_adj = fdr_bh(p);
end

function p_adj = fdr_bh(p)
p = p(:); n=numel(p); [ps,ix]=sort(p);
adj = nan(size(p)); adj(ix) = min( cummin( (n./(1:n)').*ps ), 1);
p_adj = adj;
end

function plotTwoBox(Xa,Xp,labels,ylab)
data=[]; grp=strings(0,1);
if ~isempty(Xa), data=[data; Xa(:)]; grp=[grp; repmat(string(labels{1}),numel(Xa),1)]; end
if ~isempty(Xp), data=[data; Xp(:)]; grp=[grp; repmat(string(labels{2}),numel(Xp),1)]; end
if isempty(data), axis off; return; end
boxplot(data, grp, 'Symbol','k+','Colors','k'); ylabel(ylab); grid on; box on
end

function addPBar(ax,x1,x2,yTop,p,pairLabel)
axes(ax); hold(ax,'on'); yl=ylim(ax); dy=0.05*range(yl); y=yTop+dy;
plot([x1 x1 x2 x2],[y y+dy y+dy y],'k-','LineWidth',1.1);
if ~isfinite(p), txt='p=NaN';
else, txt=sprintf('p=%.3g%s',p, iff(p<0.05,[' ' sigStar(p)],' n.s.')); end
text(mean([x1 x2]), y+dy*1.05, txt,'HorizontalAlignment','center','FontWeight','bold');
if nargin>5 && ~isempty(pairLabel)
    text(mean([x1 x2]), y+dy*2.0, pairLabel,'HorizontalAlignment','center','FontSize',9,'Color',[.3 .3 .3]);
end
end

function addSigBar(ax,x1,x2,yTop,p,pairLabel)
if ~isfinite(p) || p>=0.05, return; end
axes(ax); hold(ax,'on'); yl=ylim(ax); dy=0.05*range(yl); y=yTop+dy;
plot([x1 x1 x2 x2],[y y+dy y+dy y],'k-','LineWidth',1.2);
text(mean([x1 x2]), y+dy*1.05, sigStar(p), 'HorizontalAlignment','center','FontWeight','bold');
if nargin>5 && ~isempty(pairLabel)
    text(mean([x1 x2]), y+dy*2.0, pairLabel,'HorizontalAlignment','center','FontSize',9,'Color',[.25 .25 .25]);
end
end

function s=sigStar(p)
if p<1e-3, s='***'; elseif p<1e-2, s='**'; elseif p<0.05, s='*'; else, s=''; end
end
function s=iff(c,a,b), if c, s=a; else, s=b; end, end

function savepng(fh,fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution',150); catch, print(fh, fn, '-dpng','-r150'); end
end

function D = outerJoinInto(D, G, keys, srcName, dstName)
if isempty(G), return; end
if ~ismember(dstName,D.Properties.VariableNames), D.(dstName)=nan(height(D),1); end
keyD = joinKey(D, keys); keyG = joinKey(G, keys);
[tf,loc] = ismember(keyD, keyG);
vals = G.(srcName); D.(dstName)(tf) = vals(loc(tf));
end
function k=joinKey(T,keys)
k=repmat("",height(T),1);
for i=1:numel(keys)
    v=T.(keys{i}); if ~isstring(v), v=string(v); end
    if i==1, k=v; else, k=k+"|"+v; end
end
end
