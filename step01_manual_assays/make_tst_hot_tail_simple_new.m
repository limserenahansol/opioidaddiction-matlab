function make_tst_hot_tail_simple_new()
% Period-based analysis for Tail Immersion, TST, Hotplate using NEW cohort mapping.
%
% What this does:
% - Uses cage+color -> sex/group(A/P)/pair_id mapping (NO use of isPassive)
% - Ignores day1-2 implicitly (period bins start at day3)
% - Defines periods:
%       pre        (D3-5)
%       during     (D6-10)
%       post       (D11-13)
%       withdrawal (D14-16)
%       reexposure (D17-18)
% - Two versions:
%       ALLDAYS: use all days
%       CLEAN: excludes less reliable "switch" days (default [4 6 11 14 17])
% - Per period:
%       Unpaired ranksum (Active vs Passive)
%       Paired-by-pair signrank (Active vs median(passives in same pair))
% - Across periods within cohorts (ALL / ACTIVE / PASSIVE):
%       ANOVA + Tukey; fallback KW + BH-FDR
%
% Outputs: figures + CSVs in latest run folder under TST_HOT_TAIL_NEWCOHORT

%% ------------------- USER SETTINGS -------------------
% Less reliable days to exclude for CLEAN version:
% You mentioned day4/day6/day11/day14 are less reliable; reexposure starts day17.
% Keep 17 here if your first reexposure day is also a switch day for you.
lessReliableDays = [4 6 11 14 17];

%% -------- locate latest run + read CSV --------
tryPath = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
assert(exist(tryPath,'dir')>0,'Cannot find longitudinal_outputs');

d = dir(fullfile(tryPath,'run_*'));
assert(~isempty(d),'No run_* folders found.');
[~,idx] = max([d.datenum]);
runDir = fullfile(d(idx).folder, d(idx).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('[TST/HOT/TAIL NEWCOHORT] Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% basic types
T = ensureStringVar(T,'mouse_key');
T = ensureStringVar(T,'day_name');
if hasVar(T,'day_index'),   T.day_index   = toDouble(T.day_index);   end

%% output dir
outDir = fullfile(runDir,'TST_HOT_TAIL_NEWCOHORT');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('[TST/HOT/TAIL NEWCOHORT] Saving to: %s\n', outDir);

%% exact columns (prefer these; you can add alternates if needed)
colTail = pickVar(T, {'Immersion_Latency_s'});
colTST  = pickVar(T, {'TST_Frames_Non_moving'});
colHOT  = pickVar(T, {'HOT_Frames_Non_moving'});

if isempty(colTail) && isempty(colTST) && isempty(colHOT)
    error('None of Immersion_Latency_s, TST_Frames_Non_moving, HOT_Frames_Non_moving found.');
end

%% -------- build per-mouse-per-day medians --------
keys = {'mouse_key','day_index','day_name'};

metricVars = {};
if ~isempty(colTail), metricVars{end+1} = colTail; end %#ok<AGROW>
if ~isempty(colTST),  metricVars{end+1} = colTST;  end %#ok<AGROW>
if ~isempty(colHOT),  metricVars{end+1} = colHOT;  end %#ok<AGROW>

% keep only rows with at least one metric finite
keep = false(height(T),1);
for i=1:numel(metricVars)
    v = metricVars{i};
    keep = keep | isfinite(toDouble(T.(v)));
end
T2 = T(keep, [keys metricVars]);

% build daily medians per metric
D = unique(T2(:,keys),'rows','stable');
for i=1:numel(metricVars)
    v = metricVars{i};
    Gi = groupsummary(T2(:,[keys {v}]), keys, 'median', v);
    srcName = ['median_' v];
    % groupsummary name can vary if v not valid identifier; use last var in Gi
    srcName = Gi.Properties.VariableNames{end};
    D.(v) = nan(height(D),1);
    D = outerJoinInto(D, Gi, keys, srcName, v);
end

% Ensure numeric day_index
if ~isnumeric(D.day_index), D.day_index = toDouble(D.day_index); end

%% -------- attach NEW cohort info --------
cohort = buildNewCohortTable(); % your mapping
D = attachCohortInfo(D, cohort);

% drop rows not mapped to A/P (parsing failure)
mapped = ~isundefined(D.groupAP);
if any(~mapped)
    fprintf('[TST/HOT/TAIL NEWCOHORT] Excluding %d rows with unmapped mouse_key (cannot parse cage/color).\n', nnz(~mapped));
end
D = D(mapped,:);

%% -------- assign 5 periods --------
D.period5 = dayPeriod5(D.day_index);

% drop days outside 3-18 (day1-2 habituation, anything else)
inRange = ~isundefined(D.period5);
D = D(inRange,:);

%% -------- datasets: ALLDAYS vs CLEAN --------
D_all = D;
D_clean = D;
rm = ismember(double(D_clean.day_index), double(lessReliableDays));
D_clean(rm,:) = [];

datasets = { ...
    struct('tag','ALLDAYS','D',D_all), ...
    struct('tag','CLEAN',  'D',D_clean)};

periodNames = {'pre','during','post','withdrawal','reexposure'};

%% -------- analyze each metric --------
if ~isempty(colTail)
    analyzeMetricByPeriod(datasets, colTail, 'Tail immersion latency (s)', outDir, 'TAIL', periodNames);
end
if ~isempty(colTST)
    analyzeMetricByPeriod(datasets, colTST, 'TST non-moving (frames)', outDir, 'TST', periodNames);
end
if ~isempty(colHOT)
    analyzeMetricByPeriod(datasets, colHOT, 'Hotplate non-moving (frames)', outDir, 'HOT', periodNames);
end

fprintf('[TST/HOT/TAIL NEWCOHORT] Done.\n');
end

%% ======================= ANALYSIS CORE =======================

function analyzeMetricByPeriod(datasets, yvar, ylab, outDir, tag, periodNames)

for ds = 1:numel(datasets)
    Duse = datasets{ds}.D;
    dtag = datasets{ds}.tag;

    sub = Duse(:, {'mouse_key','day_index','period5','groupAP','pair_id',yvar});
    sub.y = toDouble(sub.(yvar));
    sub = sub(isfinite(sub.y) & ~isundefined(sub.period5) & ~isundefined(sub.groupAP), :);
    if isempty(sub)
        warning('%s (%s): no data', yvar, dtag);
        continue;
    end

    %% -------- 1) Active vs Passive per PERIOD (unpaired + paired-by-pair) --------
    fig = figure('Color','w','Position',[80 80 1680 520]);
    tlo = tiledlayout(1,numel(periodNames),'TileSpacing','compact','Padding','compact');

    csvRows = {}; % period, dataset, nA, nP, medA, medP, p_unpaired, nPairs, p_paired

    for b=1:numel(periodNames)
        nexttile; hold on
        pname = periodNames{b};

        % rows in this period
        sP = sub(string(sub.period5)==pname, :);
        if isempty(sP)
            axis off; title(pname); continue;
        end

        % per-mouse median within this period (each mouse contributes once per period)
        M = groupsummary(sP(:,{'mouse_key','groupAP','pair_id','y'}), {'mouse_key','groupAP','pair_id'}, 'median', 'y');
        M.Properties.VariableNames(end) = {'y_med'};

        Xa = M.y_med(M.groupAP=="A"); Xa = Xa(isfinite(Xa));
        Xp = M.y_med(M.groupAP=="P"); Xp = Xp(isfinite(Xp));

        plotTwoBox(Xa, Xp, {'Active','Passive'}, ylab);
        title(sprintf('%s (%s)', pname, dtag));
        yl = ylim(gca);

        % unpaired ranksum (only if both groups have data)
        p_unpaired = NaN;
        if numel(Xa)>=2 && numel(Xp)>=2
            try
                p_unpaired = ranksum(Xa, Xp);
            catch
                [~,p_unpaired] = ttest2(Xa, Xp, 'Vartype','unequal');
            end
        end
        addPBar(gca, 1, 2, yl(2), p_unpaired, 'Unpaired: Active vs Passive');

        % paired-by-pair signrank: Active vs median(passives) within each pair_id
        p_paired = NaN; nPairs = 0;
        pids = unique(M.pair_id);
        actVals = []; pasVals = [];
        for kk=1:numel(pids)
            pid = pids(kk);
            if isnan(pid) || pid==0, continue; end

            a = M.y_med(M.pair_id==pid & M.groupAP=="A");
            p = M.y_med(M.pair_id==pid & M.groupAP=="P");
            a = a(isfinite(a)); p = p(isfinite(p));
            if isempty(a) || isempty(p), continue; end

            actVals(end+1,1) = median(a,'omitnan'); %#ok<AGROW>
            pasVals(end+1,1) = median(p,'omitnan'); %#ok<AGROW>
        end

        nPairs = numel(actVals);
        if nPairs>=2
            try
                p_paired = signrank(actVals, pasVals);
            catch
                p_paired = NaN;
            end
        end

        if isfinite(p_paired)
            txt = sprintf('Paired-by-pair: p=%.3g (nPairs=%d)', p_paired, nPairs);
        else
            txt = sprintf('Paired-by-pair: p=NaN (nPairs=%d)', nPairs);
        end
        text(1.5, yl(1)+0.05*range(yl), txt, 'HorizontalAlignment','center', ...
            'FontSize',9,'Color',[0.2 0.2 0.2], 'Interpreter','none');

        csvRows(end+1,:) = {pname, dtag, numel(Xa), numel(Xp), nanmedian(Xa), nanmedian(Xp), p_unpaired, nPairs, p_paired}; %#ok<AGROW>
    end

    title(tlo, sprintf('%s — Active vs Passive per period (%s)', ylab, dtag));
    savepng(fig, fullfile(outDir, sprintf('%s_GROUP_Period_Active_vs_Passive_%s.png', tag, dtag)));
    close(fig);

    Tcsv = cell2table(csvRows, 'VariableNames', ...
        {'period','dataset','n_active','n_passive','median_active','median_passive','p_ranksum_unpaired','nPairs','p_signrank_paired'});
    try
        writetable(Tcsv, fullfile(outDir, sprintf('%s_GROUP_Period_Active_vs_Passive_%s.csv', tag, dtag)));
    catch
    end

    %% -------- 2) Across PERIODS within cohorts (ALL / ACTIVE / PASSIVE) --------
    cohorts = { ...
        struct('name','ALL',    'mask', true(height(sub),1)), ...
        struct('name','ACTIVE', 'mask', sub.groupAP=="A"), ...
        struct('name','PASSIVE','mask', sub.groupAP=="P")};

    for ci=1:numel(cohorts)
        C = cohorts{ci};
        sC = sub(C.mask, {'mouse_key','period5','y'});
        sC = sC(~isundefined(sC.period5) & isfinite(sC.y), :);
        if isempty(sC), continue; end

        % per-mouse per-period medians
        M2 = groupsummary(sC, {'mouse_key','period5'}, 'median', 'y');
        M2.Properties.VariableNames(end) = {'y_med'};
        M2 = M2(~isundefined(M2.period5) & isfinite(M2.y_med), :);
        if isempty(M2), continue; end

        X = toDouble(M2.y_med);
        G = cellstr(string(M2.period5));

        usedKW = false;
        methodName = 'ANOVA';
        posthoc = 'Tukey HSD';
        pairs = []; padj = []; grpNames = unique(G,'stable');
        p_overall = NaN;

        try
            [p_overall,~,st] = anova1(X, G, 'off');
            try
                [mcTbl,~,~,gn] = multcompare(st, 'CType','hsd', 'Display','off');
                pairs = mcTbl(:,1:2);
                padj  = mcTbl(:,6);
                grpNames = gn;
            catch
                [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
                posthoc = 'BH';
            end
        catch
            usedKW = true;
            methodName = 'KW';
            posthoc = 'BH';
            try
                p_overall = kruskalwallis(X, G, 'off');
                [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
            catch
                p_overall = NaN;
                pairs=[]; padj=[]; grpNames = unique(G,'stable');
            end
        end

        fig = figure('Color','w','Position',[100 100 1050 520]); hold on
        boxplot(X, G, 'Symbol','k+','Colors','k'); grid on; box on
        ylabel(ylab); xlabel('Period');
        title(sprintf('%s — %s across periods (%s, %s)  p=%.3g', ...
            ylab, methodName, C.name, dtag, p_overall));

        yl = ylim; yBase = yl(2); dy = 0.05*range(yl); level = 0;
        for r=1:size(pairs,1)
            a = pairs(r,1); b = pairs(r,2);
            yTop = yBase + (level+1)*dy;
            txt = sprintf('%s vs %s (%s p=%.3g)', grpNames{a}, grpNames{b}, posthoc, padj(r));
            addSigBar(gca, min(a,b), max(a,b), yTop, padj(r), txt, true);
            level = level + 1;
        end
        if size(pairs,1)>0, ylim([yl(1), yBase + (level+3)*dy]); end

        savepng(fig, fullfile(outDir, sprintf('%s_PERIOD_%s_%s_%s.png', tag, C.name, dtag, methodName)));
        close(fig);

        % CSV outputs
        pairStr = strings(numel(padj),1);
        for r=1:numel(padj)
            pairStr(r) = sprintf('%s vs %s', grpNames{pairs(r,1)}, grpNames{pairs(r,2)});
        end
        Tcsv_pairs = table(repmat(string(C.name),numel(padj),1), repmat(string(dtag),numel(padj),1), pairStr, padj, ...
            'VariableNames', {'cohort','dataset','pair','p_adj'});
        Tcsv_over  = table(string(C.name), string(dtag), usedKW, p_overall, ...
            'VariableNames', {'cohort','dataset','used_KW','p_overall'});

        try
            writetable(Tcsv_over,  fullfile(outDir, sprintf('%s_PERIOD_%s_%s_overall.csv', tag, C.name, dtag)));
            writetable(Tcsv_pairs, fullfile(outDir, sprintf('%s_PERIOD_%s_%s_pairs.csv',   tag, C.name, dtag)));
        catch
        end
    end
end
end

%% ======================= BASIC HELPERS =======================

function tf = hasVar(T,var)
tf = ismember(var, T.Properties.VariableNames);
end

function T = ensureStringVar(T, vname)
if hasVar(T,vname) && ~isstring(T.(vname))
    T.(vname) = string(T.(vname));
end
end

function x = toDouble(v)
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif iscategorical(v)
    x = double(v);
elseif isstring(v)
    x = str2double(v);
elseif iscell(v)
    if isempty(v)
        x = nan(size(v));
    elseif isnumeric(v{1})
        x = cell2mat(v);
    else
        x = str2double(string(v));
    end
else
    try
        x = double(v);
    catch
        x = str2double(string(v));
    end
end
end

function varname = pickVar(T, candidates)
varname='';
for i=1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        varname = candidates{i};
        return;
    end
end
end

function D = outerJoinInto(D, G, keys, srcName, dstName)
if isempty(G), return; end
if ~ismember(dstName, D.Properties.VariableNames)
    D.(dstName) = nan(height(D),1);
end
keyD = buildKey(D, keys);
keyG = buildKey(G, keys);
[tf, loc] = ismember(keyD, keyG);
vals = toDouble(G.(srcName));
D.(dstName)(tf) = vals(loc(tf));
end

function k = buildKey(T, keys)
k = repmat("", height(T), 1);
for i = 1:numel(keys)
    v = T.(keys{i}); if ~isstring(v), v = string(v); end
    if i==1, k=v; else, k = k + "|" + v; end
end
end

function c = dayPeriod5(day_index)
di = toDouble(day_index);
lab = strings(size(di));
lab(di>=3  & di<=5 )  = "pre";
lab(di>=6  & di<=10)  = "during";
lab(di>=11 & di<=13)  = "post";
lab(di>=14 & di<=16)  = "withdrawal";
lab(di>=17 & di<=18)  = "reexposure";
c = categorical(lab, ["pre","during","post","withdrawal","reexposure"]);
end

function plotTwoBox(Xa, Xp, xlabels, ylab)
data = []; grpNames = strings(0,1);
if ~isempty(Xa)
    data=[data; Xa(:)];
    grpNames=[grpNames; repmat(string(xlabels{1}),numel(Xa),1)];
end
if ~isempty(Xp)
    data=[data; Xp(:)];
    grpNames=[grpNames; repmat(string(xlabels{2}),numel(Xp),1)];
end
if isempty(data), axis off; title('No data'); return; end
boxplot(data, grpNames, 'Symbol','k+','Colors','k');
ylabel(ylab); grid on; box on
end

function addSigBar(ax, x1, x2, yTop, p, pairLabel, showAlways)
if nargin<7, showAlways=false; end
axes(ax); hold(ax,'on');
yl = ylim(ax); dy = 0.04*range(yl);
y  = yTop + dy;
plot([x1 x1 x2 x2], [y y+dy y+dy y], 'k-','LineWidth',1.1);
stars = sigStar(p);
if showAlways || (isfinite(p) && p<0.05)
    if isempty(stars), stars = 'n.s.'; end
    txt = sprintf('%s  p=%.3g', stars, p);
    text(mean([x1 x2]), y+dy*1.05, txt, 'HorizontalAlignment','center', ...
        'FontWeight','bold','Color','k','Interpreter','none');
    if nargin>=6 && ~isempty(pairLabel)
        text(mean([x1 x2]), y+dy*2.0, pairLabel, 'HorizontalAlignment','center', ...
            'FontSize',9,'Color',[0.2 0.2 0.2],'Interpreter','none');
    end
end
end

function addPBar(ax, x1, x2, yTop, p, pairLabel)
if nargin<6, pairLabel = ''; end
axes(ax); hold(ax,'on');
yl = ylim(ax); dy = 0.04*range(yl); y = yTop + dy;
plot([x1 x1 x2 x2], [y y+dy y+dy y], 'k-','LineWidth',1.2);

if isfinite(p)
    if p < 0.05
        txt = sprintf('p=%.3g %s', p, sigStar(p));
    else
        txt = sprintf('p=%.3g n.s.', p);
    end
else
    txt = 'p=NaN';
end

text(mean([x1 x2]), y+dy*1.05, txt, ...
    'HorizontalAlignment','center','FontWeight','bold','Color','k');

if ~isempty(pairLabel)
    text(mean([x1 x2]), y+dy*2.0, pairLabel, ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.2 0.2 0.2]);
end

ylim([yl(1), y+3*dy]);
end

function s = sigStar(p)
if ~isfinite(p), s = ''; return; end
if p < 1e-3, s = '***';
elseif p < 1e-2, s = '**';
elseif p < 0.05, s = '*';
else, s = '';
end
end

function savepng(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',150);
catch
    print(fh, fn, '-dpng','-r150');
end
end

function [pairs, p_adj, names] = fallbackPairsFDR(X, G)
names = unique(G,'stable'); k = numel(names);
pairs = nchoosek(1:k,2); p = nan(size(pairs,1),1);
for i=1:size(pairs,1)
    a = pairs(i,1); b = pairs(i,2);
    xa = X(strcmp(G, names{a})); xb = X(strcmp(G, names{b}));
    xa = xa(isfinite(xa)); xb = xb(isfinite(xb));
    if numel(xa)>=2 && numel(xb)>=2
        try
            p(i) = ranksum(xa, xb);
        catch
            [~,p(i)] = ttest2(xa, xb, 'Vartype','unequal');
        end
    else
        p(i) = NaN;
    end
end
p_adj = fdr_bh(p);
end

function p_adj = fdr_bh(p)
p = p(:); n = numel(p);
[ps,idx] = sort(p);
adj = nan(size(p));
adj(idx) = min(cummin((n./(1:n)').*ps), 1);
p_adj = adj;
end

%% ======================= COHORT (NEW) HELPERS =======================

function cohort = buildNewCohortTable()
% Your NEW cohort mapping (cage, color, sex, groupAP, pair_id)
% groupAP: A=Active, P=Passive

cage  = strings(0,1);
color = strings(0,1);
sex   = strings(0,1);
grp   = strings(0,1);
pair  = zeros(0,1);

% Pair 1: 6100 orange+red (P) vs 6100 black (A)
add('6100','orange','F','P',1);
add('6100','red',   'F','P',1);
add('6100','black', 'F','A',1);

% Pair 2: 0911 red (A) vs 0911 orange (P)
add('0911','red',   'F','A',2);
add('0911','orange','F','P',2);

% Pair 3: 0911 white (A) vs 0911 black (P)
add('0911','white','F','A',3);
add('0911','black','F','P',3);

% Pair 4: 0910 black (A) vs 0910 orange+red (P)
add('0910','black', 'M','A',4);
add('0910','orange','M','P',4);
add('0910','red',   'M','P',4);

% Pair 5: 6099 orange (A) vs 6099 red (P)
add('6099','orange','M','A',5);
add('6099','red',   'M','P',5);

% Pair 6: 6099 black (A) vs 6099 white (P) [white died day1-13; OK]
add('6099','black','M','A',6);
add('6099','white','M','P',6);

cohort = table(cage, color, sex, grp, pair, ...
    'VariableNames', {'cage','color','sex','groupAP','pair_id'});

    function add(cg, col, sx, gp, pid)
        cage(end+1,1)  = string(cg);
        color(end+1,1) = lower(string(col));
        sex(end+1,1)   = upper(string(sx));
        grp(end+1,1)   = upper(string(gp));
        pair(end+1,1)  = double(pid);
    end
end

function D = attachCohortInfo(D, cohort)
% Parse D.mouse_key into cage+color and join cohort info.
% Adds:
%   D.cage, D.color, D.sex, D.groupAP, D.pair_id

n = height(D);
cage = strings(n,1);
color = strings(n,1);

for i=1:n
    mk = lower(string(D.mouse_key(i)));
    [cg, col] = parseMouseKey(mk);
    cage(i) = cg;
    color(i) = col;
end

D.cage  = cage;
D.color = color;

keyD = D.cage + "|" + D.color;
keyC = cohort.cage + "|" + cohort.color;

sex = strings(n,1);
grp = strings(n,1);
pid = nan(n,1);

[tf, loc] = ismember(keyD, keyC);
sex(tf) = cohort.sex(loc(tf));
grp(tf) = cohort.groupAP(loc(tf));
pid(tf) = cohort.pair_id(loc(tf));

D.sex     = categorical(sex, ["F","M"]);
D.groupAP = categorical(grp, ["A","P"]);
D.pair_id = pid;
end

function [cage, color] = parseMouseKey(mk)
% Robust parsing for typical mouse_key formats.
% Works for:
%   "6100_red", "6100red", "6100/red", "cage6100_black", "0910 orange", "6099 red_m_p", etc.

cage = "";
color = "";

m = regexp(mk, '(\d{4})', 'tokens', 'once');
if ~isempty(m), cage = string(m{1}); end

colors = ["red","orange","black","white"];
for c = colors
    if contains(mk, c)
        color = c;
        break;
    end
end
end
