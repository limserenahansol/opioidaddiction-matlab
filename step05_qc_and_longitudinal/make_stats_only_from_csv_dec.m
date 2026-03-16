function make_stats_only_from_csv_dec()
% Updated stats-only script for NEW cohort + NEW timeline.
%
% Changes vs original:
% 1) Group definition from cohort mapping (Active vs Passive), NOT "had-passive ever".
% 2) Time bins (5 periods):
%       pre (D3-5), during (D6-10), post (D11-13), withdrawal (D14-16), reexposure (D17-18)
% 3) Runs two versions:
%       ALLDAYS: uses all days in each period
%       CLEAN: excludes less reliable days (default [4 6 11 14 17])
% 4) Group comparisons per period:
%       - Unpaired ranksum (Active vs Passive)
%       - Paired signrank by pair_id (Active vs median(passives in pair))
% 5) Period comparisons (one-way ANOVA + Tukey; fallback KW + BH-FDR) within cohorts:
%       ALL / ACTIVE / PASSIVE
%
% Notes:
% - If a mouse_key can’t be parsed into cage+color, groupAP becomes <undefined> and is excluded.
% - 6099 white died: missing days are naturally handled.

%% ------------------- USER SETTINGS -------------------
lessReliableDays = [4 6 11 14 17];

%% -------- locate latest run and read CSV --------
tryPath = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(tryPath,'dir')
    here = pwd; cand = here;
    for up=1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), tryPath = p; break; end
        cand = fileparts(cand);
    end
end
if ~exist(tryPath,'dir'), error('Cannot find longitudinal_outputs folder.'); end

d = dir(fullfile(tryPath,'run_*')); assert(~isempty(d),'No run_* folders found.');
[~,idx] = max([d.datenum]); runDir = fullfile(d(idx).folder, d(idx).name);

csvPath = fullfile(runDir, 'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0, 'Missing %s', csvPath);

fprintf('\n[stats-only] Reading: %s\n', csvPath);
T = readtable(csvPath, 'VariableNamingRule','preserve');

% Ensure types
T = ensureStringVar(T,'day_name');
T = ensureStringVar(T,'mouse_key');
if hasVar(T,'Session_Paradigm'), T = ensureStringVar(T,'Session_Paradigm'); end

% Robust numeric conversion (handles string/cell/categorical)
if hasVar(T,'day_index'),   T.day_index   = toDouble(T.day_index);   end
if hasVar(T,'session_idx'), T.session_idx = toDouble(T.session_idx); end

% Create output folder
figDir = fullfile(runDir, 'figs'); if ~exist(figDir,'dir'), mkdir(figDir); end
outDir = fullfile(figDir, 'stats_only_NEWCOHORT'); if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('[stats-only] Saving outputs to: %s\n', outDir);

%% -------- drop frames from TTL-only "extra" trials --------
if hasVar(T,'TrialRequirement')
    trr = T.TrialRequirement; if ~isstring(trr), trr = string(trr); end
    if hasVar(T,'Trial')
        extraMask = (toDouble(T.Trial)>0) & strcmpi(strtrim(trr),"extra");
    else
        extraMask = strcmpi(strtrim(trr),"extra");
    end
    if any(extraMask)
        fprintf('[stats-only] Excluding %d frames with TrialRequirement="extra" (%.2f%% rows)\n', ...
            nnz(extraMask), 100*nnz(extraMask)/height(T));
        T(extraMask,:) = [];
    end
end

%% -------- build per-session, then per-day table --------
tb = pickTimebase(T);
hasLick   = hasVar(T,'Lick_TTL');
hasReward = hasVar(T,'Injector_TTL');
hasPupil  = hasVar(T,'Diameter_px');

% TTL to logical
if hasLick
    T.Lick_TTL = toDouble(T.Lick_TTL);
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
    T.Lick_TTL = T.Lick_TTL > 0.5;
end
if hasReward
    T.Injector_TTL = toDouble(T.Injector_TTL);
    T.Injector_TTL(isnan(T.Injector_TTL)) = 0;
    T.Injector_TTL = T.Injector_TTL > 0.5;
end

% session keys
sessKeys = unique(T(:,{'mouse_key','day_index','day_name','session_idx'}),'rows','stable');

S = table();
S.mouse_key   = sessKeys.mouse_key;
S.day_index   = sessKeys.day_index;
S.day_name    = sessKeys.day_name;
S.session_idx = sessKeys.session_idx;
S.SessionMinutes  = nan(height(S),1);

% Lick metrics (per session)
S.lick_n             = nan(height(S),1);
S.lick_freq_per_min  = nan(height(S),1);
S.lick_meanDur_s     = nan(height(S),1);
S.lick_totalDur_s    = nan(height(S),1);
S.lick_medianIEI_s   = nan(height(S),1);

% Reward metrics (per session)
S.rew_n              = nan(height(S),1);
S.rew_freq_per_min   = nan(height(S),1);
S.rew_meanDur_s      = nan(height(S),1);
S.rew_totalDur_s     = nan(height(S),1);
S.rew_medianIRI_s    = nan(height(S),1);

% Pupil
S.pupil_mean         = nan(height(S),1);

% Optional: last Requirement per session
if hasVar(T,'RequirementLast'), S.RequirementLast = nan(height(S),1); end

for i=1:height(sessKeys)
    mk = sessKeys.mouse_key(i); di = sessKeys.day_index(i); ss = sessKeys.session_idx(i);
    r  = (T.mouse_key==mk) & (T.day_index==di) & (T.session_idx==ss);
    if ~any(r), continue; end

    tb_s = tb(r);
    dur_s = finiteRange(tb_s);
    if ~(isfinite(dur_s) && dur_s>0) && hasVar(T,'Frame')
        fr = toDouble(T.Frame(r));
        if any(isfinite(fr)), dur_s = (max(fr)-min(fr))/30; end
    end
    S.SessionMinutes(i) = dur_s/60;

    if hasVar(T,'RequirementLast')
        S.RequirementLast(i) = mean(toDouble(T.RequirementLast(r)), 'omitnan');
    end

    if hasLick
        [n, meanDur, totalDur, medIEI] = eventMetrics(tb_s, T.Lick_TTL(r));
        S.lick_n(i)           = n;
        S.lick_meanDur_s(i)   = meanDur;
        S.lick_totalDur_s(i)  = totalDur;
        S.lick_medianIEI_s(i) = medIEI;
        if S.SessionMinutes(i)>0, S.lick_freq_per_min(i) = n/S.SessionMinutes(i); end
    end

    if hasReward
        [n, meanDur, totalDur, medIRI] = eventMetrics(tb_s, T.Injector_TTL(r));
        S.rew_n(i)            = n;
        S.rew_meanDur_s(i)    = meanDur;
        S.rew_totalDur_s(i)   = totalDur;
        S.rew_medianIRI_s(i)  = medIRI;
        if S.SessionMinutes(i)>0, S.rew_freq_per_min(i) = n/S.SessionMinutes(i); end
    end

    if hasPupil
        S.pupil_mean(i) = mean(toDouble(T.Diameter_px(r)), 'omitnan');
    end
end

% collapse to per-day (median across that day's sessions)
dayKeys = unique(T(:,{'mouse_key','day_index','day_name'}),'rows','stable');
D = dayKeys;
fillCols = setdiff(S.Properties.VariableNames, {'mouse_key','day_index','day_name','session_idx'});
for c=1:numel(fillCols), D.(fillCols{c}) = nan(height(D),1); end

for i=1:height(D)
    mk = D.mouse_key(i); di = D.day_index(i);
    r  = (S.mouse_key==mk) & (S.day_index==di);
    for c=1:numel(fillCols)
        D.(fillCols{c})(i) = median(S.(fillCols{c})(r), 'omitnan');
    end
end

% ---------- bring behavior tests / percentages (median per day) ----------
D = addDayScalarFromT(D, T, 'Immersion_Latency_s');
D = addDayPatternFromT(D, T, 'TST_Pct_');
D = addDayPatternFromT(D, T, 'HOT_Pct_');

latencyCands = {'TST_latency_s','TST_s','TailSuspension_s', ...
                'Hotplate_latency_s','HotPlate_latency_s', ...
                'TailImmersion_latency_s','TailImmersion_s'};
for cc = 1:numel(latencyCands)
    if hasVar(T, latencyCands{cc})
        D = addDayScalarFromT(D, T, latencyCands{cc});
    end
end

%% -------- NEW: attach cohort info (sex/group/pair_id) from mapping --------
cohort = buildNewCohortTable();  % FIXED (pair_id always numeric)
D = attachCohortInfo(D, cohort);

%% -------- NEW: assign period bins (5) --------
D.period5 = dayPeriod5(D.day_index);

%% -------- NEW: make two datasets: ALLDAYS and CLEAN (exclude lessReliableDays) --------
D_all = D;
D_clean = D;
rm = ismember(double(D_clean.day_index), double(lessReliableDays));
D_clean(rm,:) = [];

datasets = { ...
    struct('tag','ALLDAYS','D',D_all), ...
    struct('tag','CLEAN',  'D',D_clean)};

%% -------- which metrics to analyze (auto-detect & label) --------
labels = containers.Map('KeyType','char','ValueType','char');

% Lick metrics
if hasVar(D,'lick_freq_per_min'),  labels('lick_freq_per_min')  = 'Lick freq (/min)'; end
if hasVar(D,'lick_meanDur_s'),     labels('lick_meanDur_s')     = 'Lick mean dur (s)'; end
if hasVar(D,'lick_totalDur_s'),    labels('lick_totalDur_s')    = 'Lick total dur (s)'; end
if hasVar(D,'lick_medianIEI_s'),   labels('lick_medianIEI_s')   = 'Lick median IEI (s)'; end

% Reward metrics
if hasVar(D,'rew_freq_per_min'),   labels('rew_freq_per_min')   = 'Reward freq (/min)'; end
if hasVar(D,'rew_meanDur_s'),      labels('rew_meanDur_s')      = 'Reward mean dur (s)'; end
if hasVar(D,'rew_totalDur_s'),     labels('rew_totalDur_s')     = 'Reward total dur (s)'; end
if hasVar(D,'rew_medianIRI_s'),    labels('rew_medianIRI_s')    = 'Reward median IRI (s)'; end
if hasVar(D,'rew_medianIEI_s'),    labels('rew_medianIEI_s')    = 'Reward median IRI (s)'; end

% Pupil / Requirement
if hasVar(D,'pupil_mean'),         labels('pupil_mean')         = 'Pupil mean (px)'; end
if hasVar(D,'RequirementLast'),    labels('RequirementLast')    = 'Requirement'; end

% Behaviour percentages
for v = string(D.Properties.VariableNames)
    if startsWith(v,'TST_Pct_'), labels(char(v)) = ['TST % ' stripPrefix(char(v))]; end
    if startsWith(v,'HOT_Pct_'), labels(char(v)) = ['HOT % ' stripPrefix(char(v))]; end
end

% Behaviour latencies
if hasVar(D,'Immersion_Latency_s'),     labels('Immersion_Latency_s') = 'Tail immersion latency (s)'; end
if hasVar(D,'Hotplate_latency_s'),      labels('Hotplate_latency_s')  = 'Hotplate latency (s)'; end
if hasVar(D,'HotPlate_latency_s'),      labels('HotPlate_latency_s')  = 'Hotplate latency (s)'; end
if hasVar(D,'TailImmersion_latency_s'), labels('TailImmersion_latency_s') = 'Tail-immersion latency (s)'; end
if hasVar(D,'TailImmersion_s'),         labels('TailImmersion_s')      = 'Tail-immersion latency (s)'; end
if hasVar(D,'TST_latency_s'),           labels('TST_latency_s')        = 'TST latency (s)'; end
if hasVar(D,'TST_s'),                   labels('TST_s')               = 'TST latency (s)'; end
if hasVar(D,'TailSuspension_s'),        labels('TailSuspension_s')     = 'TST latency (s)'; end

metrics = labels.keys;
assert(~isempty(metrics), 'No analyzable metrics found in the CSV.');

periodNames = {'pre','during','post','withdrawal','reexposure'};

%% =======================================================================
% 1) GROUP comparisons per period: Active vs Passive (unpaired + paired-by-pair)
% =======================================================================
for ds = 1:numel(datasets)
    tag  = datasets{ds}.tag;
    Duse = datasets{ds}.D;

    for i=1:numel(metrics)
        y   = metrics{i};
        lbl = labels(y);

        fig = figure('Color','w','Position',[80 80 1680 520]);
        tlo = tiledlayout(1,numel(periodNames),'TileSpacing','compact','Padding','compact');

        csvRows = {}; % period, dataset, nA, nP, medA, medP, p_ranksum, nPairs, p_signrank

        for b=1:numel(periodNames)
            nexttile; hold on

            sub = Duse(string(Duse.period5)==periodNames{b}, {'mouse_key','groupAP','pair_id',y});
            sub = sub(isfinite(sub.(y)), :);
            sub = sub(~isundefined(sub.groupAP), :);

            if isempty(sub)
                axis off; title(periodNames{b}); continue;
            end

            % per-mouse median within this period
            M = groupsummary(sub, {'mouse_key','groupAP','pair_id'}, 'median', y);
            M.Properties.VariableNames(end) = {'y_med'};

            Xa = M.y_med(M.groupAP=="A"); Xa = Xa(isfinite(Xa));
            Xp = M.y_med(M.groupAP=="P"); Xp = Xp(isfinite(Xp));

            % unpaired boxplot
            plotTwoBox(Xa, Xp, {'Active','Passive'}, lbl);
            title(sprintf('%s (%s)', periodNames{b}, tag));
            yl = ylim(gca);

            % unpaired Wilcoxon ranksum
            p_unpaired = NaN;
            if numel(Xa)>=2 && numel(Xp)>=2
                p_unpaired = ranksum(Xa, Xp);
            end
            addPBar(gca, 1, 2, yl(2), p_unpaired, 'Unpaired: Active vs Passive');

            % paired-by-pair signrank:
            p_paired = NaN; nPairs = 0;
            if hasVar(M,'pair_id')
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
                    p_paired = signrank(actVals, pasVals);
                end
            end

            if isfinite(p_paired)
                txt = sprintf('Paired-by-pair: p=%.3g (nPairs=%d)', p_paired, nPairs);
            else
                txt = sprintf('Paired-by-pair: p=NaN (nPairs=%d)', nPairs);
            end
            text(1.5, yl(1)+0.05*range(yl), txt, 'HorizontalAlignment','center', ...
                'FontSize',9,'Color',[0.2 0.2 0.2], 'Interpreter','none');

            csvRows(end+1,:) = {periodNames{b}, tag, numel(Xa), numel(Xp), nanmedian(Xa), nanmedian(Xp), p_unpaired, nPairs, p_paired}; %#ok<AGROW>
        end

        title(tlo, sprintf('%s — Active vs Passive per period (%s)', lbl, tag));
        savepng(fig, fullfile(outDir, sprintf('GROUP_Active_vs_Passive_%s_%s.png', y, tag)));
        close(fig);

        Tcsv = cell2table(csvRows, 'VariableNames', ...
            {'period','dataset','n_active','n_passive','median_active','median_passive','p_ranksum_unpaired','nPairs','p_signrank_paired'});
        try, writetable(Tcsv, fullfile(outDir, sprintf('GROUP_Active_vs_Passive_%s_%s.csv', y, tag))); catch, end
    end
end

%% =======================================================================
% 2) PERIOD comparison across 5 periods (ANOVA/Tukey; fallback KW + BH-FDR)
% =======================================================================
for ds = 1:numel(datasets)
    tag  = datasets{ds}.tag;
    Duse = datasets{ds}.D;

    cohorts = { ...
        struct('name','ALL',    'mask', true(height(Duse),1)), ...
        struct('name','ACTIVE', 'mask', Duse.groupAP=="A"), ...
        struct('name','PASSIVE','mask', Duse.groupAP=="P")};

    for ci = 1:numel(cohorts)
        Cc = cohorts{ci};

        for i=1:numel(metrics)
            y   = metrics{i};
            lbl = labels(y);

            sub = Duse(Cc.mask, {'mouse_key','period5',y});
            sub = sub(~isundefined(sub.period5) & isfinite(sub.(y)), :);
            if isempty(sub), continue; end

            % per-mouse, per-period medians
            M = groupsummary(sub, {'mouse_key','period5'}, 'median', y);
            M.Properties.VariableNames(end) = {'y_med'};
            M = M(~isundefined(M.period5) & isfinite(M.y_med), :);
            if isempty(M), continue; end

            X = double(M.y_med);
            G = cellstr(string(M.period5));

            usedKW     = false;
            methodName = 'ANOVA';
            posthoc    = 'Tukey HSD';
            pairs = []; padj = []; grpNames = unique(G,'stable');

            try
                [p_overall, ~, st] = anova1(X, G, 'off');
                try
                    [mcTbl, ~, ~, gn] = multcompare(st, 'CType','hsd', 'Display','off');
                    pairs    = mcTbl(:,1:2);
                    padj     = mcTbl(:,6);
                    grpNames = gn;
                catch
                    [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
                    posthoc = 'BH';
                end
            catch
                usedKW     = true;
                methodName = 'KW';
                posthoc    = 'BH';
                p_overall  = kruskalwallis(X, G, 'off');
                [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
            end

            fig = figure('Color','w','Position',[100 100 1050 520]); hold on
            boxplot(X, G, 'Symbol','k+','Colors','k'); grid on; box on
            ylabel(lbl); xlabel('Period');
            title(sprintf('%s — %s across periods (%s, %s)  p=%.3g', ...
                lbl, methodName, Cc.name, tag, p_overall));

            yl = ylim; yBase = yl(2); dy = 0.05*range(yl); level = 0;
            for r=1:size(pairs,1)
                a = pairs(r,1); b = pairs(r,2);
                yTop = yBase + (level+1)*dy;
                txt  = sprintf('%s vs %s (%s p=%.3g)', grpNames{a}, grpNames{b}, posthoc, padj(r));
                addSigBar(gca, min(a,b), max(a,b), yTop, padj(r), txt, true);
                level = level + 1;
            end
            if size(pairs,1)>0, ylim([yl(1), yBase + (level+3)*dy]); end

            savepng(fig, fullfile(outDir, sprintf('PERIOD_%s_%s_%s_%s.png', Cc.name, y, tag, methodName)));
            close(fig);

            pairStr = strings(numel(padj),1);
            for r=1:numel(padj)
                pairStr(r) = sprintf('%s vs %s', grpNames{pairs(r,1)}, grpNames{pairs(r,2)});
            end

            Tcsv_pairs = table(repmat(string(Cc.name),numel(padj),1), repmat(string(tag),numel(padj),1), pairStr, padj, ...
                'VariableNames', {'cohort','dataset','pair','p_adj'});
            Tcsv_over  = table(string(Cc.name), string(tag), usedKW, p_overall, ...
                'VariableNames', {'cohort','dataset','used_KW','p_overall'});

            try
                writetable(Tcsv_over,  fullfile(outDir, sprintf('PERIOD_%s_%s_%s_overall.csv', Cc.name, y, tag)));
                writetable(Tcsv_pairs, fullfile(outDir, sprintf('PERIOD_%s_%s_%s_pairs.csv',   Cc.name, y, tag)));
            catch
            end
        end
    end
end

fprintf('[stats-only] Done.\n');
end

%% ====================== helpers ======================

function tf = hasVar(T,var)
tf = ismember(var, T.Properties.VariableNames);
end

function T = ensureStringVar(T, vname)
if hasVar(T,vname) && ~isstring(T.(vname))
    T.(vname) = string(T.(vname));
end
end

function x = toDouble(v)
% Robust conversion to double from numeric / logical / categorical / string / cell
if isnumeric(v)
    x = double(v);
elseif islogical(v)
    x = double(v);
elseif iscategorical(v)
    x = double(v);  % category indices
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

function s = stripPrefix(nm)
s = regexprep(nm,'^(TST|HOT)_Pct_','');
s = strrep(s,'_',' ');
end

function tb = pickTimebase(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
tb = nan(height(T),1);
for i=1:numel(cands)
    if hasVar(T,cands{i})
        v = toDouble(T.(cands{i}));
        if any(isfinite(v)), tb = v; return; end
    end
end
if hasVar(T,'Frame'), tb = toDouble(T.Frame)/30; end
end

function r = finiteRange(x)
x = toDouble(x(:)); x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x)-min(x); end
end

function [n, meanDur, totalDur, medianIEI] = eventMetrics(tb_s, ttl)
ttl = logical(ttl(:)); tb_s = toDouble(tb_s(:));
if numel(tb_s) < 2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
d = diff([false; ttl; false]); onIdx = find(d==1); offIdx = find(d==-1)-1;
n = numel(onIdx); if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
dt = [diff(tb_s); median(diff(tb_s),'omitnan')]; segDur = zeros(n,1);
for k=1:n
    ii = onIdx(k):offIdx(k);
    segDur(k) = sum(dt(ii));
end
meanDur  = mean(segDur,'omitnan');
totalDur = sum(segDur,'omitnan');
onTimes = tb_s(onIdx);
if numel(onTimes)>=2, medianIEI = median(diff(onTimes),'omitnan'); else, medianIEI = NaN; end
end

function D = addDayScalarFromT(D, T, col)
if ~hasVar(T,col), return; end
G = groupsummary(T, {'mouse_key','day_index'}, 'median', col);
D = outerJoinInto(D, G, {'mouse_key','day_index'}, ['median_' col], col);
end

function D = addDayPatternFromT(D, T, prefix)
namesT = T.Properties.VariableNames;
hit = namesT(contains(namesT,prefix));
for i=1:numel(hit)
    G = groupsummary(T, {'mouse_key','day_index'}, 'median', hit{i});
    D = outerJoinInto(D, G, {'mouse_key','day_index'}, ['median_' hit{i}], hit{i});
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
        p(i) = ranksum(xa, xb);
    else
        p(i) = NaN;
    end
end
p_adj = fdr_bh(p);
end

function p_adj = fdr_bh(p)
p = p(:); n = numel(p);
[ps,idx] = sort(p); adj = nan(size(p));
adj(idx) = min( cummin( (n./(1:n)').*ps ), 1);
p_adj = adj;
end

%% ====================== NEW cohort helpers ======================

function cohort = buildNewCohortTable()
% Build mapping.
% Columns: cage (string), color (string), sex (F/M), groupAP (A/P), pair_id (double)
%
% FIX: pair_id is stored as NUMERIC vector (not cell), so no double(cell) crash.

cage  = strings(0,1);
color = strings(0,1);
sex   = strings(0,1);
grp   = strings(0,1);
pair  = zeros(0,1);

% 6100
add('6100','red','F','P',1);
add('6100','orange','F','P',1);
add('6100','black','F','A',1);

% 0911
add('0911','red','F','A',2);
add('0911','orange','F','P',2);

add('0911','black','F','P',3);
add('0911','white','F','A',3);

% 0910
add('0910','red','M','P',4);
add('0910','orange','M','P',4);
add('0910','black','M','A',4);

% 6099
add('6099','red','M','P',5);
add('6099','orange','M','A',5);

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
% Adds: D.cage, D.color, D.sex, D.groupAP, D.pair_id

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

D.sex    = categorical(sex, ["F","M"]);
D.groupAP = categorical(grp, ["A","P"]);
D.pair_id = pid;
end

function [cage, color] = parseMouseKey(mk)
% Robust parsing for common mouse_key formats.
% Expected: 4-digit cage + a color token (red/orange/black/white)

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
