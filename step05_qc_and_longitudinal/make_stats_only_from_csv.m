function make_stats_only_from_csv()
% Stand-alone stats-only script:
% 1) Group comparison (Active-only vs Had-passive) per time bin (box/whisker)
%    -> nonparam Wilcoxon ranksum; p is drawn and SHOWN even if n.s.
% 2) Period comparison across 5 bins (one-way ANOVA + Tukey HSD)
%    -> fallback to Kruskal–Wallis + BH-FDR when ANOVA/multcompare unavailable
% 3) Writes CSVs with all p-values and saves figures.
%
% Definitions:
%   "Had-passive" = mice that EVER had any passive session (isPassive==1 on any day).
%   For group plots we use ALL days from each mouse (not only passive days).

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
T = ensureString(T,'day_name');
T = ensureString(T,'mouse_key');
if hasVar(T,'Session_Paradigm'), T = ensureString(T,'Session_Paradigm'); end
if hasVar(T,'isPassive') && ~isnumeric(T.isPassive), T.isPassive = double(T.isPassive); end
if hasVar(T,'day_index')   && ~isnumeric(T.day_index),   T.day_index   = double(T.day_index);   end
if hasVar(T,'session_idx') && ~isnumeric(T.session_idx), T.session_idx = double(T.session_idx); end

% Create output folder
figDir = fullfile(runDir, 'figs'); if ~exist(figDir,'dir'), mkdir(figDir); end
outDir = fullfile(figDir, 'stats_only'); if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('[stats-only] Saving outputs to: %s\n', outDir);

%% -------- drop frames from TTL-only "extra" trials --------
if hasVar(T,'TrialRequirement')
    trr = T.TrialRequirement; if ~isstring(trr), trr = string(trr); end
    extraMask = false(height(T),1);
    if hasVar(T,'Trial'), extraMask = (T.Trial>0) & strcmpi(strtrim(trr),"extra");
    else, extraMask = strcmpi(strtrim(trr),"extra");
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
if hasLick,   T.Lick_TTL(isnan(T.Lick_TTL))         = 0; T.Lick_TTL       = T.Lick_TTL > 0.5; end
if hasReward, T.Injector_TTL(isnan(T.Injector_TTL)) = 0; T.Injector_TTL   = T.Injector_TTL > 0.5; end

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
S.rew_medianIRI_s    = nan(height(S),1);  % (IRI) intervals between rewards

% Pupil
S.pupil_mean         = nan(height(S),1);

% Optional: last Requirement per session
if hasVar(T,'RequirementLast'), S.RequirementLast = nan(height(S),1); end

for i=1:height(sessKeys)
    mk = sessKeys.mouse_key(i); di = sessKeys.day_index(i); ss = sessKeys.session_idx(i);
    r  = T.mouse_key==mk & T.day_index==di & T.session_idx==ss;
    if ~any(r), continue; end

    tb_s = tb(r);
    dur_s = finiteRange(tb_s);
    if ~(isfinite(dur_s) && dur_s>0) && hasVar(T,'Frame')
        fr = double(T.Frame(r)); if any(isfinite(fr)), dur_s = (max(fr)-min(fr))/30; end
    end
    S.SessionMinutes(i) = dur_s/60;

    if hasVar(T,'RequirementLast')
        S.RequirementLast(i) = mean(double(T.RequirementLast(r)), 'omitnan');
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
        S.pupil_mean(i) = mean(double(T.Diameter_px(r)), 'omitnan');
    end
end

% collapse to per-day (median across that day's sessions)
dayKeys = unique(T(:,{'mouse_key','day_index','day_name'}),'rows','stable');
D = dayKeys;
fillCols = setdiff(S.Properties.VariableNames, {'mouse_key','day_index','day_name','session_idx'});
for c=1:numel(fillCols), D.(fillCols{c}) = nan(height(D),1); end
for i=1:height(D)
    mk = D.mouse_key(i); di = D.day_index(i);
    r  = S.mouse_key==mk & S.day_index==di;
    for c=1:numel(fillCols)
        D.(fillCols{c})(i) = median(S.(fillCols{c})(r), 'omitnan');
    end
end

% ---------- bring behavior tests / percentages (median per day) ----------
% Immersion latency
D = addDayScalarFromT(D, T, 'Immersion_Latency_s');
% TST / HOT % (any columns with these prefixes)
D = addDayPatternFromT(D, T, 'TST_Pct_');
D = addDayPatternFromT(D, T, 'HOT_Pct_');
% Common latency field name variants
latencyCands = {'TST_latency_s','TST_s','TailSuspension_s', ...
                'Hotplate_latency_s','HotPlate_latency_s', ...
                'TailImmersion_latency_s','TailImmersion_s'};
for cc = 1:numel(latencyCands)
    if hasVar(T, latencyCands{cc})
        D = addDayScalarFromT(D, T, latencyCands{cc});
    end
end

% ---------- mouse-level group (EVER passive?) ----------
if ~hasVar(T,'isPassive')
    warning('isPassive missing in CSV; cannot form Active-only vs Had-passive groups.');
    D.mouse_hasPassive = NaN(height(D),1);
else
    hasPassiveByMouse = groupsummary(T(:,{'mouse_key','isPassive'}),'mouse_key',@(x)any(x==1),'isPassive');
    hasPassiveByMouse.Properties.VariableNames(end) = {'hasPassive'};
    D = outerJoinInto(D, hasPassiveByMouse, {'mouse_key'}, 'hasPassive', 'mouse_hasPassive');
    D.mouse_hasPassive = double(D.mouse_hasPassive); % 1 = had-passive, 0 = active-only, NaN otherwise
end

% ---------- add 5 time bins ----------
D.bin5 = dayBin5(D.day_index);

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
% Reward median interval may appear as IRI or (rarely) IEI; map both
if hasVar(D,'rew_medianIRI_s'),    labels('rew_medianIRI_s')    = 'Reward median IRI (s)'; end
if hasVar(D,'rew_medianIEI_s'),    labels('rew_medianIEI_s')    = 'Reward median IRI (s)'; end

% Pupil / Requirement
if hasVar(D,'pupil_mean'),         labels('pupil_mean')         = 'Pupil mean (px)'; end
if hasVar(D,'RequirementLast'),    labels('RequirementLast')    = 'Requirement'; end

% Behaviour percentages: any TST_Pct_* / HOT_Pct_*
for v = string(D.Properties.VariableNames)
    if startsWith(v,'TST_Pct_'), labels(char(v)) = ['TST % ' stripPrefix(char(v))]; end
    if startsWith(v,'HOT_Pct_'), labels(char(v)) = ['HOT % ' stripPrefix(char(v))]; end
end

% Behaviour latencies
if hasVar(D,'Immersion_Latency_s'),           labels('Immersion_Latency_s') = 'Tail immersion latency (s)'; end
if hasVar(D,'Hotplate_latency_s'),            labels('Hotplate_latency (s)') = 'Hotplate latency (s)'; end
if hasVar(D,'HotPlate_latency_s'),            labels('HotPlate_latency_s')   = 'Hotplate latency (s)'; end
if hasVar(D,'TailImmersion_latency_s'),       labels('TailImmersion_latency (s)') = 'Tail-immersion latency (s)'; end
if hasVar(D,'TailImmersion_s'),               labels('TailImmersion_s')      = 'Tail-immersion latency (s)'; end
if hasVar(D,'TST_latency_s'),                 labels('TST_latency_s')        = 'TST latency (s)'; end
if hasVar(D,'TST_s'),                          labels('TST_s')               = 'TST latency (s)'; end
if hasVar(D,'TailSuspension_s'),               labels('TailSuspension_s')    = 'TST latency (s)'; end

metrics = labels.keys;
assert(~isempty(metrics), 'No analyzable metrics found in the CSV.');

%% -------- 1) GROUP comparison per bin (Active-only vs Had-passive) --------
binNames = {'D3-5','D6-8','D9-11','D12-14','D15-16'};

for i=1:numel(metrics)
    y   = metrics{i};
    lbl = labels(y);

    fig = figure('Color','w','Position',[80 80 1480 520]);
    tlo = tiledlayout(1,numel(binNames),'TileSpacing','compact','Padding','compact');
    csvRows = {};

    for b=1:numel(binNames)
        nexttile; hold on

        % per-mouse median within this bin; groups by "had passive ever?"
        sub = D(string(D.bin5)==binNames{b}, {'mouse_key','mouse_hasPassive',y});
        if isempty(sub) || ~any(isfinite(sub.(y)))
            axis off; title(binNames{b}); continue;
        end
        Mb = groupsummary(sub, {'mouse_key','mouse_hasPassive'}, 'median', y);
        Mb.Properties.VariableNames(end) = {'y_med'};

        Xa = Mb.y_med(Mb.mouse_hasPassive==0); Xa = Xa(isfinite(Xa));
        Xp = Mb.y_med(Mb.mouse_hasPassive==1); Xp = Xp(isfinite(Xp));

        % boxplot
        plotTwoBox(Xa, Xp, {'Active-only','Had-passive'}, lbl);
        title(binNames{b});
        yl = ylim(gca);

        % nonparam two-sample (Wilcoxon ranksum) — ALWAYS show p (even if n.s.)
        p = NaN;
        if ~isempty(Xa) && ~isempty(Xp)
            try, p = ranksum(Xa, Xp); catch, [~,p] = ttest2(Xa, Xp, 'Vartype','unequal'); end
        end
        addPBar(gca, 1, 2, yl(2), p, 'Active-only vs Had-passive'); % always prints p

        % record row
        csvRows(end+1,:) = {binNames{b}, numel(Xa), numel(Xp), nanmedian(Xa), nanmedian(Xp), p}; %#ok<AGROW>
    end

    title(tlo, sprintf('%s — Active-only vs Had-passive (per bin)', lbl));

    % save fig + csv
    savepng(fig, fullfile(outDir, sprintf('GROUP_active_vs_passive_%s.png', y))); close(fig);
    Tcsv = cell2table(csvRows, 'VariableNames', ...
        {'bin','n_activeOnly','n_hadPassive','median_activeOnly','median_hadPassive','p_ranksum'});
    try, writetable(Tcsv, fullfile(outDir, sprintf('GROUP_active_vs_passive_%s.csv', y))); catch, end
end

%% -------- 2) PERIOD comparison (one-way ANOVA across bins) --------
cohorts = { ...
    struct('name','ALL',         'mask', true(height(D),1)), ...
    struct('name','ACTIVEONLY',  'mask', D.mouse_hasPassive==0), ...
    struct('name','HADPASSIVE',  'mask', D.mouse_hasPassive==1)};

for ci = 1:numel(cohorts)
    C = cohorts{ci};
    for i=1:numel(metrics)
        y   = metrics{i};
        lbl = labels(y);

        sub = D(C.mask, {'mouse_key','bin5',y}); % use *all days* from this cohort
        if isempty(sub) || ~any(isfinite(sub.(y))), continue; end

        % per-mouse, per-bin medians so each mouse contributes once per bin
        M = groupsummary(sub, {'mouse_key','bin5'}, 'median', y);
        M.Properties.VariableNames(end) = {'y_med'};
        M = M(~isundefined(M.bin5) & isfinite(M.y_med), :);
        if isempty(M), continue; end

        % vectors for the one-way test
        X = double(M.y_med);
        G = cellstr(string(M.bin5));   % group labels (5 bins)

        % ----- run one-way test (ANOVA; fallback to Kruskal–Wallis) -----
        usedKW     = false;
        methodName = 'ANOVA';
        posthoc    = 'Tukey HSD';
        pairs = []; padj = []; grpNames = unique(G,'stable');

        try
            % ANOVA with stats for multcompare
            [p_anova, ~, st] = anova1(X, G, 'off');

            % Tukey-HSD post-hoc (preferred)
            try
                [mcTbl, ~, ~, gn] = multcompare(st, 'CType','hsd', 'Display','off');
                pairs    = mcTbl(:,1:2);     % [i j] group indices
                padj     = mcTbl(:,6);       % adjusted p-values
                grpNames = gn;               % cellstr group names
            catch
                % if multcompare unavailable, do rank-sum + BH-FDR
                [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
                posthoc = 'BH';
            end
        catch
            % fallback: Kruskal–Wallis + BH-FDR for pairwise
            usedKW    = true;
            methodName = 'KW';
            posthoc   = 'BH';
            p_anova   = kruskalwallis(X, G, 'off');
            [pairs, padj, grpNames] = fallbackPairsFDR(X, G);
        end

        % ----- plot -----
        fig = figure('Color','w','Position',[100 100 980 520]); hold on
        boxplot(X, G, 'Symbol','k+','Colors','k'); grid on; box on
        ylabel(lbl); xlabel('Time bin');
        title(sprintf('%s — one-way %s across bins (%s)  p=%.3g', ...
                      lbl, methodName, C.name, p_anova));

        % annotate ALL pairwise p (padj) with stars (and text), always shown
        yl = ylim; yBase = yl(2); dy = 0.05*range(yl); level = 0;
        for r=1:size(pairs,1)
            a = pairs(r,1); b = pairs(r,2);
            yTop = yBase + (level+1)*dy;
            txt  = sprintf('%s vs %s (%s p=%.3g)', grpNames{a}, grpNames{b}, posthoc, padj(r));
            addSigBar(gca, min(a,b), max(a,b), yTop, padj(r), txt, true);  % ALWAYS show p; stars if p<.05
            level = level + 1;
        end
        if size(pairs,1)>0, ylim([yl(1), yBase + (level+3)*dy]); end

        savepng(fig, fullfile(outDir, sprintf('ANOVA_bins_%s_%s.png', C.name, y))); close(fig);

        % write csvs
        pairStr = strings(numel(padj),1);
        for r=1:numel(padj)
            pairStr(r) = sprintf('%s vs %s', grpNames{pairs(r,1)}, grpNames{pairs(r,2)});
        end
        Tcsv = table(repmat(string(C.name),numel(padj),1), pairStr, padj, ...
                     'VariableNames', {'cohort','pair','p_adj'});
        Tcsv_an = table(string(C.name), usedKW, p_anova, ...
                        'VariableNames', {'cohort','used_KW','p_overall'});
        try
            writetable(Tcsv_an, fullfile(outDir, sprintf('ANOVA_bins_%s_%s_overall.csv', C.name, y)));
            writetable(Tcsv,    fullfile(outDir, sprintf('ANOVA_bins_%s_%s.csv',       C.name, y)));
        catch
        end
    end
end

fprintf('[stats-only] Done.\n');
end

%% ====================== helpers ======================

function tf = hasVar(T,var)
tf = ismember(var, T.Properties.VariableNames);
end

function T = ensureString(T, vname)
if hasVar(T,vname) && ~isstring(T.(vname)), T.(vname) = string(T.(vname)); end
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
        v = double(T.(cands{i}));
        if any(isfinite(v)), tb = v; return; end
    end
end
if hasVar(T,'Frame'), tb = double(T.Frame)/30; end
end

function r = finiteRange(x)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x)-min(x); end
end

function [n, meanDur, totalDur, medianIEI] = eventMetrics(tb_s, ttl)
ttl = logical(ttl(:)); tb_s = double(tb_s(:));
if numel(tb_s) < 2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
d = diff([false; ttl; false]); onIdx = find(d==1); offIdx = find(d==-1)-1;
n = numel(onIdx); if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
dt = [diff(tb_s); median(diff(tb_s),'omitnan')]; segDur = zeros(n,1);
for k=1:n, ii = onIdx(k):offIdx(k); segDur(k) = sum(dt(ii)); end
meanDur  = mean(segDur,'omitnan'); totalDur = sum(segDur,'omitnan');
onTimes = tb_s(onIdx); if numel(onTimes)>=2, medianIEI = median(diff(onTimes),'omitnan'); else, medianIEI = NaN; end
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
vals = G.(srcName);
D.(dstName)(tf) = vals(loc(tf));
end

function k = buildKey(T, keys)
k = repmat("", height(T), 1);
for i = 1:numel(keys)
    v = T.(keys{i}); if ~isstring(v), v = string(v); end
    if i==1, k=v; else, k = k + "|" + v; end
end
end

function c = dayBin5(day_index)
di = double(day_index);
lab = strings(size(di));
lab(di>=3  & di<=5 )  = "D3-5";
lab(di>=6  & di<=8 )  = "D6-8";
lab(di>=9  & di<=11)  = "D9-11";
lab(di>=12 & di<=14)  = "D12-14";
lab(di>=15 & di<=16)  = "D15-16";
c = categorical(lab, ["D3-5","D6-8","D9-11","D12-14","D15-16"]);
end

function plotTwoBox(Xa, Xp, xlabels, ylab)
% Robust 0/1/2-group box/whisker. Uses group names as the grouping variable.
data = []; grpNames = strings(0,1);
if ~isempty(Xa), data=[data; Xa(:)]; grpNames=[grpNames; repmat(string(xlabels{1}),numel(Xa),1)]; end
if ~isempty(Xp), data=[data; Xp(:)]; grpNames=[grpNames; repmat(string(xlabels{2}),numel(Xp),1)]; end
if isempty(data), axis off; title('No data'); return; end
boxplot(data, grpNames, 'Symbol','k+','Colors','k'); % no fixed 'Labels'
ylabel(ylab); grid on; box on
end

function addSigBar(ax, x1, x2, yTop, p, pairLabel, showAlways)
% Draw bracket + stars; ALWAYS show p if showAlways==true.
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
% Draw a bracket and ALWAYS print the p-value (adds stars if significant).
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

ylim([yl(1), y+3*dy]);  % leave space for labels
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
try, exportgraphics(fh, fn, 'Resolution',150); catch, print(fh, fn, '-dpng','-r150'); end
end

function [pairs, p_adj, names] = fallbackPairsFDR(X, G)
% If multcompare is unavailable, do all unpaired tests + BH-FDR.
names = unique(G,'stable'); k = numel(names);
pairs = nchoosek(1:k,2); p = nan(size(pairs,1),1);
for i=1:size(pairs,1)
    a = pairs(i,1); b = pairs(i,2);
    xa = X(strcmp(G, names{a})); xb = X(strcmp(G, names{b}));
    xa = xa(isfinite(xa)); xb = xb(isfinite(xb));
    try, p(i) = ranksum(xa, xb); catch, [~,p(i)] = ttest2(xa, xb, 'Vartype','unequal'); end
end
p_adj = fdr_bh(p);
end

function p_adj = fdr_bh(p)
p = p(:); n = numel(p);
[ps,idx] = sort(p); adj = nan(size(p));
adj(idx) = min( cummin( (n./(1:n)').*ps ), 1);
p_adj = adj;
end
