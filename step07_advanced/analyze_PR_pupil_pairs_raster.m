function analyze_PR_pupil_pairs_raster()
% One-stop script for:
%  - PR spaghetti (active mice, days 1–10)
%  - Pupil & immersion spaghetti (all mice, plus active vs passive panels)
%  - Day5 vs Day9 TST/HOT comparisons (All + Active + Passive; paired t-tests)
%  - 4-way ANOVA (D5A vs D5P vs D9A vs D9P) for all TST/HOT metrics
%  - Pair-wise Active vs Passive zoom (pupil, immersion, TST/HOT key metrics)
%  - Active-only licking spaghetti across days (incl. lick freq /10 s)
%  - Reward & lick rasters (days 4,5,9,10; active only, y-axis = mouse IDs)
%  - Per-day and early vs late mean cumulative licking curves (active only)
%
% Uses ALL_mice_longitudinal.csv produced by the aggregator.

%% ---------- locate latest run + load ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder,d(ix).name);
csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% Normalise key columns
T = ensureString(T,'mouse_key');
T = ensureString(T,'day_name');
T = ensureString(T,'Session_Paradigm');
if ~ismember('isPassive',T.Properties.VariableNames)
    T.isPassive = nan(height(T),1);
elseif ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% TTL numerics
if ismember('Lick_TTL',T.Properties.VariableNames)
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
end
if ismember('Injector_TTL',T.Properties.VariableNames)
    T.Injector_TTL(isnan(T.Injector_TTL)) = 0;
end

%% ---------- per-session metrics, then per-day medians ----------
[S, D] = fast_session_day_metrics_basic(T, runDir);

% classify had-passive / active (mouse-level)
D.HadPassive = classifyHadPassive(S, D);
D.Group      = categorical(D.HadPassive, [false true], {'Active','Passive'});

daysAll = 1:10;

%% ---------- output dir ----------
outDir = fullfile(runDir,'figs','PR_pupil_pairs_raster');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% 1. PR spaghetti (RequirementLast) for active mice (days 1–10)
fprintf('Plotting PR spaghetti for Active mice...\n');
plot_spaghetti_metric(D, 'RequirementLast', 'PR / RequirementLast', ...
    'PR (RequirementLast) – Active mice', outDir, 'ActiveOnly', daysAll);

% PR: Active (top) vs Passive (bottom) 2-panel spaghetti
if ismember('RequirementLast', D.Properties.VariableNames)
    plot_spaghetti_metric_groups(D, 'RequirementLast', 'PR / RequirementLast', ...
        'PR (RequirementLast): Active (top) vs Passive (bottom)', ...
        fullfile(outDir,'spaghetti_PR_active_vs_passive.png'), ...
        daysAll, false);
end

%% 2. Pupil spaghetti & immersion spaghetti
fprintf('Plotting pupil & immersion spaghetti...\n');

% 2A. Pupil
if ismember('pupil_mean', D.Properties.VariableNames)
    % all mice
    plot_spaghetti_metric(D, 'pupil_mean', 'Pupil diameter (px)', ...
        'Pupil diameter (px) – all mice', outDir, 'All', daysAll);

    % Active vs Passive panels + mouse labels
    plot_spaghetti_metric_groups(D, 'pupil_mean', 'Pupil diameter (px)', ...
        'Pupil diameter (px): Active (top) vs Passive (bottom)', ...
        fullfile(outDir,'spaghetti_pupil_active_vs_passive.png'), ...
        daysAll, true);
end

% 2B. Tail immersion
% 2B. Tail immersion
% 2B. Tail immersion
if ismember('Immersion_Latency_s', D.Properties.VariableNames)
    % all mice
    plot_spaghetti_metric(D, 'Immersion_Latency_s', 'Immersion latency (s)', ...
        'Immersion latency (s) – all mice', outDir, 'All', daysAll);

    % Active vs Passive panels
    plot_spaghetti_metric_groups(D, 'Immersion_Latency_s', 'Immersion latency (s)', ...
        'Immersion latency (s): Active (top) vs Passive (bottom)', ...
        fullfile(outDir,'spaghetti_immersion_active_vs_passive.png'), ...
        daysAll, false);

    % ---- NEW: per-day Active vs Passive stats with stars ----
    daysImm = 5:10;
    compare_AP_means_per_day_with_stars( ...
        D, 'Immersion_Latency_s', daysImm, outDir);

    % ---- NEW: within-group Friedman + post-hoc (optional; keep if you want) ----
    within_group_friedman_and_posthoc( ...
        D, 'Immersion_Latency_s', daysImm, outDir);
end

%% 3. Day5 vs Day9 TST/HOT: paired t-tests (All + Active + Passive)
fprintf('Day5 vs Day9 comparisons for all TST/HOT metrics...\n');
dayA = 5; dayB = 9;

tstCols = D.Properties.VariableNames(startsWith(D.Properties.VariableNames,'TST_'));
hotCols = D.Properties.VariableNames(startsWith(D.Properties.VariableNames,'HOT_'));
allCols = [tstCols, hotCols];

for i = 1:numel(allCols)
    mName = allCols{i};
    if endsWith(mName,'_File'), continue; end
    if ~isnumeric(D.(mName)),   continue; end

    % All mice
    whisker_paired_day_compare(D, mName, dayA, dayB, outDir, 'All');

    % Active only
    if any(D.Group=="Active")
        whisker_paired_day_compare(D(D.Group=="Active",:), ...
            mName, dayA, dayB, outDir, 'Active');
    end

    % Passive only
    if any(D.Group=="Passive")
        whisker_paired_day_compare(D(D.Group=="Passive",:), ...
            mName, dayA, dayB, outDir, 'Passive');
    end
end

%% 4. 4-way ANOVA: D5A vs D5P vs D9A vs D9P (TST/HOT)
fprintf('Running 4-group ANOVA (D5A/D5P/D9A/D9P) for TST/HOT metrics...\n');
for i = 1:numel(allCols)
    mName = allCols{i};
    if endsWith(mName,'_File'), continue; end
    if ~isnumeric(D.(mName)),   continue; end
    anova4_day5_day9_active_passive(D, mName, dayA, dayB, outDir);
end
%% 4b. Non-parametric unpaired comparisons (Active vs Passive)
%  Metrics: tail immersion, HOT non-moving, TST non-moving
metricsNP = {'Immersion_Latency_s','HOT_Pct_Non_moving','TST_Pct_Non_moving'};
for i = 1:numel(metricsNP)
    mName = metricsNP{i};
    if ~ismember(mName, D.Properties.VariableNames), continue; end
    % Day 5 Active vs Passive
    unpaired_day_compare_nonparam(D, mName, 5, outDir);
    % Day 9 Active vs Passive
    unpaired_day_compare_nonparam(D, mName, 9, outDir);
end

%% 5. Pair-wise Active vs Passive zoom (selected mouse pairs)
fprintf('Generating pair-wise zoom plots...\n');

pairs = { ...
    {'2687_blue','2687_normal'}, ...
    {'8055_blue','8055_normal'}, ...
    {'5942_normal','5942_blue','2705_normal'}, ... % one active + two passive
    {'2705_black','2705_blue'}, ...
    {'8053_blue','8053_normal'} ...
    };
pairLabels = { ...
    '2687 pair (blue active vs normal passive)', ...
    '8055 pair (blue active vs normal passive)', ...
    '5942/2705 triple (5942_normal active; others passive)', ...
    '2705 pair (black active vs blue passive)', ...
    '8053 pair (blue active vs normal passive)' ...
    };

pairMetrics = { ...
    'pupil_mean', ...
    'Immersion_Latency_s', ...
    'TST_Pct_Non_moving', ...
    'HOT_Pct_Licking', ...
    'HOT_Pct_Flinching' ...
    };

for p = 1:numel(pairs)
    mkList = pairs{p};
    for mi = 1:numel(pairMetrics)
        mName = pairMetrics{mi};
        if ~ismember(mName, D.Properties.VariableNames), continue; end
        pair_zoom_plot(D, mkList, pairLabels{p}, mName, outDir);
    end
end

%% 6. Active-only licking spaghetti across days (incl. /10s)
fprintf('Plotting active-only licking spaghetti across days...\n');

plot_active_lick_spaghetti_from_S(S, D, daysAll, outDir);

%% 7. Reward & lick rasters + cumulative licking curves (active only)
fprintf('Extracting lick/reward events for active mice...\n');

[timeGridMin, meanCumPerDay, meanCumEarly, meanCumLate, sessInfo] = ...
    compute_cumulative_and_events_active(T, D, daysAll);

% Reward & lick rasters for days 4,5,9,10
rasterDays = [4 5 9 10];
for dDay = rasterDays
    raster_day_plot(sessInfo, dDay, 'reward', outDir);
    raster_day_plot(sessInfo, dDay, 'lick',   outDir);
end

% Per-day mean cumulative licking curves (active only; one line per day)
plot_cumLicks_allDays(timeGridMin, meanCumPerDay, daysAll, outDir);

% Early (1–5) vs late (6–10) mean cumulative curves
plot_cumLicks_early_late(timeGridMin, meanCumEarly, meanCumLate, outDir);

fprintf('Done. All figures are under:\n  %s\n', outDir);
end

%% ================= helpers =================

function T2 = ensureString(T2, nm)
    if ~ismember(nm, T2.Properties.VariableNames)
        T2.(nm) = repmat("",height(T2),1);
        return;
    end
    if ~isstring(T2.(nm))
        T2.(nm) = string(T2.(nm));
    end
end

function within_group_friedman_and_posthoc(D, mName, daysRange, outDir)
% Within-group time effects for Active / Passive separately.
% - Friedman test over days (repeated measures, non-parametric)
% - Post-hoc paired Wilcoxon sign-rank for Day 5 vs Day 9
%   (FDR-corrected within each group)
%
% Outputs:
%   friedman_<metric>_active_passive.csv
%   friedman_posthoc_<metric>_active.csv
%   friedman_posthoc_<metric>_passive.csv

    if ~ismember('Group', D.Properties.VariableNames)
        fprintf('within_group_friedman_and_posthoc: no Group column, skipping.\n');
        return;
    end
    if ~ismember(mName, D.Properties.VariableNames)
        fprintf('within_group_friedman_and_posthoc: %s not in D, skipping.\n', mName);
        return;
    end

    groups = ["Active","Passive"];
    friedRes = table();

    for gi = 1:numel(groups)
        gName = groups(gi);
        G = D(D.Group==gName & ismember(D.day_index,daysRange), ...
              {'mouse_key','day_index',mName});
        G = G(isfinite(G.(mName)),:);
        if isempty(G)
            fprintf('Friedman %s %s: no data.\n', mName, gName);
            continue;
        end

        mice = unique(G.mouse_key,'stable');
        nMice = numel(mice);
        nDays = numel(daysRange);

        % Build subject x day matrix
        X = nan(nMice, nDays);
        for i = 1:nMice
            for j = 1:nDays
                d  = daysRange(j);
                v  = G.(mName)(G.mouse_key==mice(i) & G.day_index==d);
                if ~isempty(v)
                    X(i,j) = mean(v,'omitnan');
                end
            end
        end

        % Require complete cases across all days
        goodMice = all(isfinite(X),2);
        X2 = X(goodMice,:);
        miceGood = mice(goodMice);

        if size(X2,1) < 3
            fprintf('Friedman %s %s: not enough complete mice (N=%d), skipping.\n', ...
                mName, gName, size(X2,1));
            continue;
        end

        % Friedman (rows = subjects, columns = days)
        [pF, ~, ~] = friedman(X2, 1, 'off');
        fprintf('Friedman %s %s: p=%.3g (N=%d, days=%d)\n', ...
            mName, gName, size(X2,1), size(X2,2));

        newRow = table(string(gName), pF, size(X2,1), ...
            'VariableNames',{'Group','p_Friedman','N_mice'});
        friedRes = [friedRes; newRow]; %#ok<AGROW>

        % ----- Post-hoc: Day5 vs Day9 (you can add more pairs here) -----
        pairs = [5 9];   % [Day1 Day2]; modify if you want more pairs
        postTbl = table();

        for k = 1:size(pairs,1)
            d1 = pairs(k,1); d2 = pairs(k,2);
            if ~ismember(d1, daysRange) || ~ismember(d2, daysRange)
                continue;
            end
            j1 = find(daysRange==d1); j2 = find(daysRange==d2);
            x1 = X2(:,j1);
            x2 = X2(:,j2);

            good = isfinite(x1) & isfinite(x2);
            x1 = x1(good); x2 = x2(good);
            if numel(x1) < 3
                continue;
            end
            pRaw = signrank(x1,x2);
            row  = table(string(gName), d1, d2, numel(x1), pRaw, ...
                'VariableNames',{'Group','Day1','Day2','N_mice','p_raw'});
            postTbl = [postTbl; row]; %#ok<AGROW>
        end

        if ~isempty(postTbl)
            % FDR within this group's pairs
            p = postTbl.p_raw;
            [p_sorted, ord] = sort(p);
            m = numel(p);
            bh = p_sorted .* (m ./ (1:m)');
            bh = min(bh,1);
            bh = cummin(bh(end:-1:1));
            bh = bh(end:-1:1);
            p_adj = nan(size(p));
            p_adj(ord) = bh;

            postTbl.p_adj = p_adj;

            star = strings(size(p));
            for ii = 1:numel(p)
                if      p_adj(ii) < 0.001, star(ii) = "***";
                elseif  p_adj(ii) < 0.01,  star(ii) = "**";
                elseif  p_adj(ii) < 0.05,  star(ii) = "*";
                else,                       star(ii) = "n.s.";
                end
            end
            postTbl.Star = star;

            fnPost = fullfile(outDir, sprintf('friedman_posthoc_%s_%s.csv', ...
                safeName(mName), lower(char(gName))));
            writetable(postTbl, fnPost);
        end
    end

    if ~isempty(friedRes)
        fnF = fullfile(outDir, sprintf('friedman_%s_active_passive.csv', ...
            safeName(mName)));
        writetable(friedRes, fnF);
    end
end

%% ---------- epoch / passive classification helpers ----------
function HadPassive = classifyHadPassive(S, D)
    mice = unique(S.mouse_key,'stable');
    HadPassive = false(height(D),1);
    havePar = ismember('Session_Paradigm', S.Properties.VariableNames);
    for i = 1:numel(mice)
        rS  = S.mouse_key==mice(i);
        ip  = S.isPassive(rS);
        di  = S.day_index(rS);
        par = strings(nnz(rS),1);
        if havePar
            p = S.Session_Paradigm(rS);
            if iscategorical(p), p = string(p); end
            par = string(p);
        end
        flag = false;
        if any(ip==1,'all')
            flag = true;
        elseif havePar && any(contains(lower(par),'passive'),'all')
            flag = true;
        elseif any(ismember(double(di),6:10) & (ip>=1 | isnan(ip)),'all')
            flag = true;
        end
        Dmask = D.mouse_key==mice(i);
        HadPassive(Dmask) = flag;
    end
end

%% ---------- generic spaghetti plotting ----------
function compare_AP_means_per_day_with_stars(D, mName, daysAll, outDir)
% Per-day Active vs Passive comparison across days.
% - Wilcoxon rank-sum (Mann–Whitney) per day
% - Benjamini–Hochberg FDR across days
% - Mean ± SEM line figure with text labels:
%       "*, p=0.012, q=0.030"  or  "n.s., p=0.53, q=0.93"
% - CSV with all numbers

    fprintf('Running compare_AP_means_per_day_with_stars for %s\n', mName);

    if ~ismember('Group', D.Properties.VariableNames)
        fprintf('  No Group column, skipping.\n');
        return;
    end
    if ~ismember(mName, D.Properties.VariableNames)
        fprintf('  Metric %s not found, skipping.\n', mName);
        return;
    end

    nDays = numel(daysAll);
    meanA = nan(nDays,1);
    meanP = nan(nDays,1);
    semA  = nan(nDays,1);
    semP  = nan(nDays,1);
    nA    = nan(nDays,1);
    nP    = nan(nDays,1);
    p_raw = nan(nDays,1);

    % ---------- per-day ranksum (unpaired) + mean/SEM ----------
    for i = 1:nDays
        d  = daysAll(i);
        xa = D.(mName)(D.day_index==d & D.Group=="Active");
        xp = D.(mName)(D.day_index==d & D.Group=="Passive");
        xa = xa(isfinite(xa));
        xp = xp(isfinite(xp));

        if numel(xa) >= 2 && numel(xp) >= 2
            nA(i)    = numel(xa);
            nP(i)    = numel(xp);
            meanA(i) = mean(xa,'omitnan');
            meanP(i) = mean(xp,'omitnan');
            semA(i)  = std(xa,'omitnan') / sqrt(nA(i));
            semP(i)  = std(xp,'omitnan') / sqrt(nP(i));
            p_raw(i) = ranksum(xa,xp);      % Mann–Whitney
        else
            fprintf('  Day %d: not enough data (NA=%d, NP=%d).\n', ...
                d, numel(xa), numel(xp));
        end
    end

    % ---------- BH-FDR across days (on raw p) ----------
    p_adj = nan(size(p_raw));
    valid = find(isfinite(p_raw) & p_raw>0);
    if ~isempty(valid)
        [p_sorted, ord] = sort(p_raw(valid));
        m = numel(p_sorted);
        q = p_sorted .* (m ./ (1:m)');   % BH
        q(q>1) = 1;
        % monotone
        for k = m-1:-1:1
            q(k) = min(q(k), q(k+1));
        end
        p_adj(valid(ord)) = q;
    end

    % ---------- star codes from q (FDR) ----------
    star = strings(nDays,1);
    for i = 1:nDays
        if ~isfinite(p_adj(i)), continue; end
        if      p_adj(i) < 0.001, star(i) = "***";
        elseif  p_adj(i) < 0.01,  star(i) = "**";
        elseif  p_adj(i) < 0.05,  star(i) = "*";
        else,                     star(i) = "n.s.";
        end
    end

    % ---------- CSV summary ----------
    resTbl = table(daysAll(:), nA, meanA, semA, ...
                   nP, meanP, semP, p_raw, p_adj, star, ...
        'VariableNames', {'day_index','N_Active','Mean_Active','SEM_Active', ...
                          'N_Passive','Mean_Passive','SEM_Passive', ...
                          'p_raw','p_adj','Star'});

    fnCsv = fullfile(outDir, ...
        sprintf('pvals_active_vs_passive_%s_ranksum.csv', safeName(mName)));
    writetable(resTbl, fnCsv);

    % ---------- FIGURE: mean ± SEM lines + p / q labels ----------
    fh1 = figure('Color','w','Position',[120 120 700 500]); hold on;

    % Active (blue-ish) / Passive (red-ish) with error bars
    ha = errorbar(daysAll, meanA, semA, '-o', ...
        'LineWidth',2, 'MarkerSize',5, ...
        'DisplayName','Active');
    hp = errorbar(daysAll, meanP, semP, '-o', ...
        'LineWidth',2, 'MarkerSize',5, ...
        'DisplayName','Passive');

    % optional: set colors if you want fixed colors
    set(ha,'Color',[0.00 0.45 0.74]);
    set(hp,'Color',[0.85 0.33 0.10]);

    xlabel('Day index');
    ylabel(strrep(mName,'_',' '), 'Interpreter','none');
    title(sprintf('%s: Active vs Passive means per day (ranksum, BH-FDR)', mName), ...
          'Interpreter','none');
    legend('Location','best');
    grid on; box on;

    yAll = [meanA-semA; meanA+semA; meanP-semP; meanP+semP];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        yMin = 0; yMax = 1;
    else
        yMin = min(yAll); yMax = max(yAll);
    end
    yRange = max(yMax - yMin, eps);
    ylim([yMin - 0.05*yRange, yMax + 0.35*yRange]);

    for i = 1:nDays
        if ~isfinite(p_raw(i)) || ~isfinite(p_adj(i)), continue; end

        % p 텍스트 만들기
        if p_raw(i) < 0.001
            pStr = 'p<0.001';
        else
            pStr = sprintf('p=%.3f', p_raw(i));
        end

        if p_adj(i) < 0.001
            qStr = 'q<0.001';
        else
            qStr = sprintf('q=%.3f', p_adj(i));
        end

        starStr = star(i);   % '***', '**', '*', or 'n.s.'

        % label: 한 day당 2줄 텍스트 (예:  "*\np=0.012, q=0.030")
        labelStr = sprintf('%s\n%s, %s', starStr, pStr, qStr);

        yStar = max([meanA(i)+semA(i), meanP(i)+semP(i)]);
        if ~isfinite(yStar), continue; end
        yStar = yStar + 0.08*yRange;

        text(daysAll(i), yStar, labelStr, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',9, ...
            'FontWeight', iff(p_adj(i)<0.05,'bold','normal'));
    end

    fn1 = fullfile(outDir, ...
        sprintf('means_active_vs_passive_%s_ranksum.png', safeName(mName)));
    printpng(fh1, fn1); close(fh1);
end


function w = iff(cond, a, b)
    if cond, w = a; else, w = b; end
end


function plot_spaghetti_metric(D, ycol, ylab, ttl, outDir, modeStr, restrictDays)
    % modeStr: 'All' or 'ActiveOnly'
    if ~ismember(ycol, D.Properties.VariableNames), return; end

    D = D(ismember(D.day_index, restrictDays), :);
    if strcmpi(modeStr,'ActiveOnly') && ismember('Group', D.Properties.VariableNames)
        D = D(D.Group=="Active", :);
    end
    if isempty(D), return; end

    mice = unique(D.mouse_key,'stable');

    fh = figure('Color','w','Position',[80 80 800 500]); hold on;
    for i = 1:numel(mice)
        di = D(D.mouse_key==mice(i),:);
        di = di(isfinite(di.(ycol)),:);
        if isempty(di), continue; end
        [~,ord] = sort(di.day_index);
        plot(di.day_index(ord), di.(ycol)(ord),'-o', ...
            'LineWidth',1.0,'MarkerSize',4);
    end

    % mean across mice per day
    G = groupsummary(D(:,{'day_index',ycol}), 'day_index','mean', ycol);
    plot(G.day_index, G.(['mean_' ycol]), 'k-','LineWidth',3);

    xlabel('Day index');
    ylabel(ylab);
    title(ttl,'Interpreter','none');
    grid on; box on;
    xlim([min(restrictDays)-0.5, max(restrictDays)+0.5]);

    fn = fullfile(outDir, sprintf('spaghetti_%s_%s.png', modeStr, safeName(ycol)));
    printpng(fh, fn); close(fh);
end

function plot_spaghetti_metric_groups(D, ycol, ylab, ttl, outFile, restrictDays, labelLines)
    if ~ismember(ycol, D.Properties.VariableNames), return; end
    if nargin < 7, labelLines = false; end

    D = D(ismember(D.day_index,restrictDays),:);

    fh = figure('Color','w','Position',[80 80 800 900]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    groups = ["Active","Passive"];
    for gi = 1:2
        nexttile; hold on;
        if ~ismember('Group', D.Properties.VariableNames)
            text(0.5,0.5,'No Group column','Units','normalized');
            axis off; continue;
        end

        G = D(D.Group==groups(gi),:);
        if isempty(G)
            text(0.5,0.5,sprintf('No %s mice',groups(gi)),'Units','normalized');
            axis off; continue;
        end

        mice = unique(G.mouse_key,'stable');
        miceStr = string(mice);
        labX = nan(numel(mice),1);
        labY = nan(numel(mice),1);

        for i = 1:numel(mice)
            di = G(G.mouse_key==mice(i),:);
            di = di(isfinite(di.(ycol)),:);
            if isempty(di), continue; end
            [~,ord] = sort(di.day_index);
            x = di.day_index(ord);
            y = di.(ycol)(ord);
            plot(x, y, '-o', 'LineWidth',1.0,'MarkerSize',4);
            if labelLines
                lastIdx = find(isfinite(y),1,'last');
                if ~isempty(lastIdx)
                    labX(i) = x(lastIdx);
                    labY(i) = y(lastIdx);
                end
            end
        end

        % mean line
        Sg = groupsummary(G(:,{'day_index',ycol}),'day_index','mean',ycol);
        plot(Sg.day_index, Sg.(['mean_' ycol]), 'k-','LineWidth',3);

        % mouse_key labels
        if labelLines
            for i = 1:numel(mice)
                if ~isfinite(labX(i)), continue; end
                text(labX(i)+0.1, labY(i), char(miceStr(i)), ...
                     'FontSize',8, 'Interpreter','none');
            end
        end

        xlabel('Day index');
        ylabel(ylab);
        title(sprintf('%s – %s mice', ylab, groups(gi)),'Interpreter','none');
        grid on; box on;
        xlim([min(restrictDays)-0.5, max(restrictDays)+0.5]);
    end

    sgtitle(ttl,'FontWeight','bold');
    printpng(fh, outFile); close(fh);
end

%% ---------- Day5 vs Day9 paired whisker + lines ----------

function whisker_paired_day_compare(D, mName, dayA, dayB, outDir, groupLabel)

    if nargin < 6 || isempty(groupLabel), groupLabel = 'All'; end

    if ~ismember(mName, D.Properties.VariableNames)
        fprintf('Metric %s not found in D, skipping.\n', mName);
        return;
    end

    Ta = D(D.day_index == dayA, {'mouse_key', mName});
    Tb = D(D.day_index == dayB, {'mouse_key', mName});

    Ta = Ta(isfinite(Ta.(mName)), :);
    Tb = Tb(isfinite(Tb.(mName)), :);

    if isempty(Ta) || isempty(Tb)
        fprintf('No data for %s on day %d or %d (%s), skipping.\n', ...
            mName, dayA, dayB, groupLabel);
        return;
    end

    try
        M = outerjoin(Ta, Tb, ...
            'Keys', 'mouse_key', ...
            'MergeKeys', true, ...
            'Type', 'inner', ...
            'Suffixes', {'_A','_B'});
    catch
        M = outerjoin(Ta, Tb, ...
            'Keys', 'mouse_key', ...
            'MergeKeys', true, ...
            'Type', 'inner');
    end

    vnames = string(M.Properties.VariableNames);
    cand = vnames(startsWith(vnames, mName + "_"));
    if numel(cand) >= 2
        leftCol  = char(cand(1));
        rightCol = char(cand(2));
    else
        if ismember(mName + "_A", vnames) && ismember(mName + "_B", vnames)
            leftCol  = char(mName + "_A");
            rightCol = char(mName + "_B");
        else
            error('Could not find paired columns for metric %s in joined table.', mName);
        end
    end

    xa = M.(leftCol);
    xb = M.(rightCol);

    good = isfinite(xa) & isfinite(xb);
    xa = xa(good);
    xb = xb(good);

    if numel(xa) < 3
        fprintf('Not enough paired mice for %s (%s, N=%d), skipping.\n', ...
            mName, groupLabel, numel(xa));
        return;
    end

    [~, p] = ttest(xa, xb);

    fh = figure('Color','w','Position',[120 120 500 420]); hold on;
    X = [ones(size(xa)); 2*ones(size(xb))];
    Y = [xa; xb];
    boxchart(X, Y, 'BoxWidth', 0.3,'BoxFaceColor',[0.7 0.85 1]);

    for i = 1:numel(xa)
        plot([1 2], [xa(i) xb(i)], '-o', ...
            'Color', [0.6 0.6 0.6], ...
            'MarkerFaceColor', [0.3 0.3 0.3], ...
            'MarkerSize', 4);
    end

    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], ...
             'XTickLabel', {sprintf('Day %d', dayA), sprintf('Day %d', dayB)});
    ylabel(strrep(mName,'_',' '), 'Interpreter','none');
    title(sprintf('%s (%s): Day%d vs Day%d (N=%d, p=%.3g)', ...
          mName, groupLabel, dayA, dayB, numel(xa), p), ...
          'Interpreter','none');
    grid on; box on;

    fn = fullfile(outDir, ...
        sprintf('whisker_paired_%s_day%d_vs_day%d_%s.png', ...
            safeName(mName), dayA, dayB, safeName(groupLabel)));
    printpng(fh, fn); close(fh);
end

%% ---------- 4-group ANOVA helper ----------

function anova4_day5_day9_active_passive(D, mName, dayA, dayB, outDir)

    if ~ismember(mName, D.Properties.VariableNames)
        return;
    end
    if ~ismember('Group', D.Properties.VariableNames)
        return;
    end

    useDays  = [dayA dayB];
    useGroup = ["Active","Passive"];

    Tsub = D(ismember(D.day_index,useDays) & ismember(D.Group,useGroup), ...
             {'mouse_key','day_index','Group',mName});

    Tsub = Tsub(isfinite(Tsub.(mName)),:);
    if height(Tsub) < 4
        fprintf('ANOVA: not enough data for %s, skipping.\n', mName);
        return;
    end

    cond4 = strings(height(Tsub),1);
    for i = 1:height(Tsub)
        d  = Tsub.day_index(i);
        gg = string(Tsub.Group(i));
        cond4(i) = sprintf('Day%d_%s', d, gg);
    end
    cond4 = categorical(cond4);

    y = Tsub.(mName);

    [p, tbl, stats] = anova1(y, cond4, 'off'); %#ok<ASGLU>
    if isnan(p)
        fprintf('ANOVA failed for %s (NaN p).\n', mName);
        return;
    end

    % post-hoc with Tukey-Kramer
    post = multcompare(stats, 'CType','tukey-kramer','Display','off');

    % figure: 4-group boxplot + jitter
    fh = figure('Color','w','Position',[120 120 560 420]); hold on;
    boxchart(cond4, y, 'BoxFaceColor',[0.8 0.9 1]);
    gCats = categories(cond4);
    for i = 1:numel(gCats)
        idx = cond4 == gCats{i};
        xj = i + (rand(sum(idx),1)-0.5)*0.15;
        plot(xj, y(idx), 'ko','MarkerFaceColor',[0.2 0.2 0.2],'MarkerSize',4);
    end
    ylabel(strrep(mName,'_',' '), 'Interpreter','none');
    title(sprintf('%s: 4-group ANOVA (p=%.3g)', mName, p),'Interpreter','none');
    grid on; box on;

    fnFig = fullfile(outDir, ...
        sprintf('anova4_%s_D5A_D5P_D9A_D9P.png', safeName(mName)));
    printpng(fh, fnFig); close(fh);

    % save post-hoc table
    postTbl = array2table(post, ...
        'VariableNames',{'Group1','Group2','LowerCI','Diff','UpperCI','p_adj'});
    for k = 1:height(postTbl)
        postTbl.Group1(k) = str2double(string(stats.gnames(post(k,1))));
        postTbl.Group2(k) = str2double(string(stats.gnames(post(k,2))));
    end
    % store names as strings as well
    postTbl.Group1Name = stats.gnames(post(:,1));
    postTbl.Group2Name = stats.gnames(post(:,2));

    fnCsv = fullfile(outDir, ...
        sprintf('anova4_%s_D5A_D5P_D9A_D9P_posthoc.csv', safeName(mName)));
    writetable(postTbl, fnCsv);
end

%% ---------- Pair-wise zoom ----------

function pair_zoom_plot(D, mkList, pairLabel, mName, outDir)
    Dp = D(ismember(D.mouse_key, mkList), :);
    if isempty(Dp) || ~ismember(mName, Dp.Properties.VariableNames), return; end

    colors = lines(numel(mkList));
    fh = figure('Color','w','Position',[80 80 800 500]); hold on;

    for i = 1:numel(mkList)
        mk = mkList{i};
        di = Dp(Dp.mouse_key==mk,:);
        di = di(isfinite(di.(mName)),:);
        if isempty(di), continue; end
        [~,ord] = sort(di.day_index);
        lab = char(mk);
        if ismember('Group',di.Properties.VariableNames)
            gstr = char(unique(string(di.Group)));
            lab = sprintf('%s (%s)', mk, gstr);
        end
        plot(di.day_index(ord), di.(mName)(ord),'-o', ...
            'Color',colors(i,:), 'LineWidth',1.5,'MarkerSize',5, ...
            'DisplayName',lab);
    end

    xlabel('Day index');
    ylabel(strrep(mName,'_','\_'));
    title({pairLabel; strrep(mName,'_','\_')},'Interpreter','none');
    grid on; box on;
    legend('Location','bestoutside');

    fn = fullfile(outDir, sprintf('pair_%s_%s.png', safeName(pairLabel), safeName(mName)));
    printpng(fh, fn); close(fh);
end

%% ---------- Active-only licking spaghetti from S ----------

function plot_active_lick_spaghetti_from_S(S, D, daysAll, outDir)

    % define active mice from D
    if ~ismember('Group', D.Properties.VariableNames)
        fprintf('No Group column in D; skipping active-only lick spaghetti.\n');
        return;
    end
    activeMice = unique(D.mouse_key(D.Group=="Active"),'stable');
    if isempty(activeMice)
        fprintf('No Active mice in D; skipping active-only lick spaghetti.\n');
        return;
    end

    % restrict S to active mice and requested days
    S2 = S(ismember(S.mouse_key,activeMice) & ismember(S.day_index,daysAll), :);
    if isempty(S2)
        fprintf('No sessions for Active mice in specified days; skipping.\n');
        return;
    end

    % build per-mouse/day means from S (lick metrics)
    lickVars = {'lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
                'lick_totalDur_s','lick_n','lick_freq_per_10s'};
    [g, mk, di] = findgroups(S2.mouse_key, S2.day_index);
    L = table(removecats(mk), double(di), ...
        'VariableNames',{'mouse_key','day_index'});
    for k = 1:numel(lickVars)
        v = lickVars{k};
        if ismember(v, S2.Properties.VariableNames)
            L.(v) = splitapply(@(x) mean(x,'omitnan'), S2.(v), g);
        else
            L.(v) = nan(numel(mk),1);
        end
    end

    % add cumulative across days using lick_n
    L = add_cumulative_across_days(L);

    lickMetrics = { ...
        'lick_freq_per_min', ...
        'lick_meanDur_s', ...
        'lick_medianIEI_s', ...
        'lick_totalDur_s', ...
        'lick_n', ...
        'lick_cumAcrossDays', ...
        'lick_freq_per_10s' ...
        };
    lickLabels  = { ...
        'Lick frequency (/min)', ...
        'Lick mean duration (s)', ...
        'Lick median IEI (s)', ...
        'Total lick duration (s)', ...
        'Total lick count (per day)', ...
        'Cumulative lick count (across days)', ...
        'Lick frequency (/10 s)' ...
        };

    for k = 1:numel(lickMetrics)
        ycol = lickMetrics{k};
        ylab = lickLabels{k};
        if ~ismember(ycol, L.Properties.VariableNames), continue; end
        plot_spaghetti_metric(L, ycol, ylab, ...
            sprintf('%s — Active mice', ylab), outDir, 'ActiveOnly', daysAll);
    end
end


%% ---------- Active-only cumulative across days (helper) ----------

function Dact2 = add_cumulative_across_days(Dact)
    Dact2 = Dact;
    if ~ismember('lick_n',Dact2.Properties.VariableNames)
        Dact2.lick_cumAcrossDays = nan(height(Dact2),1);
        return;
    end
    Dact2.lick_cumAcrossDays = nan(height(Dact2),1);
    mice = unique(Dact2.mouse_key,'stable');
    for i = 1:numel(mice)
        idx = find(Dact2.mouse_key==mice(i));
        di  = Dact2.day_index(idx);
        [~,ord] = sort(di);
        lickN = Dact2.lick_n(idx(ord));
        cumN  = cumsum(lickN);
        Dact2.lick_cumAcrossDays(idx(ord)) = cumN;
    end
end

%% ---------- Cumulative licks and event extraction ----------

function [timeGridMin, meanCumPerDay, meanCumEarly, meanCumLate, sessInfo] = ...
    compute_cumulative_and_events_active(T, D, daysAll)

    if ismember('Group', D.Properties.VariableNames)
        activeMask = (D.Group == "Active");
    else
        activeMask = false(height(D),1);
    end
    activeMice = unique(string(D.mouse_key(activeMask)), 'stable');

    mkT   = string(T.mouse_key);
    dayIx = T.day_index;

    T = T( ismember(mkT, activeMice) & ismember(dayIx, daysAll), :);

    if isempty(T)
        timeGridMin   = [];
        meanCumPerDay = [];
        meanCumEarly  = [];
        meanCumLate   = [];
        sessInfo      = struct([]);
        return;
    end

    if ismember('Lick_TTL',T.Properties.VariableNames)
        T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
    else
        T.Lick_TTL = zeros(height(T),1);
    end
    if ismember('Injector_TTL',T.Properties.VariableNames)
        T.Injector_TTL(isnan(T.Injector_TTL)) = 0;
    else
        T.Injector_TTL = zeros(height(T),1);
    end
    tb_all = pickTimebase_fast(T);

    [g, mk, di, si] = findgroups(string(T.mouse_key), T.day_index, T.session_idx);
    nG = max(g);

    sessInfo = struct('mouse_key',[], 'day_index',[], 'session_idx',[], ...
                      'lickTimesMin',[], 'rewardTimesMin',[], 'durMin',[]);
    sessInfo(nG).mouse_key = [];

    maxDurMin = 0;
    for k = 1:nG
        idx = (g==k);
        t    = tb_all(idx);
        t    = double(t(:));
        good = isfinite(t);
        if ~any(good), continue; end
        t = t(good);
        t = t - min(t);
        durMin = (max(t) - min(t))/60;
        maxDurMin = max(maxDurMin, durMin);

        L  = logical(T.Lick_TTL(idx));     L = L(good);
        R  = logical(T.Injector_TTL(idx)); R = R(good);

        lickTimesMin   = t(L)/60;
        rewardTimesMin = t(R)/60;

        sessInfo(k).mouse_key      = mk(k);
        sessInfo(k).day_index      = di(k);
        sessInfo(k).session_idx    = si(k);
        sessInfo(k).lickTimesMin   = lickTimesMin(:)';
        sessInfo(k).rewardTimesMin = rewardTimesMin(:)';
        sessInfo(k).durMin         = durMin;
    end

    if maxDurMin <= 0
        timeGridMin   = [];
        meanCumPerDay = [];
        meanCumEarly  = [];
        meanCumLate   = [];
        return;
    end

    maxDurMin   = ceil(maxDurMin);
    timeGridMin = 0:1:maxDurMin;

    nDays         = max(daysAll);
    meanCumPerDay = nan(nDays, numel(timeGridMin));
    sumCum        = zeros(nDays, numel(timeGridMin));
    nSess         = zeros(nDays,1);

    for k = 1:numel(sessInfo)
        d = sessInfo(k).day_index;
        if d<1 || d>nDays || isempty(sessInfo(k).lickTimesMin), continue; end
        cvec = arrayfun(@(tt) sum(sessInfo(k).lickTimesMin <= tt), timeGridMin);
        sumCum(d,:) = sumCum(d,:) + cvec;
        nSess(d)    = nSess(d) + 1;
    end

    for d = 1:nDays
        if nSess(d) > 0
            meanCumPerDay(d,:) = sumCum(d,:) / nSess(d);
        end
    end

    earlyDays = 1:5;
    lateDays  = 6:10;
    meanCumEarly = mean(meanCumPerDay(earlyDays, :), 1, 'omitnan');
    meanCumLate  = mean(meanCumPerDay(lateDays,  :), 1, 'omitnan');
end

function raster_day_plot(sessInfo, dayIdx, whichType, outDir)
    mask = [sessInfo.day_index] == dayIdx;
    sess = sessInfo(mask);
    if isempty(sess), return; end

    fh = figure('Color','w','Position',[80 80 900 500]); hold on;

    labels = cell(1, numel(sess));
    for i = 1:numel(sess)
        switch lower(whichType)
            case 'reward', t = sess(i).rewardTimesMin;
            otherwise,     t = sess(i).lickTimesMin;
        end
        if isempty(t), continue; end
        y = i*ones(size(t));
        scatter(t, y, 12, 'k','filled');
        labels{i} = char(sess(i).mouse_key);
    end

    xlabel('Session time (min)');
    ylabel('Mouse (session)');
    title(sprintf('%s raster – Active mice – Day %d', ...
        upper(whichType(1))+whichType(2:end), dayIdx));
    grid on; box on;

    set(gca,'YTick',1:numel(sess), 'YTickLabel',labels);
    ylim([0.5, numel(sess)+0.5]);

    fn = fullfile(outDir, sprintf('raster_%s_day%d.png', whichType, dayIdx));
    printpng(fh, fn); close(fh);
end

function plot_cumLicks_allDays(timeGridMin, meanCumPerDay, daysAll, outDir)
    if isempty(timeGridMin) || isempty(meanCumPerDay), return; end

    fh = figure('Color','w','Position',[80 80 900 550]); hold on;
    cols = lines(numel(daysAll));

    for i = 1:numel(daysAll)
        d = daysAll(i);
        if d > size(meanCumPerDay,1), continue; end
        y = meanCumPerDay(d,:);
        if all(~isfinite(y)), continue; end

        plot(timeGridMin, y, 'Color', cols(i,:), 'LineWidth',2, ...
            'DisplayName', sprintf('Day %d',d));

        lastIdx = find(isfinite(y),1,'last');
        if ~isempty(lastIdx)
            xLast = timeGridMin(lastIdx);
            yLast = y(lastIdx);
            text(xLast+0.3, yLast, sprintf('Day %d',d), ...
                'FontSize',8, 'HorizontalAlignment','left', ...
                'VerticalAlignment','middle');
        end
    end

    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title('Cumulative licking across session (Active mice)');
    grid on; box on;
    xlim([min(timeGridMin), max(timeGridMin)+1]);

    fn = fullfile(outDir,'cumLicks_allDays_active.png');
    printpng(fh, fn); close(fh);
end

function plot_cumLicks_early_late(timeGridMin, meanCumEarly, meanCumLate, outDir)
    if isempty(timeGridMin), return; end
    fh = figure('Color','w','Position',[80 80 900 550]); hold on;
    plot(timeGridMin, meanCumEarly, 'b-','LineWidth',2.5,'DisplayName','Days 1–5');
    plot(timeGridMin, meanCumLate,  'r-','LineWidth',2.5,'DisplayName','Days 6–10');
    xlabel('Session time (min)');
    ylabel('Mean cumulative licks');
    title('Early (days 1–5) vs Late (days 6–10) sessions (Active mice)');
    grid on; box on;
    legend('Location','southeast');
    fn = fullfile(outDir,'cumLicks_early_vs_late_active.png');
    printpng(fh, fn); close(fh);
end

%% ---------- fast per-session → per-day metrics (basic) ----------

function [S, D] = fast_session_day_metrics_basic(T, runDir)
    cacheMat = fullfile(runDir, 'S_D_cache_basic.mat');
    if exist(cacheMat,'file')
        L = load(cacheMat,'S','D','cache_hash');
        if isfield(L,'cache_hash') && isequal(L.cache_hash, local_hash(T))
            fprintf('Loaded S/D from cache: %s\n', cacheMat);
            S = L.S; D = L.D; return;
        end
    end

    tic;
    need = {'mouse_key','day_index','day_name','session_idx','Diameter_px', ...
            'Lick_TTL','Injector_TTL','CamTime_rel_s','PupilTimestamp_s', ...
            'CamTime_s','PlotTime_s_30fps','Session_Paradigm','isPassive', ...
            'RequirementLast','Immersion_Latency_s'};
    V = T.Properties.VariableNames;
keepExtras = V( contains(V,'TST_','IgnoreCase',true) | ...
                contains(V,'HOT_','IgnoreCase',true) | ...
                strcmpi(V,'Immersion_Latency_s') );

    need = unique([need, keepExtras], 'stable');
    need = intersect(need, V, 'stable');
    T = T(:, need);

    if ~isstring(T.mouse_key), T.mouse_key = string(T.mouse_key); end
    if ~isstring(T.day_name),  T.day_name  = string(T.day_name);  end
    if ismember('Session_Paradigm',T.Properties.VariableNames) && ~isstring(T.Session_Paradigm)
        T.Session_Paradigm = string(T.Session_Paradigm);
    end
    T.mouse_key  = categorical(T.mouse_key);
    T.day_name   = categorical(T.day_name);
    if ismember('Session_Paradigm',T.Properties.VariableNames)
        T.Session_Paradigm = categorical(T.Session_Paradigm);
    end
    if ismember('isPassive',T.Properties.VariableNames) && ~isnumeric(T.isPassive)
        T.isPassive = double(T.isPassive);
    end

    if ismember('Lick_TTL',T.Properties.VariableNames)
        T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;
    end
    if ismember('Injector_TTL',T.Properties.VariableNames)
        T.Injector_TTL(isnan(T.Injector_TTL)) = 0; T.Injector_TTL = T.Injector_TTL > 0.5;
    end

    [g, keys_mouse, keys_day, keys_dayname, keys_sess] = findgroups( ...
        T.mouse_key, T.day_index, T.day_name, T.session_idx);
    nG = max(g);
    fprintf('Computing per-session metrics for %d sessions...\n', nG);

    S = table();
    S.mouse_key        = removecats(keys_mouse);
    S.day_index        = double(keys_day);
    S.day_name         = removecats(keys_dayname);
    S.session_idx      = double(keys_sess);
    S.RequirementLast  = nan(nG,1);
    S.isPassive        = nan(nG,1);
    S.SessionMinutes   = nan(nG,1);
    S.Session_Paradigm = strings(nG,1);
    vars = {'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s', ...
            'lick_medianIEI_s','lick_freq_per_10s', ...
            'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s', ...
            'rew_medianIRI_s','pupil_mean'};
    for v = vars, S.(v{1}) = nan(nG,1); end
    for j=1:numel(keepExtras), S.(keepExtras{j}) = nan(nG,1); end

    tb_all = pickTimebase_fast(T);
    lastTick = tic; step = max(1, floor(nG/40));
    for k = 1:nG
        idx = (g == k);
        tb  = tb_all(idx);
        dur = finiteRange_fast(tb);
        S.SessionMinutes(k) = dur/60;

        if ismember('RequirementLast',T.Properties.VariableNames)
            S.RequirementLast(k) = mean(double(T.RequirementLast(idx)),'omitnan');
        end
        if ismember('isPassive',T.Properties.VariableNames)
            ip = double(T.isPassive(idx)); ip = ip(isfinite(ip));
            if ~isempty(ip), S.isPassive(k) = mode(round(ip)); end
        end
        if ismember('Session_Paradigm',T.Properties.VariableNames)
            sp = T.Session_Paradigm(idx);
            if iscategorical(sp), sp = string(mode(sp(~ismissing(sp))));
            else
                sp = string(sp); sp = sp(~ismissing(sp));
                if isempty(sp), sp = ""; else, sp = sp(1); end
            end
            S.Session_Paradigm(k) = sp;
        end

        if ismember('Diameter_px',T.Properties.VariableNames)
            S.pupil_mean(k) = mean(double(T.Diameter_px(idx)),'omitnan');
        end
        if ismember('Lick_TTL',T.Properties.VariableNames)
            [n,md,td,iei] = eventMetrics_fast(tb, logical(T.Lick_TTL(idx)));
            S.lick_n(k)=n; S.lick_meanDur_s(k)=md; S.lick_totalDur_s(k)=td; S.lick_medianIEI_s(k)=iei;
            if S.SessionMinutes(k)>0
                S.lick_freq_per_min(k)  = n / S.SessionMinutes(k);
                S.lick_freq_per_10s(k)  = n / (S.SessionMinutes(k)*6); % per 10 s
            end
        end
        if ismember('Injector_TTL',T.Properties.VariableNames)
            [n,md,td,iri] = eventMetrics_fast(tb, logical(T.Injector_TTL(idx)));
            S.rew_n(k)=n; S.rew_meanDur_s(k)=md; S.rew_totalDur_s(k)=td; S.rew_medianIRI_s(k)=iri;
            if S.SessionMinutes(k)>0
                S.rew_freq_per_min(k)=n/S.SessionMinutes(k);
            end
        end
        for j=1:numel(keepExtras)
            col = keepExtras{j};
            if ismember(col,T.Properties.VariableNames)
                S.(col)(k) = mean(double(T.(col)(idx)),'omitnan');
            end
        end
        if mod(k,step)==0 && toc(lastTick)>0.5
            fprintf('  %d/%d (%.0f%%)\n', k, nG, 100*k/nG);
            lastTick = tic;
        end
    end
    fprintf('Per-session metrics done in %.1f s.\n', toc);

    fprintf('Collapsing to per-day medians...\n');
    [g2, mk2, di2, dn2] = findgroups(S.mouse_key, S.day_index, S.day_name);
    D = table(removecats(mk2), double(di2), removecats(dn2), ...
              'VariableNames',{'mouse_key','day_index','day_name'});
    baseList  = [{'RequirementLast','Immersion_Latency_s','isPassive','SessionMinutes'}, vars];

    extraList = intersect(keepExtras, S.Properties.VariableNames, 'stable');
    list = unique([baseList, extraList], 'stable');
    for v = list
        D.(v{1}) = splitapply(@(x) median(x,'omitnan'), S.(v{1}), g2);
    end

    cache_hash = local_hash(T(:,{'mouse_key','day_index','session_idx'}));
    try
        save(cacheMat,'S','D','cache_hash','-v7.3');
        fprintf('Cached S/D to: %s\n', cacheMat);
    catch
    end
end

% --- tiny fast helpers reused above ---
function tb = pickTimebase_fast(T)
    cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
    tb = nan(height(T),1);
    for i = 1:numel(cands)
        if ismember(cands{i},T.Properties.VariableNames)
            v = double(T.(cands{i}));
            if any(isfinite(v)), tb = v; return; end
        end
    end
end

function r = finiteRange_fast(x)
    x = double(x(:)); x = x(isfinite(x));
    if isempty(x), r = 0; else, r = max(x)-min(x); end
end

function [n, meanDur, totalDur, medianIEI] = eventMetrics_fast(t, ttl)
    t   = double(t(:)); ttl = logical(ttl(:));
    good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);
    if numel(t)<2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
    dt  = diff(t); md = median(dt(isfinite(dt))); if ~isfinite(md), md=1/30; end
    t(2:end) = max(t(2:end), t(1:end-1)+md*0.5);
    d   = diff([false; ttl; false]); on = find(d==1); off = find(d==-1)-1;
    n   = numel(on);
    if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
    edges  = [t; t(end)+md];
    segDur = sum(edges(off+1) - edges(on), 2, 'omitnan');
    meanDur  = mean(segDur,'omitnan');
    totalDur = sum(segDur,'omitnan');
    if n>=2, medianIEI = median(diff(t(on)),'omitnan'); else, medianIEI = NaN; end
end

function h = local_hash(Tkeys)
    try
        raw = [uint32(double(grp2idx(categorical(Tkeys.mouse_key)))), ...
               uint32(double(Tkeys.day_index)), uint32(double(Tkeys.session_idx))];
        h = uint64(sum(uint64(raw(:)).*1664525 + 1013904223));
    catch
        h = now;
    end
end
function unpaired_day_compare_nonparam(D, mName, dayX, outDir)
% 비모수 unpaired 비교: 특정 day에서 Active vs Passive
% ranksum (Wilcoxon rank-sum) 사용

    if ~ismember(mName, D.Properties.VariableNames)
        fprintf('Nonparam: metric %s not found, skipping.\n', mName);
        return;
    end
    if ~ismember('Group', D.Properties.VariableNames)
        fprintf('Nonparam: no Group column for %s, skipping.\n', mName);
        return;
    end

    % Day + Group 별 값 추출
    idxA = D.day_index==dayX & D.Group=="Active";
    idxP = D.day_index==dayX & D.Group=="Passive";

    xa = D.(mName)(idxA); xa = xa(isfinite(xa));
    xp = D.(mName)(idxP); xp = xp(isfinite(xp));

    if numel(xa) < 2 || numel(xp) < 2
        fprintf('Nonparam %s Day %d: not enough data (NA=%d, NP=%d), skipping.\n', ...
            mName, dayX, numel(xa), numel(xp));
        return;
    end

    % Wilcoxon rank-sum test (Mann–Whitney)
    p = ranksum(xa, xp);

    % 박스플롯 + 점 찍기
    fh = figure('Color','w','Position',[120 120 500 420]); hold on;
    grp = [ones(size(xa)); 2*ones(size(xp))];
    val = [xa; xp];

    boxchart(grp, val, 'BoxWidth',0.3, 'BoxFaceColor',[0.8 0.9 1]);

    % jittered individual points
    jitter = 0.08;
    plot(1 + (rand(size(xa))-0.5)*2*jitter, xa, 'ko', ...
        'MarkerFaceColor',[0.2 0.2 0.2],'MarkerSize',4);
    plot(2 + (rand(size(xp))-0.5)*2*jitter, xp, 'ko', ...
        'MarkerFaceColor',[0.2 0.2 0.2],'MarkerSize',4);

    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2], 'XTickLabel',{'Active','Passive'});
    ylabel(strrep(mName,'_',' '), 'Interpreter','none');
    title(sprintf('%s: Day %d Active vs Passive (N_A=%d, N_P=%d, p=%.3g, ranksum)', ...
          mName, dayX, numel(xa), numel(xp), p), ...
          'Interpreter','none');
    grid on; box on;

    fn = fullfile(outDir, sprintf('unpaired_%s_day%d_active_vs_passive.png', ...
                safeName(mName), dayX));
    printpng(fh, fn); close(fh);
end

%% ---------- misc small helpers ----------

function s = safeName(nm)
    s = regexprep(nm,'[^a-zA-Z0-9]+','_');
end

function printpng(fh, fn)
    set(fh,'PaperPositionMode','auto');
    try
        exportgraphics(fh, fn, 'Resolution',180);
    catch
        print(fh, fn, '-dpng','-r180');
    end
end
