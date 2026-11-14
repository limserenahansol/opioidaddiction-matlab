function run_pr_trial_licking_trajectory
% run_pr_trial_licking_trajectory — within- vs between-trial lick rate by PR level
%
% Set USE_DEMO = false and CUSTOM_CSV to your trial-level table when ready.
% Expected columns: see README.md in this folder.
%
% Pipeline position: downstream of Step 04 lick / PR exports, parallel to Step 07 EFA.

%% ========== user config ==========
USE_DEMO  = false;  % use run_pr_licking_full_analysis.m for all-data + phase faceting
CUSTOM_CSV = ''; % e.g. 'C:\data\my_pr_trial_licking.csv'

here = fileparts(mfilename('fullpath'));
if USE_DEMO || isempty(CUSTOM_CSV)
    csvPath = fullfile(here, 'demo_data', 'demo_pr_trial_table.csv');
else
    csvPath = CUSTOM_CSV;
end
outDir = fullfile(here, 'output_demo');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ========== load ==========
T = readtable(csvPath);
req = {'mouse_key', 'day_index', 'session_idx', 'Trial', 'WithinTrialLickRate'};
assert(all(ismember(req, T.Properties.VariableNames)), ...
    'Table must contain: %s', strjoin(req, ', '));

T.Trial = double(T.Trial);
T.WithinTrialLickRate = double(T.WithinTrialLickRate);
if ismember('BetweenTrialLickRate', T.Properties.VariableNames)
    T.BetweenTrialLickRate = double(T.BetweenTrialLickRate);
    hasBetween = true;
else
    hasBetween = false;
    warning('BetweenTrialLickRate not found: only within-trial plots will be made.');
end
if ~ismember('Group', T.Properties.VariableNames)
    T.Group = repmat({'All'}, height(T), 1);
else
    T.Group = cellstr(T.Group);
end

%% ========== per mouse × trial: mean over day/session replicates ==========
G = findgroups(T.mouse_key, T.Trial);
summaryMouse = table;
summaryMouse.mouse_key = splitapply(@(x) x(1), T.mouse_key, G);
summaryMouse.Trial = splitapply(@(x) x(1), T.Trial, G);
summaryMouse.n = splitapply(@numel, T.WithinTrialLickRate, G);
summaryMouse.MeanWithinTrialLickRate = splitapply(@mean, T.WithinTrialLickRate, G);
if hasBetween
    summaryMouse.MeanBetweenTrialLickRate = splitapply(@(x) mean(x, 'omitnan'), T.BetweenTrialLickRate, G);
else
    summaryMouse.MeanBetweenTrialLickRate = nan(height(summaryMouse), 1);
end
summaryMouse.Group = splitapply(@(x) x(1), T.Group, G);
summaryMouse.mouse_key = cellstr(summaryMouse.mouse_key);
summaryMouse.Group = cellstr(summaryMouse.Group);

%% ========== plot: within-trial (group mean ± SEM + faint mouse lines) ==========
f1 = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on
mice = cellstr(unique(summaryMouse.mouse_key));
ug = unique(summaryMouse.Group, 'stable');
for im = 1:numel(mice)
    r = strcmp(summaryMouse.mouse_key, mice{im});
    plot(summaryMouse.Trial(r), summaryMouse.MeanWithinTrialLickRate(r), '-', ...
        'Color', [0.82, 0.82, 0.82], 'LineWidth', 0.6, 'HandleVisibility', 'off');
end
for ig = 1:numel(ug)
    g = ug{ig};
    sub = summaryMouse(strcmp(summaryMouse.Group, g), :);
    tr = unique(sub.Trial);
    m = nan(size(tr));
    se = nan(size(tr));
    for it = 1:numel(tr)
        w = sub.MeanWithinTrialLickRate(sub.Trial == tr(it));
        m(it) = mean(w, 'omitnan');
        se(it) = std(w, 0, 'omitnan') / max(1, sqrt(sum(~isnan(w))));
    end
    [tr, ord] = sort(tr);
    m = m(ord);
    se = se(ord);
    errorbar(tr, m, se, 'LineWidth', 2, 'CapSize', 0, 'DisplayName', char(g));
end
xlabel('PR trial / breakpoint level')
ylabel('Within-trial licking rate (1/s)')
title('Within-trial lick rate across PR levels')
grid on
box off
if numel(ug) > 1, legend('Location', 'best'); end
try
    saveas(f1, fullfile(outDir, 'WithinTrial_LickRate_by_PRLevel.png'));
    saveas(f1, fullfile(outDir, 'WithinTrial_LickRate_by_PRLevel.fig'));
catch
    saveas(f1, fullfile(outDir, 'WithinTrial_LickRate_by_PRLevel.png'));
end

%% ========== plot: between-trial ==========
if hasBetween && any(~isnan(summaryMouse.MeanBetweenTrialLickRate))
    f2 = figure('Color', 'w', 'Position', [120, 120, 900, 600]);
    hold on
    for im = 1:numel(mice)
        r = strcmp(summaryMouse.mouse_key, mice{im});
        plot(summaryMouse.Trial(r), summaryMouse.MeanBetweenTrialLickRate(r), '-', ...
            'Color', [0.82, 0.82, 0.82], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
    for ig = 1:numel(ug)
        g = ug{ig};
        sub = summaryMouse(strcmp(summaryMouse.Group, g), :);
        tr = unique(sub.Trial);
        m = nan(size(tr));
        se = nan(size(tr));
        for it = 1:numel(tr)
            w = sub.MeanBetweenTrialLickRate(sub.Trial == tr(it));
            m(it) = mean(w, 'omitnan');
            se(it) = std(w, 0, 'omitnan') / max(1, sqrt(sum(~isnan(w))));
        end
        [tr, ord] = sort(tr);
        m = m(ord);
        se = se(ord);
        okp = ~isnan(m) & tr > 0.5;
        if any(okp)
            errorbar(tr(okp), m(okp), se(okp), 'LineWidth', 2, 'CapSize', 0, 'DisplayName', char(g));
        end
    end
    xlabel('PR trial (interval ends at this trial start)')
    ylabel('Between-trial licking rate (1/s)')
    title('Between-trial lick rate across PR levels')
    grid on
    box off
    if numel(ug) > 1, legend('Location', 'best'); end
    try
        saveas(f2, fullfile(outDir, 'BetweenTrial_LickRate_by_PRLevel.png'));
        saveas(f2, fullfile(outDir, 'BetweenTrial_LickRate_by_PRLevel.fig'));
    catch
        saveas(f2, fullfile(outDir, 'BetweenTrial_LickRate_by_PRLevel.png'));
    end
end

%% ========== group-level summary (all) ==========
trials = unique(summaryMouse.Trial);
nt = numel(trials);
tw = nan(nt, 1);
semw = nan(nt, 1);
for it = 1:nt
    w = summaryMouse.MeanWithinTrialLickRate(summaryMouse.Trial == trials(it));
    tw(it) = mean(w, 'omitnan');
    semw(it) = std(w, 0, 'omitnan') / max(1, sqrt(sum(~isnan(w))));
end
summaryAll = table(trials, tw, semw, 'VariableNames', {'Trial', 'GroupMeanWithin', 'GroupSEMWithin'});

%% ========== optional PR pair (edit) ==========
PRa = 5;
PRb = 10;
a = nan(numel(mice), 1);
b = nan(numel(mice), 1);
for im = 1:numel(mice)
    r1 = strcmp(summaryMouse.mouse_key, mice{im}) & summaryMouse.Trial == PRa;
    r2 = strcmp(summaryMouse.mouse_key, mice{im}) & summaryMouse.Trial == PRb;
    if any(r1), a(im) = summaryMouse.MeanWithinTrialLickRate(find(r1, 1)); end
    if any(r2), b(im) = summaryMouse.MeanWithinTrialLickRate(find(r2, 1)); end
end
v = ~isnan(a) & ~isnan(b);
if sum(v) >= 2
    [~, p, ~, st] = ttest(a(v), b(v));
    fprintf('Paired t-test: within-trial rate PR%d vs PR%d: t(%d)=%.3f, p=%.4f (n=%d)\n', ...
        PRa, PRb, st.df, st.tstat, p, sum(v));
else
    fprintf('PR pair test skipped: need at least 2 mice with both PR%d and PR%d (check trial range).\n', PRa, PRb);
end

%% ========== save ==========
writetable(T, fullfile(outDir, 'input_trial_table_echo.csv'));
writetable(summaryMouse, fullfile(outDir, 'MouseMean_LickRate_by_PRLevel.csv'));
writetable(summaryAll, fullfile(outDir, 'GroupMean_LickRate_by_PRLevel.csv'));
save(fullfile(outDir, 'PR_LickTrajectory_workspace.mat'), 'T', 'summaryMouse', 'summaryAll');
fprintf('Step 08 PR lick trajectory done. Output: %s\n', outDir);
end
