function run_pr_trajectory_phenotype_cluster
% run_pr_trajectory_phenotype_cluster — k-means on z-scored PR-lick trajectories
% Simplified analogue of group-based trajectory modeling (e.g. GBTM / PROC TRAJ).
% Use k=2 or 3 when N is small (e.g. <~20 mice).

K = 2; % set 2–3 for small N
RATE_KIND = 'within'; % 'within' | 'between' | 'both'
rng(1);

here = fileparts(mfilename('fullpath'));
% Prefer most recent run output, else demo
outDemo = fullfile(here, 'output_demo', 'MouseMean_LickRate_by_PRLevel.csv');
if exist(outDemo, 'file') == 2
    S = readtable(outDemo);
    S.mouse_key = cellstr(S.mouse_key);
    fprintf('Using: %s\n', outDemo);
else
    csvPath = fullfile(here, 'demo_data', 'demo_pr_trial_table.csv');
    assert(exist(csvPath, 'file') == 2, 'Missing %s', csvPath);
    T0 = readtable(csvPath);
    T0.mouse_key = cellstr(T0.mouse_key);
    G = findgroups(T0.mouse_key, T0.Trial);
    S = table;
    S.mouse_key = splitapply(@(x) x(1), T0.mouse_key, G);
    S.Trial = splitapply(@(x) x(1), T0.Trial, G);
    S.MeanWithinTrialLickRate = splitapply(@mean, T0.WithinTrialLickRate, G);
    if ismember('BetweenTrialLickRate', T0.Properties.VariableNames)
        S.MeanBetweenTrialLickRate = splitapply(@(x) mean(x, 'omitnan'), T0.BetweenTrialLickRate, G);
    else
        S.MeanBetweenTrialLickRate = nan(height(S), 1);
    end
    fprintf('Built table from demo (run run_pr_trial_licking_trajectory to refresh output_demo).\n');
end

mice = cellstr(unique(S.mouse_key));
tr = unique(S.Trial);
nM = numel(mice);
nT = numel(tr);
Xw = nan(nM, nT);
Xb = nan(nM, nT);
for i = 1:nM
    for j = 1:nT
        r = strcmp(S.mouse_key, mice{i}) & S.Trial == tr(j);
        if any(r) && ismember('MeanWithinTrialLickRate', S.Properties.VariableNames)
            Xw(i, j) = S.MeanWithinTrialLickRate(find(r, 1, 'first'));
        end
        if any(r) && ismember('MeanBetweenTrialLickRate', S.Properties.VariableNames)
            Xb(i, j) = S.MeanBetweenTrialLickRate(find(r, 1, 'first'));
        end
    end
end

if strcmpi(RATE_KIND, 'within')
    X = Xw;
elseif strcmpi(RATE_KIND, 'between')
    X = Xb;
    X = X(:, 2:end); % drop col 1 if all NaN
else
    X = [Xw, Xb(:, 2:end)]; % concatenate
end
% impute column means for kmeans
mu = mean(X, 1, 'omitnan');
for c = 1:size(X, 2)
    ix = isnan(X(:, c));
    if any(ix)
        X(ix, c) = mu(c);
    end
end
Z = (X - mean(X, 2, 'omitnan')) ./ (std(X, 0, 2, 'omitnan') + eps);

try
    idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 1000, 'Display', 'off', 'Start', 'plus');
catch
    % older MATLAB
    idx = kmeans(Z, K, 'Replicates', 30, 'MaxIter', 1000, 'Display', 'off');
end
Out = table(mice(:), idx, 'VariableNames', {'mouse_key', 'trajectory_cluster'});
outDir = fullfile(here, 'output_demo');
if ~exist(outDir, 'dir'), mkdir(outDir); end
writetable(Out, fullfile(outDir, 'pr_trajectory_cluster_assignments.csv'));
fprintf('K-means k=%d on %s rates. Assignments: %s\n', K, RATE_KIND, fullfile(outDir, 'pr_trajectory_cluster_assignments.csv'));
disp(Out);
end
