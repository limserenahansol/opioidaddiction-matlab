function analyze_lick_patterns_MASTER()
% PCA + kmeans + clearly labeled mouse IDs on the PCA scatter.
% Also saves a pca_point_map.csv (PC1/PC2 + IDs) and the day/trial summaries.

%% ---------- Locate latest run ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up = 1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry = p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*')); assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]); runDir = fullfile(D(ix).folder, D(ix).name);

dataDir = fullfile(runDir, 'figs','lick_patterns_MASTER');
if ~exist(dataDir,'dir'), warning('figs folder not found, using runDir'); dataDir = runDir; end
outDir  = fullfile(dataDir,'plots_matlab'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---------- Settings ----------
K = 3;
rng(1);

% Which mice to label on the PCA plot (put {} to label ALL; or {'7597_black'} to label one)
MICE_TO_LABEL = {'7597_black'};   % <-- change this list as you like

featureList = {'lick_per_min','iei_cv','cv2_median','rhythm_index', ...
               'burst_fraction','bout_rate_per_min','iei_median','iei_mean'};

%% ---------- Load session features ----------
f_ses = fullfile(dataDir,'per_session_features.csv');
assert(exist(f_ses,'file')==2, 'Missing per_session_features.csv at %s', f_ses);
Tses  = readtable(f_ses);

need = {'mouse_key','day_index','session_idx'};
assert(all(ismember(need, Tses.Properties.VariableNames)), 'per_session_features.csv must have mouse_key/day_index/session_idx');

useFeat = featureList(ismember(featureList, Tses.Properties.VariableNames));
X = zscore(Tses{:, useFeat}, 0, 1);

% PCA
[coeff, score, ~, ~, expl] = pca(X);
cumVar = cumsum(expl); nPC = find(cumVar>=80,1); if isempty(nPC), nPC = min(2,size(score,2)); end
Z = score(:,1:nPC);

% K-means
idx = kmeans(Z, K, 'Replicates', 50, 'MaxIter', 1000, 'Display','off');

% Map table
P = table(Tses.mouse_key, Tses.day_index, Tses.session_idx, idx, score(:,1), score(:,2), ...
    'VariableNames', {'mouse_key','day_index','session_idx','cluster','PC1','PC2'});
writetable(P, fullfile(outDir,'pca_point_map.csv'));

%% ---------- PCA figure with on-plot labels ----------
C = lines(K); markers = {'o','s','^','d','>','<'};
f=figure('Color','w','Position',[100 100 760 760]); hold on
for k=1:K
    r = P.cluster==k;
    scatter(P.PC1(r), P.PC2(r), 36, 'Marker', markers{min(k,numel(markers))}, ...
        'MarkerEdgeColor','k','MarkerFaceColor',C(k,:), ...
        'DisplayName', sprintf('Cluster %d',k));
end
xlabel(sprintf('PC1 (%.1f%%)', expl(1))); ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
title('Sessions: PCA colored by k-means cluster'); legend('Location','best'); grid on

% Decide which points to label
if isempty(MICE_TO_LABEL)
    labMask = true(height(P),1);             % label ALL points
else
    labMask = ismember(P.mouse_key, MICE_TO_LABEL);
end

% Text labels (mouse_key, optionally add day/session if you want)
for i = find(labMask)'
    txt = sprintf('%s', P.mouse_key{i});     % or: sprintf('%s d%g s%g', P.mouse_key{i}, P.day_index(i), P.session_idx(i));
    text(P.PC1(i), P.PC2(i), [' ' txt], 'FontSize',8, 'FontWeight','bold', ...
        'Color',[0 0 0], 'BackgroundColor','w', 'Margin',0.1, 'Clipping','on');
    % also outline the point to pop it out
    plot(P.PC1(i), P.PC2(i), 'ko', 'MarkerSize',7, 'LineWidth',1.4);
end

saveas(f, fullfile(outDir,'pca_kmeans_clusters_LABELED.png'));

%% ---------- Clusters: mean ± SEM across day ----------
Tses.Cluster = P.cluster;
G = findgroups(Tses.day_index, Tses.Cluster);
S = table;
S.day_index = splitapply(@unique, Tses.day_index, G);
S.Cluster   = splitapply(@unique, Tses.Cluster, G);
S.n         = splitapply(@numel, Tses.lick_per_min, G);
S.mean_lick = splitapply(@mean,  Tses.lick_per_min, G);
S.sem_lick  = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), Tses.lick_per_min, G);

f2=figure('Color','w','Position',[100 100 950 520]); hold on
for k=1:K
    r = S.Cluster==k; [d, ord] = sort(S.day_index(r));
    mu=S.mean_lick(r); mu=mu(ord); se=S.sem_lick(r); se=se(ord);
    plot(d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, 'Marker', markers{min(k,numel(markers))}, ...
        'MarkerFaceColor', C(k,:), 'DisplayName', sprintf('Cluster %d',k));
    for j=1:numel(d), line([d(j) d(j)], [mu(j)-se(j) mu(j)+se(j)], 'Color', C(k,:)); end
end
xlabel('Day'); ylabel('Licks / min'); title('Clusters: mean ± SEM across day'); grid on; legend('Location','best');
saveas(f2, fullfile(outDir,'clusters_across_day_meanSEM.png'));

%% ---------- Clusters: mean ± SEM across trial (optional if file exists) ----------
f_trial = fullfile(dataDir,'per_trial_rates.csv');
if exist(f_trial,'file')==2
    Ttrial = readtable(f_trial);
    needT = {'mouse_key','day_index','session_idx'};
    if all(ismember(needT, Ttrial.Properties.VariableNames))
        Ktrial = join(Ttrial, Tses(:,[needT {'Cluster'}]), 'Keys', needT);  % older MATLAB-safe
        G2 = findgroups(Ktrial.Trial, Ktrial.Cluster);
        S2 = table;
        S2.Trial   = splitapply(@unique, Ktrial.Trial, G2);
        S2.Cluster = splitapply(@unique, Ktrial.Cluster, G2);
        S2.n       = splitapply(@numel, Ktrial.lick_rate, G2);
        S2.mean_lr = splitapply(@mean,  Ktrial.lick_rate, G2);
        S2.sem_lr  = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), Ktrial.lick_rate, G2);

        f3=figure('Color','w','Position',[100 100 1000 520]); hold on
        for k=1:K
            r = S2.Cluster==k; [t, ord] = sort(S2.Trial(r));
            mu=S2.mean_lr(r); mu=mu(ord); se=S2.sem_lr(r); se=se(ord);
            plot(t, mu, '-', 'Color', C(k,:), 'LineWidth', 1.5, 'Marker', markers{min(k,numel(markers))}, ...
                'MarkerSize', 3, 'MarkerFaceColor', C(k,:), 'DisplayName', sprintf('Cluster %d',k));
            step = max(1, round(numel(t)/50)); ii = 1:step:numel(t);
            errorbar(t(ii), mu(ii), se(ii), 'LineStyle','none', 'Color', C(k,:));
        end
        xlabel('Trial'); ylabel('Lick rate (Hz)'); title('Clusters: mean ± SEM across trial');
        grid on; legend('Location','best');
        saveas(f3, fullfile(outDir,'clusters_across_trial_meanSEM.png'));
    else
        warning('per_trial_rates.csv lacks mouse/day/session keys, skipping across-trial plot.');
    end
else
    warning('per_trial_rates.csv not found; skipping across-trial plot.');
end

fprintf('Done. Look for IDs on: %s\n', fullfile(outDir,'pca_kmeans_clusters_LABELED.png'));
fprintf('Point-to-ID map CSV: %s\n', fullfile(outDir,'pca_point_map.csv'));
end
