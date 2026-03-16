function rebuild_lick_patterns_with_ID()
% Rebuild lick-pattern plots WITH mouse IDs even if the CSVs lack them.
% - Locates latest run_* under longitudinal_outputs
% - Loads per_session_features / per_trial_rates / combined trend (if present)
% - Infers mouse IDs from available columns or file paths
% - PCA + kmeans (K=3 by default), labels dots with mouse IDs
% - Recreates "mean ± SEM across day" and "across trial" (if joinable)
% - Writes a point-to-ID map CSV
%
% Outputs: saved in <runDir>\figs\lick_patterns_MASTER\plots_matlab\

%% ---------- locate latest run_* ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up=1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry = p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*')); assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]); runDir = fullfile(D(ix).folder, D(ix).name);

dataDir = fullfile(runDir, 'figs','lick_patterns_MASTER');
if ~exist(dataDir,'dir'), warning('figs folder missing; using runDir'); dataDir = runDir; end
outDir  = fullfile(dataDir,'plots_matlab'); if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Using runDir: %s\n', runDir);

%% ---------- settings ----------
K = 3; rng(1);
LABEL_ALL = true;                 % set false to only label selected mice below
MICE_TO_LABEL = {'7597_black'};   % used if LABEL_ALL==false
featureList = {'lick_per_min','iei_cv','cv2_median','rhythm_index', ...
               'burst_fraction','bout_rate_per_min','iei_median','iei_mean'};

%% ---------- load CSVs ----------
Tses  = readCSV(dataDir, runDir, 'per_session_features.csv', true);
Ttrl  = readCSV(dataDir, runDir, 'per_trial_rates.csv', false);
Tcomb = readCSV(dataDir, runDir, 'class_combined_day_and_trial.csv', false);

%% ---------- ensure IDs exist (or infer) ----------
% preferred ID fields
idFields = {'mouse_key','mouse_id','mouse','animal','subject'};
Tses = ensureMouseID(Tses, idFields, runDir);   % adds/standardizes 'mouse_key'
assert(ismember('mouse_key', Tses.Properties.VariableNames), ...
    'Could not infer mouse_key in per_session_features.');

% also make sure 'day_index' and 'session_idx' exist (fallbacks if missing)
if ~ismember('day_index', Tses.Properties.VariableNames)
    Tses.day_index = guessOrdinal(Tses, 'day'); end
if ~ismember('session_idx', Tses.Properties.VariableNames)
    Tses.session_idx = guessOrdinal(Tses, 'session'); end

%% ---------- PCA + kmeans ----------
useFeat = featureList(ismember(featureList, Tses.Properties.VariableNames));
X = zscore(Tses{:, useFeat}, 0, 1);
[coeff,score,~,~,expl] = pca(X);
nPC = find(cumsum(expl)>=80,1,'first'); if isempty(nPC), nPC = min(2,size(score,2)); end
Z = score(:,1:nPC);
idx = kmeans(Z, K, 'Replicates', 50, 'MaxIter', 1000, 'Display','off');

P = table(Tses.mouse_key, Tses.day_index, Tses.session_idx, idx, score(:,1), score(:,2), ...
    'VariableNames', {'mouse_key','day_index','session_idx','cluster','PC1','PC2'});
writetable(P, fullfile(outDir,'pca_point_map.csv'));

%% ---------- PCA figure (IDs on dots) ----------
C = lines(K); marks = {'o','s','^','d','>','<'};
f=figure('Color','w','Position',[100 100 760 760]); hold on
for k=1:K
    r = P.cluster==k;
    scatter(P.PC1(r), P.PC2(r), 36, 'Marker', marks{min(k,numel(marks))}, ...
        'MarkerEdgeColor','k','MarkerFaceColor',C(k,:), ...
        'DisplayName', sprintf('Cluster %d',k));
end
xlabel(sprintf('PC1 (%.1f%%)',expl(1))); ylabel(sprintf('PC2 (%.1f%%)',expl(2)));
title('Sessions: PCA colored by k-means cluster'); legend('Location','best'); grid on

if LABEL_ALL
    labMask = true(height(P),1);
else
    labMask = ismember(P.mouse_key, MICE_TO_LABEL);
end
for i=find(labMask)'
    % CHANGE label text style here if you want day/session appended:
    % txt = sprintf('%s d%g s%g', P.mouse_key{i}, P.day_index(i), P.session_idx(i));
    txt = P.mouse_key{i};
    text(P.PC1(i), P.PC2(i), [' ' txt], 'FontSize',8, 'FontWeight','bold', ...
         'Color',[0 0 0],'BackgroundColor','w','Margin',0.1,'Clipping','on');
    plot(P.PC1(i), P.PC2(i), 'ko', 'MarkerSize',7, 'LineWidth',1.4);
end
saveas(f, fullfile(outDir,'pca_kmeans_clusters_LABELED.png'));

%% ---------- clusters: mean ± SEM across day ----------
Tses.Cluster = idx;
G = findgroups(Tses.day_index, Tses.Cluster);
S = table;
S.day_index = splitapply(@unique, Tses.day_index, G);
S.Cluster   = splitapply(@unique, Tses.Cluster, G);
S.n         = splitapply(@numel, Tses.lick_per_min, G);
S.mean_lick = splitapply(@mean,  Tses.lick_per_min, G);
S.sem_lick  = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), Tses.lick_per_min, G);

f2=figure('Color','w','Position',[100 100 950 520]); hold on
for k=1:K
    r = S.Cluster==k; [d,ord] = sort(S.day_index(r));
    mu=S.mean_lick(r); mu=mu(ord); se=S.sem_lick(r); se=se(ord);
    plot(d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, 'Marker', marks{min(k,numel(marks))}, ...
         'MarkerFaceColor', C(k,:), 'DisplayName', sprintf('Cluster %d',k));
    for j=1:numel(d), line([d(j) d(j)], [mu(j)-se(j) mu(j)+se(j)], 'Color',C(k,:)); end
end
xlabel('Day'); ylabel('Licks / min'); title('Clusters: mean ± SEM across day'); grid on; legend('Location','best');
saveas(f2, fullfile(outDir,'clusters_across_day_meanSEM.png'));

%% ---------- clusters: mean ± SEM across trial (if joinable) ----------
okTrial = ~isempty(Ttrl);
if okTrial
    % try to create/match IDs in trial table the same way
    Ttrl = ensureMouseID(Ttrl, idFields, runDir);
    need = {'mouse_key','day_index','session_idx'};
    if ~all(ismember(need, Ttrl.Properties.VariableNames))
        if ~ismember('day_index', Ttrl.Properties.VariableNames),  Ttrl.day_index = guessOrdinal(Ttrl,'day'); end
        if ~ismember('session_idx', Ttrl.Properties.VariableNames),Ttrl.session_idx = guessOrdinal(Ttrl,'session'); end
    end
    canJoin = all(ismember(need, Ttrl.Properties.VariableNames)) && all(ismember(need, Tses.Properties.VariableNames));
    if canJoin
        Ktrl = join(Ttrl, Tses(:,[need {'Cluster'}]), 'Keys', need);
        G2 = findgroups(Ktrl.Trial, Ktrl.Cluster);
        S2 = table;
        S2.Trial   = splitapply(@unique, Ktrl.Trial, G2);
        S2.Cluster = splitapply(@unique, Ktrl.Cluster, G2);
        S2.n       = splitapply(@numel, Ktrl.lick_rate, G2);
        S2.mean_lr = splitapply(@mean,  Ktrl.lick_rate, G2);
        S2.sem_lr  = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), Ktrl.lick_rate, G2);

        f3=figure('Color','w','Position',[100 100 1000 520]); hold on
        for k=1:K
            r = S2.Cluster==k; [t,ord] = sort(S2.Trial(r));
            mu=S2.mean_lr(r); mu=mu(ord); se=S2.sem_lr(r); se=se(ord);
            plot(t, mu, '-', 'Color', C(k,:), 'LineWidth', 1.5, 'Marker', marks{min(k,numel(marks))}, ...
                 'MarkerSize', 3, 'MarkerFaceColor', C(k,:), 'DisplayName', sprintf('Cluster %d',k));
            step = max(1, round(numel(t)/50)); ii = 1:step:numel(t);
            errorbar(t(ii), mu(ii), se(ii), 'LineStyle','none', 'Color', C(k,:));
        end
        xlabel('Trial'); ylabel('Lick rate (Hz)'); title('Clusters: mean ± SEM across trial');
        grid on; legend('Location','best');
        saveas(f3, fullfile(outDir,'clusters_across_trial_meanSEM.png'));
    else
        warning('Could not join trials to sessions (missing keys). Skipping across-trial plot.');
    end
end

%% ---------- optional: trend map if combined table present ----------
if ~isempty(Tcomb) && all(ismember({'mouse_key','rho_days','rho_trials'}, Tcomb.Properties.VariableNames))
    f4=figure('Color','w','Position',[100 100 720 720]); hold on
    scatter(Tcomb.rho_days, Tcomb.rho_trials, 60, 'o', ...
            'MarkerFaceColor',[0.3 0.5 0.9], 'MarkerEdgeColor','k');
    for i=1:height(Tcomb)
        text(Tcomb.rho_days(i)+0.01, Tcomb.rho_trials(i), [' ' Tcomb.mouse_key{i}], 'FontSize',8);
    end
    xlabel('Spearman \rho (day vs lick/min)'); ylabel('Spearman \rho (trial vs lick rate)');
    title('Mouse trend map (day-trend vs trial-trend)'); grid on
    saveas(f4, fullfile(outDir,'trend_map_day_vs_trial.png'));
end

fprintf('Done. Figures + pca_point_map.csv saved in:\n  %s\n', outDir);
end

%% ================= helpers =================
function T = readCSV(dataDir, runDir, name, must)
% Find CSV in figs folder, runDir, or recursively. Return [] if not found and must==false.
T = [];
paths = { fullfile(dataDir,name), fullfile(runDir,name) };
for i=1:numel(paths)
    if exist(paths{i},'file')==2, T = readtable(paths{i}); return; end
end
% recursive fallback
try
    F = dir(fullfile(runDir,'**',name));
    if ~isempty(F), T = readtable(fullfile(F(1).folder, F(1).name)); return; end
catch
    % no-op
end
if must, error('Required file not found: %s', name); end
end

function T = ensureMouseID(T, idFields, runDir)
% Ensure T has 'mouse_key'. If not, try to:
% 1) copy from any of idFields (case-insensitive),
% 2) parse from any filename/path columns,
% 3) parse from runDir subfolders embedded in any 'file/path/source' fields.
if isempty(T), return; end

% already present?
if ismember('mouse_key', T.Properties.VariableNames)
    T.mouse_key = toStr(T.mouse_key); return;
end

% 1) copy from alternate ID field
for f = idFields
    nm = f{1};
    hit = find(strcmpi(T.Properties.VariableNames, nm),1);
    if ~isempty(hit)
        T.mouse_key = toStr(T.(T.Properties.VariableNames{hit}));
        return;
    end
end

% 2) parse from likely filename/path columns
candCols = {'file','filename','path','source','session_file','session_path','raw_file','session_name'};
candCols = candCols(ismember(candCols, T.Properties.VariableNames));
if ~isempty(candCols)
    for c = candCols
        key = parseMouseFromStrings(toStr(T.(c{1})));
        if any(~cellfun(@isempty,key))
            T.mouse_key = key; return;
        end
    end
end

% 3) last resort: look for any string column and try parsing
strCols = T.Properties.VariableNames( varfun(@ischarcol, T, 'OutputFormat','uniform') | ...
                                      varfun(@isstringcol,T,'OutputFormat','uniform'));
for c = strCols
    key = parseMouseFromStrings(toStr(T.(c{1})));
    if any(~cellfun(@isempty,key))
        T.mouse_key = key; return;
    end
end

% If still missing, fabricate stable IDs per row (not ideal, but unblocks plotting)
warning('Could not find/parse mouse IDs; fabricating "mouse_%04d".');
T.mouse_key = arrayfun(@(i)sprintf('mouse_%04d',i), (1:height(T))', 'UniformOutput', false);
end

function out = parseMouseFromStrings(s)
% Try common patterns: "7597_black", "m7597_black", ".../7597_black/...", "mouse_7597_black"
% returns cellstr (empty if not found)
out = repmat({''}, numel(s),1);
pat = '([0-9]{3,6}_[A-Za-z]+)';  % e.g., 7597_black
for i=1:numel(s)
    txt = s{i};
    tok = regexp(txt, pat, 'tokens','once');
    if ~isempty(tok), out{i} = tok{1}; continue; end
    % alternative: "mouse_7597_black"
    tok = regexp(txt, 'mouse[_-]?([0-9]{3,6}_[A-Za-z]+)', 'tokens','once');
    if ~isempty(tok), out{i} = tok{1}; continue; end
end
end

function v = toStr(x)
if iscellstr(x) || (iscell(x) && all(cellfun(@(y)ischar(y)||isstring(y),x)))
    v = cellfun(@char, x, 'UniformOutput', false);
elseif isstring(x)
    v = cellstr(x);
elseif iscategorical(x)
    v = cellstr(x);
else
    v = cellstr(string(x));
end
end

function tf = ischarcol(x), tf = iscellstr(x); end
function tf = isstringcol(x), tf = isstring(x); end

function idx = guessOrdinal(T, base)
% Fabricate sequential indices if day/session not present.
if height(T)==0, idx = zeros(0,1); return; end
% group by mouse if available
if ismember('mouse_key', T.Properties.VariableNames)
    G = findgroups(T.mouse_key);
    idx = zeros(height(T),1);
    for g=1:max(G)
        r = G==g; idx(r) = 1:nnz(r);
    end
else
    idx = (1:height(T))';
end
warning('Guessing %s indices (not found in table).', base);
end
