function arranged_plots_all_ids()
% Clean layout:
% - Normalized session-progress curves (no legend)
% - IDs per cluster in separate right-hand tiles, multi-column
% - Across-day curves (no legend) + IDs tile
% - Autosave with timestamps

%% ----- SETUP -----
NBINS = 100; rng(1);
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here=pwd; cand=here;
    for up=1:5, p=fullfile(cand,'longitudinal_outputs'); if exist(p,'dir'), rootTry=p; break; end, cand=fileparts(cand); end
end
D = dir(fullfile(rootTry,'run_*')); assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]); runDir = fullfile(D(ix).folder, D(ix).name);
dataDir = fullfile(runDir,'figs','lick_patterns_MASTER'); if ~exist(dataDir,'dir'), dataDir = runDir; end
ts = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outDir = fullfile(dataDir,'plots_matlab'); if ~exist(outDir,'dir'), mkdir(outDir); end
keys = {'mouse_key','day_index','session_idx'};

Tses   = readtable(fullfile(dataDir,'per_session_features.csv'));
Ttrial = readtable(fullfile(dataDir,'per_trial_rates.csv'));
assert(all(ismember(keys, Tses.Properties.VariableNames)),'per_session_features missing keys');
assert(all(ismember(keys, Ttrial.Properties.VariableNames)),'per_trial_rates missing keys');

% Cluster labels if missing
if ~ismember('Cluster', Tses.Properties.VariableNames)
    feat = {'lick_per_min','iei_cv','cv2_median','rhythm_index','burst_fraction','bout_rate_per_min','iei_median','iei_mean'};
    use  = feat(ismember(feat, Tses.Properties.VariableNames));
    X    = zscore(Tses{:,use},0,1);
    Tses.Cluster = kmeans(X,3,'Replicates',50,'MaxIter',1000,'Display','off');
end
K = max(Tses.Cluster); C = lines(K); mk = {'o','s','^','d','>','<'};

%% ==== Mouse-level cluster assignment (majority) ====
K = max(Tses.Cluster);
% count sessions per mouse per cluster
Gm = findgroups(Tses.mouse_key);
mouse_ids = splitapply(@unique, Tses.mouse_key, Gm);
counts = zeros(numel(mouse_ids), K);
for k = 1:K
    counts(:,k) = splitapply(@(x)sum(x==k), Tses.Cluster, Gm);
end
[bestCount, bestK] = max(counts, [], 2);
totalCount = sum(counts,2);
purity = bestCount ./ max(1,totalCount);

MouseAssign = table(mouse_ids, totalCount, bestK, purity, ...
    'VariableNames', {'mouse_key','n_sessions','majority_cluster','purity'});

% save a readable report
writetable(MouseAssign, fullfile(outDir, ['mouse_majority_clusters_' ts '.csv']));

% print tidy ID lists: "pure" (purity==1) and "mixed" (purity<1)
for k = 1:K
    pureIDs  = MouseAssign.mouse_key(MouseAssign.majority_cluster==k & MouseAssign.purity==1);
    mixedIDs = MouseAssign(MouseAssign.majority_cluster==k & MouseAssign.purity<1, :);
    fprintf('\nCluster %d — pure mice (%d): %s\n', k, numel(pureIDs), strjoin(sort(pureIDs'), ', '));
    if ~isempty(mixedIDs)
        fprintf('Cluster %d — mixed mice (ID [purity, counts per cluster]):\n', k);
        for r = 1:height(mixedIDs)
            i = find(strcmp(mouse_ids, mixedIDs.mouse_key{r}));
            fprintf('  %s  [%.2f, %s]\n', mixedIDs.mouse_key{r}, mixedIDs.purity(r), mat2str(counts(i,:)));
        end
    end
end

%% ==== (Optional) Replot across-day using mouse-majority clusters ====
% Map each session to its mouse's majority cluster
Tmap = MouseAssign(:,{'mouse_key','majority_cluster'});
Tmap.Properties.VariableNames{2} = 'ClusterMajority';
TsesMaj = join(Tses, Tmap, 'Keys','mouse_key');

Gd = findgroups(TsesMaj.day_index, TsesMaj.ClusterMajority);
Dsum = table;
Dsum.day_index = splitapply(@unique, TsesMaj.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, TsesMaj.ClusterMajority, Gd);
Dsum.mu        = splitapply(@mean,  TsesMaj.lick_per_min, Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), TsesMaj.lick_per_min, Gd);

% plot (no legend)
C = lines(K); mk = {'o','s','^','d','>','<'};
f=figure('Color','w','Position',[60 60 900 520]); hold on
for k=1:K
    r = Dsum.Cluster==k; [d,ord] = sort(Dsum.day_index(r));
    mu=Dsum.mu(r); mu=mu(ord); se=Dsum.sem(r); se=se(ord);
    plot(d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7); set(eb,'HandleVisibility','off');
end
xlabel('Day'); ylabel('Licks / min'); title('Across day (mouse-level majority clusters)');
grid on; legend('off');
exportgraphics(gca, fullfile(outDir, ['across_day_majority_' ts '.png']), 'Resolution', 300);


%% ----- Build normalized progress summary (equal weight per session) -----
Ktrl = join(Ttrial, Tses(:,[keys {'Cluster'}]), 'Keys', keys);
Gsess = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx);
nT    = splitapply(@max, Ktrl.Trial, Gsess);
prog  = (Ktrl.Trial - 0.5) ./ nT(Gsess);
binEdges = linspace(0,1,NBINS+1);
[~,bin] = histc(prog, binEdges); bin(bin<1)=1; bin(bin>NBINS)=NBINS;

% per-session per-bin
Gsb = findgroups(Ktrl.mouse_key, Ktrl.day_index, Ktrl.session_idx, bin);
Ssb = table;
Ssb.mouse_key   = splitapply(@unique, Ktrl.mouse_key,   Gsb);
Ssb.day_index   = splitapply(@unique, Ktrl.day_index,   Gsb);
Ssb.session_idx = splitapply(@unique, Ktrl.session_idx, Gsb);
Ssb.Cluster     = splitapply(@unique, Ktrl.Cluster,     Gsb);
Ssb.bin         = splitapply(@unique, bin,              Gsb);
Ssb.meanLR      = splitapply(@mean,  Ktrl.lick_rate,    Gsb);

% cluster mean±SEM across sessions
Gcb = findgroups(Ssb.Cluster, Ssb.bin);
Csum = table;
Csum.Cluster = splitapply(@unique, Ssb.Cluster, Gcb);
Csum.bin     = splitapply(@unique, Ssb.bin,     Gcb);
Csum.nSess   = splitapply(@numel,  Ssb.meanLR,  Gcb);
Csum.mu      = splitapply(@mean,   Ssb.meanLR,  Gcb);
Csum.sem     = splitapply(@(x)nanstd(x,0)/sqrt(numel(x)), Ssb.meanLR, Gcb);
xPct = (Csum.bin - 0.5) / NBINS * 100;

%% ===== Figure A: Normalized progress (left) + ID grids (right) =====
fA = figure('Color','w','Position',[60 60 1400 600]);
tlo = tiledlayout(fA,1,2,'TileSpacing','compact','Padding','compact');

% Left: curves
ax1 = nexttile(tlo,1); hold(ax1,'on');
for k=1:K
    r = (Csum.Cluster==k); [xx,ord] = sort(xPct(r));
    mu=Csum.mu(r); mu=mu(ord); se=Csum.sem(r); se=se(ord);
    plot(ax1, xx, mu, '-', 'Color', C(k,:), 'LineWidth', 1.6, ...
         'Marker', mk{min(k,numel(mk))}, 'MarkerSize',3, 'MarkerFaceColor',C(k,:));
    ii=1:max(1,round(numel(xx)/50)):numel(xx);
    eb = errorbar(ax1, xx(ii), mu(ii), se(ii), 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off'); % avoids legend entries
end
xlabel(ax1,'Session progress (%)'); ylabel(ax1,'Lick rate (Hz)');
title(ax1, sprintf('Clusters: mean \\pm SEM across normalized session progress (%d bins)', NBINS));
grid(ax1,'on'); legend(ax1,'off');

% Right: IDs per cluster, multi-column grids
axIDs = nexttile(tlo,2); axis(axIDs,'off'); title(axIDs,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;  % adjust maxRows to make columns shorter/taller
for k=1:K
    ids = sort(unique(Ssb.mouse_key(Ssb.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);  % cellstr of column text
    ncol = numel(txtBlocks);
    x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end
outA = fullfile(outDir, ['progress_meanSEM_IDGRID_' ts '.png']);
exportgraphics(fA, outA, 'Resolution', 300); disp(['saved: ' outA]);

%% ===== Figure B: Across-day mean±SEM (no legend) + ID grids =====
% aggregate by day × cluster
Gd = findgroups(Tses.day_index, Tses.Cluster);
Dsum = table;
Dsum.day_index = splitapply(@unique, Tses.day_index, Gd);
Dsum.Cluster   = splitapply(@unique, Tses.Cluster,   Gd);
Dsum.mu        = splitapply(@mean,  Tses.lick_per_min, Gd);
Dsum.sem       = splitapply(@(x) nanstd(x,0)/sqrt(numel(x)), Tses.lick_per_min, Gd);

fB = figure('Color','w','Position',[60 60 1400 600]);
t2 = tiledlayout(fB,1,2,'TileSpacing','compact','Padding','compact');

ax2 = nexttile(t2,1); hold(ax2,'on');
for k=1:K
    r = Dsum.Cluster==k; [d,ord] = sort(Dsum.day_index(r));
    mu=Dsum.mu(r); mu=mu(ord); se=Dsum.sem(r); se=se(ord);
    p = plot(ax2, d, mu, '-', 'Color', C(k,:), 'LineWidth', 1.8, ...
             'Marker', mk{min(k,numel(mk))}, 'MarkerFaceColor', C(k,:));
    eb = errorbar(ax2, d, mu, se, 'LineStyle','none', 'Color', C(k,:)*0.7);
    set(eb,'HandleVisibility','off'); %#ok<*NASGU>
end
xlabel(ax2,'Day'); ylabel(ax2,'Licks / min'); title(ax2,'Clusters: mean \pm SEM across day');
grid(ax2,'on'); legend(ax2,'off');

axIDs2 = nexttile(t2,2); axis(axIDs2,'off'); title(axIDs2,'Mouse IDs by cluster');
y0 = 0.95; gap = 0.30; maxRows = 18;
for k=1:K
    ids = sort(unique(Tses.mouse_key(Tses.Cluster==k)));
    txtBlocks = idColumns(ids, maxRows);
    ncol = numel(txtBlocks); x0 = 0.02; colW = 0.95/max(1,ncol);
    for c=1:ncol
        text(axIDs2, x0 + (c-1)*colW, y0 - (k-1)*gap, sprintf('Cluster %d:\n%s',k, txtBlocks{c}), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8, ...
            'Interpreter','none','BackgroundColor','w','Margin',2,'EdgeColor',[.8 .8 .8]);
    end
end
outB = fullfile(outDir, ['across_day_meanSEM_IDGRID_' ts '.png']);
exportgraphics(fB, outB, 'Resolution', 300); disp(['saved: ' outB]);

%% ===== Also write plain-text & CSV ID lists =====
fid = fopen(fullfile(outDir, ['mouse_ids_by_cluster_' ts '.txt']),'w');
for k=1:K
    ids = sort(unique(Tses.mouse_key(Tses.Cluster==k)));
    fprintf(fid,'Cluster %d (%d mice):\n', k, numel(ids));
    fprintf(fid,'  %s\n\n', strjoin(ids, ', '));
end
fclose(fid);
writetable(groupsummary(Tses(:,{'mouse_key','Cluster'}),'Cluster'), ...
    fullfile(outDir, ['mouse_counts_per_cluster_' ts '.csv']));

disp('All done.');
end

% ---------- helpers ----------
function blocks = idColumns(ids, maxRows)
% Split a long list of IDs into columns of ~maxRows lines.
n = numel(ids); ncol = max(1, ceil(n/maxRows)); nrow = ceil(n/ncol);
blocks = cell(1,ncol);
for c = 1:ncol
    i1 = (c-1)*nrow + 1; i2 = min(c*nrow, n);
    blocks{c} = strjoin(ids(i1:i2), newline);
end
end
