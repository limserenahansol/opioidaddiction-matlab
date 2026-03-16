function analyze_advanced_pipeline()
% analyze_advanced_pipeline
% Implements Tasks 5.1-13.2 from the Behavior Pipeline Implementation Plan.
%
% Uses preprocessed data from analyze_passive_active_dashboard_dec2 (S_D_cache.mat).
% Loads raw CSV only for analyses requiring frame-level data (Tasks 5.3, 13.1).
%
% Tasks:
%   5.1  Motivation Across Phases (LME)
%   5.2  Licking Behavior Across Phases (distributions + LME)
%   5.3  Within-Session Temporal Dynamics (5-min bins)
%   6.1  Active vs Passive (During Phase)
%   6.2  Post-During Recovery
%   6.3  Reward Type Effects
%   7.1  Baseline Characterization (K-means clustering)
%   7.2  Individual Response to Morphine
%  10.1  Baseline Pupil Across Phases
%  10.3  Pupil-Behavior Correlations
%  11.1  TI (Tail Immersion) Analysis
%  11.2  TST/HP Analysis
%  11.3  Straub Tail Analysis
%  12.1  Composite Indices (PCA)
%  12.2  Predictive Modeling (TreeBagger)
%  13.1  Within-Session State Detection (changepoint)
%  13.2  Learning Curves (exponential fit)
%
% Requires: Statistics and Machine Learning Toolbox

%% ========== 1. LOCATE DATA ==========
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d), 'No run_* under %s', rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder, d(ix).name);

cacheMat = fullfile(runDir, 'S_D_cache.mat');
csvPath  = fullfile(runDir, 'ALL_mice_longitudinal.csv');
assert(exist(cacheMat,'file')>0, 'Missing cache: %s', cacheMat);
assert(exist(csvPath,'file')>0,  'Missing CSV: %s',   csvPath);

%% ========== 2. LOAD PREPROCESSED DATA ==========
fprintf('Loading S/D from cache: %s\n', cacheMat);
L = load(cacheMat, 'S', 'D');
S = L.S;  D = L.D;

% --- GROUP ASSIGNMENT using robust cohort roster (from dashboard_dec2) ---
% DO NOT rely on S.isPassive or D.GroupType from cache.
% Instead, parse cage+color from mouse_key and map via hard-coded roster.
cohort = buildNewCohortRoster();

% DEBUG: print raw mouse_key samples BEFORE group assignment
fprintf('\n--- DEBUG: mouse_key samples from D ---\n');
mkSamples = unique(string(D.mouse_key), 'stable');
for dbg_i = 1:min(numel(mkSamples), 14)
    mk_dbg = mkSamples(dbg_i);
    mk_norm_dbg = lower(regexprep(strtrim(mk_dbg), '[_\-/]+', ' '));
    tok_dbg = regexp(mk_norm_dbg, '(\d{4})', 'tokens', 'once');
    cage_dbg = ""; if ~isempty(tok_dbg), cage_dbg = string(tok_dbg{1}); end
    col_dbg = "";
    for cc = ["black","red","orange","white"]
        if contains(mk_norm_dbg, cc), col_dbg = cc; break; end
    end
    fprintf('  raw="%s"  norm="%s"  cage=%s  color=%s\n', mk_dbg, mk_norm_dbg, cage_dbg, col_dbg);
end
fprintf('--- END DEBUG ---\n\n');

D = attachCohortGroupAndPair(D, cohort);
S = attachCohortGroupAndPair(S, cohort);

% Diagnostic: print unique mice and their group assignment
u = unique(D(:,{'mouse_key','Group'}), 'rows', 'stable');
nA = nnz(u.Group == "Active");
nP = nnz(u.Group == "Passive");
nU = nnz(isundefined(u.Group) | ismissing(u.Group));
fprintf('\n========================================\n');
fprintf('[GROUP CHECK] Active=%d, Passive=%d, Unmapped=%d (total=%d mice)\n', nA, nP, nU, height(u));
fprintf('Expected: Active=6, Passive=8\n');
fprintf('========================================\n');
disp(u);

if nA == 0 || nP == 0
    warning(['Group assignment FAILED (one group is empty). ', ...
             'This means mouse_key format does not match cohort roster. ', ...
             'Check DEBUG output above for mouse_key format.']);
end
if (nA + nP) == 14 && nU == 0
    fprintf('[GROUP CHECK] SUCCESS: all 14 mice mapped correctly.\n\n');
elseif nU > 0
    fprintf('[GROUP CHECK] WARNING: %d mice could not be mapped.\n\n', nU);
end

% Remove rows where Group could not be assigned
D = D(~isundefined(D.Group) & ~ismissing(D.Group), :);

% Exclude habituation days (1-2)
if ismember('day_index', D.Properties.VariableNames)
    D = D(double(D.day_index) >= 3, :);
end

% Passive mice: set PR = NaN on days 6-10 (passive forced period, no real PR)
if ismember('RequirementLast', D.Properties.VariableNames)
    mask = (D.Group == "Passive") & (double(D.day_index) >= 6) & (double(D.day_index) <= 10);
    if any(mask)
        fprintf('Setting RequirementLast=NaN for Passive during days 6-10 (%d rows).\n', nnz(mask));
        D.RequirementLast(mask) = NaN;
    end
end

% Passive mice: set LICKING-related metrics = NaN on days 6-10
% (no licking data for Passive during this period)
maskLick = (D.Group == "Passive") & (double(D.day_index) >= 6) & (double(D.day_index) <= 10);
if any(maskLick)
    vnames = D.Properties.VariableNames;
    lickCols = vnames(startsWith(vnames, "lick_") | startsWith(vnames, "bout_"));
    if ~isempty(lickCols)
        for i = 1:numel(lickCols)
            col = lickCols{i};
            if isnumeric(D.(col))
                D.(col)(maskLick) = NaN;
            end
        end
        fprintf('Setting licking/bout metrics=NaN for Passive during days 6-10 (%d rows).\n', nnz(maskLick));
    end
end

% Assign period and reward type
D.Period     = periodOfDay(double(D.day_index));
D.RewardType = rewardTypeOfDay(double(D.day_index));
D = D(~isundefined(D.Period), :);

% Normalize pupil column name for downstream analyses
if ~ismember('pupil_mean', D.Properties.VariableNames)
    pupilCol = pickPupilCol(D);
    if ~isempty(pupilCol)
        D.pupil_mean = double(D.(pupilCol));
        fprintf('Using %s as pupil_mean for analyses.\n', pupilCol);
    end
end

%% ========== 2b. FINAL GROUP VERIFICATION ==========
% This block runs AFTER all filtering. If Active/Passive counts are wrong,
% the problem is in group assignment above.
fprintf('\n============ FINAL DATA SUMMARY (after all filtering) ============\n');
fprintf('Total rows in D: %d\n', height(D));
fprintf('Active rows:  %d\n', nnz(D.Group == "Active"));
fprintf('Passive rows: %d\n', nnz(D.Group == "Passive"));
uFinal = unique(D(:,{'mouse_key','Group'}), 'rows', 'stable');
fprintf('Unique mice: %d  (Active=%d, Passive=%d)\n', height(uFinal), ...
    nnz(uFinal.Group=="Active"), nnz(uFinal.Group=="Passive"));
fprintf('Active mice:  '); disp(string(uFinal.mouse_key(uFinal.Group=="Active"))');
fprintf('Passive mice: '); disp(string(uFinal.mouse_key(uFinal.Group=="Passive"))');
fprintf('==================================================================\n\n');

% HARD ASSERT: if separation failed, stop before wasting time on plots
assert(nnz(uFinal.Group=="Active") >= 1 && nnz(uFinal.Group=="Passive") >= 1, ...
    'FATAL: Group separation failed! Check DEBUG output above. Active=%d, Passive=%d', ...
    nnz(uFinal.Group=="Active"), nnz(uFinal.Group=="Passive"));

%% ========== 3. OUTPUT DIRECTORIES ==========
baseOut = fullfile(runDir, 'figs', 'advanced_pipeline');
dirs = struct();
dirs.task5_1 = fullfile(baseOut, 'task_5_1_motivation_lme');
dirs.task5_2 = fullfile(baseOut, 'task_5_2_licking_behavior');
dirs.task5_3 = fullfile(baseOut, 'task_5_3_temporal_dynamics');
dirs.task6_1 = fullfile(baseOut, 'task_6_1_group_during');
dirs.task6_2 = fullfile(baseOut, 'task_6_2_recovery');
dirs.task6_3 = fullfile(baseOut, 'task_6_3_reward_type');
dirs.task7_1 = fullfile(baseOut, 'task_7_1_clustering');
dirs.task7_2 = fullfile(baseOut, 'task_7_2_individual_morphine');
dirs.task10   = fullfile(baseOut, 'task_10_pupil');
dirs.task11   = fullfile(baseOut, 'task_11_pharmacology');
dirs.task12_1 = fullfile(baseOut, 'task_12_1_composite_indices');
dirs.task12_2 = fullfile(baseOut, 'task_12_2_predictive_modeling');
dirs.task13_1 = fullfile(baseOut, 'task_13_1_state_detection');
dirs.task13_2 = fullfile(baseOut, 'task_13_2_learning_curves');
dirs.stats    = fullfile(baseOut, 'stats');

fn = fieldnames(dirs);
for i = 1:numel(fn), if ~exist(dirs.(fn{i}),'dir'), mkdir(dirs.(fn{i})); end; end

%% ========== 4. RUN ALL TASKS ==========
fprintf('\n====== TASK 5.1: Motivation Across Phases (LME) ======\n');
try, task_5_1_motivation_lme(D, dirs.task5_1, dirs.stats);
catch ME, fprintf('  [WARN] Task 5.1 failed: %s\n', ME.message); end

fprintf('\n====== TASK 5.2: Licking Behavior Across Phases ======\n');
try, task_5_2_licking_behavior(D, dirs.task5_2, dirs.stats);
catch ME, fprintf('  [WARN] Task 5.2 failed: %s\n', ME.message); end

fprintf('\n====== TASK 5.3: Within-Session Temporal Dynamics ======\n');
try, task_5_3_temporal_dynamics(csvPath, S, D, dirs.task5_3, dirs.stats);
catch ME, fprintf('  [WARN] Task 5.3 failed: %s\n', ME.message); end

fprintf('\n====== TASK 6.1: Active vs Passive (During Phase) ======\n');
try, task_6_1_group_during(D, dirs.task6_1, dirs.stats);
catch ME, fprintf('  [WARN] Task 6.1 failed: %s\n', ME.message); end

fprintf('\n====== TASK 6.2: Post-During Recovery ======\n');
try, task_6_2_recovery(D, dirs.task6_2, dirs.stats);
catch ME, fprintf('  [WARN] Task 6.2 failed: %s\n', ME.message); end

fprintf('\n====== TASK 6.3: Reward Type Effects ======\n');
try, task_6_3_reward_type(D, dirs.task6_3, dirs.stats);
catch ME, fprintf('  [WARN] Task 6.3 failed: %s\n', ME.message); end

fprintf('\n====== TASK 7.1: Baseline Characterization (Clustering) ======\n');
try, task_7_1_clustering(D, dirs.task7_1, dirs.stats);
catch ME, fprintf('  [WARN] Task 7.1 failed: %s\n', ME.message); end

fprintf('\n====== TASK 7.2: Individual Response to Morphine ======\n');
try, task_7_2_morphine_response(D, dirs.task7_2, dirs.stats);
catch ME, fprintf('  [WARN] Task 7.2 failed: %s\n', ME.message); end

fprintf('\n====== TASK 10: Pupil Analysis ======\n');
try, task_10_pupil(D, csvPath, dirs.task10, dirs.stats);
catch ME, fprintf('  [WARN] Task 10 failed: %s\n', ME.message); end

fprintf('\n====== TASK 11: Pharmacology (TI / TST / HP / Straub) ======\n');
try, task_11_pharmacology(D, dirs.task11, dirs.stats);
catch ME, fprintf('  [WARN] Task 11 failed: %s\n', ME.message); end

fprintf('\n====== TASK 12.1: Composite Indices ======\n');
try, task_12_1_composite(D, dirs.task12_1, dirs.stats);
catch ME, fprintf('  [WARN] Task 12.1 failed: %s\n', ME.message); end

fprintf('\n====== TASK 12.2: Predictive Modeling ======\n');
try, task_12_2_predictive(D, dirs.task12_2, dirs.stats);
catch ME, fprintf('  [WARN] Task 12.2 failed: %s\n', ME.message); end

fprintf('\n====== TASK 13.1: Within-Session State Detection ======\n');
try, task_13_1_state_detection(csvPath, S, D, dirs.task13_1, dirs.stats);
catch ME, fprintf('  [WARN] Task 13.1 failed: %s\n', ME.message); end

fprintf('\n====== TASK 13.2: Learning Curves ======\n');
try, task_13_2_learning_curves(D, dirs.task13_2, dirs.stats);
catch ME, fprintf('  [WARN] Task 13.2 failed: %s\n', ME.message); end

fprintf('\n===================================================\n');
fprintf('Advanced pipeline complete. Outputs in:\n  %s\n', baseOut);
end

%% ########################################################################
%  TASK 5.1 — Motivation Across Phases (Linear Mixed-Effects Model)
%  ########################################################################
function task_5_1_motivation_lme(D, outDir, statsDir)
% Mixed-effects:  RequirementLast ~ Phase * Group + (1|mouse_key)
% Delta scores + pairwise comparisons
% Includes: model-predicted means + CI, observed-vs-fitted, plain-language report

ycol = 'RequirementLast';
if ~ismember(ycol, D.Properties.VariableNames)
    fprintf('  Skipping 5.1: %s not found.\n', ycol); return;
end

% --- LME ---
T = D(:, {'mouse_key','day_index','Period','Group', ycol});
T = T(isfinite(T.(ycol)), :);
T.Phase = removecats(categorical(string(T.Period)));
T.Grp   = removecats(categorical(string(T.Group)));
T.mouse = categorical(T.mouse_key);
T.y     = double(T.(ycol));

if height(T) < 10
    fprintf('  Skipping 5.1: too few observations (N=%d).\n', height(T)); return;
end

% EXCLUDE "During" entirely for Task 5.1
Tuse = T(string(T.Phase) ~= "During", :);
Tuse.Phase = removecats(Tuse.Phase);
if height(Tuse) < 10
    error('Task 5.1 has too few observations after excluding During.');
end

% Fit INTERACTION model on During-excluded data (this is mandatory)
lmeUse = fitlme(Tuse, 'y ~ Phase * Grp + (1|mouse)');
lmeFormulaUse = 'Phase * Grp (During excl.)';
isAdditiveUse = false;
fprintf('  Interaction model fitted successfully on During-excluded data.\n');

% --- ANOVA: use DEFAULT method first (this gives proper p-values) ---
aovDefault = anova(lmeUse);
fprintf('  ANOVA (default DF method):\n');
disp(aovDefault);

% Also try Satterthwaite — use it only if it gives FINITE p-values
aovUse = aovDefault;  % start with default
try
    aovSat = anova(lmeUse, 'DFMethod','satterthwaite');
    pTest = extractAnovP(aovSat, 'Phase');
    if isfinite(pTest)
        aovUse = aovSat;
        fprintf('  Satterthwaite ANOVA also available (Phase p=%.4g).\n', pTest);
    else
        fprintf('  Satterthwaite returned NaN p-values; using default ANOVA.\n');
    end
catch
    fprintf('  Satterthwaite not available; using default ANOVA.\n');
end

% Phases used for plots
periodsUse = ["Pre","Post","Withdrawal","Re-exposure"];

% Report which model is used
fprintf('  Using During-excluded INTERACTION model as PRIMARY result.\n');
disp(lmeUse);
fprintf('ANOVA table (formula: y ~ %s + (1|mouse)):\n', lmeFormulaUse);
disp(aovUse);

% Save LME tables
safeWriteLmeResults(lmeUse, aovUse, statsDir, 'task5_1');

% --- Extract p-values (multiple strategies) ---
pPh  = extractAnovP(aovUse, 'Phase');
pGrp = extractAnovP(aovUse, 'Grp');
pInt = extractAnovP(aovUse, 'Phase:Grp');
fprintf('  [Strategy 1] extractAnovP: Phase=%.4g, Grp=%.4g, Int=%.4g\n', pPh, pGrp, pInt);

% Strategy 2: If still NaN, compute from F-stat and DF using fcdf
if ~isfinite(pPh) || ~isfinite(pGrp)
    fprintf('  Strategy 1 returned NaN — trying fcdf fallback on default ANOVA.\n');
    [pPh, pGrp, pInt] = computePfromFcdf(aovDefault, {'Phase','Grp','Phase:Grp'});
    fprintf('  [Strategy 2] fcdf: Phase=%.4g, Grp=%.4g, Int=%.4g\n', pPh, pGrp, pInt);
end

% Strategy 3: absolute last resort — try coefTest
if ~isfinite(pPh) || ~isfinite(pGrp)
    fprintf('  Strategy 2 also NaN — trying coefTest.\n');
    try
        [~, pGrp] = coefTest(lmeUse, [0 0 0 0 0 1 0 0 0]);
    catch, end
    fprintf('  [Strategy 3] coefTest Grp: %.4g\n', pGrp);
end

fprintf('  === FINAL p-values for plots: Phase=%.4g, Grp=%.4g, Phase:Grp=%.4g ===\n', pPh, pGrp, pInt);

COL = groupColors();

% =====================================================================
% FIGURE 1: Raw data (mean +/- SEM) by phase and group
% =====================================================================
fh = figure('Color','w','Position',[80 80 900 540]); hold on;
for gi = ["Active","Passive"]
    Tg = Tuse(Tuse.Grp == gi, :);
    if isempty(Tg), continue; end
    pIdx = nan(1, numel(periodsUse));
    mu   = nan(1, numel(periodsUse));
    se   = nan(1, numel(periodsUse));
    for pi = 1:numel(periodsUse)
        vals = Tg.y(string(Tg.Phase) == periodsUse(pi));
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        pIdx(pi) = pi;
        mu(pi)   = mean(vals);
        se(pi)   = std(vals) / sqrt(numel(vals));
    end
    good = isfinite(mu);
    c = COL.(char(gi));
    fill([pIdx(good), fliplr(pIdx(good))], [mu(good)-se(good), fliplr(mu(good)+se(good))], ...
        c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(pIdx(good), mu(good), '-o', 'Color', c, 'LineWidth', 2, 'MarkerFaceColor', c, ...
        'DisplayName', char(gi));
end

% Per-phase between-group ranksum (Holm-corrected)
rawPtraj = nan(1, numel(periodsUse));
for pi = 1:numel(periodsUse)
    vA = Tuse.y(Tuse.Grp=="Active"  & string(Tuse.Phase)==periodsUse(pi));
    vP = Tuse.y(Tuse.Grp=="Passive" & string(Tuse.Phase)==periodsUse(pi));
    rawPtraj(pi) = safeRanksum(vA, vP);
end
adjPtraj = holmCorrect(rawPtraj);
yl = ylim; topY = yl(2); pad = 0.06*diff(yl);
for pi = 1:numel(periodsUse)
    if isfinite(adjPtraj(pi)) && adjPtraj(pi) < 0.05
        text(pi, topY + pad*0.2, starStr(adjPtraj(pi)), ...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',10,'Color',[0.8 0 0]);
    end
end
ylim([yl(1), topY + pad*0.8]);

set(gca,'XTick',1:numel(periodsUse),'XTickLabel',periodsUse);
xlabel('Phase'); ylabel('PR Breakpoint (RequirementLast)');
titleParts = {};
titleParts{end+1} = sprintf('Phase p=%s %s', fmtP(pPh), starStr(pPh));
titleParts{end+1} = sprintf('Grp p=%s %s', fmtP(pGrp), starStr(pGrp));
if isfinite(pInt)
    titleParts{end+1} = sprintf('Phase:Grp p=%s %s', fmtP(pInt), starStr(pInt));
end
title({sprintf('Motivation (LME: %s, During excl.)', lmeFormulaUse), ...
       strjoin(titleParts, ', ')});
legend('Location','eastoutside'); grid off; box off;

% Annotation box with per-phase tests
statsLines = {};
for pi = 1:numel(periodsUse)
    if isfinite(rawPtraj(pi))
        nA = nnz(Tuse.Grp=="Active"  & string(Tuse.Phase)==periodsUse(pi) & isfinite(Tuse.y));
        nP = nnz(Tuse.Grp=="Passive" & string(Tuse.Phase)==periodsUse(pi) & isfinite(Tuse.y));
        statsLines{end+1} = sprintf('%s: Ranksum A(n=%d) vs P(n=%d) p=%s (Holm adj=%s) %s', ...
            periodsUse(pi), nA, nP, fmtP(rawPtraj(pi)), fmtP(adjPtraj(pi)), starStr(adjPtraj(pi))); %#ok<AGROW>
    end
end
if ~isempty(statsLines), addStatsTextbox(gca, statsLines); end
printpng(fh, fullfile(outDir, 'motivation_trajectory_lme.png')); close(fh);

% =====================================================================
% FIGURE 2: Model-predicted means + 95% CI  (the LME fit visualized)
% =====================================================================
% Build a "new data" table with every Phase x Group combination
% and use predict() to get model estimates with CI
predPhases = categories(Tuse.Phase);
predGrps   = categories(Tuse.Grp);
newRows = table();
for pi = 1:numel(predPhases)
    for gi = 1:numel(predGrps)
        row = table();
        row.Phase = categorical(predPhases(pi), predPhases);
        row.Grp   = categorical(predGrps(gi), predGrps);
        row.mouse = categorical("dummy_mouse");  % needed for random effect syntax
        newRows = [newRows; row]; %#ok<AGROW>
    end
end

% Predict with CI (conditional=false => marginal / fixed-effects only)
try
    [yPred, yCI] = predict(lmeUse, newRows, 'Conditional', false, 'Alpha', 0.05);
    newRows.yPred = yPred;
    newRows.yLo   = yCI(:,1);
    newRows.yHi   = yCI(:,2);
    hasPredictions = true;
catch ME
    fprintf('  Could not generate LME predictions: %s\n', ME.message);
    hasPredictions = false;
end

if hasPredictions
    fh = figure('Color','w','Position',[80 80 950 520]); hold on;
    predActive = nan(1, numel(periodsUse));
    predPassive = nan(1, numel(periodsUse));

    for gi = ["Active","Passive"]
        c = COL.(char(gi));
        xPlot = nan(1, numel(periodsUse));
        yPlot = nan(1, numel(periodsUse));
        yLo   = nan(1, numel(periodsUse));
        yHi   = nan(1, numel(periodsUse));
        for pi = 1:numel(periodsUse)
            idx = (string(newRows.Phase)==periodsUse(pi)) & (string(newRows.Grp)==gi);
            if ~any(idx), continue; end
            xPlot(pi) = pi;
            yPlot(pi) = newRows.yPred(idx);
            yLo(pi)   = newRows.yLo(idx);
            yHi(pi)   = newRows.yHi(idx);
        end
        good = isfinite(yPlot);
        if ~any(good), continue; end

        % CI shading
        fill([xPlot(good), fliplr(xPlot(good))], ...
             [yLo(good), fliplr(yHi(good))], ...
             c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
        % Predicted mean line
        plot(xPlot(good), yPlot(good), '-s', 'Color', c, 'LineWidth', 2.5, ...
            'MarkerFaceColor', c, 'MarkerSize', 9, 'DisplayName', sprintf('%s (model)', char(gi)));

        % Also overlay raw means as open circles (per-mouse median to avoid bias)
        Tg = Tuse(Tuse.Grp == gi, :);
        [gpm, mkp, php] = findgroups(Tg.mouse, Tg.Phase);
        pm = table(mkp, php, splitapply(@(x) median(x,'omitnan'), Tg.y, gpm), ...
            'VariableNames', {'mouse','Phase','y'});
        for pi2 = 1:numel(periodsUse)
            vals = pm.y(string(pm.Phase) == periodsUse(pi2));
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            plot(pi2 + 0.08*(gi=="Passive"), mean(vals), 'o', 'Color', c, ...
                'MarkerSize', 7, 'LineWidth', 1.5, 'HandleVisibility','off');
        end

        if gi == "Active"
            predActive = yPlot;
        else
            predPassive = yPlot;
        end
    end

    % Dummy handles for legend
    hRaw = plot(nan, nan, 'ko', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'Raw mean');

    set(gca,'XTick',1:numel(periodsUse),'XTickLabel',periodsUse);
    xlabel('Phase'); ylabel('PR Breakpoint (RequirementLast)');
    % Build informative title with actual p-values
    titleLine1 = sprintf('LME Predicted Means (During excl.)  [%s]', lmeFormulaUse);
    titleLine2Parts = {};
    titleLine2Parts{end+1} = sprintf('Phase p=%s %s', fmtP(pPh), starStr(pPh));
    titleLine2Parts{end+1} = sprintf('Grp p=%s %s', fmtP(pGrp), starStr(pGrp));
    if isfinite(pInt)
        titleLine2Parts{end+1} = sprintf('Phase:Grp p=%s %s', fmtP(pInt), starStr(pInt));
    end
    title({titleLine1, strjoin(titleLine2Parts, ', ')});
    legend('Location','eastoutside'); grid off; box off;

    % Add annotation box with key stats
    statTxt = {};
    statTxt{end+1} = sprintf('Phase: F=%.2f, p=%s %s', safeAnovF(aovUse,'Phase'), fmtP(pPh), starStr(pPh));
    statTxt{end+1} = sprintf('Group: F=%.2f, p=%s %s', safeAnovF(aovUse,'Grp'), fmtP(pGrp), starStr(pGrp));
    if isfinite(pInt)
        statTxt{end+1} = sprintf('Phase x Group: F=%.2f, p=%s %s', safeAnovF(aovUse,'Phase:Grp'), fmtP(pInt), starStr(pInt));
    end
    statTxt{end+1} = sprintf('N = %d obs, %d mice | squares = model, circles = raw mean', height(Tuse), numel(unique(Tuse.mouse)));
    annotation('textbox', [0.14 0.02 0.55 0.15], 'String', strjoin(statTxt, '\n'), ...
        'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], ...
        'FontSize', 9, 'Interpreter','none');

    printpng(fh, fullfile(outDir, 'motivation_LME_predicted_means.png')); close(fh);
end

% (Figure 2c removed — During-excluded model is now the primary model)

% =====================================================================
% FIGURE 2b: Separate LME per group (Active-only / Passive-only)
% =====================================================================
fh = figure('Color','w','Position',[80 80 950 520]); hold on;
for gi = ["Active","Passive"]
    Tg = Tuse(Tuse.Grp == gi, :);
    fprintf('  Fig2b: %s group has %d rows, %d mice\n', gi, height(Tg), numel(unique(Tg.mouse)));
    if height(Tg) < 6
        fprintf('  Fig2b: Skipping %s — too few rows.\n', gi);
        continue;
    end
    try
        lmeG = fitlme(Tg, 'y ~ Phase + (1|mouse)');
    catch ME
        fprintf('  Fig2b: LME failed for %s: %s\n', gi, ME.message);
        continue;
    end
    aovG = anova(lmeG);
    pPhG = extractAnovP(aovG, 'Phase');
    if ~isfinite(pPhG)
        pPhG = computePfromFcdf(aovG, {'Phase'});
    end

    % Predict per-phase means with CI
    predPhasesG = categories(Tg.Phase);
    newG = table();
    for pi = 1:numel(predPhasesG)
        row = table();
        row.Phase = categorical(predPhasesG(pi), predPhasesG);
        row.mouse = categorical("dummy_mouse");
        newG = [newG; row]; %#ok<AGROW>
    end
    try
        [yPredG, yCIG] = predict(lmeG, newG, 'Conditional', false, 'Alpha', 0.05);
    catch
        continue;
    end

    % Align predictions to global phase order
    xPlot = nan(1, numel(periodsUse));
    yPlot = nan(1, numel(periodsUse));
    yLo   = nan(1, numel(periodsUse));
    yHi   = nan(1, numel(periodsUse));
    for pi = 1:numel(periodsUse)
        idx = string(newG.Phase) == periodsUse(pi);
        if any(idx)
            xPlot(pi) = pi;
            yPlot(pi) = yPredG(idx);
            yLo(pi)   = yCIG(idx,1);
            yHi(pi)   = yCIG(idx,2);
        end
    end
    good = isfinite(yPlot);
    c = COL.(char(gi));
    fill([xPlot(good), fliplr(xPlot(good))], [yLo(good), fliplr(yHi(good))], ...
        c, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(xPlot(good), yPlot(good), '-s', 'Color', c, 'LineWidth', 2.5, ...
        'MarkerFaceColor', c, 'MarkerSize', 9, ...
        'DisplayName', sprintf('%s model (Phase p=%s %s)', char(gi), fmtP(pPhG), starStr(pPhG)));

    % Overlay raw per-mouse medians as open circles
    [gpm, mkp, php] = findgroups(Tg.mouse, Tg.Phase);
    pm = table(mkp, php, splitapply(@(x) median(x,'omitnan'), Tg.y, gpm), ...
        'VariableNames', {'mouse','Phase','y'});
    for pi2 = 1:numel(periodsUse)
        vals = pm.y(string(pm.Phase) == periodsUse(pi2));
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        plot(pi2 + 0.08*(gi=="Passive"), mean(vals), 'o', 'Color', c, ...
            'MarkerSize', 7, 'LineWidth', 1.5, 'HandleVisibility','off');
    end
end
set(gca,'XTick',1:numel(periodsUse),'XTickLabel',periodsUse);
xlabel('Phase'); ylabel('PR Breakpoint (RequirementLast)');
title({'Separate LME per Group (model vs raw per-mouse median)', ...
       '(squares = model prediction, circles = raw median)'}); 
legend('Location','eastoutside'); grid off; box off;
printpng(fh, fullfile(outDir, 'motivation_LME_predicted_means_separate_groups.png')); close(fh);

% =====================================================================
% FIGURE 3: Observed vs Fitted per mouse (spaghetti + model fit)
% =====================================================================
Tuse.yFitted = fitted(lmeUse);
fh = figure('Color','w','Position',[80 80 1050 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel A: observed data per mouse (spaghetti) + group means
ax1 = nexttile; hold(ax1,'on');
title(ax1, sprintf('A) Observed Data  |  Phase p=%s, Grp p=%s', fmtP(pPh), fmtP(pGrp)));
for gi = ["Active","Passive"]
    c = COL.(char(gi));
    Tg = Tuse(Tuse.Grp == gi, :);
    mk = unique(Tg.mouse,'stable');
    % Mouse-level lines (light)
    for mi = 1:numel(mk)
        Tm = Tg(Tg.mouse == mk(mi), :);
        xm = nan(1, numel(periodsUse));
        ym = nan(1, numel(periodsUse));
        for pi = 1:numel(periodsUse)
            vals = Tm.y(string(Tm.Phase) == periodsUse(pi));
            if isempty(vals), continue; end
            xm(pi) = pi; ym(pi) = mean(vals);
        end
        good = isfinite(ym);
        if sum(good) >= 2
            plot(ax1, xm(good), ym(good), '-', 'Color', [c 0.25], 'LineWidth', 0.7, ...
                'HandleVisibility','off');
        end
        scatter(ax1, xm(good), ym(good), 15, c, 'filled', 'MarkerFaceAlpha', 0.4, ...
            'HandleVisibility','off');
    end
    % Group mean (thick)
    muG = nan(1, numel(periodsUse));
    for pi = 1:numel(periodsUse)
        v = Tg.y(string(Tg.Phase)==periodsUse(pi));
        v = v(isfinite(v));
        if ~isempty(v), muG(pi) = mean(v); end
    end
    good = isfinite(muG);
    plot(ax1, find(good), muG(good), '-o', 'Color', c, 'LineWidth', 2.5, ...
        'MarkerFaceColor', c, 'MarkerSize', 8, 'DisplayName', char(gi));
end
set(ax1,'XTick',1:numel(periodsUse),'XTickLabel',periodsUse);
ylabel(ax1,'RequirementLast'); legend(ax1,'Location','eastoutside');
grid(ax1,'off'); box(ax1,'off');

% Panel B: model fitted values per mouse (spaghetti) + model means
ax2 = nexttile; hold(ax2,'on');
title(ax2, sprintf('B) LME Fitted (interaction)  |  Phase:Grp p=%s %s', fmtP(pInt), starStr(pInt)));
for gi = ["Active","Passive"]
    c = COL.(char(gi));
    Tg = Tuse(Tuse.Grp == gi, :);
    mk = unique(Tg.mouse,'stable');
    for mi = 1:numel(mk)
        Tm = Tg(Tg.mouse == mk(mi), :);
        xm = nan(1, numel(periodsUse));
        ym = nan(1, numel(periodsUse));
        for pi = 1:numel(periodsUse)
            vals = Tm.yFitted(string(Tm.Phase) == periodsUse(pi));
            if isempty(vals), continue; end
            xm(pi) = pi; ym(pi) = mean(vals);
        end
        good = isfinite(ym);
        if sum(good) >= 2
            plot(ax2, xm(good), ym(good), '-', 'Color', [c 0.25], 'LineWidth', 0.7, ...
                'HandleVisibility','off');
        end
        scatter(ax2, xm(good), ym(good), 15, c, 'filled', 'MarkerFaceAlpha', 0.4, ...
            'HandleVisibility','off');
    end
    % Model predicted group mean (from newRows if available)
    if hasPredictions
        xp = nan(1,numel(periodsUse)); yp = nan(1,numel(periodsUse));
        for pi = 1:numel(periodsUse)
            idx = (string(newRows.Phase)==periodsUse(pi)) & (string(newRows.Grp)==gi);
            if any(idx), xp(pi) = pi; yp(pi) = newRows.yPred(idx); end
        end
        good = isfinite(yp);
        plot(ax2, xp(good), yp(good), '-s', 'Color', c, 'LineWidth', 2.5, ...
            'MarkerFaceColor', c, 'MarkerSize', 9, 'DisplayName', sprintf('%s (model)', char(gi)));
    end
end
set(ax2,'XTick',1:numel(periodsUse),'XTickLabel',periodsUse);
ylabel(ax2,'Fitted RequirementLast'); legend(ax2,'Location','eastoutside');
grid(ax2,'off'); box(ax2,'off');

printpng(fh, fullfile(outDir, 'motivation_observed_vs_fitted.png')); close(fh);

% =====================================================================
% FIGURE 4: Residual diagnostic plots
% =====================================================================
Tuse.residuals = residuals(lmeUse);
fh = figure('Color','w','Position',[80 80 950 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Residuals vs Fitted
ax1 = nexttile; hold(ax1,'on');
scatter(ax1, Tuse.yFitted, Tuse.residuals, 25, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.5);
yline(ax1, 0, 'k:');
xlabel(ax1,'Fitted values'); ylabel(ax1,'Residuals');
title(ax1,'Residuals vs Fitted'); grid(ax1,'off'); box(ax1,'off');

% QQ plot of residuals
ax2 = nexttile;
try
    qqplot(ax2, Tuse.residuals);
catch
    % fallback: manual QQ
    r = sort(Tuse.residuals);
    n = numel(r);
    theoretical = norminv(((1:n)' - 0.5) / n);
    plot(ax2, theoretical, r, 'o', 'MarkerSize', 4);
    hold(ax2,'on');
    plot(ax2, [-3 3], [-3 3]*std(r), 'r-');
end
title(ax2,'QQ Plot of Residuals'); grid(ax2,'off'); box(ax2,'off');

printpng(fh, fullfile(outDir, 'motivation_LME_diagnostics.png')); close(fh);

% =====================================================================
% FIGURE 5: Delta from Pre per mouse
% =====================================================================
mice = unique(Tuse.mouse,'stable');
preMed = nan(numel(mice),1);
for i = 1:numel(mice)
    v = Tuse.y(Tuse.mouse==mice(i) & string(Tuse.Phase)=="Pre");
    preMed(i) = median(v,'omitnan');
end
preMap = containers.Map(cellstr(string(mice)), num2cell(preMed));

Tuse.delta = nan(height(Tuse),1);
for r = 1:height(Tuse)
    mk = char(string(Tuse.mouse(r)));
    if isKey(preMap, mk) && isfinite(preMap(mk))
        Tuse.delta(r) = Tuse.y(r) - preMap(mk);
    end
end

% Use per-mouse median delta per phase (avoid repeated points per mouse)
[gmd, mkd, phd, grd] = findgroups(Tuse.mouse, Tuse.Phase, Tuse.Grp);
Td = table(mkd, phd, grd, splitapply(@(x) median(x,'omitnan'), Tuse.delta, gmd), ...
    'VariableNames', {'mouse','Phase','Grp','delta'});

% Sanity check: each mouse must map to exactly one Group
[gm1, mkc1] = findgroups(Td.mouse);
bad1 = false(numel(mkc1),1);
for ii = 1:numel(mkc1)
    ug = unique(string(Td.Grp(gm1==ii)));
    if numel(ug) > 1, bad1(ii) = true; end
end
if any(bad1)
    badMice = string(mkc1(bad1));
    error('Group mapping conflict in Task 5.1 delta: mouse appears in both groups: %s', ...
        strjoin(badMice, ', '));
end

fh = figure('Color','w','Position',[80 80 900 480]); hold on;
plotDeltaBoxByPhaseGroup(Td, periodsUse, COL);
title('RequirementLast: Change from Pre');
ylabel('RequirementLast - Pre median');
printpng(fh, fullfile(outDir, 'motivation_delta_from_pre.png')); close(fh);

% =====================================================================
% PLAIN-LANGUAGE INTERPRETATION REPORT
% =====================================================================
writeLmeInterpretation(lmeUse, aovUse, lmeFormulaUse, isAdditiveUse, Tuse, outDir, statsDir);

fprintf('  Task 5.1 done.\n');
end

%% --- Helper: Write plain-language LME interpretation ---
function writeLmeInterpretation(lme, aov, lmeFormula, isAdditive, T, outDir, statsDir)
fid = fopen(fullfile(statsDir, 'task5_1_LME_interpretation.txt'), 'w');
if fid < 0, return; end

nObs  = lme.NumObservations;
nMice = numel(unique(T.mouse));
pPh   = extractAnovP(aov, 'Phase');
pGrp  = extractAnovP(aov, 'Grp');
pInt  = extractAnovP(aov, 'Phase:Grp');
% fcdf fallback if extractAnovP fails
if ~isfinite(pPh) || ~isfinite(pGrp)
    [pPh, pGrp, pInt] = computePfromFcdf(aov, {'Phase','Grp','Phase:Grp'});
end

% --- Header ---
fprintf(fid, '=============================================================\n');
fprintf(fid, ' TASK 5.1 — LME Analysis: Plain-Language Interpretation\n');
fprintf(fid, '=============================================================\n\n');

% --- What is LME ---
fprintf(fid, '--- WHAT IS AN LME MODEL? ---\n\n');
fprintf(fid, 'LME = Linear Mixed-Effects model.\n');
fprintf(fid, 'It is a statistical method for analyzing data where you have\n');
fprintf(fid, 'REPEATED MEASUREMENTS from the same subjects (mice).\n\n');
fprintf(fid, 'A regular ANOVA or t-test treats every data point as independent.\n');
fprintf(fid, 'But your data has multiple days per mouse. An LME accounts for\n');
fprintf(fid, 'the fact that measurements from the same mouse are correlated.\n\n');
fprintf(fid, 'Think of it as: "What is the effect of Phase and Group on PR\n');
fprintf(fid, 'breakpoint, after accounting for the fact that each mouse has\n');
fprintf(fid, 'its own baseline level?"\n\n');

% --- Model formula ---
fprintf(fid, '--- YOUR MODEL ---\n\n');
fprintf(fid, 'Formula: RequirementLast ~ %s + (1 | mouse)\n\n', lmeFormula);
fprintf(fid, 'This means:\n');
fprintf(fid, '  RequirementLast = the outcome variable (PR breakpoint)\n');
fprintf(fid, '  Phase           = which experimental phase (Pre, During, Post, etc.)\n');
fprintf(fid, '  Grp             = Active vs Passive group\n');
fprintf(fid, '  (1 | mouse)     = "each mouse can have its own baseline"\n\n');

if contains(lmeFormula, 'During excluded')
    fprintf(fid, 'NOTE: During (days 6-10) is EXCLUDED from the model because\n');
    fprintf(fid, 'Passive mice have no PR values there. This allows a valid\n');
    fprintf(fid, 'Phase x Group interaction model to be fit.\n\n');
elseif isAdditive
    fprintf(fid, 'NOTE: The interaction model (Phase * Grp) was tried first, but\n');
    fprintf(fid, 'it was "rank-deficient" — this means some Phase x Group\n');
    fprintf(fid, 'combinations have no data (Passive mice have no PR score\n');
    fprintf(fid, 'during days 6-10). So the simpler additive model was used.\n');
    fprintf(fid, 'The additive model assumes the Group difference is the SAME\n');
    fprintf(fid, 'size across all phases.\n\n');
end

% --- Fixed vs Random effects ---
fprintf(fid, '--- FIXED vs RANDOM EFFECTS ---\n\n');
fprintf(fid, 'FIXED EFFECTS (population-level):\n');
fprintf(fid, '  These are the effects we CARE about scientifically.\n');
fprintf(fid, '  - Phase: Does PR breakpoint change across experimental phases?\n');
fprintf(fid, '  - Group: Do Active and Passive mice differ?\n');
fprintf(fid, '  These represent the average effect across ALL mice.\n\n');
fprintf(fid, 'RANDOM EFFECTS (subject-level):\n');
fprintf(fid, '  - (1 | mouse): Each mouse gets its own "baseline offset"\n');
fprintf(fid, '    Some mice naturally press more, some less.\n');
fprintf(fid, '    The random effect captures this individual variation\n');
fprintf(fid, '    so it does not contaminate the Phase/Group estimates.\n\n');

% --- Key results ---
fprintf(fid, '--- KEY RESULTS (ANOVA TABLE) ---\n\n');
fprintf(fid, 'The ANOVA table tests whether each factor has a\n');
fprintf(fid, 'statistically significant effect on PR breakpoint.\n\n');

fprintf(fid, '  Phase effect:  F = %.2f, p = %s %s\n', ...
    safeAnovF(aov,'Phase'), fmtP(pPh), starStr(pPh));
if pPh < 0.05
    fprintf(fid, '    --> SIGNIFICANT: PR breakpoint changes across phases.\n');
    fprintf(fid, '        The mice do not perform the same in all phases.\n\n');
else
    fprintf(fid, '    --> Not significant: No strong evidence that phases differ.\n\n');
end

fprintf(fid, '  Group effect:  F = %.2f, p = %s %s\n', ...
    safeAnovF(aov,'Grp'), fmtP(pGrp), starStr(pGrp));
if pGrp < 0.05
    fprintf(fid, '    --> SIGNIFICANT: Active and Passive mice differ overall.\n\n');
else
    fprintf(fid, '    --> Not significant: No strong evidence groups differ.\n\n');
end

if isfinite(pInt)
    fprintf(fid, '  Phase x Group: F = %.2f, p = %s %s\n', ...
        safeAnovF(aov,'Phase:Grp'), fmtP(pInt), starStr(pInt));
    if pInt < 0.05
        fprintf(fid, '    --> SIGNIFICANT: The Phase effect is DIFFERENT for Active\n');
        fprintf(fid, '        vs Passive mice (the lines are not parallel).\n\n');
    else
        fprintf(fid, '    --> Not significant: Both groups show similar phase patterns.\n\n');
    end
end

% --- Interpreting coefficients ---
fprintf(fid, '--- INTERPRETING THE COEFFICIENTS ---\n\n');
fprintf(fid, 'The fixed effects coefficients tell you:\n\n');

try
    coefTbl = lme.Coefficients;
    if isa(coefTbl, 'dataset'), coefTbl = dataset2table(coefTbl); end
    
    % Get reference levels
    fprintf(fid, '(Intercept) = %.1f\n', lme.fixedEffects(1));
    fprintf(fid, '  This is the PREDICTED mean PR breakpoint for the\n');
    fprintf(fid, '  REFERENCE group. The reference is typically the first\n');
    fprintf(fid, '  alphabetical level of each factor.\n');
    fprintf(fid, '  In your case: Phase = "During", Group = "Active"\n');
    fprintf(fid, '  So: Active mice in the During phase have a predicted\n');
    fprintf(fid, '  mean PR breakpoint of ~%.1f.\n\n', lme.fixedEffects(1));

    coefNames = lme.CoefficientNames;
    coefVals  = lme.fixedEffects;
    for ci = 2:numel(coefNames)
        nm = coefNames{ci};
        val = coefVals(ci);
        fprintf(fid, '%s = %.2f\n', nm, val);
        if contains(nm, 'Phase_')
            phaseName = strrep(nm, 'Phase_', '');
            if val > 0
                fprintf(fid, '  -> %s phase is %.1f HIGHER than the reference phase.\n\n', phaseName, val);
            elseif val < 0
                fprintf(fid, '  -> %s phase is %.1f LOWER than the reference phase.\n\n', phaseName, abs(val));
            else
                fprintf(fid, '  -> %s phase is the same as reference.\n\n', phaseName);
            end
        elseif contains(nm, 'Grp_')
            grpName = strrep(nm, 'Grp_', '');
            if val > 0
                fprintf(fid, '  -> %s mice score %.1f HIGHER than Active (on average).\n\n', grpName, val);
            elseif val < 0
                fprintf(fid, '  -> %s mice score %.1f LOWER than Active (on average).\n\n', grpName, abs(val));
            end
        end
    end
catch
    fprintf(fid, '  (Could not extract coefficient details automatically.)\n\n');
end

% --- Random effects interpretation ---
fprintf(fid, '--- RANDOM EFFECTS ---\n\n');
try
    reVar = lme.MSE;  % residual variance
    fprintf(fid, 'Residual standard deviation = %.2f\n', sqrt(reVar));
    fprintf(fid, '  This is the typical spread of individual observations\n');
    fprintf(fid, '  around the fitted line. A PR breakpoint observation is\n');
    fprintf(fid, '  typically within +/- %.1f of the model prediction.\n\n', 2*sqrt(reVar));
catch
end

fprintf(fid, 'Mouse random intercept:\n');
fprintf(fid, '  If this is very small (near zero), it means that after\n');
fprintf(fid, '  accounting for Phase and Group, there is little\n');
fprintf(fid, '  mouse-to-mouse baseline variation. The Phase and Group\n');
fprintf(fid, '  effects already explain the between-mouse differences.\n\n');

% --- Model fit ---
fprintf(fid, '--- MODEL FIT STATISTICS ---\n\n');
fprintf(fid, 'Number of observations: %d\n', nObs);
fprintf(fid, 'Number of mice: %d\n', nMice);
fprintf(fid, 'AIC = %.1f  (lower = better fit, used to compare models)\n', lme.ModelCriterion.AIC);
fprintf(fid, 'BIC = %.1f  (like AIC but penalizes complexity more)\n\n', lme.ModelCriterion.BIC);

% --- Significance stars ---
fprintf(fid, '--- WHAT DO THE STARS MEAN? ---\n\n');
fprintf(fid, '  ***  p < 0.001  (very strong evidence)\n');
fprintf(fid, '  **   p < 0.01   (strong evidence)\n');
fprintf(fid, '  *    p < 0.05   (evidence, but moderate)\n');
fprintf(fid, '  n.s. p >= 0.05  (not significant / no strong evidence)\n\n');
fprintf(fid, 'p-value = the probability of seeing this result (or more extreme)\n');
fprintf(fid, 'if there were truly NO effect. A small p-value means it is\n');
fprintf(fid, 'unlikely the data would look like this by random chance.\n\n');

% --- Summary ---
fprintf(fid, '--- PLAIN SUMMARY ---\n\n');
fprintf(fid, 'Using %d observations from %d mice:\n', nObs, nMice);
if pPh < 0.05
    fprintf(fid, '- PR breakpoint CHANGES significantly across phases (p=%s).\n', fmtP(pPh));
else
    fprintf(fid, '- PR breakpoint does NOT change significantly across phases (p=%s).\n', fmtP(pPh));
end
if pGrp < 0.05
    fprintf(fid, '- Active mice have SIGNIFICANTLY different PR than Passive (p=%s).\n', fmtP(pGrp));
    grpCoef = NaN;
    try
        idx = find(contains(lme.CoefficientNames, 'Grp_Passive'));
        if ~isempty(idx), grpCoef = lme.fixedEffects(idx(1)); end
    catch
    end
    if isfinite(grpCoef)
        if grpCoef < 0
            fprintf(fid, '  Passive mice score ~%.1f LOWER than Active on average.\n', abs(grpCoef));
        else
            fprintf(fid, '  Passive mice score ~%.1f HIGHER than Active on average.\n', grpCoef);
        end
    end
else
    fprintf(fid, '- Active and Passive do NOT differ significantly (p=%s).\n', fmtP(pGrp));
end

fprintf(fid, '\n--- FIGURES PRODUCED ---\n\n');
fprintf(fid, '1. motivation_trajectory_lme.png\n');
fprintf(fid, '   Raw data: mean +/- SEM of PR breakpoint by phase and group.\n\n');
fprintf(fid, '2. motivation_LME_predicted_means.png\n');
fprintf(fid, '   Model predictions: what the LME estimates the mean PR to be\n');
fprintf(fid, '   for each Phase x Group combination, with 95%% confidence interval.\n');
fprintf(fid, '   Squares = model estimate. Circles = raw data mean.\n');
fprintf(fid, '   If they overlap closely, the model fits the data well.\n\n');
fprintf(fid, '3. motivation_observed_vs_fitted.png\n');
fprintf(fid, '   Left: actual data per mouse (thin lines) + group mean (thick).\n');
fprintf(fid, '   Right: model-fitted values per mouse + model predicted mean.\n');
fprintf(fid, '   This shows how well the model captures individual trajectories.\n\n');
fprintf(fid, '4. motivation_LME_diagnostics.png\n');
fprintf(fid, '   Left: residuals vs fitted — should show no pattern (random cloud).\n');
fprintf(fid, '   Right: QQ plot — residuals should fall on the diagonal line\n');
fprintf(fid, '   if the normality assumption is met.\n\n');
fprintf(fid, '5. motivation_delta_from_pre.png\n');
fprintf(fid, '   Change from Pre-phase baseline per mouse, by phase and group.\n\n');

fprintf(fid, '=============================================================\n');
fclose(fid);
fprintf('  Plain-language interpretation saved to: task5_1_LME_interpretation.txt\n');
end

function F = safeAnovF(aov, termName)
% Extract F-statistic from anova output for a given term.
F = NaN;
try
    if isa(aov, 'dataset'), aov = dataset2table(aov); end
    if ~istable(aov), return; end
    
    % Get term names
    terms = string.empty;
    if ismember('Term', aov.Properties.VariableNames)
        terms = string(aov.Term);
    end
    if isempty(terms) && ~isempty(aov.Properties.RowNames)
        terms = string(aov.Properties.RowNames);
    end
    if isempty(terms), return; end
    
    idx = find(strcmpi(terms, termName), 1, 'first');
    if isempty(idx)
        idx = find(contains(lower(terms), lower(string(termName))), 1, 'first');
    end
    if isempty(idx), return; end
    
    vnames = aov.Properties.VariableNames;
    fCandidates = {'FStat','F','fStat','Fstat','F_value'};
    for ci = 1:numel(fCandidates)
        match = find(strcmpi(vnames, fCandidates{ci}), 1, 'first');
        if ~isempty(match)
            colData = aov.(vnames{match});
            if isnumeric(colData) && numel(colData) >= idx
                F = double(colData(idx));
                if isfinite(F), return; end
            end
        end
    end
catch
end
end

%% ########################################################################
%  TASK 5.2 — Licking Behavior Across Phases
%  ########################################################################
function task_5_2_licking_behavior(D, outDir, statsDir)
% Analyze lick patterns by phase: lick rate, ILI, bout metrics
% Box/violin plots + LME for each metric

metrics = {'lick_freq_per_min','lick_medianIEI_s','lick_meanDur_s', ...
           'bout_n','bout_meanDur_s','bout_freq_per_min'};
labels  = {'Lick rate (/min)','Median ILI (s)','Mean lick dur (s)', ...
           'Bout count','Bout mean dur (s)','Bout freq (/min)'};

allStats = cell(0,6);
periods  = ["Pre","During","Post","Withdrawal","Re-exposure"];

for mi = 1:numel(metrics)
    ycol = metrics{mi};
    ylab = labels{mi};
    if ~ismember(ycol, D.Properties.VariableNames)
        fprintf('  Skipping %s (not found).\n', ycol); continue;
    end

    T = D(:, {'mouse_key','Period','Group', ycol});
    T = T(isfinite(T.(ycol)), :);
    T.Phase = removecats(categorical(string(T.Period)));
    T.Grp   = removecats(categorical(string(T.Group)));
    T.mouse = categorical(T.mouse_key);
    T.y     = double(T.(ycol));
    if height(T) < 10, continue; end

    % Sanity check: each mouse must map to exactly one Group
    [gm, mkc] = findgroups(T.mouse);
    bad = false(numel(mkc),1);
    for ii = 1:numel(mkc)
        ug = unique(string(T.Grp(gm==ii)));
        if numel(ug) > 1
            bad(ii) = true;
        end
    end
    if any(bad)
        badMice = string(mkc(bad));
        error('Group mapping conflict in Task 5.2: mouse appears in both groups: %s', ...
            strjoin(badMice, ', '));
    end

    % LME (try interaction; fall back to additive)
    try
        lme = fitlme(T, 'y ~ Phase * Grp + (1|mouse)');
    catch
        lme = fitlme(T, 'y ~ Phase + Grp + (1|mouse)');
    end
    aov = anova(lme);
    pPhase = extractAnovP(aov, 'Phase');
    pInter = extractAnovP(aov, 'Phase:Grp');

    allStats(end+1,:) = {ycol, 'LME_Phase',       'Phase',        pPhase, height(T), ''}; %#ok<AGROW>
    allStats(end+1,:) = {ycol, 'LME_PhaseXGroup',  'Phase:Grp',   pInter, height(T), ''}; %#ok<AGROW>

    % Per-mouse median (avoid pseudoreplication in plots and per-period stats)
    [gpm, mkp, php, grpp] = findgroups(T.mouse, T.Phase, T.Grp);
    pm = table(mkp, php, grpp, splitapply(@(x) median(x,'omitnan'), T.y, gpm), ...
        'VariableNames', {'mouse','Phase','Grp','y'});

    % Between-group per period (ranksum on per-mouse medians)
    for pi = 1:numel(periods)
        if periods(pi) == "During"
            % No Passive licking data during 6-10
            continue;
        end
        a = pm.y(pm.Grp=="Active"  & string(pm.Phase)==periods(pi));
        p = pm.y(pm.Grp=="Passive" & string(pm.Phase)==periods(pi));
        a = a(isfinite(a)); p = p(isfinite(p));
        pv = safeRanksum(a, p);
        allStats(end+1,:) = {ycol, 'Ranksum_AvP_mouseMedian', char(periods(pi)), pv, numel(a)+numel(p), ''}; %#ok<AGROW>
    end

    % --- Friedman test (within-group across phases, RM) ---
    for gi = ["Passive","Active"]
        mice_gi = unique(pm.mouse(pm.Grp==gi), 'stable');
        validPhFr = [];
        for pi = 1:numel(periods)
            if any(ismember(mice_gi, pm.mouse(pm.Grp==gi & string(pm.Phase)==periods(pi))))
                validPhFr(end+1) = pi; %#ok<AGROW>
            end
        end
        if numel(validPhFr) >= 2 && numel(mice_gi) >= 3
            frMat = nan(numel(mice_gi), numel(validPhFr));
            for ci = 1:numel(validPhFr)
                for mi = 1:numel(mice_gi)
                    v = pm.y(pm.mouse==mice_gi(mi) & string(pm.Phase)==periods(validPhFr(ci)));
                    if ~isempty(v) && isfinite(v(1)), frMat(mi,ci) = v(1); end
                end
            end
            okR = all(isfinite(frMat),2);
            frP = NaN;
            if sum(okR) >= 3
                try, frP = friedman(frMat(okR,:), 1, 'off'); catch, end
            end
            allStats(end+1,:) = {ycol, sprintf('Friedman_%s', gi), 'across_phases', frP, sum(okR), ''}; %#ok<AGROW>
        end
    end

    % --- Distribution plot (box + scatter, split by group) ---
    COL = groupColors();
    fh = figure('Color','w','Position',[60 60 1050 500]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    for gi = ["Passive","Active"]
        ax = nexttile; hold(ax,'on');
        Tg = pm(pm.Grp == gi, :);
        mice_g = unique(Tg.mouse, 'stable');
        for pi = 1:numel(periods)
            vals = Tg.y(string(Tg.Phase)==periods(pi));
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            b = boxchart(ax, repmat(pi, numel(vals),1), vals, 'BoxWidth',0.35, ...
                'MarkerStyle','none','BoxFaceAlpha',0.2);
            b.BoxFaceColor = COL.(char(gi));
            scatter(ax, repmat(pi, numel(vals),1), vals, 20, 'k', 'filled', 'MarkerFaceAlpha',0.6);
            text(ax, pi, min(vals)-0.05*range(Tg.y(isfinite(Tg.y))), sprintf('n=%d',numel(vals)), ...
                'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
        end
        set(ax,'XTick',1:numel(periods),'XTickLabel',periods);

        % Friedman test (non-parametric RM ANOVA) across phases
        % Requires: one row per mouse, one column per phase (balanced)
        friedP = NaN;
        validPh = [];
        for pi = 1:numel(periods)
            if all(ismember(mice_g, Tg.mouse(string(Tg.Phase)==periods(pi))))
                validPh(end+1) = pi; %#ok<AGROW>
            end
        end
        if numel(validPh) >= 2 && numel(mice_g) >= 3
            friedMat = nan(numel(mice_g), numel(validPh));
            for ci = 1:numel(validPh)
                for mi = 1:numel(mice_g)
                    v = Tg.y(Tg.mouse==mice_g(mi) & string(Tg.Phase)==periods(validPh(ci)));
                    if ~isempty(v) && isfinite(v(1)), friedMat(mi,ci) = v(1); end
                end
            end
            okRows = all(isfinite(friedMat),2);
            if sum(okRows) >= 3 && size(friedMat,2) >= 2
                try, friedP = friedman(friedMat(okRows,:), 1, 'off'); catch, end
            end
        end

        % Paired post-hoc: Signrank (Wilcoxon signed-rank) vs Pre, Holm-corrected
        postHocStr = '';
        prePh = find(periods == "Pre");
        if ~isempty(prePh)
            phRaw = nan(1, numel(periods));
            for pi2 = 1:numel(periods)
                if pi2 == prePh, continue; end
                % Paired: match by mouse
                paired_pre = nan(numel(mice_g),1);
                paired_oth = nan(numel(mice_g),1);
                for mi = 1:numel(mice_g)
                    vp = Tg.y(Tg.mouse==mice_g(mi) & string(Tg.Phase)==periods(prePh));
                    vo = Tg.y(Tg.mouse==mice_g(mi) & string(Tg.Phase)==periods(pi2));
                    if ~isempty(vp) && isfinite(vp(1)), paired_pre(mi) = vp(1); end
                    if ~isempty(vo) && isfinite(vo(1)), paired_oth(mi) = vo(1); end
                end
                ok_pair = isfinite(paired_pre) & isfinite(paired_oth);
                if sum(ok_pair) >= 3
                    try, phRaw(pi2) = signrank(paired_pre(ok_pair), paired_oth(ok_pair)); catch, end
                end
            end
            phAdj = holmCorrect(phRaw);
            sigPh = {};
            for pi2 = 1:numel(periods)
                if isfinite(phAdj(pi2)) && phAdj(pi2) < 0.05
                    sigPh{end+1} = sprintf('Pre vs %s p=%s%s', periods(pi2), fmtP(phAdj(pi2)), starStr(phAdj(pi2))); %#ok<AGROW>
                end
            end
            if ~isempty(sigPh), postHocStr = strjoin(sigPh, '; '); end
        end

        titleStr = sprintf('%s - %s (per-mouse median)  Friedman p=%s %s', gi, ylab, fmtP(friedP), starStr(friedP));
        if ~isempty(postHocStr)
            titleStr = sprintf('%s\n%s', titleStr, postHocStr);
        end
        title(ax, titleStr, 'FontSize', 8);
        ylabel(ax, ylab); grid(ax,'off'); box(ax,'off');
    end
    printpng(fh, fullfile(outDir, sprintf('distribution_%s.png', safeName(ycol)))); close(fh);
end

% Save stats
if ~isempty(allStats)
    ST = cell2table(allStats, 'VariableNames', {'metric','test','effect','p','N','note'});
    writetable(ST, fullfile(statsDir, 'task5_2_licking_stats.csv'));
end
fprintf('  Task 5.2 done.\n');
end

%% ########################################################################
%  TASK 5.3 — Within-Session Temporal Dynamics (5-min bins)
%  ########################################################################
function task_5_3_temporal_dynamics(csvPath, S, D, outDir, statsDir) %#ok<INUSL>
% Load raw CSV, bin lick TTL into 5-min windows, plot by phase

binEdges = [0, 300, 600, 900];   % seconds
binLabels = {'0-5 min','5-10 min','10-15 min'};
nBins = numel(binLabels);

% Load raw frame-level data
fprintf('  Loading raw CSV for time-binned analysis...\n');
T = readtable(csvPath, 'VariableNamingRule','preserve');
T = ensureString(T, 'mouse_key');
T.mouse_key = categorical(T.mouse_key);

needCols = {'mouse_key','day_index','session_idx','Lick_TTL'};
tbCol = pickTimebaseCol(T);
if isempty(tbCol)
    fprintf('  Skipping 5.3: no valid time column.\n'); return;
end

% Threshold Lick_TTL
if ismember('Lick_TTL', T.Properties.VariableNames)
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
    T.Lick_TTL = T.Lick_TTL > 0.5;
else
    fprintf('  Skipping 5.3: Lick_TTL not found.\n'); return;
end

% Group sessions
[g, km, kd, ks] = findgroups(T.mouse_key, T.day_index, T.session_idx);
nSes = max(g);
fprintf('  Computing time-binned lick rates for %d sessions...\n', nSes);

maxRows = nSes * nBins;
arr_mk   = strings(maxRows, 1);
arr_day  = nan(maxRows, 1);
arr_ses  = nan(maxRows, 1);
arr_bin  = nan(maxRows, 1);
arr_lick = nan(maxRows, 1);
arr_rate = nan(maxRows, 1);
rowIdx = 0;

for si = 1:nSes
    idx = (g == si);
    tb  = double(T.(tbCol)(idx));
    ttl = logical(T.Lick_TTL(idx));
    good = isfinite(tb);
    tb  = tb(good);  ttl = ttl(good);
    if numel(tb) < 10, continue; end

    % Shift to start at 0
    tb = tb - min(tb);
    edg = diff([false; ttl(:); false]);
    onsets = find(edg == 1);
    if isempty(onsets), onsetTimes = []; else
        onsetTimes = tb(min(onsets, numel(tb)));
    end

    for bi = 1:nBins
        inBin = onsetTimes >= binEdges(bi) & onsetTimes < binEdges(bi+1);
        nLicks = sum(inBin);
        dur_s  = binEdges(bi+1) - binEdges(bi);
        rate   = nLicks / (dur_s / 60);  % licks per minute
        rowIdx = rowIdx + 1;
        arr_mk(rowIdx)   = string(km(si));
        arr_day(rowIdx)  = double(kd(si));
        arr_ses(rowIdx)  = double(ks(si));
        arr_bin(rowIdx)  = bi;
        arr_lick(rowIdx) = nLicks;
        arr_rate(rowIdx) = rate;
    end
end

if rowIdx == 0
    fprintf('  Skipping 5.3: no valid sessions.\n'); return;
end

B = table(categorical(arr_mk(1:rowIdx)), arr_day(1:rowIdx), arr_ses(1:rowIdx), ...
    arr_bin(1:rowIdx), arr_lick(1:rowIdx), arr_rate(1:rowIdx), ...
    'VariableNames', {'mouse_key','day_index','session_idx','bin_idx','lick_count','lick_rate_per_min'});

% Assign period and group
B.Period = periodOfDay(B.day_index);
B = B(~isundefined(B.Period), :);

% Attach group from D
mouseGrp = unique(D(:,{'mouse_key','Group'}), 'rows', 'stable');
B.Group = categorical(repmat("",height(B),1), {'Active','Passive'});
for r = 1:height(mouseGrp)
    mk = mouseGrp.mouse_key(r);
    m  = B.mouse_key == mk;
    B.Group(m) = mouseGrp.Group(r);
end
B = B(~isundefined(B.Group), :);

% Collapse to day-level median per bin
[gB, mk2, di2, bi2, per2, grp2] = findgroups(B.mouse_key, B.day_index, B.bin_idx, B.Period, B.Group);
dayBin = table(removecats(mk2), double(di2), double(bi2), removecats(per2), removecats(grp2), ...
    splitapply(@(x) median(x,'omitnan'), B.lick_rate_per_min, gB), ...
    'VariableNames',{'mouse_key','day_index','bin_idx','Period','Group','rate'});

% --- Plot: lick rate by time bin, separate panels per phase ---
COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];
fh = figure('Color','w','Position',[60 60 1200 450]);
tiledlayout(1, numel(periods), 'TileSpacing','compact','Padding','compact');

for pi = 1:numel(periods)
    ax = nexttile; hold(ax,'on');
    groupRates = struct();  % store per-bin rates for stats
    for gi = ["Active","Passive"]
        Tg = dayBin(dayBin.Period==periods(pi) & dayBin.Group==gi, :);
        if isempty(Tg), continue; end
        mu = nan(1, nBins);
        se = nan(1, nBins);
        binVals = cell(1, nBins);
        for bi = 1:nBins
            v = Tg.rate(Tg.bin_idx == bi);
            v = v(isfinite(v));
            mu(bi) = mean(v);
            se(bi) = std(v) / sqrt(max(numel(v),1));
            binVals{bi} = v;
        end
        c = COL.(char(gi));
        errorbar(ax, (1:nBins)+0.05*(gi=="Passive"), mu, se, '-o', 'Color', c, ...
            'LineWidth', 1.5, 'MarkerFaceColor', c, 'CapSize', 6, 'DisplayName', char(gi));
        groupRates.(char(gi)) = binVals;
    end

    % Per-bin between-group ranksum (skip During: no Passive licking data)
    binStats = {};
    if periods(pi) ~= "During" && isfield(groupRates,'Active') && isfield(groupRates,'Passive')
        rawBinP = nan(1, nBins);
        for bi = 1:nBins
            vA = groupRates.Active{bi}; vP = groupRates.Passive{bi};
            rawBinP(bi) = safeRanksum(vA, vP);
        end
        adjBinP = holmCorrect(rawBinP);
        for bi = 1:nBins
            if isfinite(adjBinP(bi)) && adjBinP(bi) < 0.05
                binStats{end+1} = sprintf('bin%d %s', bi, starStr(adjBinP(bi))); %#ok<AGROW>
            end
        end
    end

    set(ax,'XTick',1:nBins,'XTickLabel',binLabels);
    if isempty(binStats)
        if periods(pi) == "During"
            title(ax, sprintf('%s (Active only)', char(periods(pi))));
        else
            title(ax, char(periods(pi)));
        end
    else
        title(ax, sprintf('%s\nA vs P: %s', char(periods(pi)), strjoin(binStats, ', ')), 'FontSize', 8);
    end
    if pi==1, ylabel(ax,'Lick rate (/min)'); end
    grid(ax,'off'); box(ax,'off');
    if pi == numel(periods), legend(ax,'Location','eastoutside'); end
end
printpng(fh, fullfile(outDir, 'temporal_dynamics_by_phase.png')); close(fh);

% --- Plot: early (bin1) vs late (bin3) by group across phases ---
fh = figure('Color','w','Position',[80 80 800 500]); hold on;
earlyLate = table();
for pi = 1:numel(periods)
    for gi = ["Active","Passive"]
        Tg = dayBin(dayBin.Period==periods(pi) & dayBin.Group==gi, :);
        mice = unique(Tg.mouse_key,'stable');
        for mi = 1:numel(mice)
            Tm = Tg(Tg.mouse_key==mice(mi),:);
            early = median(Tm.rate(Tm.bin_idx==1),'omitnan');
            late  = median(Tm.rate(Tm.bin_idx==nBins),'omitnan');
            earlyLate = [earlyLate; table(mice(mi), periods(pi), gi, early, late, late-early, ...
                'VariableNames',{'mouse_key','Period','Group','early','late','diff'})]; %#ok<AGROW>
        end
    end
end
if ~isempty(earlyLate)
    earlyLate.Period = categorical(earlyLate.Period, cellstr(periods), 'Ordinal', true);
    earlyLate.Group  = categorical(earlyLate.Group, {'Active','Passive'});
    plotDualBox(earlyLate, 'diff', 'Late - Early lick rate (/min)', periods, COL);
    title('Within-Session Engagement Decline (Late - Early bin)');
    writetable(earlyLate, fullfile(statsDir, 'task5_3_early_late_engagement.csv'));
end
printpng(fh, fullfile(outDir, 'early_vs_late_engagement.png')); close(fh);

fprintf('  Task 5.3 done.\n');
end

%% ########################################################################
%  TASK 6.1 — Active vs Passive (During Phase)
%  ########################################################################
function task_6_1_group_during(D, outDir, statsDir)

metrics = {'lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
           'bout_n','bout_meanDur_s','rew_freq_per_min','rew_totalDur_s','pupil_mean'};
labels  = {'Lick rate','Lick mean dur','Lick total dur','ILI median', ...
           'Bout count','Bout mean dur','Reward freq','Reward total dur','Pupil mean'};

Dd = D(string(D.Period)=="During", :);
if isempty(Dd), fprintf('  No During-phase data.\n'); return; end

allStats = cell(0,6);
COL = groupColors();

for mi = 1:numel(metrics)
    ycol = metrics{mi}; ylab = labels{mi};
    if ~ismember(ycol, Dd.Properties.VariableNames), continue; end

    A = Dd.(ycol)(Dd.Group=="Active");  A = A(isfinite(A));
    P = Dd.(ycol)(Dd.Group=="Passive"); P = P(isfinite(P));

    isLickMetric = startsWith(ycol, "lick_") || startsWith(ycol, "bout_");
    passiveEmpty = isempty(P) || all(~isfinite(P));

    % For lick metrics: still plot Active-only even if Passive has no data
    if isLickMetric && passiveEmpty
        fprintf('  During %s: Passive has no data, showing Active only.\n', ycol);
        if isempty(A), continue; end
        allStats(end+1,:) = {ycol, 'Active_only', 'During', NaN, numel(A), ...
                             sprintf('mean=%.3g (Passive N/A)', mean(A,'omitnan'))}; %#ok<AGROW>
        fh = figure('Color','w','Position',[100 100 400 500]); hold on;
        boxchart(ones(numel(A),1), A, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', COL.Active);
        scatter(ones(numel(A),1), A, 28, COL.Active, 'filled', 'MarkerFaceAlpha', 0.7);
        set(gca,'XTick',1,'XTickLabel',{sprintf('Active (n=%d)\nmean=%.2f', numel(A), mean(A))});
        ylabel(ylab);
        title(sprintf('During Phase: %s (Active only)\n(Passive has no licking data)', ylab));
        grid off; box off;
        printpng(fh, fullfile(outDir, sprintf('during_%s.png', safeName(ycol)))); close(fh);
        continue;
    end

    if isempty(A) && isempty(P), continue; end

    pv = safeRanksum(A, P);
    effectD = (mean(A,'omitnan') - mean(P,'omitnan'));
    pooledSD = sqrt(((numel(A)-1)*var(A) + (numel(P)-1)*var(P)) / max(1, numel(A)+numel(P)-2));
    cohensD = effectD / max(eps, pooledSD);
    allStats(end+1,:) = {ycol, 'Ranksum_AvP', 'During', pv, numel(A)+numel(P), ...
                         sprintf('diff=%.3g, d=%.2f', effectD, cohensD)}; %#ok<AGROW>

    fh = figure('Color','w','Position',[100 100 480 500]); hold on;
    allVals = [A(:); P(:)];
    if ~isempty(A)
        boxchart(ones(numel(A),1), A, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', COL.Active);
        scatter(ones(numel(A),1), A, 28, COL.Active, 'filled', 'MarkerFaceAlpha', 0.7);
        text(1, min(allVals)-0.08*range(allVals), sprintf('n=%d', numel(A)), ...
            'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
    end
    if ~isempty(P)
        boxchart(2*ones(numel(P),1), P, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', COL.Passive);
        scatter(2*ones(numel(P),1), P, 28, COL.Passive, 'filled', 'MarkerFaceAlpha', 0.7);
        text(2, min(allVals)-0.08*range(allVals), sprintf('n=%d', numel(P)), ...
            'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
    end
    % Significance bracket
    yMax = max(allVals); yBr = yMax * 1.1;
    plot([1 1 2 2], [yBr*0.97, yBr, yBr, yBr*0.97], 'k-', 'LineWidth', 1);
    text(1.5, yBr*1.03, sprintf('Ranksum p=%s %s', fmtP(pv), starStr(pv)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',10);

    set(gca,'XTick',[1 2],'XTickLabel', ...
        {sprintf('Active\n(%.1f +/- %.1f)', mean(A,'omitnan'), std(A,'omitnan')/sqrt(max(numel(A),1))), ...
         sprintf('Passive\n(%.1f +/- %.1f)', mean(P,'omitnan'), std(P,'omitnan')/sqrt(max(numel(P),1)))});
    ylabel(ylab);
    title(sprintf('During Phase: %s\nCohen''s d = %.2f', ylab, cohensD));
    grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('during_%s.png', safeName(ycol)))); close(fh);
end

if ~isempty(allStats)
    ST = cell2table(allStats, 'VariableNames', {'metric','test','phase','p','N','note'});
    writetable(ST, fullfile(statsDir, 'task6_1_during_comparison.csv'));
end
fprintf('  Task 6.1 done.\n');
end

%% ########################################################################
%  TASK 6.2 — Post-During Recovery
%  ########################################################################
function task_6_2_recovery(D, outDir, statsDir)
% Track recovery after yoked period. Slope of RequirementLast across days.

ycol = 'RequirementLast';
if ~ismember(ycol, D.Properties.VariableNames)
    fprintf('  Skipping 6.2: %s not found.\n', ycol); return;
end

COL = groupColors();
mice = unique(D.mouse_key, 'stable');

% Compute per-mouse slope in Post+Withdrawal+Re-exposure
recPhases = ["Post","Withdrawal","Re-exposure"];
rows = {};
for i = 1:numel(mice)
    mk = mice(i);
    Dm = D(D.mouse_key==mk & ismember(string(D.Period), recPhases), :);
    y  = double(Dm.(ycol));
    x  = double(Dm.day_index);
    good = isfinite(x) & isfinite(y);
    x = x(good); y = y(good);
    if numel(x) < 2, continue; end
    p = polyfit(x, y, 1);
    g = string(Dm.Group(1));
    rows(end+1,:) = {char(string(mk)), g, p(1), p(2), numel(x)}; %#ok<AGROW>
end

if isempty(rows)
    fprintf('  Skipping 6.2: insufficient data.\n'); return;
end

recTbl = cell2table(rows, 'VariableNames', {'mouse_key','Group','slope','intercept','N'});
recTbl.slope = double(recTbl.slope);
recTbl.Group = categorical(recTbl.Group, {'Active','Passive'});

% Compare slopes between groups
sA = recTbl.slope(recTbl.Group=="Active");
sP = recTbl.slope(recTbl.Group=="Passive");
pv = safeRanksum(sA, sP);

% Cohen's d for slope comparison
pooledSD_s = sqrt(((numel(sA)-1)*var(sA) + (numel(sP)-1)*var(sP)) / max(1, numel(sA)+numel(sP)-2));
cohenD_s = (mean(sA,'omitnan') - mean(sP,'omitnan')) / max(eps, pooledSD_s);

% Plot: recovery slope comparison with full stats
fh = figure('Color','w','Position',[100 100 480 500]); hold on;
if ~isempty(sA)
    b1 = boxchart(ones(numel(sA),1), sA, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b1.BoxFaceColor = COL.Active;
    scatter(ones(numel(sA),1), sA, 30, COL.Active, 'filled', 'MarkerFaceAlpha',0.7);
    text(1, min(sA)-0.1*range([sA(:); sP(:)]), sprintf('n=%d', numel(sA)), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
end
if ~isempty(sP)
    b2 = boxchart(2*ones(numel(sP),1), sP, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b2.BoxFaceColor = COL.Passive;
    scatter(2*ones(numel(sP),1), sP, 30, COL.Passive, 'filled', 'MarkerFaceAlpha',0.7);
    text(2, min(sP)-0.1*range([sA(:); sP(:)]), sprintf('n=%d', numel(sP)), ...
        'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
end
% Significance bracket
yMax = max([sA(:); sP(:)]);
yBr  = yMax * 1.15;
plot([1 1 2 2], [yBr*0.97, yBr, yBr, yBr*0.97], 'k-', 'LineWidth', 1);
text(1.5, yBr*1.03, sprintf('Ranksum p=%s %s', fmtP(pv), starStr(pv)), ...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
yline(0,'k:');
set(gca,'XTick',[1 2],'XTickLabel', ...
    {sprintf('Active\n(%.3f +/- %.3f)', mean(sA), std(sA)/sqrt(numel(sA))), ...
     sprintf('Passive\n(%.3f +/- %.3f)', mean(sP), std(sP)/sqrt(numel(sP)))});
ylabel('Slope (RequirementLast / day)');
title(sprintf('Recovery Slope (Post-Reexposure)\nCohen''s d = %.2f', cohenD_s));
grid off; box off;
printpng(fh, fullfile(outDir, 'recovery_slope_comparison.png')); close(fh);

% Plot: spaghetti trajectory for Post-Reexposure with per-day stats
fh = figure('Color','w','Position',[80 80 900 540]); hold on;
Drec = D(double(D.day_index)>=11, :);
for gi = ["Active","Passive"]
    Dg = Drec(Drec.Group==gi, :);
    mk = unique(Dg.mouse_key,'stable');
    c = COL.(char(gi));
    for i = 1:numel(mk)
        r = Dg.mouse_key==mk(i);
        dd = double(Dg.day_index(r)); yy = double(Dg.(ycol)(r));
        [dd,ord] = sort(dd); yy = yy(ord);
        plot(dd, yy, '-', 'Color', [c 0.35], 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    % group mean +/- SEM
    dAll = unique(double(Dg.day_index)); dAll = sort(dAll);
    mu = nan(size(dAll)); se = nan(size(dAll));
    for j = 1:numel(dAll)
        v = double(Dg.(ycol)(double(Dg.day_index)==dAll(j)));
        v = v(isfinite(v));
        mu(j) = mean(v); se(j) = std(v)/sqrt(max(numel(v),1));
    end
    errorbar(dAll, mu, se, '-o', 'Color', c, 'LineWidth', 2.5, 'MarkerFaceColor', c, ...
        'CapSize', 4, 'DisplayName', char(gi));
end

% Per-day between-group ranksum (Holm-corrected)
allDaysRec = unique(double(Drec.day_index)); allDaysRec = sort(allDaysRec);
rawPr = nan(1, numel(allDaysRec));
for di = 1:numel(allDaysRec)
    vA = double(Drec.(ycol)(Drec.Group=="Active"  & double(Drec.day_index)==allDaysRec(di)));
    vP = double(Drec.(ycol)(Drec.Group=="Passive" & double(Drec.day_index)==allDaysRec(di)));
    rawPr(di) = safeRanksum(vA, vP);
end
adjPr = holmCorrect(rawPr);
yl = ylim; topY = yl(2);
for di = 1:numel(allDaysRec)
    if isfinite(adjPr(di)) && adjPr(di) < 0.05
        text(allDaysRec(di), topY, starStr(adjPr(di)), ...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Color',[0.8 0 0]);
    end
end

xlabel('Day'); ylabel('RequirementLast');
title('Recovery Trajectory (Post, Withdrawal, Re-exposure)  stars = Holm-adj A vs P');
legend('Location','eastoutside'); grid off; box off;
xline(14,'k--','Withdrawal','LabelVerticalAlignment','top');
xline(17,'k--','Re-exposure','LabelVerticalAlignment','top');
printpng(fh, fullfile(outDir, 'recovery_trajectory.png')); close(fh);

writetable(recTbl, fullfile(statsDir, 'task6_2_recovery_slopes.csv'));
fprintf('  Task 6.2 done.\n');
end

%% ########################################################################
%  TASK 6.3 — Reward Type Effects
%  ########################################################################
function task_6_3_reward_type(D, outDir, statsDir)
% Compare phases with different reward types:
% Pre(water) vs During(morphine), Post(morphine) vs Withdrawal(water),
% Withdrawal(water) vs Re-exposure(morphine), Re-exposure vs During (sensitization?)

ycol = 'RequirementLast';
if ~ismember(ycol, D.Properties.VariableNames)
    fprintf('  Skipping 6.3: %s not found.\n', ycol); return;
end

COL = groupColors();
comparisons = {
    "Pre",        "During",       "Water->Morphine"
    "Post",       "Withdrawal",   "Morphine->Water"
    "Withdrawal", "Re-exposure",  "Water->Morphine (return)"
    "During",     "Re-exposure",  "Sensitization?"
};

allStats = cell(0,6);

for ci = 1:size(comparisons,1)
    ph1 = comparisons{ci,1}; ph2 = comparisons{ci,2}; label = comparisons{ci,3};

    for gi = ["Active","Passive"]
        Dg = D(D.Group==gi, :);
        mice = unique(Dg.mouse_key,'stable');

        med1 = nan(numel(mice),1); med2 = nan(numel(mice),1);
        for mi = 1:numel(mice)
            v1 = Dg.(ycol)(Dg.mouse_key==mice(mi) & string(Dg.Period)==ph1);
            v2 = Dg.(ycol)(Dg.mouse_key==mice(mi) & string(Dg.Period)==ph2);
            med1(mi) = median(double(v1),'omitnan');
            med2(mi) = median(double(v2),'omitnan');
        end
        good = isfinite(med1) & isfinite(med2);
        pv = NaN;
        if sum(good) >= 3
            pv = signrank(med1(good), med2(good));
        end
        allStats(end+1,:) = {ycol, sprintf('Signrank_%s', gi), label, pv, sum(good), ...
                             sprintf('%s vs %s', ph1, ph2)}; %#ok<AGROW>
    end
end

% Plot: paired bar for each comparison — BOTH groups
for gi = ["Active","Passive"]
    Dg = D(D.Group==gi, :);
    mice = unique(Dg.mouse_key,'stable');
    if isempty(mice), continue; end

    fh = figure('Color','w','Position',[60 60 1100 500]);
    tiledlayout(1, size(comparisons,1), 'TileSpacing','compact','Padding','compact');

    for ci = 1:size(comparisons,1)
        ax = nexttile; hold(ax,'on');
        ph1 = comparisons{ci,1}; ph2 = comparisons{ci,2}; label = comparisons{ci,3};

        med1 = nan(numel(mice),1); med2 = nan(numel(mice),1);
        for mi = 1:numel(mice)
            v1 = Dg.(ycol)(Dg.mouse_key==mice(mi) & string(Dg.Period)==ph1);
            v2 = Dg.(ycol)(Dg.mouse_key==mice(mi) & string(Dg.Period)==ph2);
            med1(mi) = median(double(v1),'omitnan');
            med2(mi) = median(double(v2),'omitnan');
        end
        good = isfinite(med1) & isfinite(med2);

        c = COL.(char(gi));
        for mi = find(good(:))'
            plot(ax, [1 2], [med1(mi) med2(mi)], '-', 'LineWidth', 0.6, 'Color', [c 0.3]);
        end
        scatter(ax, ones(sum(good),1), med1(good), 30, c, 'filled', 'MarkerFaceAlpha',0.7);
        scatter(ax, 2*ones(sum(good),1), med2(good), 30, c, 'filled', 'MarkerFaceAlpha',0.7);

        pv = NaN;
        if sum(good) >= 3, pv = signrank(med1(good), med2(good)); end

        % Significance bracket
        allVals = [med1(good); med2(good)];
        if isempty(allVals), continue; end
        yMax = max(allVals);
        yBr  = yMax * 1.15;
        plot(ax, [1 1 2 2], [yBr*0.97, yBr, yBr, yBr*0.97], 'k-', 'LineWidth', 1);
        text(ax, 1.5, yBr*1.04, sprintf('p=%s %s', fmtP(pv), starStr(pv)), ...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',8);
        text(ax, 1.5, min(allVals)-0.08*range(allVals), ...
            sprintf('n=%d pairs, Signrank', sum(good)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);

        set(ax,'XTick',[1 2],'XTickLabel',{char(ph1), char(ph2)});
        title(ax, sprintf('%s\np=%s %s', label, fmtP(pv), starStr(pv)), 'FontSize',9);
        ylabel(ax,'RequirementLast');
        grid(ax,'off'); box(ax,'off');
    end
    sgtitle(sprintf('Reward Type Effects (%s group, paired signed-rank test)', gi));
    printpng(fh, fullfile(outDir, sprintf('reward_type_effects_%s.png', lower(char(gi))))); close(fh);
end

if ~isempty(allStats)
    ST = cell2table(allStats, 'VariableNames',{'metric','test','comparison','p','N','note'});
    writetable(ST, fullfile(statsDir, 'task6_3_reward_type_stats.csv'));
end
fprintf('  Task 6.3 done.\n');
end

%% ########################################################################
%  TASK 7.1 — Baseline Characterization (K-means clustering)
%  ########################################################################
function task_7_1_clustering(D, outDir, statsDir)
% Cluster mice by Pre phase behavior using K-means
% Includes behavioral, pupil, AND pharmacological features

behFeats = {'RequirementLast','lick_freq_per_min','lick_medianIEI_s','bout_n', ...
            'bout_meanDur_s','rew_freq_per_min','pupil_mean'};
% Add pharmacological features if available
pharmaSearch = {'TI_Latency','TI_latency','tail_immersion','TST_Pct','TST_pct', ...
                'HOT_Pct','HP_Pct','hot_plate','STRAUB','straub','Straub_score'};
pharmaFeats = {};
for psi = 1:numel(pharmaSearch)
    col = findCol(D, pharmaSearch{psi});
    if ~isempty(col) && ~any(strcmp(pharmaFeats, col)) && ~any(strcmp(behFeats, col))
        pharmaFeats{end+1} = col; %#ok<AGROW>
    end
end
if ~isempty(pharmaFeats)
    fprintf('  7.1: Including pharmacological features: %s\n', strjoin(pharmaFeats, ', '));
end
feats = [behFeats, pharmaFeats];
feats = feats(ismember(feats, D.Properties.VariableNames));

if numel(feats) < 3
    fprintf('  Skipping 7.1: fewer than 3 features available.\n'); return;
end
fprintf('  7.1: Using %d features: %s\n', numel(feats), strjoin(feats, ', '));

Dpre = D(string(D.Period)=="Pre", :);
mice = unique(Dpre.mouse_key, 'stable');
if numel(mice) < 4
    fprintf('  Skipping 7.1: too few mice in Pre phase (%d).\n', numel(mice)); return;
end

% Mouse-level medians for Pre phase
Xraw = nan(numel(mice), numel(feats));
grp  = strings(numel(mice),1);
for i = 1:numel(mice)
    r = Dpre.mouse_key == mice(i);
    for j = 1:numel(feats)
        Xraw(i,j) = median(double(Dpre.(feats{j})(r)), 'omitnan');
    end
    g = string(Dpre.Group(find(r,1,'first')));
    grp(i) = g;
end

% Drop mice with any NaN feature
ok = all(isfinite(Xraw), 2);
Xraw = Xraw(ok,:); mice = mice(ok); grp = grp(ok);

if size(Xraw,1) < 4
    fprintf('  Skipping 7.1: too few complete mice (%d).\n', size(Xraw,1)); return;
end

% Z-score
mu = mean(Xraw,1); sd = std(Xraw,0,1); sd(sd==0) = 1;
Xz = (Xraw - mu) ./ sd;

% PCA for visualization
[~,score,~,~,expl] = pca(Xz);

COL = groupColors();

% --- Run K-means for K=2 AND K=3 ---
for K = [2 3]
    if size(Xz,1) <= K, continue; end
    rng(42);
    [cidx, ~] = kmeans(Xz, K, 'Replicates', 100, 'MaxIter', 500);
    clrMap = lines(K);

    % Plot 1: PCA colored by cluster + group labels
    fh = figure('Color','w','Position',[80 80 750 550]); hold on;
    for ki = 1:K
        m = cidx == ki;
        scatter(score(m,1), score(m,2), 60, clrMap(ki,:), 'filled', ...
            'MarkerEdgeColor','k', 'DisplayName', sprintf('Cluster %d (n=%d)', ki, sum(m)));
    end
    % Add group + mouse labels
    for i = 1:numel(mice)
        text(score(i,1)+0.1, score(i,2)+0.1, ...
            sprintf('%s (%s)', char(string(mice(i))), char(grp(i))), 'FontSize', 6);
    end
    xlabel(sprintf('PC1 (%.1f%%)', expl(1)));
    ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
    title(sprintf('Pre-Phase Clustering K=%d (K-means)', K));
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('baseline_clusters_K%d_pca.png', K))); close(fh);

    % Plot 2: cluster profiles (bar)
    fh = figure('Color','w','Position',[80 80 800 450]);
    tiledlayout(1,K,'TileSpacing','compact','Padding','compact');
    for ki = 1:K
        ax = nexttile; hold(ax,'on');
        m = cidx == ki;
        muC = mean(Xz(m,:), 1);
        bar(ax, 1:numel(feats), muC, 'FaceColor', clrMap(ki,:), 'FaceAlpha', 0.6);
        set(ax,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',30);
        % Show group composition
        nA = sum(grp(m) == "Active"); nP = sum(grp(m) == "Passive");
        title(ax, sprintf('Cluster %d (N=%d: A=%d,P=%d)', ki, sum(m), nA, nP));
        ylabel(ax,'Z-score');
        grid(ax,'off'); box(ax,'off');
    end
    printpng(fh, fullfile(outDir, sprintf('cluster_profiles_K%d.png', K))); close(fh);

    % Save cluster assignments for this K
    clTbl = table(string(mice), grp, cidx, 'VariableNames', {'mouse_key','Group','Cluster'});
    writetable(clTbl, fullfile(statsDir, sprintf('task7_1_cluster_K%d.csv', K)));
end

% === PER-PHASE PCA (separate PCA for each phase, not just Pre) ===
allPeriods = ["Pre","During","Post","Withdrawal","Re-exposure"];
allMice = unique(D.mouse_key,'stable');
grpAll = strings(numel(allMice),1);
for i = 1:numel(allMice)
    ri = find(D.mouse_key==allMice(i), 1, 'first');
    if ~isempty(ri), grpAll(i) = string(D.Group(ri)); end
end

for ppi = 1:numel(allPeriods)
    pName = allPeriods(ppi);
    Dph = D(string(D.Period)==pName, :);
    miceP = unique(Dph.mouse_key,'stable');
    if numel(miceP) < 4, continue; end

    Xph = nan(numel(miceP), numel(feats));
    grpPh = strings(numel(miceP),1);
    for i = 1:numel(miceP)
        r = Dph.mouse_key == miceP(i);
        for j = 1:numel(feats)
            Xph(i,j) = median(double(Dph.(feats{j})(r)), 'omitnan');
        end
        grpPh(i) = string(Dph.Group(find(r,1,'first')));
    end
    okPh = all(isfinite(Xph),2);
    if sum(okPh) < 4, continue; end

    XzPh = (Xph(okPh,:) - mean(Xph(okPh,:))) ./ max(std(Xph(okPh,:)),eps);
    [~,scorePh,~,~,explPh] = pca(XzPh);

    fh = figure('Color','w','Position',[80 80 700 530]); hold on;
    for gi = ["Active","Passive"]
        m = grpPh(okPh)==gi;
        if ~any(m), continue; end
        scatter(scorePh(m,1), scorePh(m,2), 60, COL.(char(gi)), 'filled', ...
            'MarkerEdgeColor','k','DisplayName',char(gi));
    end
    miceOkPh = miceP(okPh);
    for i = 1:numel(miceOkPh)
        text(scorePh(i,1)+0.1, scorePh(i,2)+0.1, char(string(miceOkPh(i))), 'FontSize',6);
    end
    xlabel(sprintf('PC1 (%.1f%%)', explPh(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explPh(2)));
    title(sprintf('%s Phase PCA (colored by Group)', pName));
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('pca_%s_phase_by_group.png', lower(char(pName))))); close(fh);
end

% === COMBINED PHASE-DEPENDENT PCA: each observation = (mouse, phase) ===
% This shows both group AND phase effects on PCA simultaneously
phaseColMap = struct('Pre',[0.3 0.7 0.3], 'During',[1.0 0.5 0.0], ...
    'Post',[0.5 0.0 0.5], 'Withdrawal',[0.0 0.5 1.0], 'Re_exposure',[0.6 0.2 0.0]);
markerMap = struct('Active','o', 'Passive','s');  % circle=Active, square=Passive

combRows = [];
combPhase = {};
combGroup = {};
combMouse = {};
nF = numel(feats);
for ppi = 1:numel(allPeriods)
    pName = allPeriods(ppi);
    for i = 1:numel(allMice)
        r = D.mouse_key==allMice(i) & string(D.Period)==pName;
        if ~any(r), continue; end
        row = nan(1, nF);
        for j = 1:nF
            row(j) = median(double(D.(feats{j})(r)), 'omitnan');
        end
        if all(isfinite(row))
            combRows(end+1,:) = row; %#ok<AGROW>
            combPhase{end+1} = char(pName); %#ok<AGROW>
            combGroup{end+1} = char(grpAll(i)); %#ok<AGROW>
            combMouse{end+1} = char(string(allMice(i))); %#ok<AGROW>
        end
    end
end

if size(combRows,1) >= 6
    XzComb = (combRows - mean(combRows)) ./ max(std(combRows),eps);
    [~,scoreComb,~,~,explComb] = pca(XzComb);

    % Plot: colored by PHASE, shaped by GROUP
    fh = figure('Color','w','Position',[80 80 850 600]); hold on;
    legendH = gobjects(0); legendL = {};
    for ppi = 1:numel(allPeriods)
        pName = allPeriods(ppi);
        pKey = strrep(char(pName),'-','_');
        if isfield(phaseColMap, pKey)
            pc = phaseColMap.(pKey);
        else
            pc = [0.5 0.5 0.5];
        end
        for gi = ["Active","Passive"]
            mKey = char(gi);
            mk = markerMap.(mKey);
            idx = strcmp(combPhase, char(pName)) & strcmp(combGroup, char(gi));
            if ~any(idx), continue; end
            h = scatter(scoreComb(idx,1), scoreComb(idx,2), 50, pc, mk, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.5);
            legendH(end+1) = h; %#ok<AGROW>
            legendL{end+1} = sprintf('%s-%s (n=%d)', pName, gi, sum(idx)); %#ok<AGROW>
        end
    end
    xlabel(sprintf('PC1 (%.1f%%)', explComb(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explComb(2)));
    title({'Phase x Group PCA (all observations)', 'circle=Active, square=Passive, color=Phase'});
    legend(legendH, legendL, 'Location','eastoutside','FontSize',7);
    grid off; box off;
    printpng(fh, fullfile(outDir, 'pca_phase_x_group_combined.png')); close(fh);

    % Plot: colored by GROUP only (to see separation) with phase labels
    fh = figure('Color','w','Position',[80 80 800 580]); hold on;
    for gi = ["Active","Passive"]
        idx = strcmp(combGroup, char(gi));
        if ~any(idx), continue; end
        scatter(scoreComb(idx,1), scoreComb(idx,2), 50, COL.(char(gi)), 'filled', ...
            'MarkerEdgeColor','k','DisplayName',char(gi));
    end
    for i = 1:size(combRows,1)
        text(scoreComb(i,1)+0.08, scoreComb(i,2)+0.08, ...
            sprintf('%s', combPhase{i}(1:min(3,end)), combMouse{i}), 'FontSize',5, 'Color',[0.4 0.4 0.4]);
    end
    xlabel(sprintf('PC1 (%.1f%%)', explComb(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explComb(2)));
    title('Combined PCA colored by Group (all phases pooled)');
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, 'pca_combined_by_group_all_phases.png')); close(fh);

    % Plot: colored by PHASE only (to see phase effect regardless of group)
    fh = figure('Color','w','Position',[80 80 800 580]); hold on;
    hPh = gobjects(0); lPh = {};
    for ppi = 1:numel(allPeriods)
        pName = allPeriods(ppi);
        pKey = strrep(char(pName),'-','_');
        if isfield(phaseColMap, pKey)
            pc = phaseColMap.(pKey);
        else
            pc = [0.5 0.5 0.5];
        end
        idx = strcmp(combPhase, char(pName));
        if ~any(idx), continue; end
        hPh(end+1) = scatter(scoreComb(idx,1), scoreComb(idx,2), 50, pc, 'filled', ...
            'MarkerEdgeColor','k'); %#ok<AGROW>
        lPh{end+1} = sprintf('%s (n=%d)', pName, sum(idx)); %#ok<AGROW>
    end
    xlabel(sprintf('PC1 (%.1f%%)', explComb(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explComb(2)));
    title('Combined PCA colored by Phase (both groups pooled)');
    legend(hPh, lPh, 'Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, 'pca_combined_by_phase_only.png')); close(fh);
end

% --- Clustering across ALL phases (concatenated features) ---
periods = ["Pre","Post","Withdrawal","Re-exposure"];
nPh = numel(periods);
Xall = nan(numel(allMice), nF * nPh);
for i = 1:numel(allMice)
    for pi = 1:nPh
        r = D.mouse_key==allMice(i) & string(D.Period)==periods(pi);
        for j = 1:nF
            if ismember(feats{j}, D.Properties.VariableNames)
                Xall(i, (pi-1)*nF + j) = median(double(D.(feats{j})(r)), 'omitnan');
            end
        end
    end
end
okAll = all(isfinite(Xall),2);
if sum(okAll) >= 4
    XzAll = (Xall(okAll,:) - mean(Xall(okAll,:))) ./ max(std(Xall(okAll,:)),eps);
    [~, scoreAll,~,~,explAll] = pca(XzAll);

    fh = figure('Color','w','Position',[80 80 750 550]); hold on;
    for gi = ["Active","Passive"]
        m = grpAll(okAll) == gi;
        if ~any(m), continue; end
        scatter(scoreAll(m,1), scoreAll(m,2), 60, COL.(char(gi)), 'filled', ...
            'MarkerEdgeColor','k', 'DisplayName', char(gi));
    end
    miceOk = allMice(okAll);
    for i = 1:numel(miceOk)
        text(scoreAll(i,1)+0.1, scoreAll(i,2)+0.1, char(string(miceOk(i))), 'FontSize', 6);
    end
    xlabel(sprintf('PC1 (%.1f%%)', explAll(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explAll(2)));
    title('All-Phase Concatenated PCA (colored by Group)');
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, 'allphase_pca_by_group.png')); close(fh);
end

fprintf('  Task 7.1 done.\n');
end

%% ########################################################################
%  TASK 7.2 — Individual Response to Morphine
%  ########################################################################
function task_7_2_morphine_response(D, outDir, statsDir)
% Morphine response: compare Pre vs Post AND Pre vs Re-exposure (NOT During,
% since Passive mice lack PR data during morphine-yoked days 6-10).
% Morphine phases: During (Active only), Post, Re-exposure
% Water phases:    Pre, Withdrawal

ycol = 'RequirementLast';
if ~ismember(ycol, D.Properties.VariableNames)
    fprintf('  Skipping 7.2: %s not found.\n', ycol); return;
end

mice = unique(D.mouse_key, 'stable');
phases = ["Pre","During","Post","Withdrawal","Re-exposure"];
phaseMed = nan(numel(mice), numel(phases));
grp = strings(numel(mice),1);

for i = 1:numel(mice)
    r = D.mouse_key == mice(i);
    grp(i) = string(D.Group(find(r,1,'first')));
    for pi = 1:numel(phases)
        phaseMed(i,pi) = median(double(D.(ycol)(r & string(D.Period)==phases(pi))), 'omitnan');
    end
end

preMed = phaseMed(:,1);   % Pre
postMed = phaseMed(:,3);   % Post (both groups have this)
reexpMed = phaseMed(:,5);  % Re-exposure (both groups have this)
COL = groupColors();

% === Plot 1: Baseline vs Morphine effect — POST - Pre (both groups) ===
morphEffect_Post = postMed - preMed;   % Primary: available for all mice
morphEffect_Reexp = reexpMed - preMed;

fh = figure('Color','w','Position',[80 80 1100 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for ei = 1:2
    if ei==1, effect = morphEffect_Post; effLabel = 'Post - Pre'; else, effect = morphEffect_Reexp; effLabel = 'Re-exposure - Pre'; end
    ax = nexttile; hold(ax,'on');
    for gi = ["Active","Passive"]
        m = grp == gi;
        if ~any(m), continue; end
        ok = isfinite(preMed(m)) & isfinite(effect(m));
        xv = preMed(m); yv = effect(m);
        scatter(ax, xv(ok), yv(ok), 60, COL.(char(gi)), 'filled', 'MarkerEdgeColor','k', ...
            'DisplayName', char(gi));
        mIdx = find(m);
        for j = find(ok(:))'
            text(ax, xv(j)+0.2, yv(j)+0.2, char(string(mice(mIdx(j)))), 'FontSize',6);
        end
    end
    yline(ax, 0, 'k:', 'No change');
    xlabel(ax, 'Baseline PR (Pre median)');
    ylabel(ax, sprintf('Morphine Effect (%s)', effLabel));
    % Between-group test
    eA = effect(grp=="Active"); eP = effect(grp=="Passive");
    pv = safeRanksum(eA(isfinite(eA)), eP(isfinite(eP)));
    title(ax, sprintf('%s (Ranksum A vs P: p=%s %s)', effLabel, fmtP(pv), starStr(pv)), 'FontSize',9);
    legend(ax, 'Location','eastoutside'); grid(ax,'off'); box(ax,'off');
end
printpng(fh, fullfile(outDir, 'morphine_response_scatter.png')); close(fh);

% === Plot 1b: Baseline PR vs morphine effect with Spearman correlation + regression line ===
fh = figure('Color','w','Position',[80 80 1100 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for ei = 1:2
    if ei==1, effect = morphEffect_Post; effLabel = 'Post - Pre'; else, effect = morphEffect_Reexp; effLabel = 'Reexp - Pre'; end
    ax = nexttile; hold(ax,'on');
    okC = isfinite(preMed) & isfinite(effect);
    % Per-group correlation
    corrLines = {};
    for gi = ["Active","Passive"]
        m = grp == gi & okC;
        if sum(m) >= 3
            [rG, pG] = corr(preMed(m), effect(m), 'type','Spearman');
            corrLines{end+1} = sprintf('%s: Spearman r=%.2f, p=%s, n=%d', gi, rG, fmtP(pG), sum(m)); %#ok<AGROW>
        end
        scatter(ax, preMed(m), effect(m), 50, COL.(char(gi)), 'filled', 'MarkerEdgeColor','k', 'DisplayName', char(gi));
    end
    % Overall correlation
    if sum(okC) >= 4
        [rAll, pAll] = corr(preMed(okC), effect(okC), 'type','Spearman');
        % Add regression line
        p = polyfit(preMed(okC), effect(okC), 1);
        xr = linspace(min(preMed(okC)), max(preMed(okC)), 50);
        plot(ax, xr, polyval(p, xr), 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
        corrLines{end+1} = sprintf('Overall: r=%.2f, p=%s, n=%d', rAll, fmtP(pAll), sum(okC)); %#ok<AGROW>
    end
    yline(ax, 0, 'k:', 'No change');
    xlabel(ax, 'Baseline PR (Pre median)');
    ylabel(ax, sprintf('Morphine Effect (%s)', effLabel));
    title(ax, sprintf('Baseline vs %s\n(Spearman, dashed=linear fit)', effLabel), 'FontSize',9);
    legend(ax, 'Location','eastoutside'); grid(ax,'off'); box(ax,'off');
    if ~isempty(corrLines), addStatsTextbox(ax, corrLines); end
end
printpng(fh, fullfile(outDir, 'baseline_vs_morphine_effect_correlation.png')); close(fh);

% === Plot 2: Delta from Pre for EACH phase (box+scatter, both groups) ===
% Skip "During" for between-group comparison (Passive has no data)
plotPhases = ["Post","Withdrawal","Re-exposure"];  % phases both groups have
plotPhaseIdx = [3, 4, 5];  % index into phaseMed columns

fh = figure('Color','w','Position',[80 80 850 500]); hold on;
off = [-0.12 0.12];
for gi_idx = 1:2
    gi = ["Active","Passive"]; gi = gi(gi_idx);
    m = grp == gi;
    c = COL.(char(gi));
    for pp = 1:numel(plotPhases)
        pi = plotPhaseIdx(pp);
        delta = phaseMed(m, pi) - preMed(m);
        delta = delta(isfinite(delta));
        if isempty(delta), continue; end
        x = pp + off(gi_idx);
        boxchart(repmat(x, numel(delta),1), delta, 'BoxWidth',0.18, ...
            'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
        scatter(repmat(x, numel(delta),1), delta, 20, c, 'filled', 'MarkerFaceAlpha',0.65);
    end
end
% Between-group brackets per phase
rawPdelta = nan(1, numel(plotPhases));
for pp = 1:numel(plotPhases)
    pi = plotPhaseIdx(pp);
    dA = phaseMed(grp=="Active", pi) - preMed(grp=="Active");
    dP = phaseMed(grp=="Passive", pi) - preMed(grp=="Passive");
    rawPdelta(pp) = safeRanksum(dA(isfinite(dA)), dP(isfinite(dP)));
end
adjPdelta = holmCorrect(rawPdelta);
yl = ylim; baseY = yl(2); pad = 0.08*diff(yl);
for pp = 1:numel(plotPhases)
    pv = adjPdelta(pp);
    if isfinite(pv) && pv < 0.05
        y = baseY + pad*0.3;
        x1 = pp+off(1); x2 = pp+off(2);
        plot([x1 x1 x2 x2], [y-pad*0.15, y, y, y-pad*0.15], 'k-', 'LineWidth',0.9);
        text(pp, y+pad*0.1, sprintf('%s p=%s', starStr(pv), fmtP(pv)), ...
            'HorizontalAlignment','center','FontSize',8);
        baseY = max(baseY, y+pad*0.4);
    end
end
ylim([yl(1), baseY+pad*0.3]);
set(gca,'XTick',1:numel(plotPhases),'XTickLabel',plotPhases);
yline(0,'k:'); ylabel('Delta from Pre (median)'); xlabel('Phase');
title({'Change from Pre Phase by Group', '(During excluded: Passive has no PR during yoked days)'});
hA = scatter(nan,nan,30,COL.Active,'filled','DisplayName','Active');
hP = scatter(nan,nan,30,COL.Passive,'filled','DisplayName','Passive');
legend([hA hP],'Location','eastoutside'); grid off; box off;
statsLines = {};
for pp = 1:numel(plotPhases)
    nA = sum(isfinite(phaseMed(grp=="Active", plotPhaseIdx(pp)) - preMed(grp=="Active")));
    nP = sum(isfinite(phaseMed(grp=="Passive", plotPhaseIdx(pp)) - preMed(grp=="Passive")));
    statsLines{end+1} = sprintf('%s: Ranksum A(n=%d) vs P(n=%d) adj.p=%s %s', ...
        plotPhases(pp), nA, nP, fmtP(adjPdelta(pp)), starStr(adjPdelta(pp))); %#ok<AGROW>
end
addStatsTextbox(gca, statsLines);
printpng(fh, fullfile(outDir, 'morphine_response_delta_by_phase.png')); close(fh);

% === Plot 3: Water vs Morphine comparison ===
% Water phases: Pre, Withdrawal  |  Morphine phases: Post, Re-exposure (both groups have)
fh = figure('Color','w','Position',[80 80 900 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gi = ["Active","Passive"]
    ax = nexttile; hold(ax,'on');
    m = grp == gi;
    c = COL.(char(gi));
    % Compute per-mouse median across water phases and morphine phases
    waterMed = mean([phaseMed(m,1), phaseMed(m,4)], 2, 'omitnan');   % Pre + Withdrawal
    morphMed = mean([phaseMed(m,3), phaseMed(m,5)], 2, 'omitnan');   % Post + Re-exposure
    okWM = isfinite(waterMed) & isfinite(morphMed);
    nPairs = sum(okWM);
    % Paired plot with connecting lines
    for mi = find(okWM(:))'
        plot(ax, [1 2], [waterMed(mi) morphMed(mi)], '-', 'Color', [c 0.3], 'LineWidth', 0.8);
    end
    scatter(ax, ones(nPairs,1), waterMed(okWM), 40, c, 'filled', 'MarkerEdgeColor','k');
    scatter(ax, 2*ones(nPairs,1), morphMed(okWM), 40, c, 'filled', 'MarkerEdgeColor','k');
    set(ax,'XTick',[1 2],'XTickLabel',{'Water (Pre+With)','Morphine (Post+Reexp)'});
    xlim(ax,[0.5 2.5]);
    pv = NaN;
    if nPairs >= 3, try, pv = signrank(waterMed(okWM), morphMed(okWM)); catch, end, end
    ylabel(ax, 'PR Breakpoint (mean of phase medians)');
    title(ax, sprintf('%s: Water vs Morphine (Signrank p=%s %s, n=%d)', gi, fmtP(pv), starStr(pv), nPairs), 'FontSize',9);
    grid(ax,'off'); box(ax,'off');
end
printpng(fh, fullfile(outDir, 'water_vs_morphine_comparison.png')); close(fh);

% === Plot 4: Individual trajectories across ALL phases ===
fh = figure('Color','w','Position',[80 80 900 480]); hold on;
for gi = ["Active","Passive"]
    m = find(grp == gi);
    c = COL.(char(gi));
    for mi = 1:numel(m)
        y = phaseMed(m(mi), :);
        ok = isfinite(y);
        plot(find(ok), y(ok), '-o', 'Color', [c 0.35], 'LineWidth', 0.8, ...
            'MarkerFaceColor', c, 'MarkerSize', 4, 'HandleVisibility','off');
    end
    % Group mean
    muG = mean(phaseMed(m,:), 1, 'omitnan');
    seG = std(phaseMed(m,:), 0, 1, 'omitnan') ./ sqrt(sum(isfinite(phaseMed(m,:)),1));
    good = isfinite(muG);
    errorbar(find(good), muG(good), seG(good), '-o', 'Color', c, 'LineWidth', 2.5, ...
        'MarkerFaceColor', c, 'CapSize', 5, 'DisplayName', char(gi));
end
set(gca,'XTick',1:numel(phases),'XTickLabel',phases);
xlabel('Phase'); ylabel('RequirementLast (per-mouse median)');
title('Individual Trajectories Across All Phases');
legend('Location','eastoutside'); grid off; box off;
printpng(fh, fullfile(outDir, 'individual_trajectories_all_phases.png')); close(fh);

% === Plot 5: Active-only During effect (since only Active has During PR) ===
duringMed = phaseMed(:,2);  % During
activeIdx = grp=="Active";
effDur = duringMed(activeIdx) - preMed(activeIdx);
effPost = postMed(activeIdx) - preMed(activeIdx);
effReexp = reexpMed(activeIdx) - preMed(activeIdx);
fh = figure('Color','w','Position',[80 80 650 500]); hold on;
c = COL.Active;
effects = {effDur, effPost, effReexp};
effLabels = {"During-Pre","Post-Pre","Reexp-Pre"};
for ei = 1:3
    e = effects{ei}; e = e(isfinite(e));
    if isempty(e), continue; end
    boxchart(repmat(ei, numel(e),1), e, 'BoxWidth',0.4, ...
        'MarkerStyle','none','BoxFaceAlpha',0.2,'BoxFaceColor',c);
    scatter(repmat(ei, numel(e),1), e, 25, c, 'filled','MarkerFaceAlpha',0.7);
end
% Paired Signrank: During vs Post, During vs Reexp
pairs = [1 2; 1 3; 2 3]; rawPE = nan(1,3);
for pp = 1:3
    e1 = effects{pairs(pp,1)}; e2 = effects{pairs(pp,2)};
    okP = isfinite(e1) & isfinite(e2);
    if sum(okP) >= 3
        try, rawPE(pp) = signrank(e1(okP), e2(okP)); catch, end
    end
end
adjPE = holmCorrect(rawPE);
yl = ylim; baseY = yl(2); padE = 0.08*diff(yl);
for pp = 1:3
    if isfinite(adjPE(pp)) && adjPE(pp) < 0.05
        y = baseY + padE*(0.3 + 0.5*(pp-1));
        x1 = pairs(pp,1); x2 = pairs(pp,2);
        plot([x1 x1 x2 x2], [y-padE*0.1, y, y, y-padE*0.1], 'k-', 'LineWidth',0.9);
        text(mean([x1 x2]), y+padE*0.1, sprintf('%s p=%s', starStr(adjPE(pp)), fmtP(adjPE(pp))), ...
            'HorizontalAlignment','center','FontSize',8);
    end
end
yline(0,'k:');
set(gca,'XTick',1:3,'XTickLabel',effLabels);
ylabel('Delta from Pre'); xlabel('Morphine Phase');
title({'Active Group: Morphine Effect Across Phases', '(Signrank pairwise, Holm-adj)'});
grid off; box off;
printpng(fh, fullfile(outDir, 'active_morphine_effect_by_phase.png')); close(fh);

% Save
resTbl = table(string(mice), grp, phaseMed, ...
    'VariableNames', {'mouse_key','Group','phase_medians'});
writetable(resTbl, fullfile(statsDir, 'task7_2_individual_morphine_response.csv'));

nResp = sum(morphEffect_Post > 0 & isfinite(morphEffect_Post));
nTotal = sum(isfinite(morphEffect_Post));
fprintf('  Responders (Post>Pre): %d / %d\n', nResp, nTotal);
fprintf('  Task 7.2 done.\n');
end

%% ########################################################################
%  TASK 10 — Pupil Analysis
%  ########################################################################
function task_10_pupil(D, csvPath, outDir, statsDir)
% 10.1: Baseline pupil across phases (LME)
% 10.3: Pupil-behavior correlations

pupilCol = pickPupilCol(D);
if isempty(pupilCol)
    fprintf('  Skipping Task 10: pupil column not found.\n'); return;
end
ycol = pupilCol;

COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% --- 10.1 LME: pupil_mean ~ Phase * Group + (1|mouse) ---
T = D(:, {'mouse_key','Period','Group', ycol});
T = T(isfinite(T.(ycol)), :);
T.Phase = removecats(categorical(string(T.Period)));
T.Grp   = removecats(categorical(string(T.Group)));
T.mouse = categorical(T.mouse_key);
T.y     = double(T.(ycol));

if height(T) >= 10
    try
        lme = fitlme(T, 'y ~ Phase * Grp + (1|mouse)');
    catch
        lme = fitlme(T, 'y ~ Phase + Grp + (1|mouse)');
    end
    aov = anova(lme);
    safeWriteLmeResults(lme, aov, statsDir, 'task10_pupil');
    fprintf('  10.1 LME fitted. Phase p=%.3g, Phase:Grp p=%.3g\n', ...
        extractAnovP(aov,'Phase'), extractAnovP(aov,'Phase:Grp'));
end

% Plot: pupil trajectory by phase and group with per-phase stats
fh = figure('Color','w','Position',[80 80 900 540]); hold on;
for gi = ["Active","Passive"]
    Tg = T(T.Grp == gi, :);
    if isempty(Tg), continue; end
    mu = nan(1,numel(periods)); se = nan(1,numel(periods));
    for pi = 1:numel(periods)
        v = Tg.y(string(Tg.Phase)==periods(pi));
        v = v(isfinite(v));
        if isempty(v), continue; end
        mu(pi) = mean(v); se(pi) = std(v)/sqrt(numel(v));
    end
    c = COL.(char(gi));
    good = isfinite(mu);
    errorbar((1:numel(periods))+0.05*(gi=="Passive"), mu, se, '-o', 'Color', c, ...
        'LineWidth', 2, 'MarkerFaceColor', c, 'CapSize', 6, 'DisplayName', char(gi));
end

% Per-phase between-group ranksum (Holm-corrected)
rawPpupil = nan(1, numel(periods));
for pi = 1:numel(periods)
    vA = T.y(T.Grp=="Active"  & string(T.Phase)==periods(pi));
    vP = T.y(T.Grp=="Passive" & string(T.Phase)==periods(pi));
    rawPpupil(pi) = safeRanksum(vA, vP);
end
adjPpupil = holmCorrect(rawPpupil);

yl = ylim; baseY = yl(2); pad = 0.06*diff(yl);
for pi = 1:numel(periods)
    if isfinite(adjPpupil(pi)) && adjPpupil(pi) < 0.05
        y = baseY + pad*0.3;
        text(pi, y, starStr(adjPpupil(pi)), 'HorizontalAlignment','center', ...
            'FontWeight','bold','FontSize',10,'Color',[0.8 0 0]);
        baseY = max(baseY, y+pad);
    end
end
ylim([yl(1), baseY+pad*0.5]);

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Phase'); ylabel('Pupil Mean (px)');

% Build title with LME + per-phase stats summary
statsStr = '';
if exist('aov','var')
    pPh = extractAnovP(aov,'Phase');
    pGrp = extractAnovP(aov,'Phase:Grp');
    statsStr = sprintf('  (LME Phase p=%s %s, Interaction p=%s %s)', ...
        fmtP(pPh), starStr(pPh), fmtP(pGrp), starStr(pGrp));
end
title(sprintf('Baseline Pupil Size Across Phases%s', statsStr));
legend('Location','eastoutside'); grid off; box off;

% Annotation box with all per-phase tests
statsLines = {};
for pi = 1:numel(periods)
    if isfinite(rawPpupil(pi))
        statsLines{end+1} = sprintf('%s: Ranksum A vs P p=%s (Holm adj=%s) %s', ...
            periods(pi), fmtP(rawPpupil(pi)), fmtP(adjPpupil(pi)), starStr(adjPpupil(pi))); %#ok<AGROW>
    end
end
if ~isempty(statsLines), addStatsTextbox(gca, statsLines); end
printpng(fh, fullfile(outDir, 'pupil_trajectory_by_phase.png')); close(fh);

% --- 10.3 Pupil-Behavior Correlations ---
behavCols = {'lick_freq_per_min','RequirementLast','bout_n','rew_freq_per_min','lick_medianIEI_s'};
behavLabs = {'Lick rate','PR breakpoint','Bout count','Reward freq','ILI'};
behavCols = behavCols(ismember(behavCols, D.Properties.VariableNames));
behavLabs = behavLabs(1:numel(behavCols));

if ~isempty(behavCols)
    corrRows = cell(0,5);
    nPlots = numel(behavCols);
    fh = figure('Color','w','Position',[40 40 350*min(nPlots,4) 380*ceil(nPlots/4)]);
    tiledlayout(ceil(nPlots/4), min(nPlots,4), 'TileSpacing','compact','Padding','compact');

    for bi = 1:numel(behavCols)
        ax = nexttile; hold(ax,'on');
        bc = behavCols{bi};
        x = double(D.(ycol));
        y = double(D.(bc));
        ok = isfinite(x) & isfinite(y);
        if sum(ok) < 5, continue; end

        grpStr = string(D.Group);
        for gi = ["Active","Passive"]
            m = ok & (grpStr == gi);
            scatter(ax, x(m), y(m), 25, COL.(char(gi)), 'filled', 'MarkerFaceAlpha', 0.6, ...
                'DisplayName', char(gi));
        end

        % Overall correlation (Spearman)
        [rho, pCorr] = corr(x(ok), y(ok), 'type','Spearman');
        corrRows(end+1,:) = {ycol, bc, rho, pCorr, sum(ok)}; %#ok<AGROW>

        % Per-group correlations
        grpCorrStr = '';
        for gcorr = ["Active","Passive"]
            m2 = ok & (grpStr == gcorr);
            if sum(m2) >= 5
                [rG, pG] = corr(x(m2), y(m2), 'type','Spearman');
                grpCorrStr = sprintf('%s %s: r=%.2f,p=%s%s', grpCorrStr, ...
                    char(gcorr), rG, fmtP(pG), starStr(pG));
            end
        end

        title(ax, sprintf('%s  (n=%d)\nr=%.2f, p=%s %s\n%s', behavLabs{bi}, sum(ok), ...
            rho, fmtP(pCorr), starStr(pCorr), strtrim(grpCorrStr)), 'FontSize', 7);
        xlabel(ax, 'Pupil mean'); ylabel(ax, behavLabs{bi});
        grid(ax,'off'); box(ax,'off');
    end
    printpng(fh, fullfile(outDir, 'pupil_behavior_correlations.png')); close(fh);

    corrTbl = cell2table(corrRows, 'VariableNames',{'pupil_var','behav_var','rho','p','N'});
    writetable(corrTbl, fullfile(statsDir, 'task10_pupil_behavior_corr.csv'));
end

% --- 10.2 Event-locked pupil (reward / lick) from raw CSV ---
try
    task_10_event_locked(csvPath, D, outDir);
catch ME
    fprintf('  [WARN] Task 10.2 event-locked failed: %s\n', ME.message);
end

fprintf('  Task 10 done.\n');
end

%% ########################################################################
%  TASK 11 — Pharmacology (TI / TST / HP / Straub)
%  ########################################################################
function task_11_pharmacology(D, outDir, statsDir)
COL = groupColors();
allStats = cell(0,6);

% --- 11.1 TI (Tail Immersion) ---
tiCol = findCol(D, 'Immersion_Latency');
if ~isempty(tiCol)
    fprintf('  11.1 TI: using column %s\n', tiCol);

    % LME: TI ~ day * group + (1|mouse)
    T = D(:, {'mouse_key','day_index','Period','Group'});
    T.y = double(D.(tiCol));
    T = T(isfinite(T.y), :);
    T.day  = double(T.day_index);
    T.Grp  = removecats(categorical(string(T.Group)));
    T.mouse = categorical(T.mouse_key);

    if height(T) >= 10
        try
            lme = fitlme(T, 'y ~ day * Grp + (1|mouse)');
        catch
            lme = fitlme(T, 'y ~ day + Grp + (1|mouse)');
        end
        aov = anova(lme);
        safeWriteLmeResults(lme, aov, statsDir, 'task11_1_TI');
        fprintf('    TI LME: day p=%.3g, Grp p=%.3g\n', extractAnovP(aov,'day'), extractAnovP(aov,'Grp'));
    end

    % Trajectory plot with per-day between-group stats
    fh = figure('Color','w','Position',[80 80 900 540]); hold on;
    for gi = ["Active","Passive"]
        Tg = T(T.Grp==gi,:);
        mk = unique(Tg.mouse,'stable');
        c = COL.(char(gi));
        for i = 1:numel(mk)
            r = Tg.mouse==mk(i);
            dd = Tg.day(r); yy = Tg.y(r);
            [dd,ord] = sort(dd); yy = yy(ord);
            plot(dd, yy, '-', 'Color', [c 0.3], 'LineWidth', 0.7, 'HandleVisibility','off');
        end
        dAll = unique(Tg.day); dAll = sort(dAll);
        mu = nan(size(dAll));
        for j = 1:numel(dAll)
            mu(j) = mean(Tg.y(Tg.day==dAll(j)), 'omitnan');
        end
        plot(dAll, mu, '-o', 'Color', c, 'LineWidth', 2.2, 'MarkerFaceColor', c, 'DisplayName', char(gi));
    end

    % Per-day between-group ranksum (Holm-corrected)
    allDays = unique(T.day); allDays = sort(allDays);
    rawPd = nan(1, numel(allDays));
    for di = 1:numel(allDays)
        vA = T.y(T.Grp=="Active"  & T.day==allDays(di));
        vP = T.y(T.Grp=="Passive" & T.day==allDays(di));
        rawPd(di) = safeRanksum(vA, vP);
    end
    adjPd = holmCorrect(rawPd);
    yl = ylim; topY = yl(2);
    for di = 1:numel(allDays)
        if isfinite(adjPd(di)) && adjPd(di) < 0.05
            text(allDays(di), topY, starStr(adjPd(di)), ...
                'HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Color',[0.8 0 0]);
        end
    end

    xlabel('Day'); ylabel('Immersion Latency (s)');
    % Include LME results in title
    lmeStr = '';
    if exist('aov','var')
        lmeStr = sprintf('  (LME day p=%s, Grp p=%s)', ...
            fmtP(extractAnovP(aov,'day')), fmtP(extractAnovP(aov,'Grp')));
    end
    title(sprintf('Tail Immersion Latency Across Days%s', lmeStr));
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, 'TI_trajectory.png')); close(fh);

    % Correlation with RequirementLast (if available)
    if ismember('RequirementLast', D.Properties.VariableNames)
        x = double(D.(tiCol)); y = double(D.RequirementLast);
        ok = isfinite(x) & isfinite(y);
        if sum(ok) >= 5
            [rho, pC] = corr(x(ok), y(ok), 'type','Spearman');
            allStats(end+1,:) = {tiCol, 'Spearman_vs_PR', 'all', pC, sum(ok), sprintf('rho=%.3f', rho)}; %#ok<AGROW>

            fh = figure('Color','w','Position',[80 80 580 500]); hold on;
            grpStr = string(D.Group);
            grpCorrLines = {};
            for gi = ["Active","Passive"]
                m = ok & (grpStr==gi);
                scatter(x(m), y(m), 30, COL.(char(gi)), 'filled', 'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', char(gi));
                if sum(m) >= 5
                    [rG, pG] = corr(x(m), y(m), 'type','Spearman');
                    grpCorrLines{end+1} = sprintf('%s: r=%.2f, p=%s %s (n=%d)', ...
                        char(gi), rG, fmtP(pG), starStr(pG), sum(m)); %#ok<AGROW>
                end
            end
            xlabel('TI Latency (s)'); ylabel('RequirementLast');
            title(sprintf('TI vs PR (Spearman r=%.2f, p=%s %s, n=%d)', rho, fmtP(pC), starStr(pC), sum(ok)));
            legend('Location','eastoutside'); grid off; box off;
            if ~isempty(grpCorrLines), addStatsTextbox(gca, grpCorrLines); end
            printpng(fh, fullfile(outDir, 'TI_vs_PR_correlation.png')); close(fh);
        end
    end
else
    fprintf('  11.1 TI: column not found, skipping.\n');
end

% --- 11.2 TST / HP (with correlations to RequirementLast) ---
tstCol = findCol(D, 'TST_Pct');
hotCol = findCol(D, 'HOT_Pct');

pharmaMetrics = {};
pharmaLabels  = {};
pharmaShort   = {};
if ~isempty(tstCol), pharmaMetrics{end+1} = tstCol; pharmaLabels{end+1} = 'TST (% nonmoving)'; pharmaShort{end+1} = 'TST'; end
if ~isempty(hotCol), pharmaMetrics{end+1} = hotCol; pharmaLabels{end+1} = 'HP (% nonmoving)'; pharmaShort{end+1} = 'HP'; end

periods = ["Pre","During","Post","Withdrawal","Re-exposure"];
for pi2 = 1:numel(pharmaMetrics)
    pc = pharmaMetrics{pi2}; pl = pharmaLabels{pi2}; ps = pharmaShort{pi2};
    T2 = D(:,{'mouse_key','Period','Group'});
    T2.y = double(D.(pc));
    T2 = T2(isfinite(T2.y), :);
    if height(T2) < 5, continue; end

    % Phase box plot
    fh = figure('Color','w','Position',[80 80 850 450]); hold on;
    plotDualBox(T2, 'y', pl, periods, COL);
    title(sprintf('%s by Phase', pl));
    printpng(fh, fullfile(outDir, sprintf('%s_by_phase.png', safeName(pc)))); close(fh);

    % Correlation with RequirementLast
    if ismember('RequirementLast', D.Properties.VariableNames)
        x = double(D.(pc)); y = double(D.RequirementLast);
        ok = isfinite(x) & isfinite(y);
        if sum(ok) >= 5
            [rho, pC] = corr(x(ok), y(ok), 'type','Spearman');
            allStats(end+1,:) = {pc, 'Spearman_vs_PR', 'all', pC, sum(ok), sprintf('rho=%.3f', rho)}; %#ok<AGROW>

            fh = figure('Color','w','Position',[80 80 580 500]); hold on;
            grpStr = string(D.Group);
            grpCorrLines = {};
            for gi = ["Active","Passive"]
                m = ok & (grpStr==gi);
                scatter(x(m), y(m), 30, COL.(char(gi)), 'filled', 'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', char(gi));
                if sum(m) >= 5
                    [rG, pG] = corr(x(m), y(m), 'type','Spearman');
                    grpCorrLines{end+1} = sprintf('%s: r=%.2f, p=%s %s (n=%d)', ...
                        char(gi), rG, fmtP(pG), starStr(pG), sum(m)); %#ok<AGROW>
                end
            end
            xlabel(pl); ylabel('RequirementLast');
            title(sprintf('%s vs PR (Spearman r=%.2f, p=%s %s, n=%d)', ps, rho, fmtP(pC), starStr(pC), sum(ok)));
            legend('Location','eastoutside'); grid off; box off;
            if ~isempty(grpCorrLines), addStatsTextbox(gca, grpCorrLines); end
            printpng(fh, fullfile(outDir, sprintf('%s_vs_PR_correlation.png', safeName(pc)))); close(fh);
        end
    end
end

% --- 11.3 Straub Tail ---
straubCol = findCol(D, 'STRAUB');
if ~isempty(straubCol)
    fprintf('  11.3 Straub: using column %s\n', straubCol);
    T3 = D(:, {'mouse_key','day_index','Period','Group'});
    T3.y = double(D.(straubCol));
    T3 = T3(isfinite(T3.y), :);

    if ~isempty(T3)
        % Trajectory plot
        fh = figure('Color','w','Position',[80 80 700 500]); hold on;
        for gi = ["Active","Passive"]
            Tg = T3(T3.Group==gi, :);
            dAll = unique(double(Tg.day_index)); dAll = sort(dAll);
            mu = nan(size(dAll)); se = nan(size(dAll));
            for j = 1:numel(dAll)
                v = Tg.y(double(Tg.day_index)==dAll(j));
                v = v(isfinite(v));
                mu(j) = mean(v); se(j) = std(v)/sqrt(max(numel(v),1));
            end
            c = COL.(char(gi));
            errorbar(dAll, mu, se, '-o', 'Color', c, 'LineWidth', 2, 'MarkerFaceColor', c, ...
                'CapSize', 4, 'DisplayName', char(gi));
        end
        % Per-day between-group ranksum (Holm-corrected)
        allDays = unique(double(T3.day_index)); allDays = sort(allDays);
        rawPs = nan(1, numel(allDays));
        for di = 1:numel(allDays)
            vA = T3.y(T3.Group=="Active"  & double(T3.day_index)==allDays(di));
            vP = T3.y(T3.Group=="Passive" & double(T3.day_index)==allDays(di));
            rawPs(di) = safeRanksum(vA, vP);
        end
        adjPs = holmCorrect(rawPs);
        yl = ylim; topY = yl(2);
        for di = 1:numel(allDays)
            if isfinite(adjPs(di)) && adjPs(di) < 0.05
                text(allDays(di), topY, starStr(adjPs(di)), ...
                    'HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Color',[0.8 0 0]);
            end
        end
        xlabel('Day'); ylabel('Straub Tail Score');
        title('Straub Tail Score Across Days (stars = Holm-adj A vs P)');
        legend('Location','eastoutside'); grid off; box off;
        printpng(fh, fullfile(outDir, 'straub_trajectory.png')); close(fh);

        % Day 5 vs Day 10 vs Day 16 comparison (key timepoints, paired Signrank)
        keyDays = [5, 10, 16];
        keyLabels = {"Day5 (Pre end)", "Day10 (During end)", "Day16 (Withdrawal end)"};
        for gi = ["Active","Passive"]
            fh = figure('Color','w','Position',[80 80 550 480]); hold on;
            c = COL.(char(gi));
            Tg3 = T3(T3.Group==gi, :);
            mice_g3 = unique(Tg3.mouse_key, 'stable');
            dayVals = cell(1, numel(keyDays));
            for ki = 1:numel(keyDays)
                v = Tg3.y(double(Tg3.day_index)==keyDays(ki));
                v = v(isfinite(v));
                dayVals{ki} = v;
                if ~isempty(v)
                    boxchart(repmat(ki, numel(v),1), v, 'BoxWidth',0.4, ...
                        'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
                    scatter(repmat(ki, numel(v),1), v, 25, c, 'filled', 'MarkerFaceAlpha',0.7);
                end
            end
            % Paired Signrank pairwise brackets (all pairs, Holm-corrected)
            pairs = [1 2; 1 3; 2 3]; rawPP = nan(1,3);
            for pp = 1:3
                d1 = keyDays(pairs(pp,1)); d2 = keyDays(pairs(pp,2));
                paired1 = nan(numel(mice_g3),1); paired2 = nan(numel(mice_g3),1);
                for mi = 1:numel(mice_g3)
                    v1 = Tg3.y(Tg3.mouse_key==mice_g3(mi) & double(Tg3.day_index)==d1);
                    v2 = Tg3.y(Tg3.mouse_key==mice_g3(mi) & double(Tg3.day_index)==d2);
                    if ~isempty(v1) && isfinite(v1(1)), paired1(mi) = v1(1); end
                    if ~isempty(v2) && isfinite(v2(1)), paired2(mi) = v2(1); end
                end
                ok_p = isfinite(paired1) & isfinite(paired2);
                if sum(ok_p) >= 3
                    try, rawPP(pp) = signrank(paired1(ok_p), paired2(ok_p)); catch, end
                end
            end
            adjPP = holmCorrect(rawPP);
            yl = ylim; baseY = yl(2); pad = 0.08*diff(yl);
            for pp = 1:3
                if isfinite(adjPP(pp)) && adjPP(pp) < 0.05
                    y = baseY + pad*(0.3 + 0.5*(pp-1));
                    x1 = pairs(pp,1); x2 = pairs(pp,2);
                    plot([x1 x1 x2 x2], [y-pad*0.1, y, y, y-pad*0.1], 'k-', 'LineWidth', 0.9);
                    text(mean([x1 x2]), y+pad*0.1, sprintf('%s p=%s (Signrank)', starStr(adjPP(pp)), fmtP(adjPP(pp))), ...
                        'HorizontalAlignment','center','FontSize',8);
                end
            end
            set(gca,'XTick',1:numel(keyDays),'XTickLabel',keyLabels);
            ylabel('Straub Tail Score');
            title(sprintf('Straub: Key Timepoints (%s) [paired Signrank, Holm-adj]', gi));
            grid off; box off;
            printpng(fh, fullfile(outDir, sprintf('straub_key_days_%s.png', lower(char(gi))))); close(fh);
        end

        % Correlation with RequirementLast
        if ismember('RequirementLast', D.Properties.VariableNames)
            x = double(D.(straubCol)); y = double(D.RequirementLast);
            ok = isfinite(x) & isfinite(y);
            if sum(ok) >= 5
                [rho, pC] = corr(x(ok), y(ok), 'type','Spearman');
                allStats(end+1,:) = {straubCol, 'Spearman_vs_PR', 'all', pC, sum(ok), sprintf('rho=%.3f', rho)}; %#ok<AGROW>
                fh = figure('Color','w','Position',[80 80 580 500]); hold on;
                grpStr = string(D.Group);
                grpCorrLines = {};
                for gi = ["Active","Passive"]
                    m = ok & (grpStr==gi);
                    scatter(x(m), y(m), 30, COL.(char(gi)), 'filled', 'MarkerFaceAlpha', 0.7, ...
                        'DisplayName', char(gi));
                    if sum(m) >= 5
                        [rG, pG] = corr(x(m), y(m), 'type','Spearman');
                        grpCorrLines{end+1} = sprintf('%s: r=%.2f, p=%s %s (n=%d)', ...
                            char(gi), rG, fmtP(pG), starStr(pG), sum(m)); %#ok<AGROW>
                    end
                end
                xlabel('Straub Tail Score'); ylabel('RequirementLast');
                title(sprintf('Straub vs PR (Spearman r=%.2f, p=%s)', rho, fmtP(pC)));
                legend('Location','eastoutside'); grid off; box off;
                if ~isempty(grpCorrLines), addStatsTextbox(gca, grpCorrLines); end
                printpng(fh, fullfile(outDir, 'straub_vs_PR_correlation.png')); close(fh);
            end
        end
    end
else
    fprintf('  11.3 Straub: column not found, skipping.\n');
end

if ~isempty(allStats)
    ST = cell2table(allStats, 'VariableNames',{'metric','test','level','p','N','note'});
    writetable(ST, fullfile(statsDir, 'task11_pharmacology_stats.csv'));
end
fprintf('  Task 11 done.\n');
end

%% ########################################################################
%  TASK 12.1 — Composite Indices (PCA)
%  ########################################################################
function task_12_1_composite(D, outDir, statsDir)
% Create composite measures via PCA:
% Motivation Index, Withdrawal Index, Morphine Response Index
% NOW includes pharmacological data (TI, TST, HP, Straub) when available

COL = groupColors();

% Behavioral + pupil features
behFeats = {'RequirementLast','lick_freq_per_min','rew_freq_per_min','lick_totalDur_s', ...
            'bout_n','bout_meanDur_s','rew_totalDur_s','pupil_mean'};
% Pharmacological features — try multiple naming conventions
pharmaSearch = {'TI_Latency','TI_latency','tail_immersion','TST_Pct','TST_pct', ...
                'HOT_Pct','HP_Pct','hot_plate','STRAUB','straub','Straub_score'};
pharmaFeats = {};
for pi = 1:numel(pharmaSearch)
    col = findCol(D, pharmaSearch{pi});
    if ~isempty(col) && ~any(strcmp(pharmaFeats, col))
        pharmaFeats{end+1} = col; %#ok<AGROW>
    end
end
if ~isempty(pharmaFeats)
    fprintf('  12.1: Including pharmacological features: %s\n', strjoin(pharmaFeats, ', '));
end

feats = [behFeats, pharmaFeats];
feats = feats(ismember(feats, D.Properties.VariableNames));

if numel(feats) < 3
    fprintf('  Skipping 12.1: fewer than 3 features.\n'); return;
end

% Mouse-level medians per period
mice = unique(D.mouse_key,'stable');
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

mousePeriod = table();
for i = 1:numel(mice)
    for pi = 1:numel(periods)
        r = D.mouse_key==mice(i) & string(D.Period)==periods(pi);
        if ~any(r), continue; end
        row = table();
        row.mouse_key = mice(i);
        row.Period = periods(pi);
        row.Group  = string(D.Group(find(r,1,'first')));
        for fi = 1:numel(feats)
            row.(feats{fi}) = median(double(D.(feats{fi})(r)),'omitnan');
        end
        mousePeriod = [mousePeriod; row]; %#ok<AGROW>
    end
end

if height(mousePeriod) < 5
    fprintf('  Skipping 12.1: insufficient data.\n'); return;
end

% Build feature matrix
X = mousePeriod{:, feats};
ok = all(isfinite(X),2);
X = X(ok,:);
mousePeriod = mousePeriod(ok,:);

if size(X,1) < 5
    fprintf('  Skipping 12.1: insufficient complete rows.\n'); return;
end

% Z-score
mu = mean(X,1); sd = std(X,0,1); sd(sd==0) = 1;
Xz = (X - mu) ./ sd;

% PCA
[coeff, score, ~, ~, expl] = pca(Xz);

% Motivation Index = PC1 score
mousePeriod.MotivationIndex = score(:,1);

% Withdrawal Index = Withdrawal - Post delta of PC1 (both groups have Post)
mice2 = unique(mousePeriod.mouse_key,'stable');
withdrawalIdx = nan(numel(mice2),1);
grp2 = strings(numel(mice2),1);
for i = 1:numel(mice2)
    r = mousePeriod.mouse_key == mice2(i);
    postV  = mousePeriod.MotivationIndex(r & mousePeriod.Period=="Post");
    withdr = mousePeriod.MotivationIndex(r & mousePeriod.Period=="Withdrawal");
    if ~isempty(postV) && ~isempty(withdr)
        withdrawalIdx(i) = mean(withdr,'omitnan') - mean(postV,'omitnan');
    end
    g = mousePeriod.Group(find(r,1,'first'));
    grp2(i) = g;
end

% --- PC1, PC2, PC3 trajectory plots ---
mousePeriod.PeriodOrd = categorical(mousePeriod.Period, cellstr(periods), 'Ordinal', true);
nPC = min(3, size(score,2));
pcNames = arrayfun(@(i) sprintf('PC%d', i), 1:nPC, 'UniformOutput', false);
for pci = 1:nPC
    mousePeriod.(pcNames{pci}) = score(:, pci);
end

for pci = 1:nPC
    pcCol = pcNames{pci};
    fh = figure('Color','w','Position',[80 80 900 480]); hold on;
    for gi = ["Active","Passive"]
        Tg = mousePeriod(mousePeriod.Group==gi, :);
        if isempty(Tg), continue; end
        mk = unique(Tg.mouse_key,'stable');
        c = COL.(char(gi));
        for mi = 1:numel(mk)
            r = Tg.mouse_key == mk(mi);
            x = double(Tg.PeriodOrd(r));
            y = Tg.(pcCol)(r);
            [x,ord] = sort(x); y = y(ord);
            plot(x, y, '-', 'Color', [c 0.3], 'LineWidth', 0.8, 'HandleVisibility','off');
        end
        mu = nan(1,numel(periods)); se = nan(1,numel(periods));
        for pi = 1:numel(periods)
            v = Tg.(pcCol)(Tg.Period==periods(pi));
            v = v(isfinite(v));
            if isempty(v), continue; end
            mu(pi) = mean(v); se(pi) = std(v)/sqrt(numel(v));
        end
        good = isfinite(mu);
        errorbar((1:numel(periods))+0.05*(gi=="Passive"), mu, se, '-o', 'Color', c, ...
            'LineWidth', 2, 'MarkerFaceColor', c, 'CapSize', 6, 'DisplayName', char(gi));
    end
    % Per-phase between-group test
    rawPpc = nan(1,numel(periods));
    for pi = 1:numel(periods)
        vA = mousePeriod.(pcCol)(mousePeriod.Group=="Active" & mousePeriod.Period==periods(pi));
        vP = mousePeriod.(pcCol)(mousePeriod.Group=="Passive" & mousePeriod.Period==periods(pi));
        rawPpc(pi) = safeRanksum(vA(isfinite(vA)), vP(isfinite(vP)));
    end
    adjPpc = holmCorrect(rawPpc);
    yl = ylim; topY = yl(2);
    for pi = 1:numel(periods)
        if isfinite(adjPpc(pi)) && adjPpc(pi) < 0.05
            text(pi, topY+0.03*diff(yl), starStr(adjPpc(pi)), ...
                'HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Color',[0.8 0 0]);
        end
    end
    set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
    xlabel('Phase'); ylabel(sprintf('%s score', pcCol));
    title(sprintf('%s Trajectory (%.1f%% var)', pcCol, expl(pci)));
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('%s_trajectory.png', lower(pcCol)))); close(fh);
end

% Plot: PCA loadings
fh = figure('Color','w','Position',[80 80 600 400]);
bar(coeff(:,1:nPC));
set(gca,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',30);
ylabel('Loading'); title('PCA Loadings (PC1-PC3)');
legend(arrayfun(@(i) sprintf('PC%d (%.0f%%)',i,expl(i)), 1:nPC, 'UniformOutput',false), ...
    'Location','best');
grid off; box off;
printpng(fh, fullfile(outDir, 'pca_loadings.png')); close(fh);

% --- Withdrawal Index (During vs Withdrawal delta of PC1) ---
% Already computed above as withdrawalIdx

% --- Morphine Response Index (Post - Pre delta of PC1, both groups have Post) ---
% Also compute Re-exposure - Pre as alternative
morphRespIdx = nan(numel(mice2),1);
reexpRespIdx = nan(numel(mice2),1);
for i = 1:numel(mice2)
    r = mousePeriod.mouse_key == mice2(i);
    postV = mousePeriod.PC1(r & mousePeriod.Period=="Post");
    pre   = mousePeriod.PC1(r & mousePeriod.Period=="Pre");
    reexp = mousePeriod.PC1(r & mousePeriod.Period=="Re-exposure");
    if ~isempty(postV) && ~isempty(pre)
        morphRespIdx(i) = mean(postV,'omitnan') - mean(pre,'omitnan');
    end
    if ~isempty(reexp) && ~isempty(pre)
        reexpRespIdx(i) = mean(reexp,'omitnan') - mean(pre,'omitnan');
    end
end

% Plot: Withdrawal Index comparison
fh = figure('Color','w','Position',[80 80 500 480]); hold on;
for gi_idx = 1:2
    gi = ["Active","Passive"]; gi = gi(gi_idx);
    m = grp2 == gi;
    v = withdrawalIdx(m); v = v(isfinite(v));
    if isempty(v), continue; end
    c = COL.(char(gi));
    boxchart(repmat(gi_idx, numel(v),1), v, 'BoxWidth',0.4, 'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
    scatter(repmat(gi_idx, numel(v),1), v, 30, c, 'filled', 'MarkerFaceAlpha',0.7);
end
pWI = safeRanksum(withdrawalIdx(grp2=="Active"), withdrawalIdx(grp2=="Passive"));
yl = ylim; yBr = yl(2)+0.05*diff(yl);
plot([1 1 2 2], [yBr*0.97 yBr yBr yBr*0.97], 'k-', 'LineWidth',1);
text(1.5, yBr*1.03, sprintf('p=%s %s', fmtP(pWI), starStr(pWI)), ...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel('Withdrawal Index (Withdrawal - Post of PC1)');
title(sprintf('Withdrawal Index (Ranksum p=%s)', fmtP(pWI)));
grid off; box off;
printpng(fh, fullfile(outDir, 'withdrawal_index_comparison.png')); close(fh);

% Plot: Morphine Response Index comparison
fh = figure('Color','w','Position',[80 80 500 480]); hold on;
for gi_idx = 1:2
    gi = ["Active","Passive"]; gi = gi(gi_idx);
    m = grp2 == gi;
    v = morphRespIdx(m); v = v(isfinite(v));
    if isempty(v), continue; end
    c = COL.(char(gi));
    boxchart(repmat(gi_idx, numel(v),1), v, 'BoxWidth',0.4, 'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
    scatter(repmat(gi_idx, numel(v),1), v, 30, c, 'filled', 'MarkerFaceAlpha',0.7);
end
pMR = safeRanksum(morphRespIdx(grp2=="Active"), morphRespIdx(grp2=="Passive"));
yl = ylim; yBr = yl(2)+0.05*diff(yl);
plot([1 1 2 2], [yBr*0.97 yBr yBr yBr*0.97], 'k-', 'LineWidth',1);
text(1.5, yBr*1.03, sprintf('p=%s %s', fmtP(pMR), starStr(pMR)), ...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel('Morphine Response Index (Post - Pre of PC1)');
title(sprintf('Morphine Response Index (Ranksum p=%s)', fmtP(pMR)));
grid off; box off;
printpng(fh, fullfile(outDir, 'morphine_response_index_comparison.png')); close(fh);

% Plot: Re-exposure Response Index comparison
fh = figure('Color','w','Position',[80 80 500 480]); hold on;
for gi_idx = 1:2
    gi = ["Active","Passive"]; gi = gi(gi_idx);
    m = grp2 == gi;
    v = reexpRespIdx(m); v = v(isfinite(v));
    if isempty(v), continue; end
    c = COL.(char(gi));
    boxchart(repmat(gi_idx, numel(v),1), v, 'BoxWidth',0.4, 'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
    scatter(repmat(gi_idx, numel(v),1), v, 30, c, 'filled', 'MarkerFaceAlpha',0.7);
end
pRE = safeRanksum(reexpRespIdx(grp2=="Active"), reexpRespIdx(grp2=="Passive"));
yl = ylim; yBr = yl(2)+0.05*diff(yl);
plot([1 1 2 2], [yBr*0.97 yBr yBr yBr*0.97], 'k-', 'LineWidth',1);
text(1.5, yBr*1.03, sprintf('p=%s %s', fmtP(pRE), starStr(pRE)), ...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel('Re-exposure Response Index (Reexp - Pre of PC1)');
title(sprintf('Re-exposure Response Index (Ranksum p=%s)', fmtP(pRE)));
grid off; box off;
printpng(fh, fullfile(outDir, 'reexposure_response_index_comparison.png')); close(fh);

% Save all
mousePeriod.PeriodOrd = [];
writetable(mousePeriod, fullfile(statsDir, 'task12_1_composite_indices.csv'));
wiTbl = table(string(mice2), grp2, withdrawalIdx, morphRespIdx, reexpRespIdx, ...
    'VariableNames', {'mouse_key','Group','WithdrawalIndex','MorphineResponseIndex','ReexpResponseIndex'});
writetable(wiTbl, fullfile(statsDir, 'task12_1_all_indices.csv'));

fprintf('  Task 12.1 done.\n');
end

%% ########################################################################
%  TASK 12.2 — Predictive Modeling (TreeBagger)
%  ########################################################################
function task_12_2_predictive(D, outDir, statsDir)
% Predict morphine response from baseline (Pre) features using Random Forest
% Run SEPARATE models for each target phase:
%   Pre vs Post      (morphine — both groups)
%   Pre vs Re-exposure (morphine — both groups, often peak effect)
%   Pre vs During    (morphine — Active only, Passive has no PR)
%   Pre vs Withdrawal (WATER control — should NOT be predictable if model is morphine-specific)

feats = {'lick_freq_per_min','lick_medianIEI_s','lick_meanDur_s', ...
         'bout_n','bout_meanDur_s','rew_freq_per_min','pupil_mean'};
% Also try adding pharmacological features
pharmaSearch = {'TI_Latency','TI_latency','TST_Pct','HOT_Pct','STRAUB'};
for pi2 = 1:numel(pharmaSearch)
    col = findCol(D, pharmaSearch{pi2});
    if ~isempty(col) && ~ismember(col, feats)
        feats{end+1} = col; %#ok<AGROW>
    end
end
feats = feats(ismember(feats, D.Properties.VariableNames));

if numel(feats) < 3
    fprintf('  Skipping 12.2: fewer than 3 features.\n'); return;
end

mice = unique(D.mouse_key,'stable');
COL = groupColors();

% Pre features (mouse-level medians)
Xpre_all = nan(numel(mice), numel(feats));
grp_all  = strings(numel(mice),1);
phaseMed = struct();
targetPhases = ["During","Post","Withdrawal","Re-exposure"];

for i = 1:numel(mice)
    r = D.mouse_key == mice(i);
    g = string(D.Group(find(r,1,'first')));
    grp_all(i) = g;
    for fi = 1:numel(feats)
        Xpre_all(i,fi) = median(double(D.(feats{fi})(r & string(D.Period)=="Pre")), 'omitnan');
    end
end

% Compute per-mouse phase medians for RequirementLast
if ~ismember('RequirementLast', D.Properties.VariableNames)
    fprintf('  Skipping 12.2: RequirementLast not found.\n'); return;
end
preMed_all = nan(numel(mice),1);
for i = 1:numel(mice)
    r = D.mouse_key==mice(i);
    preMed_all(i) = median(double(D.RequirementLast(r & string(D.Period)=="Pre")), 'omitnan');
    for tp = 1:numel(targetPhases)
        ph = targetPhases(tp);
        phKey = strrep(char(ph),'-','_');
        if ~isfield(phaseMed, phKey), phaseMed.(phKey) = nan(numel(mice),1); end
        phaseMed.(phKey)(i) = median(double(D.RequirementLast(r & string(D.Period)==ph)), 'omitnan');
    end
end

% === Run separate prediction for each target phase ===
targets = struct( ...
    'name',  {'Post-Pre', 'Reexposure-Pre', 'During-Pre', 'Withdrawal-Pre'}, ...
    'phase', {'Post',     'Re_exposure',     'During',     'Withdrawal'}, ...
    'label', {'Post - Pre (morphine)', 'Reexp - Pre (morphine, often peak)', ...
              'During - Pre (morphine, Active only)', 'Withdrawal - Pre (WATER control)'}, ...
    'short', {'post_pre', 'reexp_pre', 'during_pre', 'withdrawal_pre'}, ...
    'isControl', {false, false, false, true} ...
);

allResults = {};

for ti = 1:numel(targets)
    tgt = targets(ti);
    phKey = tgt.phase;
    yDelta = phaseMed.(phKey) - preMed_all;

    % Select mice with complete features AND target
    ok = all(isfinite(Xpre_all),2) & isfinite(yDelta);

    % During only available for Active
    if strcmp(tgt.phase, 'During')
        ok = ok & (grp_all == "Active");
    end

    Xuse = Xpre_all(ok,:);
    yUse = yDelta(ok);
    mice_use = mice(ok);
    grp_use  = grp_all(ok);

    fprintf('  12.2 [%s]: %d mice (Active=%d, Passive=%d)\n', ...
        tgt.name, numel(yUse), sum(grp_use=="Active"), sum(grp_use=="Passive"));

    if numel(yUse) < 4
        fprintf('    Too few mice, skipping.\n'); continue;
    end

    % LOO cross-validation
    nM = numel(yUse);
    yPred = nan(nM,1);
    importance = zeros(1, numel(feats));

    for i = 1:nM
        trainIdx = setdiff(1:nM, i);
        try
            mdl = TreeBagger(50, Xuse(trainIdx,:), yUse(trainIdx), ...
                'Method','regression','OOBPredictorImportance','on','MinLeafSize',1);
            yPred(i) = predict(mdl, Xuse(i,:));
            importance = importance + mdl.OOBPermutedPredictorDeltaError;
        catch
            b = [ones(numel(trainIdx),1), Xuse(trainIdx,:)] \ yUse(trainIdx);
            yPred(i) = [1, Xuse(i,:)] * b;
        end
    end
    importance = importance / nM;

    ok2 = isfinite(yPred) & isfinite(yUse);
    rho = NaN; pCorr = NaN; mse = NaN;
    if sum(ok2) >= 3
        [rho, pCorr] = corr(yUse(ok2), yPred(ok2), 'type','Spearman');
        mse = mean((yUse(ok2) - yPred(ok2)).^2);
        fprintf('    LOO: Spearman r=%.3f, p=%.3g, MSE=%.3f\n', rho, pCorr, mse);
    end

    % Store result
    allResults{end+1} = struct('name',tgt.name, 'label',tgt.label, 'short',tgt.short, ...
        'rho',rho, 'p',pCorr, 'mse',mse, 'n',numel(yUse), 'isControl',tgt.isControl); %#ok<AGROW>

    % Plot: predicted vs actual
    fh = figure('Color','w','Position',[80 80 580 500]); hold on;
    for gi = ["Active","Passive"]
        m = grp_use == gi & ok2;
        if ~any(m), continue; end
        scatter(yUse(m), yPred(m), 50, COL.(char(gi)), 'filled', 'MarkerEdgeColor','k', ...
            'DisplayName', char(gi));
    end
    if any(ok2)
        mn = min([yUse(ok2); yPred(ok2)]); mx = max([yUse(ok2); yPred(ok2)]);
        plot([mn mx], [mn mx], 'k--', 'HandleVisibility','off');
    end
    xlabel(sprintf('Actual: %s', tgt.name));
    ylabel('Predicted (LOO)');
    ctrlStr = '';
    if tgt.isControl, ctrlStr = ' [WATER CONTROL]'; end
    title(sprintf('Prediction: %s%s\nr=%.2f, p=%s, n=%d', tgt.label, ctrlStr, rho, fmtP(pCorr), numel(yUse)), 'FontSize',9);
    legend('Location','eastoutside'); grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('prediction_%s.png', tgt.short))); close(fh);

    % Plot: feature importance
    fh = figure('Color','w','Position',[80 80 600 400]);
    [impSort, impOrd] = sort(importance, 'descend');
    barh(impSort, 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.7);
    set(gca,'YTick',1:numel(feats),'YTickLabel',feats(impOrd));
    xlabel('Permutation Importance');
    title(sprintf('Feature Importance: %s', tgt.name));
    grid off; box off;
    printpng(fh, fullfile(outDir, sprintf('importance_%s.png', tgt.short))); close(fh);

    % Save per-target
    predTbl = table(string(mice_use), grp_use, yUse, yPred, ...
        'VariableNames',{'mouse_key','Group','actual','predicted'});
    writetable(predTbl, fullfile(statsDir, sprintf('task12_2_%s_predictions.csv', tgt.short)));
end

% === Summary comparison: prediction accuracy across targets ===
if numel(allResults) >= 2
    fh = figure('Color','w','Position',[80 80 700 450]); hold on;
    names = {}; rhos = []; cols = {};
    for ri = 1:numel(allResults)
        R = allResults{ri};
        names{ri} = R.name; %#ok<AGROW>
        rhos(ri) = R.rho; %#ok<AGROW>
        if R.isControl
            cols{ri} = [0.6 0.6 0.6]; %#ok<AGROW>
        else
            cols{ri} = [0.3 0.5 0.8]; %#ok<AGROW>
        end
    end
    for ri = 1:numel(rhos)
        bar(ri, rhos(ri), 'FaceColor', cols{ri}, 'FaceAlpha', 0.7, 'EdgeColor','k');
    end
    set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',20);
    ylabel('Spearman r (LOO prediction)');
    title({'Prediction Accuracy: Morphine Phases vs Water Control', ...
           'If model captures morphine effect, Water control should have low r'});
    yline(0,'k:');
    % Add p-values above bars
    for ri = 1:numel(allResults)
        R = allResults{ri};
        pStr = '';
        if isfinite(R.p)
            pStr = sprintf('p=%s %s\nn=%d', fmtP(R.p), starStr(R.p), R.n);
        end
        text(ri, max(0, rhos(ri))+0.03, pStr, 'HorizontalAlignment','center','FontSize',7);
    end
    grid off; box off;
    printpng(fh, fullfile(outDir, 'prediction_summary_all_targets.png')); close(fh);
end

% Save feature importance summary
impTbl = table(feats(:), importance(:), 'VariableNames',{'feature','importance'});
writetable(impTbl, fullfile(statsDir, 'task12_2_feature_importance.csv'));

fprintf('  Task 12.2 done.\n');
end

%% ########################################################################
%  TASK 13.1 — Within-Session State Detection (Changepoint)
%  ########################################################################
function task_13_1_state_detection(csvPath, S, D, outDir, statsDir) %#ok<INUSL>
% Detect engaged vs disengaged states within sessions using changepoint detection

binSize_s = 30;  % 30-second bins for state detection

fprintf('  Loading raw CSV for state detection...\n');
T = readtable(csvPath, 'VariableNamingRule','preserve');
T = ensureString(T, 'mouse_key');
T.mouse_key = categorical(T.mouse_key);

tbCol = pickTimebaseCol(T);
if isempty(tbCol) || ~ismember('Lick_TTL', T.Properties.VariableNames)
    fprintf('  Skipping 13.1: missing time or Lick_TTL columns.\n'); return;
end
T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
T.Lick_TTL = T.Lick_TTL > 0.5;

[g, km, kd, ks] = findgroups(T.mouse_key, T.day_index, T.session_idx);
nSes = max(g);
fprintf('  Processing %d sessions for state detection...\n', nSes);

st_mk   = strings(nSes, 1);
st_day  = nan(nSes, 1);
st_ses  = nan(nSes, 1);
st_ncp  = nan(nSes, 1);
st_pct  = nan(nSes, 1);
st_thr  = nan(nSes, 1);
stIdx = 0;

for si = 1:nSes
    idx = (g == si);
    tb  = double(T.(tbCol)(idx));
    ttl = logical(T.Lick_TTL(idx));
    good = isfinite(tb);
    tb = tb(good); ttl = ttl(good);
    if numel(tb) < 100, continue; end

    tb = tb - min(tb);
    maxT = max(tb);
    edg = 0:binSize_s:maxT;
    if numel(edg) < 3, continue; end

    % Count lick onsets per bin
    d_ttl = diff([false; ttl(:); false]);
    onsets = find(d_ttl == 1);
    if isempty(onsets), onsetTimes = []; else
        onsetTimes = tb(min(onsets, numel(tb)));
    end

    binRate = nan(numel(edg)-1, 1);
    for bi = 1:numel(edg)-1
        inBin = onsetTimes >= edg(bi) & onsetTimes < edg(bi+1);
        binRate(bi) = sum(inBin) / (binSize_s/60);  % per minute
    end

    % Changepoint detection
    try
        ipt = findchangepts(binRate, 'MaxNumChanges', 2, 'Statistic', 'mean');
    catch
        % Fallback: simple threshold
        thrFB = median(binRate,'omitnan');
        ipt = find(diff(binRate > thrFB) ~= 0);
        if numel(ipt) > 2, ipt = ipt(1:2); end
    end

    % Classify states (above/below mean)
    thr = mean(binRate, 'omitnan');
    engaged = binRate > thr;
    pctEngaged = sum(engaged) / numel(engaged) * 100;

    stIdx = stIdx + 1;
    st_mk(stIdx)  = string(km(si));
    st_day(stIdx) = double(kd(si));
    st_ses(stIdx) = double(ks(si));
    st_ncp(stIdx) = numel(ipt);
    st_pct(stIdx) = pctEngaged;
    st_thr(stIdx) = thr;
end

if stIdx == 0
    fprintf('  No sessions processed.\n'); return;
end

stateTbl = table(categorical(st_mk(1:stIdx)), st_day(1:stIdx), st_ses(1:stIdx), ...
    st_ncp(1:stIdx), st_pct(1:stIdx), st_thr(1:stIdx), ...
    'VariableNames', {'mouse_key','day_index','session_idx','n_changepoints','pct_engaged','threshold_rate'});

% Assign period and group
stateTbl.Period = periodOfDay(stateTbl.day_index);

% Attach group
mouseGrp = unique(D(:,{'mouse_key','Group'}), 'rows', 'stable');
stateTbl.Group = categorical(repmat("",height(stateTbl),1), {'Active','Passive'});
for r = 1:height(mouseGrp)
    mk = mouseGrp.mouse_key(r);
    m = stateTbl.mouse_key == mk;
    stateTbl.Group(m) = mouseGrp.Group(r);
end
stateTbl = stateTbl(~isundefined(stateTbl.Group) & ~isundefined(stateTbl.Period), :);

COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% Plot 1: % engaged by phase/group
fh = figure('Color','w','Position',[80 80 900 500]); hold on;
plotDualBox(stateTbl, 'pct_engaged', '% Time Engaged', periods, COL);
title({'Within-Session Engagement (% bins with lick rate > mean)', ...
       'Higher = more consistently engaged throughout session'});
printpng(fh, fullfile(outDir, 'engagement_by_phase.png')); close(fh);

% Plot 2: # changepoints by phase/group
fh = figure('Color','w','Position',[80 80 900 500]); hold on;
plotDualBox(stateTbl, 'n_changepoints', '# State Transitions', periods, COL);
title({'State Transitions per Session (changepoint detection on lick rate)', ...
       'More transitions = more switching between engaged/disengaged states'});
printpng(fh, fullfile(outDir, 'changepoints_by_phase.png')); close(fh);

% Plot 3: engagement trajectory by day (more intuitive)
fh = figure('Color','w','Position',[80 80 900 480]); hold on;
for gi = ["Active","Passive"]
    Tg = stateTbl(stateTbl.Group==gi, :);
    if isempty(Tg), continue; end
    dAll = unique(Tg.day_index); dAll = sort(dAll);
    mu = nan(size(dAll)); se = nan(size(dAll));
    for di = 1:numel(dAll)
        v = Tg.pct_engaged(Tg.day_index==dAll(di));
        v = v(isfinite(v));
        if ~isempty(v), mu(di) = mean(v); se(di) = std(v)/sqrt(numel(v)); end
    end
    c = COL.(char(gi));
    good = isfinite(mu);
    errorbar(dAll(good), mu(good), se(good), '-o', 'Color', c, 'LineWidth', 2, ...
        'MarkerFaceColor', c, 'CapSize', 4, 'DisplayName', char(gi));
end
xlabel('Day'); ylabel('% Time Engaged');
title({'Engagement % Across Days', ...
       '(% of 30s bins with lick rate above session mean)'});
legend('Location','eastoutside'); grid off; box off;

% Add phase background shading
yl = ylim;
phases = {[3 5],[6 10],[11 13],[14 16],[17 18]};
phColors = [0.8 0.9 1.0; 1.0 0.85 0.85; 0.85 1.0 0.85; 0.9 0.9 0.7; 1.0 0.85 1.0];
phLabels = {'Pre','During','Post','Withdr.','Re-exp.'};
for phi = 1:numel(phases)
    fill([phases{phi}(1)-0.5 phases{phi}(2)+0.5 phases{phi}(2)+0.5 phases{phi}(1)-0.5], ...
         [yl(1) yl(1) yl(2) yl(2)], phColors(phi,:), 'FaceAlpha', 0.15, 'EdgeColor','none', 'HandleVisibility','off');
    text(mean(phases{phi}), yl(2)-0.02*diff(yl), phLabels{phi}, 'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4]);
end
printpng(fh, fullfile(outDir, 'engagement_trajectory_by_day.png')); close(fh);

writetable(stateTbl, fullfile(statsDir, 'task13_1_state_detection.csv'));
fprintf('  Task 13.1 done.\n');
end

%% ########################################################################
%  TASK 13.2 — Learning Curves (Exponential Fit)
%  ########################################################################
function task_13_2_learning_curves(D, outDir, statsDir)
% Fit: y(d) = A - (A - B) * exp(-k * d)
% where A=asymptote, B=initial, k=learning rate

ycol = 'RequirementLast';
if ~ismember(ycol, D.Properties.VariableNames)
    fprintf('  Skipping 13.2: %s not found.\n', ycol); return;
end

COL = groupColors();
mice = unique(D.mouse_key, 'stable');

fitRows = {};
for i = 1:numel(mice)
    mk = mice(i);
    r = D.mouse_key == mk;
    g = string(D.Group(find(r,1,'first')));
    xRaw = double(D.day_index(r));
    yRaw = double(D.(ycol)(r));
    [x, y] = dailyMedianXY(xRaw, yRaw);
    if numel(x) < 4, continue; end

    % Robust bounded exponential fit (prevents unrealistic k outliers)
    [b, yhat, R2] = fitExpLearningRobust(x, y);
    fitRows(end+1,:) = {char(string(mk)), g, b(1), b(2), b(3), R2, numel(x)}; %#ok<AGROW>
end

if isempty(fitRows)
    fprintf('  No mice fitted.\n'); return;
end

fitTbl = cell2table(fitRows, 'VariableNames', ...
    {'mouse_key','Group','asymptote','initial','k_rate','R2','N_days'});
fitTbl.asymptote = double(fitTbl.asymptote);
fitTbl.initial   = double(fitTbl.initial);
fitTbl.k_rate    = double(fitTbl.k_rate);
fitTbl.R2        = double(fitTbl.R2);

% Compare learning rate (k) between groups
fitTbl.Group = categorical(fitTbl.Group, {'Active','Passive'});
kA = fitTbl.k_rate(fitTbl.Group=="Active"  & isfinite(fitTbl.k_rate));
kP = fitTbl.k_rate(fitTbl.Group=="Passive" & isfinite(fitTbl.k_rate));
pv = safeRanksum(kA, kP);

% Plot 1: learning rate comparison with bracket
fh = figure('Color','w','Position',[100 100 480 500]); hold on;
if ~isempty(kA)
    b1 = boxchart(ones(numel(kA),1), kA, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b1.BoxFaceColor = COL.Active;
    scatter(ones(numel(kA),1), kA, 30, COL.Active, 'filled', 'MarkerFaceAlpha',0.7);
end
if ~isempty(kP)
    b2 = boxchart(2*ones(numel(kP),1), kP, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b2.BoxFaceColor = COL.Passive;
    scatter(2*ones(numel(kP),1), kP, 30, COL.Passive, 'filled', 'MarkerFaceAlpha',0.7);
end
% Bracket
allK = [kA(:); kP(:)];
if ~isempty(allK)
    yBr = max(allK)*1.15;
    plot([1 1 2 2], [yBr*0.97 yBr yBr yBr*0.97], 'k-', 'LineWidth',1);
    text(1.5, yBr*1.03, sprintf('Ranksum p=%s %s', fmtP(pv), starStr(pv)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
end
set(gca,'XTick',[1 2],'XTickLabel', ...
    {sprintf('Active (n=%d)', numel(kA)), sprintf('Passive (n=%d)', numel(kP))});
ylabel('Learning Rate (k)');
title({'Learning Rate Comparison', ...
       '(higher k = faster approach to asymptote)'});
grid off; box off;
printpng(fh, fullfile(outDir, 'learning_rate_comparison.png')); close(fh);

% Plot 2: asymptote comparison
aA = fitTbl.asymptote(fitTbl.Group=="Active"  & isfinite(fitTbl.asymptote));
aP = fitTbl.asymptote(fitTbl.Group=="Passive" & isfinite(fitTbl.asymptote));
pvA = safeRanksum(aA, aP);
fh = figure('Color','w','Position',[100 100 480 500]); hold on;
if ~isempty(aA)
    b1 = boxchart(ones(numel(aA),1), aA, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b1.BoxFaceColor = COL.Active;
    scatter(ones(numel(aA),1), aA, 30, COL.Active, 'filled', 'MarkerFaceAlpha',0.7);
end
if ~isempty(aP)
    b2 = boxchart(2*ones(numel(aP),1), aP, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2);
    b2.BoxFaceColor = COL.Passive;
    scatter(2*ones(numel(aP),1), aP, 30, COL.Passive, 'filled', 'MarkerFaceAlpha',0.7);
end
allA = [aA(:); aP(:)];
if ~isempty(allA)
    yBr = max(allA)*1.15;
    plot([1 1 2 2], [yBr*0.97 yBr yBr yBr*0.97], 'k-', 'LineWidth',1);
    text(1.5, yBr*1.03, sprintf('Ranksum p=%s %s', fmtP(pvA), starStr(pvA)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
end
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel('Asymptote (max PR level)');
title({'Asymptote Comparison', '(maximum PR breakpoint the mouse reaches)'});
grid off; box off;
printpng(fh, fullfile(outDir, 'asymptote_comparison.png')); close(fh);

% Plot 3: fitted curves overlaid on data (separate panels per group)
model = @(b,xd) b(1) - (b(1)-b(2)) .* exp(-b(3) .* (xd - min(xd)));
fh = figure('Color','w','Position',[80 80 1100 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi = ["Active","Passive"]
    ax = nexttile; hold(ax,'on');
    c = COL.(char(gi));
    Fg = fitTbl(fitTbl.Group==gi & isfinite(fitTbl.k_rate), :);
    for fi = 1:height(Fg)
        mk = categorical(string(Fg.mouse_key(fi)));
        r = D.mouse_key == mk;
        xRaw = double(D.day_index(r)); yRaw = double(D.(ycol)(r));
        [x, y] = dailyMedianXY(xRaw, yRaw);
        if numel(x) < 4, continue; end

        plot(ax, x, y, 'o', 'Color', [c 0.4], 'MarkerSize', 4, 'HandleVisibility','off');
        xf = linspace(min(x), max(x), 50);
        b = [Fg.asymptote(fi), Fg.initial(fi), Fg.k_rate(fi)];
        if all(isfinite(b))
            plot(ax, xf, model(b, xf), '-', 'Color', [c 0.6], 'LineWidth', 1.2, 'HandleVisibility','off');
        end
    end
    xlabel(ax, 'Day'); ylabel(ax, 'RequirementLast');
    title(ax, sprintf('%s (n=%d mice)\nrobust bounded: y = A - (A-B)*exp(-k*day)', gi, height(Fg)));
    grid(ax,'off'); box(ax,'off');
end
sgtitle({'Learning Curves: Exponential Fit per Mouse', ...
         '(dots = daily medians, lines = robust bounded fit)'});
printpng(fh, fullfile(outDir, 'learning_curves_fitted.png')); close(fh);

% Plot 4: R-squared distribution
fh = figure('Color','w','Position',[100 100 480 400]); hold on;
r2A = fitTbl.R2(fitTbl.Group=="Active"  & isfinite(fitTbl.R2));
r2P = fitTbl.R2(fitTbl.Group=="Passive" & isfinite(fitTbl.R2));
if ~isempty(r2A)
    boxchart(ones(numel(r2A),1), r2A, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', COL.Active);
    scatter(ones(numel(r2A),1), r2A, 30, COL.Active, 'filled', 'MarkerFaceAlpha',0.7);
end
if ~isempty(r2P)
    boxchart(2*ones(numel(r2P),1), r2P, 'BoxWidth',0.35,'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', COL.Passive);
    scatter(2*ones(numel(r2P),1), r2P, 30, COL.Passive, 'filled', 'MarkerFaceAlpha',0.7);
end
set(gca,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
ylabel('R^2 (goodness of fit)');
title({'Model Fit Quality', '(R^2: how well the exponential curve fits each mouse)'});
grid off; box off;
printpng(fh, fullfile(outDir, 'learning_curves_R2.png')); close(fh);

writetable(fitTbl, fullfile(statsDir, 'task13_2_learning_curve_params.csv'));
fprintf('  Task 13.2 done.\n');
end

%% ########################################################################
%  HELPER FUNCTIONS
%  ########################################################################

function [xDay, yDay] = dailyMedianXY(xRaw, yRaw)
% Collapse repeated samples within day to one value/day.
xRaw = double(xRaw(:)); yRaw = double(yRaw(:));
ok = isfinite(xRaw) & isfinite(yRaw);
xRaw = xRaw(ok); yRaw = yRaw(ok);
if isempty(xRaw), xDay = []; yDay = []; return; end
[uDay, ~, g] = unique(xRaw);
yDay = splitapply(@(v) median(v,'omitnan'), yRaw, g);
xDay = uDay(:);
[xDay, ord] = sort(xDay);
yDay = yDay(ord);
end

function [b, yhat, R2] = fitExpLearningRobust(x, y)
% Robust bounded fit for y = A - (A-B)*exp(-k*(x-min(x))).
% Bounds reduce unstable, biologically implausible k estimates.
x = double(x(:)); y = double(y(:));
b = [NaN NaN NaN]; yhat = nan(size(y)); R2 = NaN;
if numel(x) < 4, return; end

x0 = x - min(x);
model = @(bb, xd) bb(1) - (bb(1)-bb(2)) .* exp(-bb(3) .* xd);

lb = [min(y)-10, min(y)-10, 0];
ub = [max(y)+10, max(y)+10, 1.5];
b0 = [prctile(y,80), prctile(y,20), 0.2];
b0 = min(max(b0, lb), ub);

fitOK = false;
if exist('lsqcurvefit','file') == 2
    try
        opts = optimset('Display','off', 'MaxIter', 500, 'TolX', 1e-6, 'TolFun', 1e-6);
        bTry = lsqcurvefit(model, b0, x0, y, lb, ub, opts);
        if all(isfinite(bTry))
            b = bTry; fitOK = true;
        end
    catch
    end
end

if ~fitOK
    % Fallback without Optimization Toolbox: sigmoid transform to enforce bounds
    toBound = @(u, lo, hi) lo + (hi-lo) ./ (1 + exp(-u));
    toUnb   = @(v, lo, hi) log((v-lo) ./ max(eps, (hi-v)));
    u0 = [toUnb(b0(1),lb(1),ub(1)), toUnb(b0(2),lb(2),ub(2)), toUnb(b0(3),lb(3),ub(3))];

    robustLoss = @(r, d) (abs(r)<=d).*0.5.*r.^2 + (abs(r)>d).*(d.*(abs(r)-0.5*d));
    delta = max(1e-6, 1.5*std(y,0,'omitnan'));
    obj = @(u) sum(robustLoss(y - model([toBound(u(1),lb(1),ub(1)), ...
                                        toBound(u(2),lb(2),ub(2)), ...
                                        toBound(u(3),lb(3),ub(3))], x0), delta), 'omitnan');
    try
        u = fminsearch(obj, u0, optimset('Display','off','MaxIter',1000,'MaxFunEvals',3000));
        b = [toBound(u(1),lb(1),ub(1)), toBound(u(2),lb(2),ub(2)), toBound(u(3),lb(3),ub(3))];
        fitOK = all(isfinite(b));
    catch
    end
end

if ~fitOK, return; end
yhat = model(b, x0);
ssRes = sum((y - yhat).^2, 'omitnan');
ssTot = sum((y - mean(y,'omitnan')).^2, 'omitnan');
if ssTot > 0
    R2 = 1 - ssRes/ssTot;
end
end

% --- Data Helpers ---

function P = periodOfDay(d)
p = strings(size(d));
p(d>=3  & d<=5 )  = "Pre";
p(d>=6  & d<=10)  = "During";
p(d>=11 & d<=13)  = "Post";
p(d>=14 & d<=16)  = "Withdrawal";
p(d>=17 & d<=18)  = "Re-exposure";
P = categorical(p, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

function R = rewardTypeOfDay(d)
r = strings(size(d));
r(d>=3  & d<=5 ) = "Water";
r(d>=6  & d<=13) = "Morphine";
r(d>=14 & d<=16) = "Water";
r(d>=17 & d<=18) = "Morphine";
R = categorical(r, ["Water","Morphine"]);
end

function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end

function col = pickTimebaseCol(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
col = '';
for i = 1:numel(cands)
    if ismember(cands{i}, T.Properties.VariableNames)
        v = double(T.(cands{i}));
        if any(isfinite(v)), col = cands{i}; return; end
    end
end
end

function nm = findCol(D, key)
nm = '';
V = D.Properties.VariableNames;
hit = find(contains(lower(V), lower(key)), 1, 'first');
if ~isempty(hit), nm = V{hit}; end
end

function col = pickPupilCol(T)
% Pick a reasonable pupil summary column from day/session-level tables.
cands = {'pupil_mean','pupil_mean_px','PupilMean','Diameter_px','pupil'};
col = pickFirstVar(T, cands);
end

function col = pickPupilRawCol(T)
% Pick a pupil column from raw frame-level table.
cands = {'Diameter_px','PupilDiameter','pupil','Pupil','pupil_px'};
col = pickFirstVar(T, cands);
end

function col = pickFirstVar(T, cands)
col = '';
for i = 1:numel(cands)
    if ismember(cands{i}, T.Properties.VariableNames)
        col = cands{i};
        return;
    end
end
end

function mk = normalizeMouseKey(mkRaw)
% Normalize mouse keys to "cage_color" when possible.
mkRaw = string(mkRaw);
mk_norm = lower(regexprep(strtrim(mkRaw), '[_\-/]+', ' '));
[cage, color] = parseMouseKeyRobust(mk_norm);
mk = lower(mkRaw);
hasBoth = cage ~= "" & color ~= "";
mk(hasBoth) = lower(cage(hasBoth) + "_" + color(hasBoth));
mk = regexprep(mk, '\s+', '');
end

% --- Cohort Roster + Robust Group Assignment (from dashboard_dec2) ---

function cohort = buildNewCohortRoster()
% Hard roster: cage + color -> Group + PairID
% 6100: black=Active, orange/red=Passive  (Pair1)
% 0911: red=Active vs orange=Passive      (Pair2)
% 0911: white=Active vs black=Passive     (Pair3)
% 0910: black=Active, orange/red=Passive  (Pair4)
% 6099: orange=Active vs red=Passive      (Pair5)
% 6099: black=Active vs white=Passive     (Pair6; white died day13)

cage  = strings(0,1);
color = strings(0,1);
group = strings(0,1);
pair  = strings(0,1);
role  = strings(0,1);
sex   = strings(0,1);

    function add(cg, col, gp, pid, rl, sx)
        cage(end+1,1)  = string(cg);
        color(end+1,1) = lower(string(col));
        group(end+1,1) = string(gp);
        pair(end+1,1)  = string(pid);
        role(end+1,1)  = string(rl);
        sex(end+1,1)   = string(sx);
    end

add("6100","black","Active","P1","Active","F");
add("6100","orange","Passive","P1","Passive","F");
add("6100","red","Passive","P1","Passive","F");

add("0911","red","Active","P2","Active","F");
add("0911","orange","Passive","P2","Passive","F");

add("0911","white","Active","P3","Active","F");
add("0911","black","Passive","P3","Passive","F");

add("0910","black","Active","P4","Active","M");
add("0910","orange","Passive","P4","Passive","M");
add("0910","red","Passive","P4","Passive","M");

add("6099","orange","Active","P5","Active","M");
add("6099","red","Passive","P5","Passive","M");

add("6099","black","Active","P6","Active","M");
add("6099","white","Passive","P6","Passive","M");

cohort = table(cage, color, group, pair, role, sex, ...
    'VariableNames', {'cage','color','Group','PairID','RoleInPair','Sex'});
end

function T = attachCohortGroupAndPair(T, cohort)
% Adds Group/PairID/RoleInPair/Sex to T using robust mouse_key parsing.
% Strategy:
%   1) Convert mouse_key to string (handles categorical, cellstr, etc.)
%   2) Parse cage (4-digit) + color using simple contains() matching
%   3) Build "cage|color" key and match against cohort roster
%   4) Unmatched rows get <undefined> Group (NOT defaulted to Active)
%
% NO fragile "_a"/"_p" override -- rely ONLY on cage+color cohort map.

% Ensure mouse_key is string (handles categorical, cellstr, double, etc.)
if iscategorical(T.mouse_key)
    mk = string(T.mouse_key);
elseif iscellstr(T.mouse_key) || isstring(T.mouse_key)
    mk = string(T.mouse_key);
elseif isnumeric(T.mouse_key)
    mk = string(T.mouse_key);
else
    mk = string(T.mouse_key);
end

% Normalize: lowercase, replace separators with spaces
mk_norm = lower(regexprep(strtrim(mk), '[_\-/\\]+', ' '));

% Parse cage + color robustly
[cg, col] = parseMouseKeyRobust(mk_norm);

% Build cage|color keys for matching
keyT = cg + "|" + col;
keyC = cohort.cage + "|" + cohort.color;
[tf, loc] = ismember(keyT, keyC);

% Assign from cohort map ONLY (no override)
grp  = strings(height(T), 1);
pair = strings(height(T), 1);
role = strings(height(T), 1);
sx   = strings(height(T), 1);

grp(tf)  = string(cohort.Group(loc(tf)));
pair(tf) = string(cohort.PairID(loc(tf)));
role(tf) = string(cohort.RoleInPair(loc(tf)));
sx(tf)   = string(cohort.Sex(loc(tf)));

% Diagnostic: print mapping results
nMapped   = sum(tf);
nUnmapped = sum(~tf);
nUniqueUnmapped = numel(unique(mk(~tf)));
fprintf('  attachCohortGroupAndPair: %d/%d rows mapped, %d unmapped (%d unique keys)\n', ...
    nMapped, height(T), nUnmapped, nUniqueUnmapped);

if nUnmapped > 0
    unmappedKeys = unique(mk(~tf));
    fprintf('  [WARN] Unmapped mouse_keys (will be excluded from analyses):\n');
    for ui = 1:min(numel(unmappedKeys), 20)
        idx1 = find(mk == unmappedKeys(ui), 1);
        fprintf('    raw="%s"  norm="%s"  cage="%s"  color="%s"  key="%s"\n', ...
            unmappedKeys(ui), mk_norm(idx1), cg(idx1), col(idx1), keyT(idx1));
    end
end

% Print unique mice and their assigned group
umk = unique(mk(tf), 'stable');
fprintf('  Mapped mice:\n');
for ui = 1:numel(umk)
    idx1 = find(mk == umk(ui) & tf, 1);
    fprintf('    "%s" -> %s (%s)\n', umk(ui), grp(idx1), pair(idx1));
end

% Apply as categoricals
T.Group      = categorical(grp, ["Active","Passive"]);
T.PairID     = categorical(pair);
T.RoleInPair = categorical(role, ["Active","Passive"]);
T.Sex        = categorical(sx, ["F","M"]);
end

function [cage, color] = parseMouseKeyRobust(mk_norm)
% Parse cage (4-digit number) and color from normalized mouse_key.
% mk_norm should be lowercased with underscores/dashes replaced by spaces.
% Uses simple 'contains' for color matching (no word boundaries needed).
%
% Handles: "6100red f s", "6100 orange f p", "6099 red m p", "0910 black m a"

n = numel(mk_norm);
cage  = strings(n, 1);
color = strings(n, 1);

colors = ["black","red","orange","white"];

for i = 1:n
    s = mk_norm(i);

    % cage: first 4-digit number
    tok = regexp(s, '(\d{4})', 'tokens', 'once');
    if ~isempty(tok), cage(i) = string(tok{1}); end

    % color: simple substring match (no word boundary needed)
    for c = colors
        if contains(s, c)
            color(i) = c;
            break;
        end
    end
end
end

% --- Plotting Helpers ---

function COL = groupColors()
COL.Active  = [0.9 0.3 0.2];
COL.Passive = [0.2 0.4 0.9];
end

function printpng(fh, fn)
set(fh,'PaperPositionMode','auto');
try, exportgraphics(fh, fn, 'Resolution', 180);
catch, print(fh, fn, '-dpng', '-r180');
end
end

function s = safeName(nm), s = regexprep(nm, '[^a-zA-Z0-9]+', '_'); end

function s = starStr(p)
if ~isfinite(p), s = 'n.s.'; return; end
if p < 0.001, s = '***';
elseif p < 0.01, s = '**';
elseif p < 0.05, s = '*';
else, s = 'n.s.';
end
end

function t = fmtP(p)
if ~isfinite(p), t = 'NaN';
elseif p < 1e-4, t = sprintf('%.2e', p);
else, t = sprintf('%.4g', p);
end
end

function p = safeRanksum(a, b)
a = a(isfinite(a)); b = b(isfinite(b));
p = NaN;
if numel(a) >= 2 && numel(b) >= 2
    try, p = ranksum(a, b); catch, end
end
end

function pAdj = holmCorrect(pRaw)
% Holm-Bonferroni step-down correction for multiple comparisons.
% Input: vector of raw p-values (can contain NaN).
% Output: adjusted p-values (NaN preserved).
pAdj = nan(size(pRaw));
valid = isfinite(pRaw);
if sum(valid) == 0, return; end
pv = pRaw(valid);
m = numel(pv);
[pSorted, idx] = sort(pv);
pHolm = nan(m,1);
for i = 1:m
    pHolm(i) = min(1, pSorted(i) * (m - i + 1));
end
% Enforce monotonicity (each adjusted p >= previous)
for i = 2:m
    pHolm(i) = max(pHolm(i), pHolm(i-1));
end
% Unsort
pOut = nan(m,1);
pOut(idx) = pHolm;
pAdj(valid) = pOut;
end

function addStatsTextbox(ax, lines)
% Add a stats annotation box below the plot (no overlap).
if isempty(lines), return; end
txt = strjoin(lines, '\n');
fig = ancestor(ax, 'figure');

% Reserve bottom space once per figure
if ~isappdata(fig, 'statsBoxAdjusted')
    axList = findall(fig, 'Type','axes');
    for k = 1:numel(axList)
        pos = get(axList(k), 'Position');
        % shrink height and move up to create bottom margin
        pos(2) = pos(2) + 0.10;
        pos(4) = pos(4) - 0.10;
        if pos(4) > 0.1
            set(axList(k), 'Position', pos);
        end
    end
    setappdata(fig, 'statsBoxAdjusted', true);
end

% Bottom-centered box (outside data area)
annotation(fig, 'textbox', [0.12 0.01 0.76 0.10], ...
    'String', txt, 'FitBoxToText','off', 'BackgroundColor','w', ...
    'EdgeColor',[0.7 0.7 0.7], 'FontSize',8, 'Interpreter','none', ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left');
end

function [varargout] = computePfromFcdf(aov, termNames)
% Compute p-values from F-stat, DF1, DF2 using fcdf. Bypasses pValue column.
% Returns one p-value per requested termName.
varargout = cell(1, numel(termNames));
for ti = 1:numel(termNames)
    varargout{ti} = NaN;
end
try
    % Convert to table if needed
    if isa(aov, 'dataset'), aov = dataset2table(aov); end
    if ~istable(aov), return; end
    
    % Get row identifiers
    if ismember('Term', aov.Properties.VariableNames)
        terms = string(aov.Term);
    elseif ~isempty(aov.Properties.RowNames)
        terms = string(aov.Properties.RowNames);
    else
        return;
    end
    
    % Find F, DF1, DF2 columns
    vn = aov.Properties.VariableNames;
    fCol = find(strcmpi(vn,'FStat') | strcmpi(vn,'F'), 1);
    d1Col = find(strcmpi(vn,'DF1'), 1);
    d2Col = find(strcmpi(vn,'DF2'), 1);
    if isempty(fCol) || isempty(d1Col) || isempty(d2Col)
        fprintf('  computePfromFcdf: cannot find F/DF1/DF2 columns in: %s\n', strjoin(vn,', '));
        return;
    end
    
    for ti = 1:numel(termNames)
        idx = find(strcmpi(terms, termNames{ti}), 1);
        if isempty(idx)
            idx = find(contains(lower(terms), lower(string(termNames{ti}))), 1);
        end
        if isempty(idx), continue; end
        
        F   = double(aov{idx, fCol});
        df1 = double(aov{idx, d1Col});
        df2 = double(aov{idx, d2Col});
        if isfinite(F) && isfinite(df1) && isfinite(df2) && df1 > 0 && df2 > 0
            varargout{ti} = 1 - fcdf(F, df1, df2);
            fprintf('  fcdf: %s F=%.3f, DF1=%.0f, DF2=%.0f => p=%.4g\n', ...
                termNames{ti}, F, df1, df2, varargout{ti});
        end
    end
catch ME
    fprintf('  computePfromFcdf error: %s\n', ME.message);
end
end

function p = extractAnovP(aov, termName)
% Robust extraction of p-value from anova() output.
% Handles all MATLAB versions: table, dataset, different column names.
% Also tries computing from F,DF1,DF2 via fcdf as last resort.
p = NaN;
try
    % Convert dataset to table if needed
    if isa(aov, 'dataset'), aov = dataset2table(aov); end
    if ~istable(aov), return; end

    % Get term names (row identifiers)
    terms = string.empty;
    if ismember('Term', aov.Properties.VariableNames)
        terms = string(aov.Term);
    end
    if isempty(terms) && ~isempty(aov.Properties.RowNames)
        terms = string(aov.Properties.RowNames);
    end
    if isempty(terms), return; end

    % Find the row matching termName (exact first, then contains)
    idx = find(strcmpi(terms, termName), 1, 'first');
    if isempty(idx)
        idx = find(contains(lower(terms), lower(string(termName))), 1, 'first');
    end
    if isempty(idx), return; end

    vnames = aov.Properties.VariableNames;

    % Strategy A: try known p-value column names
    pColCandidates = {'pValue','pvalue','p_Value','p','Prob_F','ProbF','pValueLRT'};
    for ci = 1:numel(pColCandidates)
        match = find(strcmpi(vnames, pColCandidates{ci}), 1, 'first');
        if ~isempty(match)
            val = aov{idx, match};
            if iscell(val), val = val{1}; end  % handle cell extraction
            p = double(val);
            if isfinite(p), return; end
        end
    end

    % Strategy B: try accessing by variable name directly (dot notation)
    for ci = 1:numel(pColCandidates)
        match = find(strcmpi(vnames, pColCandidates{ci}), 1, 'first');
        if ~isempty(match)
            colData = aov.(vnames{match});
            if isnumeric(colData) && numel(colData) >= idx
                p = double(colData(idx));
                if isfinite(p), return; end
            end
        end
    end

    % Strategy C: compute p from F, DF1, DF2 using fcdf
    fCol = find(strcmpi(vnames,'FStat') | strcmpi(vnames,'F'), 1);
    d1Col = find(strcmpi(vnames,'DF1'), 1);
    d2Col = find(strcmpi(vnames,'DF2'), 1);
    if ~isempty(fCol) && ~isempty(d1Col) && ~isempty(d2Col)
        F   = double(aov.(vnames{fCol})(idx));
        df1 = double(aov.(vnames{d1Col})(idx));
        df2 = double(aov.(vnames{d2Col})(idx));
        if isfinite(F) && isfinite(df1) && isfinite(df2) && df1 > 0 && df2 > 0
            p = 1 - fcdf(F, df1, df2);
            return;
        end
    end

    % Strategy D: last resort — try the last numeric column
    for ci = numel(vnames):-1:1
        colData = aov.(vnames{ci});
        if isnumeric(colData) && numel(colData) >= idx
            val = double(colData(idx));
            if isfinite(val) && val >= 0 && val <= 1
                p = val;
                return;
            end
        end
    end
catch
end
end

function pKW = safeKruskalWallis(vals, groups)
% Safe wrapper for kruskalwallis. Returns NaN on failure.
pKW = NaN;
vals = vals(isfinite(vals));
if numel(vals) < 5, return; end
try
    pKW = kruskalwallis(vals, groups, 'off');
catch
end
end

function [pPairwise, compLabels] = pairwiseRanksum(vals, groups, groupLabels)
% Pairwise ranksum tests between all group pairs.
% Returns vector of p-values and cell array of comparison labels.
nG = numel(groupLabels);
nComp = nG*(nG-1)/2;
pPairwise = nan(nComp,1);
compLabels = cell(nComp,1);
ci = 0;
for i = 1:nG
    for j = (i+1):nG
        ci = ci + 1;
        a = vals(groups == groupLabels(i)); a = a(isfinite(a));
        b = vals(groups == groupLabels(j)); b = b(isfinite(b));
        pPairwise(ci) = safeRanksum(a, b);
        compLabels{ci} = sprintf('%s vs %s', groupLabels(i), groupLabels(j));
    end
end
end

function safeWriteLmeResults(lme, aov, statsDir, prefix)
% Safely write LME coefficients and ANOVA table to CSV files.
% Handles dataset vs table discrepancies across MATLAB versions.

% --- Coefficients ---
try
    coefTbl = lme.Coefficients;
    if isa(coefTbl, 'dataset'), coefTbl = dataset2table(coefTbl); end
    if ~istable(coefTbl)
        % Manual extraction
        coefTbl = table(lme.CoefficientNames(:), lme.fixedEffects, ...
            'VariableNames', {'Name','Estimate'});
    end
    writetable(coefTbl, fullfile(statsDir, sprintf('%s_lme_coefficients.csv', prefix)));
catch ME
    fprintf('  Could not save LME coefficients: %s\n', ME.message);
end

% --- ANOVA ---
try
    if istable(aov)
        writetable(aov, fullfile(statsDir, sprintf('%s_lme_anova.csv', prefix)));
    elseif isa(aov, 'dataset')
        writetable(dataset2table(aov), fullfile(statsDir, sprintf('%s_lme_anova.csv', prefix)));
    else
        % Convert anova result to a plain table
        vn = aov.Properties.VariableNames;
        rn = aov.Properties.RowNames;
        vals = zeros(numel(rn), numel(vn));
        for vi = 1:numel(vn)
            vals(:,vi) = double(aov.(vn{vi}));
        end
        aovTbl = array2table(vals, 'VariableNames', vn);
        aovTbl.Term = rn;
        aovTbl = movevars(aovTbl, 'Term', 'Before', 1);
        writetable(aovTbl, fullfile(statsDir, sprintf('%s_lme_anova.csv', prefix)));
    end
catch ME
    fprintf('  Could not save ANOVA table: %s\n', ME.message);
end
end

function plotDeltaBoxByPhaseGroup(T, periods, COL)
% Box + scatter of T.delta by Phase and Group with between-group stats
off = [-0.12 0.12];
groups = ["Passive","Active"];
for gi = 1:numel(groups)
    g = groups(gi);
    c = COL.(char(g));
    for pi = 1:numel(periods)
        vals = T.delta(string(T.Phase)==periods(pi) & string(T.Grp)==g);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        x = pi + off(gi);
        boxchart(repmat(x, numel(vals),1), vals, 'BoxWidth',0.18, ...
            'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
        scatter(repmat(x, numel(vals),1), vals, 20, c, 'filled', 'MarkerFaceAlpha',0.65);
    end
end

% Between-group significance per period (Holm-corrected)
nPer = numel(periods);
rawP = nan(1, nPer);
for pi = 1:nPer
    A = T.delta(string(T.Phase)==periods(pi) & string(T.Grp)=="Active");
    P = T.delta(string(T.Phase)==periods(pi) & string(T.Grp)=="Passive");
    rawP(pi) = safeRanksum(A, P);
end
adjP = holmCorrect(rawP);

yl = ylim; baseY = yl(2); pad = 0.08*diff(yl);
for pi = 1:nPer
    pv = adjP(pi);
    if isfinite(pv) && pv < 0.05
        x1 = pi + off(1); x2 = pi + off(2);
        y  = baseY + pad*0.3;
        plot([x1 x1 x2 x2], [y-pad*0.15, y, y, y-pad*0.15], 'k-', 'LineWidth', 0.9);
        text(pi, y+pad*0.1, starStr(pv), 'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
        baseY = max(baseY, y+pad*0.4);
    end
end
ylim([yl(1), baseY+pad*0.3]);

set(gca,'XTick',1:nPer,'XTickLabel',periods);
yline(0,'k:');
hA = scatter(nan,nan,30,COL.Active,'filled','DisplayName','Active');
hP = scatter(nan,nan,30,COL.Passive,'filled','DisplayName','Passive');
legend([hA hP],'Location','eastoutside');
grid off; box off;

% Stats box with method + N
statsLines = {};
for pi = 1:nPer
    if isfinite(rawP(pi))
        nA = nnz(string(T.Grp)=="Active"  & string(T.Phase)==periods(pi) & isfinite(T.delta));
        nP = nnz(string(T.Grp)=="Passive" & string(T.Phase)==periods(pi) & isfinite(T.delta));
        statsLines{end+1} = sprintf('%s: Ranksum A(n=%d) vs P(n=%d) p=%s (Holm adj=%s) %s', ...
            periods(pi), nA, nP, fmtP(rawP(pi)), fmtP(adjP(pi)), starStr(adjP(pi))); %#ok<AGROW>
    end
end
if ~isempty(statsLines), addStatsTextbox(gca, statsLines); end
end

function plotDualBox(T, ycol, ylab, periods, COL)
% Box + scatter of T.(ycol) by Period and Group with full statistics:
%   - Friedman test (RM non-parametric ANOVA) across periods within each group
%   - Signrank (paired Wilcoxon) for post-hoc vs Pre within each group
%   - Ranksum Active vs Passive per period (with Holm-Bonferroni correction)
%   - Significance brackets on plot for between-group per period
%   - Stats annotation box with overall tests
off = [-0.12 0.12];
groups = ["Active","Passive"];
for gi = 1:numel(groups)
    g = groups(gi);
    c = COL.(char(g));
    mask = string(T.Group) == g;
    for pi = 1:numel(periods)
        vals = double(T.(ycol)(mask & string(T.Period)==periods(pi)));
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        x = pi + off(gi);
        boxchart(repmat(x, numel(vals),1), vals, 'BoxWidth',0.18, ...
            'MarkerStyle','none','BoxFaceAlpha',0.2, 'BoxFaceColor', c);
        scatter(repmat(x, numel(vals),1), vals, 20, c, 'filled', 'MarkerFaceAlpha',0.65);
    end
end

% --- Between-group significance per period (Holm-corrected) ---
nPer = numel(periods);
rawP = nan(1, nPer);
for pi = 1:nPer
    A = double(T.(ycol)(string(T.Group)=="Active"  & string(T.Period)==periods(pi)));
    P = double(T.(ycol)(string(T.Group)=="Passive" & string(T.Period)==periods(pi)));
    rawP(pi) = safeRanksum(A, P);
end
adjP = holmCorrect(rawP);

yl = ylim; baseY = yl(2); pad = 0.08*diff(yl);
for pi = 1:nPer
    pv = adjP(pi);
    if isfinite(pv) && pv < 0.05
        x1 = pi + off(1); x2 = pi + off(2);
        y  = baseY + pad*0.3;
        plot([x1 x1 x2 x2], [y-pad*0.15, y, y, y-pad*0.15], 'k-', 'LineWidth', 0.9);
        text(pi, y+pad*0.1, starStr(pv), 'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
        baseY = max(baseY, y+pad*0.4);
    end
end
ylim([yl(1), baseY+pad*0.3]);

set(gca,'XTick',1:nPer,'XTickLabel',periods);
ylabel(ylab);
hA = scatter(nan,nan,30,COL.Active,'filled','DisplayName','Active');
hP = scatter(nan,nan,30,COL.Passive,'filled','DisplayName','Passive');
legend([hA hP],'Location','eastoutside');
grid off; box off;

% --- Annotation box with overall statistics + post-hoc details ---
statsLines = {};

% Detect mouse column for paired (repeated-measures) tests
mouseCol = '';
if ismember('mouse_key', T.Properties.VariableNames), mouseCol = 'mouse_key';
elseif ismember('mouse', T.Properties.VariableNames), mouseCol = 'mouse';
elseif ismember('MouseID', T.Properties.VariableNames), mouseCol = 'MouseID';
end

for grpLabel = ["Active","Passive"]
    mask_g = string(T.Group) == grpLabel;
    if ~any(mask_g), continue; end

    % Aggregate per mouse per period (median) if mouse col is available
    if ~isempty(mouseCol)
        mice_g = unique(T.(mouseCol)(mask_g), 'stable');

        % Build balanced [mice x phases] matrix for Friedman
        validPh = [];
        for pi = 1:nPer
            hasData = false(numel(mice_g),1);
            for mi = 1:numel(mice_g)
                v = double(T.(ycol)(mask_g & string(T.Period)==periods(pi) & T.(mouseCol)==mice_g(mi)));
                if any(isfinite(v)), hasData(mi) = true; end
            end
            if any(hasData), validPh(end+1) = pi; end %#ok<AGROW>
        end

        % Build data matrix for Friedman (rows=mice, cols=phases)
        friedMat = nan(numel(mice_g), numel(validPh));
        for ci = 1:numel(validPh)
            for mi = 1:numel(mice_g)
                v = double(T.(ycol)(mask_g & string(T.Period)==periods(validPh(ci)) & T.(mouseCol)==mice_g(mi)));
                v = v(isfinite(v));
                if ~isempty(v), friedMat(mi,ci) = median(v,'omitnan'); end
            end
        end
        okRows = all(isfinite(friedMat),2);
        friedP = NaN;
        if sum(okRows) >= 3 && size(friedMat,2) >= 2
            try, friedP = friedman(friedMat(okRows,:), 1, 'off'); catch, end
        end

        omniLine = sprintf('%s Friedman p=%s %s (n=%d mice, %d phases)', ...
            grpLabel, fmtP(friedP), starStr(friedP), sum(okRows), numel(validPh));

        % Paired post-hoc: Signrank vs Pre (Holm-corrected)
        prePh = find(periods == "Pre");
        if friedP < 0.05 && ~isempty(prePh)
            phLines = {};
            rawPH = []; compNames = {};
            for pi2 = 1:nPer
                if pi2 == prePh, continue; end
                paired_pre = nan(numel(mice_g),1);
                paired_oth = nan(numel(mice_g),1);
                for mi = 1:numel(mice_g)
                    vp = double(T.(ycol)(mask_g & string(T.Period)==periods(prePh) & T.(mouseCol)==mice_g(mi)));
                    vp = vp(isfinite(vp));
                    vo = double(T.(ycol)(mask_g & string(T.Period)==periods(pi2) & T.(mouseCol)==mice_g(mi)));
                    vo = vo(isfinite(vo));
                    if ~isempty(vp), paired_pre(mi) = median(vp,'omitnan'); end
                    if ~isempty(vo), paired_oth(mi) = median(vo,'omitnan'); end
                end
                ok_pair = isfinite(paired_pre) & isfinite(paired_oth);
                if sum(ok_pair) >= 3
                    try, rawPH(end+1) = signrank(paired_pre(ok_pair), paired_oth(ok_pair)); %#ok<AGROW>
                    catch, rawPH(end+1) = NaN; end %#ok<AGROW>
                else
                    rawPH(end+1) = NaN; %#ok<AGROW>
                end
                compNames{end+1} = sprintf('Pre vs %s', periods(pi2)); %#ok<AGROW>
            end
            adjPH = holmCorrect(rawPH);
            for ci = 1:numel(rawPH)
                phLines{end+1} = sprintf('  %s: adj.p=%s %s (Signrank,Holm)', compNames{ci}, fmtP(adjPH(ci)), starStr(adjPH(ci))); %#ok<AGROW>
            end
            omniLine = [omniLine '; post-hoc:'];
            statsLines{end+1} = omniLine; %#ok<AGROW>
            for ci = 1:numel(phLines), statsLines{end+1} = phLines{ci}; end %#ok<AGROW>
        else
            if friedP >= 0.05
                omniLine = [omniLine ' (post-hoc not performed: omnibus n.s.)'];
            end
            statsLines{end+1} = omniLine; %#ok<AGROW>
        end
    else
        % Fallback: no mouse column found — use KW (less appropriate but safe)
        vG = []; gGidx = [];
        for pi = 1:nPer
            vals = double(T.(ycol)(mask_g & string(T.Period)==periods(pi)));
            vals = vals(isfinite(vals));
            vG = [vG; vals(:)]; gGidx = [gGidx; repmat(pi, numel(vals),1)]; %#ok<AGROW>
        end
        pKW = NaN;
        if numel(unique(gGidx)) >= 2 && numel(vG) >= 5
            try, pKW = kruskalwallis(vG, gGidx, 'off'); catch, end
        end
        statsLines{end+1} = sprintf('%s KW(fallback) p=%s %s', grpLabel, fmtP(pKW), starStr(pKW)); %#ok<AGROW>
    end
end

% Per-period between-group (ranksum, Holm)
for pi = 1:nPer
    if isfinite(rawP(pi))
        nA = nnz(string(T.Group)=="Active"  & string(T.Period)==periods(pi) & isfinite(double(T.(ycol))));
        nP = nnz(string(T.Group)=="Passive" & string(T.Period)==periods(pi) & isfinite(double(T.(ycol))));
        statsLines{end+1} = sprintf('%s: Ranksum A(n=%d) vs P(n=%d) adj.p=%s %s', ...
            periods(pi), nA, nP, fmtP(adjP(pi)), starStr(adjP(pi))); %#ok<AGROW>
    end
end

if ~isempty(statsLines)
    addStatsTextbox(gca, statsLines);
end
end

%% --- Event-locked pupil analysis (Task 10.2) ---
function task_10_event_locked(csvPath, D, outDir)
preWin_s  = 2;
postWin_s = 2;
bin_s     = 0.2;
earlyDays = 3:5;
lateDays  = 6:10;
focusA    = 7:9;
focusB    = 12:13;

subDir = fullfile(outDir, 'event_locked');
if ~exist(subDir,'dir'), mkdir(subDir); end

T = readtable(csvPath, 'VariableNamingRule','preserve');
T = ensureString(T, 'mouse_key');
if ~ismember('day_index', T.Properties.VariableNames)
    fprintf('  Skipping 10.2: day_index not found in CSV.\n'); return;
end
if ~ismember('session_idx', T.Properties.VariableNames)
    T.session_idx = ones(height(T),1);
end

tbCol = pickTimebaseCol(T);
pupilCol = pickPupilRawCol(T);
if isempty(tbCol) || isempty(pupilCol)
    fprintf('  Skipping 10.2: missing time or pupil columns.\n'); return;
end

hasLick = ismember('Lick_TTL', T.Properties.VariableNames);
hasRew  = ismember('Injector_TTL', T.Properties.VariableNames);
if ~hasLick && ~hasRew
    fprintf('  Skipping 10.2: no Lick_TTL or Injector_TTL columns.\n'); return;
end

if hasLick, T.Lick_TTL(isnan(T.Lick_TTL)) = 0; end
if hasRew,  T.Injector_TTL(isnan(T.Injector_TTL)) = 0; end

% Attach Group from D using normalized mouse_key
T.mouse_key_norm = normalizeMouseKey(T.mouse_key);
map = unique(D(:,{'mouse_key','Group'}), 'rows', 'stable');
map.mouse_key_norm = normalizeMouseKey(map.mouse_key);
map.Group = string(map.Group);
T = outerjoin(T, map(:,{'mouse_key_norm','Group'}), 'Keys','mouse_key_norm', ...
    'MergeKeys', true, 'Type','left');
if ~ismember('Group', T.Properties.VariableNames)
    T.Group = repmat("Active", height(T), 1);
else
    g = string(T.Group);
    g(ismissing(g)) = "Active";
    T.Group = g;
end

keys = {'mouse_key_norm','day_index','session_idx'};
G = findgroups(T(:,keys));
tAxis = -preWin_s:bin_s:postWin_s;

E = table( ...
    'Size', [0 7], ...
    'VariableTypes', {'string','double','double','string','string','string','cell'}, ...
    'VariableNames', {'mouse_key','day_index','session_idx','Group','eventType','subType','trace'});

for gi = 1:max(G)
    r = (G == gi);
    t = double(T.(tbCol)(r));
    pup = double(T.(pupilCol)(r));
    if hasLick, lk = double(T.Lick_TTL(r)); else, lk = []; end
    if hasRew,  rw = double(T.Injector_TTL(r)); else, rw = []; end
    if isempty(t) || all(isnan(t)) || all(isnan(pup)), continue; end
    [t, ord] = sort(t);
    pup = pup(ord);
    if hasLick, lk = lk(ord); end
    if hasRew,  rw = rw(ord); end

    mk = string(T.mouse_key_norm(find(r,1,'first')));
    day = double(T.day_index(find(r,1,'first')));
    ses = double(T.session_idx(find(r,1,'first')));
    grp = string(T.Group(find(r,1,'first')));

    if hasRew
        rewardTimes = detectRisingEdges(t, rw);
        rewardTimes = enforceMinSeparation(rewardTimes, 0.5);
        if ~isempty(rewardTimes)
            ER = extractEventLocked(t, pup, rewardTimes, tAxis);
            E = addEventRow(E, mk, day, ses, grp, "reward", "all", ER);
        end
    end

    if hasLick
        lickTimes = detectRisingEdges(t, lk);
        lickTimes = enforceMinSeparation(lickTimes, 0.02);
        if ~isempty(lickTimes)
            [boutStart, boutEnd] = makeLickBouts(lickTimes, 2.0);
            isRewarded = false(size(boutStart));
            if hasRew
                rewardTimes = detectRisingEdges(t, rw);
                rewardTimes = enforceMinSeparation(rewardTimes, 0.5);
                for bi = 1:numel(boutStart)
                    isRewarded(bi) = any(rewardTimes >= boutStart(bi) & rewardTimes <= boutEnd(bi) + 1.0);
                end
            end
            EB_rs = extractEventLocked(t, pup, boutStart(isRewarded), tAxis);
            EB_ns = extractEventLocked(t, pup, boutStart(~isRewarded), tAxis);
            E = addEventRow(E, mk, day, ses, grp, "lickBoutStart", "rewarded", EB_rs);
            E = addEventRow(E, mk, day, ses, grp, "lickBoutStart", "nonreward", EB_ns);
        end
    end
end

if isempty(E)
    fprintf('  Skipping 10.2: no event-locked traces.\n'); return;
end

if hasLick
    plotEventLockedComparison(E, tAxis, subDir, "Passive", "lickBoutStart", "rewarded", focusA, focusB, ...
        'PUPIL_lickBoutStart_rewarded_PASSIVE_days7-9_vs_12-13.png');
    plotEventLockedComparison(E, tAxis, subDir, "Passive", "lickBoutStart", "nonreward", focusA, focusB, ...
        'PUPIL_lickBoutStart_nonreward_PASSIVE_days7-9_vs_12-13.png');
    % Active lick bout: day 7-9 vs 12-13
    plotEventLockedComparison(E, tAxis, subDir, "Active", "lickBoutStart", "rewarded", focusA, focusB, ...
        'PUPIL_lickBoutStart_rewarded_ACTIVE_days7-9_vs_12-13.png');
    plotEventLockedComparison(E, tAxis, subDir, "Active", "lickBoutStart", "nonreward", focusA, focusB, ...
        'PUPIL_lickBoutStart_nonreward_ACTIVE_days7-9_vs_12-13.png');
end
if hasRew
    plotRewardLockedEarlyLate(E, tAxis, subDir, earlyDays, lateDays);
    plotRewardLockedByPeriod(E, tAxis, subDir);

    % *** Day 7-9 vs 12-13: reward-locked comparison for BOTH groups ***
    plotEventLockedComparison(E, tAxis, subDir, "Active", "reward", "all", focusA, focusB, ...
        'PUPIL_reward_ACTIVE_days7-9_vs_12-13.png');
    plotEventLockedComparison(E, tAxis, subDir, "Passive", "reward", "all", focusA, focusB, ...
        'PUPIL_reward_PASSIVE_days7-9_vs_12-13.png');

    % *** Day 7-9 vs 12-13: Active vs Passive direct comparison (same days) ***
    plotRewardLockedActiveVsPassive(E, tAxis, subDir, focusA, 'days7-9');
    plotRewardLockedActiveVsPassive(E, tAxis, subDir, focusB, 'days12-13');
end
end

function E = addEventRow(E, mk, day, ses, grp, eventType, subType, trace)
if isempty(trace), return; end
rowTbl = table(mk, day, ses, grp, eventType, subType, {trace}, ...
    'VariableNames', {'mouse_key','day_index','session_idx','Group','eventType','subType','trace'});
E = [E; rowTbl];
end

function times = detectRisingEdges(t, x01)
x01 = x01(:);
x01(~isfinite(x01)) = 0;
x01 = x01 > 0.5;
dx = diff([false; x01]);
idx = find(dx == 1);
times = t(idx);
times = times(isfinite(times));
end

function times = enforceMinSeparation(times, minSep)
if isempty(times), return; end
times = sort(times(:));
keep = true(size(times));
last = -Inf;
for i = 1:numel(times)
    if times(i) - last < minSep
        keep(i) = false;
    else
        last = times(i);
    end
end
times = times(keep);
end

function [boutStart, boutEnd] = makeLickBouts(lickTimes, gap_s)
lickTimes = sort(lickTimes(:));
if isempty(lickTimes)
    boutStart = []; boutEnd = []; return;
end
boutStart = lickTimes(1);
boutEnd   = lickTimes(1);
starts = [];
ends = [];
for i = 2:numel(lickTimes)
    if lickTimes(i) - lickTimes(i-1) > gap_s
        starts(end+1,1) = boutStart; %#ok<AGROW>
        ends(end+1,1)   = boutEnd;   %#ok<AGROW>
        boutStart = lickTimes(i);
        boutEnd   = lickTimes(i);
    else
        boutEnd = lickTimes(i);
    end
end
starts(end+1,1) = boutStart;
ends(end+1,1)   = boutEnd;
boutStart = starts;
boutEnd   = ends;
end

function trace = extractEventLocked(t, pup, eventTimes, tAxis)
if isempty(eventTimes), trace = []; return; end
finiteMask = isfinite(t) & isfinite(pup);
t = t(finiteMask);
pup = pup(finiteMask);
if numel(t) < 2, trace = []; return; end
[t, ord] = sort(t);
pup = pup(ord);
[t, ia] = unique(t, 'stable');
pup = pup(ia);

M = nan(numel(eventTimes), numel(tAxis));
for i = 1:numel(eventTimes)
    te = eventTimes(i) + tAxis;
    pi = interp1(t, pup, te, 'linear', nan);
    pre = (tAxis >= min(tAxis) & tAxis < 0);
    b = mean(pi(pre), 'omitnan');
    pi = pi - b;
    M(i,:) = pi;
end
trace = mean(M, 1, 'omitnan');
end

function plotEventLockedComparison(E, tAxis, outDir, groupName, eventType, subType, daysA, daysB, fname)
rG = (string(E.Group)==groupName) & (string(E.eventType)==eventType) & (string(E.subType)==subType);
A = stackTraces(E.trace(rG & ismember(E.day_index, daysA)));
B = stackTraces(E.trace(rG & ismember(E.day_index, daysB)));
if isempty(A) || isempty(B), return; end

muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));

% --- Statistics: compare mean post-event response (t>0) ---
postIdx = tAxis >= 0;
meanPostA = mean(A(:,postIdx), 2, 'omitnan');  % one value per trial
meanPostB = mean(B(:,postIdx), 2, 'omitnan');
pvPost = safeRanksum(meanPostA, meanPostB);

% --- Statistics: compare peak response ---
peakA = max(A(:,postIdx), [], 2, 'omitnan');
peakB = max(B(:,postIdx), [], 2, 'omitnan');
pvPeak = safeRanksum(peakA, peakB);

fig = figure('Color','w','Position',[80 80 900 580]); hold on
shaded(tAxis, muA, seA, [0.2 0.4 0.9]);
hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', [0.2 0.4 0.9]);
shaded(tAxis, muB, seB, [0.9 0.3 0.2]);
hB = plot(tAxis, muB, 'LineWidth', 2.5, 'Color', [0.9 0.3 0.2]);
xline(0,'k-'); yline(0,'k:');
xlabel('Time from event (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('%s %s (%s): Days %s vs %s', groupName, eventType, subType, rangeStr(daysA), rangeStr(daysB)));
legend([hA hB], {sprintf('Days %s (n=%d)', rangeStr(daysA), size(A,1)), ...
                  sprintf('Days %s (n=%d)', rangeStr(daysB), size(B,1))}, 'Location','best');
grid off; box off

% Stats annotation
statsLines = { ...
    sprintf('Mean post-event: Ranksum p=%s %s', fmtP(pvPost), starStr(pvPost)), ...
    sprintf('Peak response: Ranksum p=%s %s', fmtP(pvPeak), starStr(pvPeak)), ...
    sprintf('n(A)=%d trials, n(B)=%d trials', size(A,1), size(B,1))};
addStatsTextbox(gca, statsLines);

exportgraphics(fig, fullfile(outDir,fname), 'Resolution', 220);
close(fig);
end

function plotRewardLockedEarlyLate(E, tAxis, outDir, earlyDays, lateDays)
eventType = "reward";
subType   = "all";
fig = figure('Color','w','Position',[80 80 1200 580]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

postIdx = tAxis >= 0;

for gi = 1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);
    ax = nexttile; hold(ax,'on');
    rG = (string(E.Group)==g) & (string(E.eventType)==eventType) & (string(E.subType)==subType);
    A = stackTraces(E.trace(rG & ismember(E.day_index, earlyDays)));
    B = stackTraces(E.trace(rG & ismember(E.day_index, lateDays)));
    hA = []; hB = [];
    nA = 0; nB = 0;
    if ~isempty(A)
        muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
        shaded(tAxis, muA, seA, [0.2 0.4 0.9]); hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', [0.2 0.4 0.9]);
        nA = size(A,1);
    end
    if ~isempty(B)
        muB = mean(B,1,'omitnan'); seB = std(B,0,1,'omitnan')/sqrt(size(B,1));
        shaded(tAxis, muB, seB, [0.9 0.3 0.2]); hB = plot(tAxis, muB, 'LineWidth', 2.5, 'Color', [0.9 0.3 0.2]);
        nB = size(B,1);
    end
    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');

    % Stats: compare post-event mean pupil
    pvPost = NaN; pvPeak = NaN;
    if ~isempty(A) && ~isempty(B)
        meanPostA = mean(A(:,postIdx), 2, 'omitnan');
        meanPostB = mean(B(:,postIdx), 2, 'omitnan');
        pvPost = safeRanksum(meanPostA, meanPostB);
        pvPeak = safeRanksum(max(A(:,postIdx),[],2,'omitnan'), max(B(:,postIdx),[],2,'omitnan'));
    end

    title(ax, sprintf('%s: Days %s vs %s\nPost-event Ranksum p=%s %s', ...
        g, rangeStr(earlyDays), rangeStr(lateDays), fmtP(pvPost), starStr(pvPost)));
    if ~isempty(hA) && ~isempty(hB)
        legend([hA hB], {sprintf('Days %s (n=%d)',rangeStr(earlyDays), nA), ...
                          sprintf('Days %s (n=%d)',rangeStr(lateDays), nB)}, 'Location','best');
    end
    grid off; box off
end
exportgraphics(fig, fullfile(outDir,'PUPIL_reward_locked_early_vs_late_active_passive.png'), 'Resolution', 220);
close(fig);
end

function plotRewardLockedByPeriod(E, tAxis, outDir)
eventType = "reward";
subType   = "all";
periods = struct( ...
    'name',  {'pre','during','post','withdrawal','reexposure'}, ...
    'days',  {3:5, 6:10, 11:13, 14:16, 17:18}, ...
    'color', {[0.0 0.6 0.0], [1.0 0.5 0.0], [0.5 0.0 0.5], [0.5 0.8 1.0], [0.4 0.2 0.0]} ...
);

postIdx = tAxis >= 0;

fig = figure('Color','w','Position',[80 80 1200 620]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for gi = 1:2
    gLabel = ["Active","Passive"];
    g = gLabel(gi);
    ax = nexttile; hold(ax,'on');
    rG = (string(E.Group)==g) & (string(E.eventType)==eventType) & (string(E.subType)==subType);
    h = gobjects(0);
    leg = {};
    perMeanPost = cell(numel(periods),1);  % store mean post-event per trial for stats
    perNames = {};
    for pi = 1:numel(periods)
        A = stackTraces(E.trace(rG & ismember(E.day_index, periods(pi).days)));
        if isempty(A), perMeanPost{pi} = []; continue; end
        mu = mean(A,1,'omitnan'); se = std(A,0,1,'omitnan')/sqrt(size(A,1));
        shaded(tAxis, mu, se, periods(pi).color);
        h(end+1) = plot(tAxis, mu, 'LineWidth', 2.5, 'Color', periods(pi).color); %#ok<AGROW>
        leg{end+1} = sprintf('%s d%s (n=%d)', periods(pi).name, rangeStr(periods(pi).days), size(A,1)); %#ok<AGROW>
        perMeanPost{pi} = mean(A(:,postIdx), 2, 'omitnan');
        perNames{end+1} = periods(pi).name; %#ok<AGROW>
    end
    xline(0,'k-'); yline(0,'k:');
    xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');

    % KW across periods for post-event mean pupil (trial-level, not mouse-level)
    allVals = []; allGrps = [];
    for pi = 1:numel(periods)
        v = perMeanPost{pi};
        if isempty(v), continue; end
        allVals = [allVals; v(:)]; allGrps = [allGrps; repmat(pi, numel(v), 1)]; %#ok<AGROW>
    end
    kwP = NaN;
    if numel(unique(allGrps)) >= 2 && numel(allVals) >= 5
        try, kwP = kruskalwallis(allVals, allGrps, 'off'); catch, end
    end

    % Pairwise post-hoc: pre vs each other period (trial-level Ranksum)
    postHocStr = '';
    if ~isempty(perMeanPost{1}) && numel(perMeanPost{1}) >= 2
        phRaw = nan(1, numel(periods));
        for pi2 = 2:numel(periods)
            if isempty(perMeanPost{pi2}), continue; end
            phRaw(pi2) = safeRanksum(perMeanPost{1}, perMeanPost{pi2});
        end
        phAdj = holmCorrect(phRaw);
        sigPh = {};
        for pi2 = 2:numel(periods)
            if isfinite(phAdj(pi2)) && phAdj(pi2) < 0.05
                sigPh{end+1} = sprintf('pre vs %s p=%s%s', periods(pi2).name, fmtP(phAdj(pi2)), starStr(phAdj(pi2))); %#ok<AGROW>
            end
        end
        if ~isempty(sigPh), postHocStr = strjoin(sigPh, '; '); end
    end

    titleStr = sprintf('Reward-locked by period (%s)  KW(trial-level) p=%s %s', g, fmtP(kwP), starStr(kwP));
    if ~isempty(postHocStr)
        titleStr = sprintf('%s\n%s', titleStr, postHocStr);
    end
    title(ax, titleStr, 'FontSize', 8);

    if ~isempty(h)
        legend(h, leg, 'Location','eastoutside');
    end
    grid off; box off
end
exportgraphics(fig, fullfile(outDir,'PUPIL_reward_locked_by_period_active_passive.png'), 'Resolution', 220);
close(fig);
end

function plotRewardLockedActiveVsPassive(E, tAxis, outDir, days, dayLabel)
% Direct Active vs Passive comparison of reward-locked pupil for given days
eventType = "reward"; subType = "all";
rA = (string(E.Group)=="Active") & (string(E.eventType)==eventType) & (string(E.subType)==subType) & ismember(E.day_index, days);
rP = (string(E.Group)=="Passive") & (string(E.eventType)==eventType) & (string(E.subType)==subType) & ismember(E.day_index, days);
A = stackTraces(E.trace(rA));
B = stackTraces(E.trace(rP));
if isempty(A) || isempty(B), return; end

postIdx = tAxis >= 0;
meanPostA = mean(A(:,postIdx), 2, 'omitnan');
meanPostP = mean(B(:,postIdx), 2, 'omitnan');
pvPost = safeRanksum(meanPostA, meanPostP);

fig = figure('Color','w','Position',[80 80 900 580]); hold on;
muA = mean(A,1,'omitnan'); seA = std(A,0,1,'omitnan')/sqrt(size(A,1));
muP = mean(B,1,'omitnan'); seP = std(B,0,1,'omitnan')/sqrt(size(B,1));
COL = groupColors();
shaded(tAxis, muA, seA, COL.Active);
hA = plot(tAxis, muA, 'LineWidth', 2.5, 'Color', COL.Active);
shaded(tAxis, muP, seP, COL.Passive);
hP = plot(tAxis, muP, 'LineWidth', 2.5, 'Color', COL.Passive);
xline(0,'k-'); yline(0,'k:');
xlabel('Time from reward (s)'); ylabel('\Delta pupil (baseline-subtracted)');
title(sprintf('Reward-locked: Active vs Passive (%s)\nPost-event Ranksum p=%s %s', ...
    dayLabel, fmtP(pvPost), starStr(pvPost)));
legend([hA hP], {sprintf('Active (n=%d trials)', size(A,1)), ...
                  sprintf('Passive (n=%d trials)', size(B,1))}, 'Location','eastoutside');
grid off; box off;
exportgraphics(fig, fullfile(outDir, sprintf('PUPIL_reward_AvP_%s.png', dayLabel)), 'Resolution', 220);
close(fig);
end

function M = stackTraces(tracesCell)
if isempty(tracesCell), M=[]; return; end
n = numel(tracesCell{1});
M = nan(numel(tracesCell), n);
for i = 1:numel(tracesCell)
    v = tracesCell{i};
    if isempty(v) || numel(v) ~= n, continue; end
    M(i,:) = v(:)';
end
M = M(any(isfinite(M),2),:);
end

function shaded(x, mu, se, color)
if nargin < 4, color = 0.9*[1 1 1]; end
x = x(:)'; mu = mu(:)'; se = se(:)';
fill([x fliplr(x)], [mu-se fliplr(mu+se)], color, 'EdgeColor','none', 'FaceAlpha',0.25);
end

function s = rangeStr(v)
v = unique(v); v = sort(v);
if isempty(v), s=''; return; end
if numel(v)==1, s=sprintf('%d',v); return; end
if all(diff(v)==1), s=sprintf('%d-%d', v(1), v(end));
else, s=strjoin(string(v),'-'); end
end
