function analyze_modules_5_to_11_new()
% analyze_modules_5_to_11 (includes Module 12: Predictive Modeling)
% Implementation of Modules 5-12 from the Behavior Pipeline Implementation Plan.
% Uses preprocessed data from analyze_passive_active_dashboard_dec2 (S_D_cache.mat).
%
% Module 05: Feature QC — outlier detection, missingness heatmap, 10s bins
% Module 06: Mixed-effects / GLMM — LME, GLMM(Poisson), reward_type, random slopes
% Module 07: PCA + Clustering + EFA — parallel analysis, phase-specific, multi-modal
% Module 08: Event-locked analyses — pupil & lick around lick/reward events
% Module 09: Cumulative curve fitting — sigmoid/exponential per session
% Module 10: Cross-modal integration — partial corr, LME with assay covariates
% Module 11: RL scaffold — Actor-Critic TD model per mouse
%
% Requires: Statistics and Machine Learning Toolbox

%% ========== 1. LOCATE AND LOAD DATA ==========
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d), 'No run_* under %s', rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder, d(ix).name);

cacheMat = fullfile(runDir, 'S_D_cache.mat');
csvPath  = fullfile(runDir, 'ALL_mice_longitudinal.csv');
assert(exist(cacheMat,'file')>0, 'Missing cache: %s', cacheMat);

fprintf('Loading S/D from cache: %s\n', cacheMat);
L = load(cacheMat, 'S', 'D');
S = L.S; D = L.D;

% Convert mouse_key to string for robust parsing
S.mouse_key = string(S.mouse_key);
D.mouse_key = string(D.mouse_key);

% --- Group assignment: same logic as analyze_passive_active_dashboard_dec2 ---
cohort = buildNewCohortRoster();
D = attachCohortGroupAndPair(D, cohort);
S = attachCohortGroupAndPair(S, cohort);

% Sanity check group counts
uD = unique(D(:,{'mouse_key','Group'}),'rows','stable');
nA = nnz(uD.Group=="Active"); nP = nnz(uD.Group=="Passive");
fprintf('[GROUP CHECK] Active=%d, Passive=%d (total=%d)\n', nA, nP, height(uD));
if nA==0 || nP==0
    warning('Group assignment failed — one group is empty! Check mouse_key parsing.');
end

% Remove rows where Group could not be assigned
D = D(~isundefined(D.Group) & ~ismissing(D.Group), :);
S = S(~isundefined(S.Group) & ~ismissing(S.Group), :);

D = D(double(D.day_index) >= 3, :);

% -----------------------------------------------------------------------
% CRITICAL:  Passive mice do NOT do the PR task during "During" (day 6-10).
% They receive forced replay rewards — NO licking, NO PR breakpoint.
% → Set ALL lick / bout / reward-delivery metrics to NaN for Passive×During
%   so that (a) stats never compare Active vs Passive on licking in "During",
%   (b) dots / means are simply absent for Passive in that phase, and
%   (c) Active "During" data is still available for within-Active analysis.
% -----------------------------------------------------------------------
lickRelatedCols = {'RequirementLast', ...
    'lick_n','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s','lick_totalDur_s', ...
    'bout_n','bout_meanDur_s', ...
    'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s'};
maskPassiveDuring = (D.Group=="Passive") & (double(D.day_index)>=6) & (double(D.day_index)<=10);
nMasked = 0;
for ci = 1:numel(lickRelatedCols)
    col = lickRelatedCols{ci};
    if ismember(col, D.Properties.VariableNames)
        D.(col)(maskPassiveDuring) = NaN;
        nMasked = nMasked + 1;
    end
end
fprintf('[PASSIVE MASK] Set %d lick-related columns to NaN for %d Passive×During rows.\n', ...
    nMasked, nnz(maskPassiveDuring));

D.Period     = periodOfDay(double(D.day_index));
D.RewardType = rewardTypeOfDay(double(D.day_index));
D = D(~isundefined(D.Period), :);

% Session-level: same mask for S
S = S(double(S.day_index) >= 3, :);
S.Period = periodOfDay(double(S.day_index));
S = S(~isundefined(S.Period) & ~ismissing(S.Group), :);

maskPassiveDuringS = (S.Group=="Passive") & (double(S.day_index)>=6) & (double(S.day_index)<=10);
for ci = 1:numel(lickRelatedCols)
    col = lickRelatedCols{ci};
    if ismember(col, S.Properties.VariableNames)
        S.(col)(maskPassiveDuringS) = NaN;
    end
end

%% ========== 2. OUTPUT DIRECTORIES ==========
baseOut = fullfile(runDir, 'figs', 'modules_5_to_11');
dirs = struct();
dirs.mod5  = fullfile(baseOut, 'module_05_feature_qc');
dirs.mod6  = fullfile(baseOut, 'module_06_glmm');
dirs.mod7  = fullfile(baseOut, 'module_07_pca_efa');
dirs.mod8  = fullfile(baseOut, 'module_08_event_locked');
dirs.mod9  = fullfile(baseOut, 'module_09_cumulative_fit');
dirs.mod10 = fullfile(baseOut, 'module_10_crossmodal');
dirs.mod11 = fullfile(baseOut, 'module_11_rl_model');
dirs.mod12 = fullfile(baseOut, 'module_12_predictive');
dirs.stats = fullfile(baseOut, 'stats');
fn = fieldnames(dirs);
for i = 1:numel(fn), if ~exist(dirs.(fn{i}),'dir'), mkdir(dirs.(fn{i})); end; end

%% ========== 3. RUN MODULES ==========
fprintf('\n====== MODULE 05: Feature QC ======\n');
try, module_05_feature_qc(D, S, dirs.mod5, dirs.stats);
catch ME, fprintf('  [WARN] Module 5 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 06: GLMM / LME ======\n');
try, module_06_glmm(D, S, dirs.mod6, dirs.stats);
catch ME, fprintf('  [WARN] Module 6 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 07: PCA + Clustering + EFA ======\n');
try, module_07_pca_efa(D, dirs.mod7, dirs.stats);
catch ME, fprintf('  [WARN] Module 7 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 08: Event-locked Analyses ======\n');
try, module_08_event_locked(csvPath, D, dirs.mod8, dirs.stats);
catch ME, fprintf('  [WARN] Module 8 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 09: Cumulative Curve Fitting ======\n');
try, module_09_cumulative_fit(csvPath, D, dirs.mod9, dirs.stats);
catch ME, fprintf('  [WARN] Module 9 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 10: Cross-Modal Integration ======\n');
try, module_10_crossmodal(D, dirs.mod10, dirs.stats);
catch ME, fprintf('  [WARN] Module 10 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 11: RL Model Scaffold ======\n');
try, module_11_rl_model(D, dirs.mod11, dirs.stats);
catch ME, fprintf('  [WARN] Module 11 failed: %s\n', ME.message); end

fprintf('\n====== MODULE 12: Predictive Modeling ======\n');
try, module_12_predictive(D, dirs.mod12, dirs.stats);
catch ME, fprintf('  [WARN] Module 12 failed: %s\n', ME.message); end

fprintf('\n====== Writing Module 11 & 12 Explanation (Korean) ======\n');
try, write_module_11_12_explanation(baseOut);
catch ME, fprintf('  [WARN] Explanation write failed: %s\n', ME.message); end

fprintf('\n===================================================\n');
fprintf('Modules 5-12 complete. Outputs in:\n  %s\n', baseOut);
end

%% ########################################################################
%  MODULE 05 — Feature QC: Outliers, Missingness, Distribution, Export
%  ########################################################################
function module_05_feature_qc(D, S, outDir, statsDir)
feats = pickNumericFeatures(D);
fprintf('  QC on %d features across %d day-level rows.\n', numel(feats), height(D));

% --- 5a. Missingness heatmap (mouse x feature) ---
%  Separate features into:  (A) daily behavioral  (B) sporadic-day assay (TST/HOT/Straub)
%  For (B), missingness is expected on most days, so only show (A) in the main heatmap.
%  Also, Passive mice have NaN lick data during "During" period by design — mark as expected.

mice = unique(D.mouse_key, 'stable');

% Classify features
isSporadicAssay = contains(feats,'TST') | contains(feats,'HOT') | contains(feats,'STRAUB') | ...
    contains(feats,'Straub','IgnoreCase',true);
dailyFeats   = feats(~isSporadicAssay);
assayFeats   = feats(isSporadicAssay);

lickRelated = {'RequirementLast','lick_n','lick_freq_per_min','lick_meanDur_s', ...
    'lick_medianIEI_s','lick_totalDur_s','bout_n','bout_meanDur_s', ...
    'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s'};

% --- Heatmap A: Daily behavioral features (excluding expected NaN) ---
missMat = nan(numel(mice), numel(dailyFeats));
for i = 1:numel(mice)
    r = D.mouse_key == mice(i);
    grp = string(D.Group(find(r,1,'first')));
    for j = 1:numel(dailyFeats)
        fn = dailyFeats{j};
        vals = double(D.(fn)(r));
        days = double(D.day_index(r));
        % Mark expected NaN: Passive lick cols during During (day 6-10)
        if grp == "Passive" && ismember(fn, lickRelated)
            expected = (days >= 6 & days <= 10);
            vals(expected) = [];  % exclude from missingness calc
        end
        if isempty(vals)
            missMat(i,j) = 0;  % all data was expected NaN → not a QC issue
        else
            missMat(i,j) = mean(~isfinite(vals)) * 100;
        end
    end
end

fh = figure('Color','w','Position',[40 40 max(600, 30*numel(dailyFeats)) 50*numel(mice)+120]);
imagesc(missMat); colorbar; colormap(flipud(hot));
set(gca,'XTick',1:numel(dailyFeats),'XTickLabel',dailyFeats,'XTickLabelRotation',40, ...
    'YTick',1:numel(mice),'YTickLabel',cellstr(string(mice)));
title('Unexpected Missingness — daily behavioral features (% NaN, expected NaN excluded)');
caxis([0 100]); grid off; box on;
addBottomStatNote(fh, 'Passive lick/PR data during During (day6-10) excluded (expected NaN). TST/HOT/Straub shown separately.');
printpng(fh, fullfile(outDir, 'missingness_heatmap.png')); close(fh);

% --- Heatmap B: Sporadic assay features (TST / HOT / Straub) ---
if ~isempty(assayFeats)
    % Show which days each assay is available (across all mice)
    assayMiss = nan(numel(mice), numel(assayFeats));
    for i = 1:numel(mice)
        r = D.mouse_key == mice(i);
        for j = 1:numel(assayFeats)
            vals = double(D.(assayFeats{j})(r));
            assayMiss(i,j) = sum(isfinite(vals));  % count of AVAILABLE days (not %)
        end
    end
    fh = figure('Color','w','Position',[40 40 max(400, 30*numel(assayFeats)) 50*numel(mice)+120]);
    imagesc(assayMiss); colorbar; colormap(parula);
    set(gca,'XTick',1:numel(assayFeats),'XTickLabel',assayFeats,'XTickLabelRotation',40, ...
        'YTick',1:numel(mice),'YTickLabel',cellstr(string(mice)));
    title('Assay data availability — # of days with data (TST/HOT/Straub)');
    grid off; box on;
    addBottomStatNote(fh, 'These assays are only collected on specific test days — not daily. Color = number of available days.');
    printpng(fh, fullfile(outDir, 'assay_availability_heatmap.png')); close(fh);
end

% --- 5b. Outlier detection (IQR method, threshold=3) ---
iqrThr = 3;
outlierCount = zeros(numel(feats),1);
outlierTbl = table();
for j = 1:numel(feats)
    x = double(D.(feats{j}));
    q1 = prctile(x, 25); q3 = prctile(x, 75);
    iqr = q3 - q1;
    lo = q1 - iqrThr*iqr; hi = q3 + iqrThr*iqr;
    isOut = isfinite(x) & (x < lo | x > hi);
    outlierCount(j) = sum(isOut);
end
outlierTbl = table(feats(:), outlierCount, 'VariableNames', {'Feature','N_outliers'});
writetable(outlierTbl, fullfile(statsDir, 'mod05_outlier_counts.csv'));

% --- 5c. Distribution plots (histogram per feature) ---
nPerPage = 6;
nPages = ceil(numel(feats)/nPerPage);
for pg = 1:nPages
    fh = figure('Color','w','Position',[40 40 1100 700]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for k = 1:nPerPage
        idx = (pg-1)*nPerPage + k;
        if idx > numel(feats), break; end
        ax = nexttile; hold(ax,'on');
        x = double(D.(feats{idx})); x = x(isfinite(x));
        if isempty(x), title(ax, feats{idx}); continue; end
        histogram(ax, x, 20, 'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.7, 'EdgeColor','w');
        title(ax, feats{idx}, 'Interpreter','none', 'FontSize',9);
        grid(ax,'off'); box(ax,'off');
    end
    printpng(fh, fullfile(outDir, sprintf('distributions_page%d.png', pg))); close(fh);
end

% --- 5d. Export features_session CSV ---
writetable(D, fullfile(statsDir, 'features_day_level.csv'));
fprintf('  Exported features_day_level.csv (%d rows, %d cols).\n', height(D), width(D));
fprintf('  Module 05 done.\n');
end

%% ########################################################################
%  MODULE 06 — Mixed-Effects / GLMM
%  ########################################################################
function module_06_glmm(D, S, outDir, statsDir)
COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% === 6a. LME for continuous outcomes with reward_type ===
contMetrics = {'RequirementLast','pupil_mean','lick_freq_per_min','lick_medianIEI_s','bout_meanDur_s'};
contLabels  = {'PR Breakpoint','Pupil Mean','Lick Rate','ILI','Bout Mean Dur'};

lmeSummaryRows = cell(0,8);   % compact comparison table (as before)
allFixedCoefs  = table();      % full coefficient table across all metrics/models
allAnovaTbls   = table();      % full ANOVA tables
allFitStats    = table();      % AIC / BIC / LogLik / R2
forestData     = cell(0,6);    % for forest plot: {metric, term, beta, lower, upper, pval}

for mi = 1:numel(contMetrics)
    ycol = contMetrics{mi};
    if ~ismember(ycol, D.Properties.VariableNames), continue; end

    T = D(:, {'mouse_key','day_index','Period','Group'});
    T.y = double(D.(ycol));
    T.RewardType = rewardTypeOfDay(double(D.day_index));
    T = T(isfinite(T.y), :);

    % Exclude "During" period for lick-related metrics:
    % Passive mice have no lick data in During, creating structural imbalance
    % that harms Phase*Grp interaction estimation.
    isLick = ismember(ycol, {'RequirementLast','lick_freq_per_min','lick_meanDur_s', ...
        'lick_medianIEI_s','lick_totalDur_s','bout_n','bout_meanDur_s','bout_totalDur_s', ...
        'rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s', ...
        'Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min'});
    if isLick
        T = T(string(T.Period) ~= "During", :);
        fprintf('  LME(%s): excluded During period (Passive has no PR data).\n', ycol);
    end

    T.Phase = removecats(categorical(string(T.Period)));
    T.Grp   = removecats(categorical(string(T.Group)));
    T.Rew   = removecats(T.RewardType);
    T.mouse = categorical(T.mouse_key);
    T.day   = double(T.day_index);
    if height(T) < 15, continue; end

    % NOTE: RewardType is perfectly collinear with Phase (deterministic from day),
    %   so 'Phase + Rew' is rank-deficient.  We use Rew*Grp as a separate model instead.
    modelSpecs = { ...
        'LME_additive',    'y ~ Phase + Grp + (1|mouse)'; ...
        'LME_interaction', 'y ~ Phase * Grp + (1|mouse)'; ...
        'LME_randSlope',   'y ~ Phase + Grp + (day|mouse)'; ...
        'LME_rewXgrp',     'y ~ Rew * Grp + (1|mouse)'};

    for mdi = 1:size(modelSpecs,1)
        mName = modelSpecs{mdi,1}; mFormula = modelSpecs{mdi,2};
        try
            lme = fitlme(T, mFormula);
            aov = anova(lme);

            % --- Fixed-effects coefficient table ---
            fe = dataset2table(lme.Coefficients);
            fe.metric = repmat(string(ycol), height(fe), 1);
            fe.model  = repmat(string(mName), height(fe), 1);
            % Normalise column names (MATLAB versions differ)
            fe = normCoefCols(fe);
            allFixedCoefs = [allFixedCoefs; fe]; %#ok<AGROW>

            % --- ANOVA table ---
            aovT = safeAnov2Table(aov);
            aovT.metric = repmat(string(ycol), height(aovT), 1);
            aovT.model  = repmat(string(mName), height(aovT), 1);
            allAnovaTbls = [allAnovaTbls; aovT]; %#ok<AGROW>

            % --- Fit statistics ---
            [~,~,feStat] = fixedEffects(lme); % used only for R2
            condR2 = lme.Rsquared.Adjusted;
            margR2 = lme.Rsquared.Ordinary;
            fitRow = table(string(ycol), string(mName), lme.LogLikelihood, ...
                lme.ModelCriterion.AIC, lme.ModelCriterion.BIC, ...
                margR2, condR2, height(T), ...
                'VariableNames', {'metric','model','LogLik','AIC','BIC','R2_marginal','R2_conditional','N'});
            allFitStats = [allFitStats; fitRow]; %#ok<AGROW>

            % --- Forest-plot data (Group effect from interaction model) ---
            if strcmp(mName, 'LME_interaction')
                feT = normCoefCols(dataset2table(lme.Coefficients));
                grpRow = feT(contains(string(feT.Name),'Grp','IgnoreCase',true) & ...
                             ~contains(string(feT.Name),':'), :);
                if height(grpRow) >= 1
                    beta = grpRow.Estimate(1); se_ = grpRow.SE(1); p_ = grpRow.pValue(1);
                    forestData(end+1,:) = {ycol, char(grpRow.Name(1)), beta, beta-1.96*se_, beta+1.96*se_, p_}; %#ok<AGROW>
                end
            end

            % --- Compact summary row ---
            pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp');
            pInt = NaN;
            if strcmp(mName,'LME_interaction'), pInt = safeAnovP(aov,'Phase:Grp');
            elseif strcmp(mName,'LME_rewXgrp')
                pPh = safeAnovP(aov,'Rew'); pInt = safeAnovP(aov,'Rew:Grp');
            end
            lmeSummaryRows(end+1,:) = {ycol, mName, mFormula, pPh, pGr, pInt, lme.LogLikelihood, height(T)}; %#ok<AGROW>

            % --- Console printout ---
            fprintf('\n======== LME: %s | %s ========\n', ycol, mName);
            fprintf('  Formula: %s\n', mFormula);
            fprintf('  N obs = %d, N mice = %d\n', height(T), numel(unique(T.mouse)));
            fprintf('  AIC = %.1f, BIC = %.1f, LogLik = %.1f\n', ...
                lme.ModelCriterion.AIC, lme.ModelCriterion.BIC, lme.LogLikelihood);
            fprintf('  R2 marginal = %.3f, R2 conditional = %.3f\n', margR2, condR2);
            fprintf('  --- Fixed Effects ---\n');
            disp(normCoefCols(dataset2table(lme.Coefficients)));
            fprintf('  --- ANOVA (Type III) ---\n');
            disp(aovT(aovT.metric==string(ycol),:));
        catch ME
            fprintf('  LME(%s, %s) failed: %s\n', ycol, mName, ME.message);
        end
    end
end

% Save comprehensive LME output tables
if ~isempty(lmeSummaryRows)
    lmeTbl = cell2table(lmeSummaryRows, 'VariableNames', ...
        {'metric','model','formula','pPhase','pGroup','pInteraction_or_Rew','LogLik','N'});
    writetable(lmeTbl, fullfile(statsDir, 'mod06_lme_model_comparison.csv'));
end
if ~isempty(allFixedCoefs)
    writetable(allFixedCoefs, fullfile(statsDir, 'mod06_lme_fixed_effects.csv'));
    fprintf('  Saved mod06_lme_fixed_effects.csv (%d rows).\n', height(allFixedCoefs));
end
if ~isempty(allAnovaTbls)
    writetable(allAnovaTbls, fullfile(statsDir, 'mod06_lme_anova_tables.csv'));
    fprintf('  Saved mod06_lme_anova_tables.csv.\n');
end
if ~isempty(allFitStats)
    writetable(allFitStats, fullfile(statsDir, 'mod06_lme_fit_statistics.csv'));
    fprintf('  Saved mod06_lme_fit_statistics.csv.\n');
    % Print quick fit-stat comparison to console
    fprintf('\n  ---- LME Fit Statistics Summary ----\n');
    disp(allFitStats);
end

% === 6a-forest. Forest plot of Group effect (from LME_interaction) ===
if ~isempty(forestData)
    nF = size(forestData,1);
    figH = max(300, 100 + 65*nF);  % enough height for labels + bottom note
    fh = figure('Color','w','Position',[80 80 800 figH]); 
    % Leave bottom margin for stat note
    ax = axes('Position', [0.22 0.15 0.70 0.75]); hold(ax, 'on');
    yPos = 1:nF;
    revIdx = nF:-1:1;
    for fi = 1:nF
        ri = revIdx(fi);
        beta = forestData{fi,3}; lo = forestData{fi,4}; hi = forestData{fi,5}; p_ = forestData{fi,6};
        plot(ax, [lo hi], [ri ri], 'k-', 'LineWidth', 1.5);
        if p_ < 0.05, mc = [0.85 0.15 0.15]; else, mc = [0.3 0.3 0.3]; end
        plot(ax, beta, ri, 'o', 'MarkerSize', 8, 'MarkerFaceColor', mc, 'MarkerEdgeColor','k');
        text(ax, hi+0.02*abs(hi-lo+1), ri, sprintf(' %s', starStr(p_)), 'FontSize', 9);
    end
    xline(ax, 0, 'k--', 'LineWidth', 0.8);
    set(ax, 'YTick', yPos, 'YTickLabel', forestData(revIdx,1), 'YLim', [0.5 nF+0.5]);
    xlabel(ax, 'Group Effect (\beta_{Active vs Passive})  [95% CI]');
    title(ax, 'LME Fixed Effect: Group (Interaction Model)');
    grid(ax,'off'); box(ax,'off');
    addBottomStatNote(fh, 'LME: y~Phase*Grp+(1|mouse). Dot=beta, line=95%CI. Red=p<.05, gray=n.s.');
    printpng(fh, fullfile(outDir, 'forest_group_effect_LME.png')); close(fh);
    fprintf('  Saved forest_group_effect_LME.png\n');
end

% === 6b. GLMM for count outcomes (Poisson) ===
% Use session-level S for integer counts, fallback to D
% Also include RequirementLast as pseudo-count (it IS a count: PR ratio)
countMetrics = {'lick_n','rew_n','bout_n','RequirementLast'};
countLabels  = {'Lick Count','Reward Count','Bout Count','PR Breakpoint (count)'};

allGlmmCoefs = table();
allGlmmAnova = table();
glmmSummary  = cell(0,8);

for mi = 1:numel(countMetrics)
    ycol = countMetrics{mi};
    % Try S first (session-level), fall back to D (day-level)
    if ismember(ycol, S.Properties.VariableNames)
        T = S(:, intersect({'mouse_key','day_index','Period','Group'}, S.Properties.VariableNames));
        T.y = round(double(S.(ycol)));
    elseif ismember(ycol, D.Properties.VariableNames)
        T = D(:, {'mouse_key','day_index','Period','Group'});
        T.y = round(double(D.(ycol)));
        fprintf('  GLMM(%s): using D (day-level) instead of S.\n', ycol);
    else
        fprintf('  GLMM(%s): column not found, skipping.\n', ycol);
        continue;
    end
    T = T(isfinite(T.y) & T.y >= 0, :);

    % Exclude During for lick-related counts (Passive has no data)
    T = T(string(T.Period) ~= "During", :);
    fprintf('  GLMM(%s): excluded During period (Passive has no lick data). N=%d\n', ycol, height(T));

    T.Phase = removecats(categorical(string(T.Period)));
    T.Grp   = removecats(categorical(string(T.Group)));
    T.mouse = categorical(T.mouse_key);
    if height(T) < 15, fprintf('  GLMM(%s): too few obs (%d), skipping.\n', ycol, height(T)); continue; end

    modelSpecs = { ...
        'GLMM_Poisson',     'y ~ Phase + Grp + (1|mouse)', 'Poisson'; ...
        'GLMM_Poisson_Int', 'y ~ Phase * Grp + (1|mouse)', 'Poisson'};

    for mdi = 1:size(modelSpecs,1)
        mName = modelSpecs{mdi,1}; mFormula = modelSpecs{mdi,2}; mDist = modelSpecs{mdi,3};
        try
            glme = fitglme(T, mFormula, 'Distribution', mDist);
            aov = anova(glme);

            % --- Fixed-effects coefficients ---
            fe = dataset2table(glme.Coefficients);
            fe.metric = repmat(string(ycol), height(fe), 1);
            fe.model  = repmat(string(mName), height(fe), 1);
            fe = normCoefCols(fe);
            allGlmmCoefs = [allGlmmCoefs; fe]; %#ok<AGROW>

            % --- ANOVA ---
            aovT = safeAnov2Table(aov);
            aovT.metric = repmat(string(ycol), height(aovT), 1);
            aovT.model  = repmat(string(mName), height(aovT), 1);
            allGlmmAnova = [allGlmmAnova; aovT]; %#ok<AGROW>

            % --- Summary row ---
            pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp');
            pInt = NaN;
            if contains(mName,'Int'), pInt = safeAnovP(aov,'Phase:Grp'); end
            glmmSummary(end+1,:) = {ycol, mName, mFormula, pPh, pGr, pInt, glme.LogLikelihood, height(T)}; %#ok<AGROW>

            % --- Console printout ---
            fprintf('\n======== GLMM: %s | %s ========\n', ycol, mName);
            fprintf('  Formula: %s  (Distribution: %s)\n', mFormula, mDist);
            fprintf('  N obs = %d, N mice = %d, LogLik = %.1f\n', height(T), numel(unique(T.mouse)), glme.LogLikelihood);
            fprintf('  --- Fixed Effects (log scale) ---\n');
            disp(normCoefCols(dataset2table(glme.Coefficients)));
            fprintf('  --- ANOVA ---\n');
            disp(aovT(aovT.metric==string(ycol),:));
        catch ME
            fprintf('  GLMM(%s, %s) failed: %s\n', ycol, mName, ME.message);
        end
    end
end

% Save comprehensive GLMM output tables
if ~isempty(glmmSummary)
    glmmTbl = cell2table(glmmSummary, 'VariableNames', ...
        {'metric','model','formula','pPhase','pGroup','pInteraction','LogLik','N'});
    writetable(glmmTbl, fullfile(statsDir, 'mod06_glmm_model_comparison.csv'));
end
if ~isempty(allGlmmCoefs)
    writetable(allGlmmCoefs, fullfile(statsDir, 'mod06_glmm_fixed_effects.csv'));
    fprintf('  Saved mod06_glmm_fixed_effects.csv (%d rows).\n', height(allGlmmCoefs));
end
if ~isempty(allGlmmAnova)
    writetable(allGlmmAnova, fullfile(statsDir, 'mod06_glmm_anova_tables.csv'));
    fprintf('  Saved mod06_glmm_anova_tables.csv.\n');
end

% GLMM effect plots for count metrics (with p-values)
for mi = 1:numel(countMetrics)
    ycol = countMetrics{mi};
    ylab = countLabels{mi};
    if ismember(ycol, D.Properties.VariableNames)
        Dg = D; useD = true;
    elseif ismember(ycol, S.Properties.VariableNames)
        Dg = S; useD = false;
    else
        continue;
    end
    if ~ismember('Period', Dg.Properties.VariableNames), continue; end

    fh = figure('Color','w','Position',[80 80 950 520]); hold on;
    allPeriods = ["Pre","Post","Withdrawal","Re-exposure"];
    muAll = struct(); seAll = struct();
    mouseMeansGLMM = struct('Active',{cell(1,numel(allPeriods))}, 'Passive',{cell(1,numel(allPeriods))});
    for gi = ["Active","Passive"]
        Wg = Dg(Dg.Group==gi,:);
        muV = nan(1,numel(allPeriods)); seV = nan(1,numel(allPeriods));
        c = COL.(char(gi)); xoff = 0.05*(gi=="Passive");
        for pi = 1:numel(allPeriods)
            sub = Wg(string(removecats(categorical(Wg.Period)))==allPeriods(pi),:);
            if isempty(sub), continue; end
            [G,~] = findgroups(sub.mouse_key);
            mMeans = splitapply(@(x) mean(x,'omitnan'), double(sub.(ycol)), G);
            mMeans = mMeans(isfinite(mMeans));
            mouseMeansGLMM.(char(gi)){pi} = mMeans;
            if isempty(mMeans), continue; end
            muV(pi) = mean(mMeans); seV(pi) = std(mMeans)/sqrt(numel(mMeans));
            jit = (rand(numel(mMeans),1)-0.5)*0.08;
            scatter(pi+xoff+jit, mMeans, 20, c, 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility','off');
        end
        muAll.(char(gi)) = muV; seAll.(char(gi)) = seV;
        errorbar((1:numel(allPeriods))+xoff, muV, seV, '-o', 'Color', c, ...
            'LineWidth', 1.8, 'MarkerFaceColor', c, 'CapSize', 5, 'DisplayName', char(gi));
    end

    % Per-phase ranksum p-values (Active vs Passive per-mouse means)
    ax = gca;
    for pi = 1:numel(allPeriods)
        vA = mouseMeansGLMM.Active{pi};
        vP = mouseMeansGLMM.Passive{pi};
        pRS = safeRanksum(vA, vP);
        addStarsBetweenGroups(ax, pi, muAll.Active(pi), muAll.Passive(pi), ...
            seAll.Active(pi), seAll.Passive(pi), pRS);
    end

    % GLMM Poisson p-values for title annotation
    pPh = NaN; pGr = NaN; pInt = NaN;
    try
        Tglmm = Dg(:, {'mouse_key','Period','Group'}); Tglmm.y = double(Dg.(ycol));
        Tglmm = Tglmm(isfinite(Tglmm.y),:);
        Tglmm.Phase = removecats(categorical(string(Tglmm.Period)));
        Tglmm.Grp   = removecats(categorical(string(Tglmm.Group)));
        Tglmm.mouse  = categorical(Tglmm.mouse_key);
        % Exclude During period
        Tglmm = Tglmm(Tglmm.Phase ~= "During",:);
        Tglmm.Phase = removecats(Tglmm.Phase);
        if height(Tglmm) > 10
            glme = fitglme(Tglmm, 'y ~ Phase * Grp + (1|mouse)', 'Distribution', 'Poisson');
            aov = anova(glme);
            pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp'); pInt = safeAnovP(aov,'Phase:Grp');
        end
    catch, end

    set(gca,'XTick',1:numel(allPeriods),'XTickLabel',allPeriods);
    xlabel('Phase'); ylabel(ylab);
    title(sprintf('GLMM: %s  (Phase %s, Grp %s, Interaction %s)', ylab, fmtP(pPh), fmtP(pGr), fmtP(pInt)));
    legend('Location','northeastoutside'); grid off; box off;
    nA = numel(unique(Dg.mouse_key(Dg.Group=="Active")));
    nP = numel(unique(Dg.mouse_key(Dg.Group=="Passive")));
    addBottomStatNote(fh, sprintf('GLMM Poisson: y~Phase*Grp+(1|mouse). During excluded. Stars=Wilcoxon rank-sum on per-mouse means. N: Active=%d, Passive=%d.', nA, nP));
    printpng(fh, fullfile(outDir, sprintf('glmm_effect_%s.png', safeName(ycol)))); close(fh);
end

% === 6c. Effect plots (Phase x Group) for key metrics + stats ===
for mi = 1:numel(contMetrics)
    ycol = contMetrics{mi}; ylab = contLabels{mi};
    if ~ismember(ycol, D.Properties.VariableNames), continue; end

    fh = figure('Color','w','Position',[80 80 950 500]); hold on;
    muAll = struct(); seAll = struct();
    mouseMeansStore = struct('Active',{cell(1,numel(periods))}, 'Passive',{cell(1,numel(periods))});
    for gi = ["Active","Passive"]
        Dg = D(D.Group==gi, :);
        mu = nan(1,numel(periods)); se = nan(1,numel(periods));
        xoff = 0.06*(gi=="Passive");
        c = COL.(char(gi));
        for pi = 1:numel(periods)
            sub = Dg(string(Dg.Period)==periods(pi), :);
            if isempty(sub), continue; end
            % Per-mouse mean (one value per mouse, averaged across days)
            [G, ~] = findgroups(sub.mouse_key);
            mMeans = splitapply(@(x) mean(x,'omitnan'), double(sub.(ycol)), G);
            mMeans = mMeans(isfinite(mMeans));
            mouseMeansStore.(char(gi)){pi} = mMeans;
            if isempty(mMeans), continue; end
            mu(pi) = mean(mMeans); se(pi) = std(mMeans)/sqrt(numel(mMeans));
            % Individual mouse dots (one dot = one mouse, jittered)
            jit = (rand(numel(mMeans),1)-0.5)*0.08;
            scatter(pi+xoff+jit, mMeans, 22, c, 'filled', 'MarkerFaceAlpha', 0.45, 'HandleVisibility','off');
        end
        muAll.(char(gi)) = mu; seAll.(char(gi)) = se;
        errorbar((1:numel(periods))+xoff, mu, se, '-o', ...
            'Color',c,'LineWidth',2,'MarkerFaceColor',c,'CapSize',6,'DisplayName',char(gi));
    end

    % Per-phase ranksum (Active vs Passive) using per-mouse means
    ax = gca;
    isLickMetric = ismember(ycol, {'RequirementLast','lick_n','lick_freq_per_min', ...
        'lick_meanDur_s','lick_medianIEI_s','lick_totalDur_s', ...
        'bout_n','bout_meanDur_s','rew_n','rew_freq_per_min','rew_meanDur_s', ...
        'rew_totalDur_s','rew_medianIRI_s'});
    for pi = 1:numel(periods)
        vA = mouseMeansStore.Active{pi};
        vP = mouseMeansStore.Passive{pi};
        % If Passive has no data for a lick metric in "During" — annotate why
        if isLickMetric && periods(pi)=="During" && isempty(vP)
            yLim = get(ax, 'YLim');
            text(ax, pi+0.06, yLim(1)+0.05*diff(yLim), 'Passive: no PR', ...
                'FontSize',7, 'Color',[0.5 0.5 0.5], 'FontAngle','italic', ...
                'HorizontalAlignment','center');
            continue;  % skip ranksum — comparison is meaningless
        end
        pRS = safeRanksum(vA, vP);
        addStarsBetweenGroups(ax, pi, muAll.Active(pi), muAll.Passive(pi), ...
            seAll.Active(pi), seAll.Passive(pi), pRS);
    end

    % Overall LME for title annotation
    Tlme = D(:, {'mouse_key','Period','Group'}); Tlme.y = double(D.(ycol));
    Tlme = Tlme(isfinite(Tlme.y),:);
    Tlme.Phase = removecats(categorical(string(Tlme.Period)));
    Tlme.Grp   = removecats(categorical(string(Tlme.Group)));
    Tlme.mouse = categorical(Tlme.mouse_key);
    pPh = NaN; pGr = NaN; pInt = NaN;
    if height(Tlme) >= 15
        try
            lme = fitlme(Tlme, 'y ~ Phase * Grp + (1|mouse)');
            aov = anova(lme); pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp'); pInt = safeAnovP(aov,'Phase:Grp');
        catch
            try lme = fitlme(Tlme, 'y ~ Phase + Grp + (1|mouse)'); aov = anova(lme); pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp'); catch, end
        end
    end
    set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
    xlabel('Phase'); ylabel(ylab);
    nA = numel(unique(D.mouse_key(D.Group=="Active")));
    nP = numel(unique(D.mouse_key(D.Group=="Passive")));
    title(sprintf('%s — Phase %s, Grp %s, Int %s', ylab, fmtP(pPh), fmtP(pGr), fmtP(pInt)));
    legend('Location','northeastoutside'); grid off; box off;
    addBottomStatNote(fh, sprintf( ...
        'LME: y~Phase*Grp+(1|mouse). Stars: Wilcoxon rank-sum on per-mouse means. nActive=%d, nPassive=%d. *p<.05 **p<.01 ***p<.001', nA, nP));
    printpng(fh, fullfile(outDir, sprintf('effect_plot_%s.png', safeName(ycol)))); close(fh);
end

fprintf('  Module 06 done.\n');
end

%% ########################################################################
%  MODULE 07 — PCA + Clustering + EFA
%  ########################################################################
function module_07_pca_efa(D, outDir, statsDir)
COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% Behavioral features
behavFeats = {'RequirementLast','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
         'lick_totalDur_s','bout_n','bout_meanDur_s','bout_totalDur_s', ...
         'rew_freq_per_min','rew_meanDur_s','rew_totalDur_s', ...
         'pupil_mean','Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min'};
behavFeats = behavFeats(ismember(behavFeats, D.Properties.VariableNames));

% Pharmacological / assay features (Immersion, TST, HOT plate, Straub tail)
pharmaFeats = {};
if ismember('Immersion_Latency_s', D.Properties.VariableNames), pharmaFeats{end+1}='Immersion_Latency_s'; end
pharmaCand = D.Properties.VariableNames(contains(D.Properties.VariableNames,'TST_Pct') | ...
    contains(D.Properties.VariableNames,'HOT_Pct') | ...
    contains(D.Properties.VariableNames,'STRAUB','IgnoreCase',true) | ...
    contains(D.Properties.VariableNames,'Straub','IgnoreCase',true));
pharmaFeats = [pharmaFeats, pharmaCand];
% Filter out sparse / zero-variance pharmacological columns
keepPh = true(1, numel(pharmaFeats));
for ci = 1:numel(pharmaFeats)
    vals = double(D.(pharmaFeats{ci}));
    nFin = sum(isfinite(vals)); valsFin = vals(isfinite(vals));
    vv = 0; if numel(valsFin) > 1, vv = var(valsFin); end
    if nFin < 15 || vv < 1e-8 || nFin < 0.1*height(D)
        fprintf('  Module 7: dropping pharma feature %s (nFinite=%d, var=%.2g).\n', pharmaFeats{ci}, nFin, vv);
        keepPh(ci) = false;
    end
end
pharmaFeats = pharmaFeats(keepPh);
fprintf('  Module 7: %d behavioral + %d pharmacological features.\n', numel(behavFeats), numel(pharmaFeats));

feats = [behavFeats, pharmaFeats];
fprintf('  Module 7: %d features total: %s\n', numel(feats), strjoin(feats,', '));
if numel(feats) < 4
    fprintf('  Skipping Module 7: too few features.\n'); return;
end

% Build mouse x period matrix (median per mouse per period)
mice = unique(D.mouse_key,'stable');
nCols = 3 + numel(feats);
rows = cell(0, nCols);  % explicit pre-allocation with correct width
for i = 1:numel(mice)
    for pi = 1:numel(periods)
        r = D.mouse_key==mice(i) & string(D.Period)==periods(pi);
        if ~any(r), continue; end
        row = cell(1, nCols);
        row{1} = char(string(mice(i)));
        row{2} = char(periods(pi));
        gVal = D.Group(find(r,1,'first'));
        if isundefined(gVal), continue; end   % skip undefined groups
        row{3} = char(string(gVal));
        for fi = 1:numel(feats)
            row{3+fi} = median(double(D.(feats{fi})(r)),'omitnan');
        end
        rows(end+1,:) = row; %#ok<AGROW>
    end
end

if isempty(rows), fprintf('  No data for Module 7.\n'); return; end
fprintf('  Module 7: built %d mouse x period rows, %d cols.\n', size(rows,1), size(rows,2));

varNames = [{'mouse_key','Period','Group'}, feats];
assert(numel(varNames) == size(rows,2), 'Column count mismatch: varNames=%d, rows=%d', numel(varNames), size(rows,2));
W = cell2table(rows, 'VariableNames', varNames);
for fi = 1:numel(feats)
    col = W.(feats{fi});
    if iscell(col)
        W.(feats{fi}) = cellfun(@double, col);
    elseif ~isnumeric(col)
        W.(feats{fi}) = double(col);
    end
end
W.Group  = categorical(string(W.Group), {'Active','Passive'});
W.Period = categorical(string(W.Period), cellstr(periods), 'Ordinal', true);

% Extract numeric matrix — column by column to avoid type mismatch
X = nan(height(W), numel(feats));
for fi = 1:numel(feats)
    X(:, fi) = double(W.(feats{fi}));
end

% Impute NaN with column median (important for pharma features which are sparse)
nImputed = 0;
for fi = 1:size(X,2)
    nanIdx = ~isfinite(X(:,fi));
    if any(nanIdx)
        colMed = median(X(~nanIdx,fi),'omitnan');
        if isfinite(colMed)
            X(nanIdx,fi) = colMed;
            nImputed = nImputed + nnz(nanIdx);
        end
    end
end
if nImputed > 0, fprintf('  Module 7: imputed %d NaN values with column median.\n', nImputed); end

% Drop rows where ALL values are still NaN (no data at all) or zero-var columns
ok = any(isfinite(X),2);
X = X(ok,:); W = W(ok,:);
% Drop any columns that are still all-NaN or zero variance after imputation
goodCol = true(1,size(X,2));
for fi = 1:size(X,2)
    v = X(:,fi); v = v(isfinite(v));
    if isempty(v) || var(v) < 1e-10, goodCol(fi) = false; end
end
if any(~goodCol)
    fprintf('  Module 7: removed %d zero-var columns after imputation.\n', nnz(~goodCol));
    X = X(:,goodCol); feats = feats(goodCol);
end
if size(X,1) < 6, fprintf('  Too few complete rows (%d).\n', size(X,1)); return; end
fprintf('  Module 7: %d rows, %d features for PCA/EFA.\n', size(X,1), size(X,2));

% Z-score (global, log1p for counts)
for j = 1:size(X,2)
    if all(X(:,j) >= 0), X(:,j) = log1p(X(:,j)); end
end
mu = mean(X,1); sd = std(X,0,1); sd(sd==0)=1;
Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;

% Pre-build Wper for reliable period comparison (ordinal categorical can mismatch with ==)
Wper = string(removecats(categorical(W.Period)));

% === 7A. PCA ===
try
[coeff, score, ~, ~, expl] = pca(Xz);
W.PC1 = score(:,1); W.PC2 = score(:,2);
if size(score,2)>=3, W.PC3 = score(:,3); end

% PCA variance explained bar chart
nPC = min(7, size(coeff,2));  % use 7 PCs
nShow = min(numel(expl), 10);
fh = figure('Color','w','Position',[80 80 600 350]); hold on;
bar(1:nShow, expl(1:nShow), 'FaceColor',[0.3 0.6 0.8],'EdgeColor','w');
plot(1:nShow, cumsum(expl(1:nShow)), 'ko-','LineWidth',1.5,'MarkerFaceColor','r');
xlabel('Principal Component'); ylabel('% Variance');
title(sprintf('PCA Variance Explained (cumulative %.0f%% in %d PCs)', sum(expl(1:min(nPC,numel(expl)))), nPC));
legend({'Individual','Cumulative'},'Location','east');
grid off; box off;
printpng(fh, fullfile(outDir, 'pca_variance_explained.png')); close(fh);
fprintf('  PCA: top %d PCs explain %.1f%% of variance.\n', nPC, sum(expl(1:nPC)));

% PCA scatter (PC1 vs PC2 colored by Group, shape by Period)
fh = figure('Color','w','Position',[80 80 750 550]); hold on;
markers = {'o','^','s','d','p'};
Wper = string(removecats(categorical(W.Period)));  % ensure plain string for comparison
for pi = 1:numel(periods)
    for gi = ["Active","Passive"]
        m = (Wper == periods(pi)) & (W.Group == gi);
        if ~any(m), continue; end
        c = COL.(char(gi));
        scatter(W.PC1(m), W.PC2(m), 50, c, 'filled', 'Marker', markers{min(pi,numel(markers))}, ...
            'MarkerEdgeColor','k', 'MarkerFaceAlpha', 0.7, 'HandleVisibility','off');
    end
end
scatter(nan,nan,50,COL.Active,'filled','DisplayName','Active');
scatter(nan,nan,50,COL.Passive,'filled','DisplayName','Passive');
for pi = 1:numel(periods)
    plot(nan,nan,'k','LineStyle','none','Marker',markers{min(pi,numel(markers))},'MarkerSize',8, ...
        'DisplayName',char(periods(pi)));
end
legend('Location','bestoutside');
xlabel(sprintf('PC1 (%.1f%%)',expl(1))); ylabel(sprintf('PC2 (%.1f%%)',expl(2)));
title('PCA of Mouse x Period Feature Vectors'); grid off; box off;
printpng(fh, fullfile(outDir, 'pca_scatter.png')); close(fh);

% PCA trajectory (PC1 across phases per group) + stats + dots
fh = figure('Color','w','Position',[80 80 850 500]); hold on;
muPC_grp = struct(); sePC_grp = struct();
for gi = ["Active","Passive"]
    Wg = W(W.Group==gi,:);
    muPC = nan(1,numel(periods)); sePC = nan(1,numel(periods));
    xoff = 0.05*(gi=="Passive"); c = COL.(char(gi));
    Wgper = string(removecats(categorical(Wg.Period)));
    for pi = 1:numel(periods)
        v = Wg.PC1(Wgper==periods(pi));
        v = v(isfinite(v));
        if isempty(v), continue; end
        muPC(pi) = mean(v); sePC(pi) = std(v)/sqrt(numel(v));
        jit = (rand(numel(v),1)-0.5)*0.08;
        scatter(pi+xoff+jit, v, 18, c, 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility','off');
    end
    muPC_grp.(char(gi)) = muPC; sePC_grp.(char(gi)) = sePC;
    errorbar((1:numel(periods))+xoff, muPC, sePC, '-o', ...
        'Color',c,'LineWidth',2,'MarkerFaceColor',c,'CapSize',6,'DisplayName',char(gi));
end
% Per-phase ranksum
ax = gca;
for pi = 1:numel(periods)
    vA = W.PC1(W.Group=="Active"  & Wper==periods(pi));
    vP = W.PC1(W.Group=="Passive" & Wper==periods(pi));
    pRS = safeRanksum(vA, vP);
    addStarsBetweenGroups(ax, pi, muPC_grp.Active(pi), muPC_grp.Passive(pi), ...
        sePC_grp.Active(pi), sePC_grp.Passive(pi), pRS);
end
% LME for title
pPh = NaN; pGr = NaN;
Tpc = W(:,{'mouse_key','Period','Group','PC1'}); Tpc.y = Tpc.PC1;
Tpc.Phase = removecats(categorical(string(Tpc.Period))); Tpc.Grp = removecats(Tpc.Group);
Tpc.mouse = categorical(Tpc.mouse_key);
if height(Tpc) >= 10
    try lme = fitlme(Tpc, 'y ~ Phase + Grp + (1|mouse)'); aov = anova(lme);
        pPh = safeAnovP(aov,'Phase'); pGr = safeAnovP(aov,'Grp'); catch, end
end
set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Phase'); ylabel('PC1 Score');
nA = nnz(W.Group=="Active"); nP = nnz(W.Group=="Passive");
title(sprintf('PC1 Trajectory — Phase %s, Grp %s', fmtP(pPh), fmtP(pGr)));
legend('Location','northeastoutside'); grid off; box off;
addBottomStatNote(fh, sprintf( ...
    'LME: PC1~Phase+Grp+(1|mouse). Stars: Wilcoxon rank-sum (per-mouse medians). nA=%d, nP=%d. *p<.05 **p<.01 ***p<.001', nA, nP));
printpng(fh, fullfile(outDir, 'pca_trajectory_pc1.png')); close(fh);

% PCA loadings heatmap (nPC already defined above)
fhW = max(700, 35*numel(feats)); fhH = max(350, 60*nPC);
fh = figure('Color','w','Position',[80 80 fhW fhH]);
imagesc(coeff(:,1:nPC)'); colorbar; caxis([-1 1]); colormap(redblue());
set(gca,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',45, ...
    'YTick',1:nPC,'YTickLabel',arrayfun(@(i) sprintf('PC%d (%.0f%%)',i,expl(i)),1:nPC,'Uni',false), ...
    'FontSize', max(6, 10-floor(numel(feats)/8)));
title(sprintf('PCA Loadings (top %d PCs)', nPC)); grid off; box on;
printpng(fh, fullfile(outDir, 'pca_loadings_heatmap.png')); close(fh);

% --- t-SNE visualization ---
if size(Xz,1) >= 10
    rng(42);
    perpVal = min(30, floor(size(Xz,1)/3));
    Y_tsne = tsne(Xz, 'Perplexity', max(5, perpVal), 'NumDimensions', 2);
    W.tSNE1 = Y_tsne(:,1); W.tSNE2 = Y_tsne(:,2);

    fh = figure('Color','w','Position',[80 80 750 550]); hold on;
    markers = {'o','^','s','d','p'};
    for pi = 1:numel(periods)
        for gi = ["Active","Passive"]
            m = (Wper == periods(pi)) & (W.Group == gi);
            if ~any(m), continue; end
            c = COL.(char(gi));
            scatter(W.tSNE1(m), W.tSNE2(m), 50, c, 'filled', 'Marker', markers{min(pi,numel(markers))}, ...
                'MarkerEdgeColor','k', 'MarkerFaceAlpha', 0.7, 'HandleVisibility','off');
        end
    end
    scatter(nan,nan,50,COL.Active,'filled','DisplayName','Active');
    scatter(nan,nan,50,COL.Passive,'filled','DisplayName','Passive');
    for pi = 1:numel(periods)
        plot(nan,nan,'k','LineStyle','none','Marker',markers{min(pi,numel(markers))},'MarkerSize',8,'DisplayName',char(periods(pi)));
    end
    legend('Location','bestoutside');
    xlabel('t-SNE 1'); ylabel('t-SNE 2');
    title('t-SNE of Mouse x Period Feature Vectors'); grid off; box off;
    printpng(fh, fullfile(outDir, 'tsne_scatter.png')); close(fh);
    fprintf('  t-SNE scatter saved.\n');
end

catch ME7a
    fprintf('  [WARN] Module 7A (PCA/t-SNE) failed: %s\n', ME7a.message);
    % Ensure W has PC columns even on failure, so downstream doesn't crash
    if ~ismember('PC1', W.Properties.VariableNames), W.PC1 = nan(height(W),1); end
    if ~ismember('PC2', W.Properties.VariableNames), W.PC2 = nan(height(W),1); end
    if ~exist('nPC','var'), nPC = min(7, numel(feats)); end
    if ~exist('coeff','var') || isempty(coeff), coeff = nan(numel(feats), numel(feats)); end
    if ~exist('expl','var') || isempty(expl), expl = nan(numel(feats),1); end
    if ~exist('Wper','var'), Wper = string(removecats(categorical(W.Period))); end
end

% === 7B. K-means + DBSCAN Clustering ===
try
Kopt = 3;  % user-specified K=3
Kmax = min(5, floor(size(Xz,1)/3));
silVals = nan(1,Kmax);
for K = 2:Kmax
    rng(42);
    [cidx,~] = kmeans(Xz, K, 'Replicates', 50, 'MaxIter', 300);
    silVals(K) = mean(silhouette(Xz, cidx));
end
rng(42);
[cidx, ~] = kmeans(Xz, Kopt, 'Replicates', 100, 'MaxIter', 500);
W.Cluster = cidx;

fh = figure('Color','w','Position',[80 80 500 400]);
bar(2:Kmax, silVals(2:Kmax), 'FaceColor', [0.3 0.6 0.8]);
hold on; bar(Kopt, silVals(Kopt), 'FaceColor', [0.9 0.3 0.2]);
xlabel('K (number of clusters)'); ylabel('Mean Silhouette');
title(sprintf('Silhouette Analysis (K=%d selected)', Kopt));
grid off; box off;
printpng(fh, fullfile(outDir, 'silhouette_analysis.png')); close(fh);

% PCA scatter colored by K-means cluster + dot identity labels
if ismember('PC1', W.Properties.VariableNames) && ismember('PC2', W.Properties.VariableNames)
    fh = figure('Color','w','Position',[80 80 1050 600]); hold on;
    cmap = lines(Kopt);
    periodMarkers = {'o','^','s','d','p'};  % Pre, During, Post, Withdrawal, Re-exposure
    periodList = ["Pre","During","Post","Withdrawal","Re-exposure"];

    % Draw cluster-colored dots with shape = period, fill = group
    for k = 1:Kopt
        % Dummy scatter for cluster legend entry
        scatter(nan, nan, 80, cmap(k,:), 'filled', 'MarkerEdgeColor','k', ...
            'DisplayName', sprintf('Cluster %d (n=%d)', k, nnz(cidx==k)));
    end
    for i = 1:height(W)
        k = cidx(i);
        per = string(Wper(i));
        grpStr = string(W.Group(i));
        pidx = find(periodList == per, 1);
        if isempty(pidx), pidx = 1; end
        mkr = periodMarkers{min(pidx, numel(periodMarkers))};
        if grpStr == "Active"
            scatter(W.PC1(i), W.PC2(i), 70, cmap(k,:), 'filled', 'Marker', mkr, ...
                'MarkerEdgeColor','k', 'LineWidth', 0.8, 'HandleVisibility','off');
        else  % Passive = open marker
            scatter(W.PC1(i), W.PC2(i), 70, cmap(k,:), 'Marker', mkr, ...
                'MarkerEdgeColor', cmap(k,:), 'LineWidth', 1.5, 'HandleVisibility','off');
        end
    end

    % Add text labels (short mouse name + period abbreviation)
    perAbbrev = ["Pr","Du","Po","Wd","Re"];
    for i = 1:height(W)
        mk = string(W.mouse_key(i));
        % Shorten mouse name: keep cage_color (e.g. "6100_blk")
        mk = regexprep(mk, '(\d{4})_?(orange)', '$1_org');
        mk = regexprep(mk, '(\d{4})_?(black)', '$1_blk');
        mk = regexprep(mk, '(\d{4})_?(white)', '$1_wht');
        mk = regexprep(mk, '(\d{4})_?(red)', '$1_red');
        per = string(Wper(i));
        pidx = find(periodList == per, 1);
        if isempty(pidx), pidx = 1; end
        lbl = sprintf('%s.%s', mk, perAbbrev(pidx));
        text(W.PC1(i)+0.12, W.PC2(i)+0.12, lbl, 'FontSize', 5, 'Color', [0.3 0.3 0.3], ...
            'Clipping','on');
    end

    % Legend entries for period shapes and group fill style
    for pi = 1:numel(periodList)
        plot(nan, nan, 'k', 'LineStyle','none', 'Marker', periodMarkers{pi}, ...
            'MarkerSize', 7, 'MarkerFaceColor','k', 'DisplayName', char(periodList(pi)));
    end
    scatter(nan, nan, 50, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor','k', 'DisplayName','Active (filled)');
    scatter(nan, nan, 50, [0.5 0.5 0.5], 'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth',1.5, 'DisplayName','Passive (open)');

    xlabel(sprintf('PC1 (%.1f%%)', expl(1))); ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
    title(sprintf('K-means Clusters on PCA (K=%d)', Kopt));
    legend('Location','bestoutside','FontSize',7); grid off; box off;
    addBottomStatNote(fh, sprintf('K-means (K=%d, silhouette=%.2f). %d obs, %d features. Shape=period, Filled=Active/Open=Passive. Labels=mouse.period.', ...
        Kopt, silVals(Kopt), size(Xz,1), size(Xz,2)));
    printpng(fh, fullfile(outDir, 'pca_kmeans_clusters.png')); close(fh);
end

% Cluster composition table (which groups/phases are in each cluster)
clusterComp = table();
for k = 1:Kopt
    m = cidx == k;
    nA = nnz(W.Group(m)=="Active"); nP = nnz(W.Group(m)=="Passive");
    for pi = 1:numel(periods)
        nPh = nnz(Wper(m)==periods(pi));
        clusterComp = [clusterComp; table(k, nA, nP, periods(pi), nPh, ...
            'VariableNames', {'Cluster','nActive','nPassive','Phase','nPhase'})]; %#ok<AGROW>
    end
end
writetable(clusterComp, fullfile(statsDir, 'mod07_cluster_composition.csv'));
fprintf('  K-means cluster composition saved.\n');

% DBSCAN clustering (density-based, no need to specify K)
try
    epsilon = median(pdist(Xz));  % data-driven epsilon
    minPts = max(3, floor(size(Xz,1)/10));
    dbIdx = dbscan(Xz, epsilon, minPts);
    W.DBSCAN_Cluster = dbIdx;
    nClusters = numel(unique(dbIdx(dbIdx > 0)));
    nNoise = nnz(dbIdx == -1);
    fprintf('  DBSCAN: %d clusters, %d noise points (eps=%.2f, minPts=%d).\n', nClusters, nNoise, epsilon, minPts);
catch ME_db
    fprintf('  DBSCAN skipped: %s\n', ME_db.message);
end

catch ME7b
    fprintf('  [WARN] Module 7B (Clustering) failed: %s\n', ME7b.message);
end

% === 7C. EFA with Parallel Analysis ===
try
pFeat = size(Xz,2);
nObs  = size(Xz,1);

% Parallel analysis: compare eigenvalues with random data
nPerm = 100;
R = corrcoef(Xz);
eigReal = sort(eig(R), 'descend');
eigRand = zeros(pFeat, nPerm);
for p = 1:nPerm
    Xrand = randn(nObs, pFeat);
    Rrand = corrcoef(Xrand);
    eigRand(:,p) = sort(eig(Rrand), 'descend');
end
eigRand95 = prctile(eigRand, 95, 2);
nFactorsPA = sum(eigReal > eigRand95);
% Use at least 4 factors (or up to pFeat-1), even if parallel analysis suggests fewer
nFactors = max(4, nFactorsPA);
nFactors = min(nFactors, pFeat-1);
nFactors = min(nFactors, floor(nObs/3));  % must have enough obs per factor

fprintf('  Parallel analysis: %d factors; using %d factors (min 4).\n', nFactorsPA, nFactors);

% Scree plot with parallel analysis
fh = figure('Color','w','Position',[80 80 600 400]); hold on;
plot(1:pFeat, eigReal, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Actual');
plot(1:pFeat, eigRand95, 'r--', 'LineWidth', 1.5, 'DisplayName', '95th pctile random');
yline(1, 'b:', 'Kaiser criterion');
xlabel('Factor'); ylabel('Eigenvalue');
title(sprintf('Scree Plot with Parallel Analysis (%d factors)', nFactors));
legend('Location','northeast'); grid off; box off;
printpng(fh, fullfile(outDir, 'scree_parallel_analysis.png')); close(fh);

% Run EFA (factoran)
try
    [Lambda, Psi, ~, stats] = factoran(Xz, nFactors, 'rotate', 'varimax', 'maxit', 2000);
catch
    try
        Xj = Xz + 1e-5*randn(size(Xz));
        [Lambda, Psi, ~, stats] = factoran(Xj, nFactors, 'rotate', 'varimax', 'maxit', 2000);
    catch ME
        fprintf('  EFA failed: %s\n', ME.message);
        Lambda = []; Psi = [];
    end
end

if ~isempty(Lambda)
    % EFA loadings heatmap
    fhW = max(700, 35*numel(feats)); fhH = max(350, 60*nFactors);
    fh = figure('Color','w','Position',[80 80 fhW fhH]);
    imagesc(Lambda'); colorbar; caxis([-1 1]); colormap(redblue());
    set(gca,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',45, ...
        'YTick',1:nFactors,'YTickLabel',arrayfun(@(i) sprintf('F%d',i),1:nFactors,'Uni',false), ...
        'FontSize', max(6, 10-floor(numel(feats)/8)));
    title(sprintf('EFA Loadings (%d factors, varimax)', nFactors)); grid off; box on;
    printpng(fh, fullfile(outDir, 'efa_loadings_heatmap.png')); close(fh);

    % Factor scores (regression method)
    invPsi = diag(1 ./ max(Psi, 1e-6));
    A = Lambda' * invPsi * Lambda + eye(nFactors);
    Wmat = A \ (Lambda' * invPsi);
    Fscores = (Wmat * Xz')';
    for k = 1:nFactors
        W.(sprintf('F%d',k)) = Fscores(:,k);
    end

    % Communality bar chart (1-Psi for each variable)
    communality = 1 - Psi;
    fh = figure('Color','w','Position',[80 80 max(600,35*numel(feats)) 350]);
    bar(communality, 'FaceColor',[0.4 0.7 0.4],'EdgeColor','w');
    set(gca,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',45, ...
        'FontSize', max(6, 10-floor(numel(feats)/8)));
    ylabel('Communality (1-uniqueness)'); title(sprintf('EFA Communality (%d factors)', nFactors));
    yline(0.5,'r--','LineWidth',1); grid off; box off;
    printpng(fh, fullfile(outDir, 'efa_communality.png')); close(fh);

    % Variance explained per factor
    varPerFactor = sum(Lambda.^2,1) / numel(feats) * 100;
    fprintf('  EFA variance per factor: %s\n', strjoin(arrayfun(@(x) sprintf('%.1f%%',x), varPerFactor, 'Uni', false), ', '));

    % EFA loadings table
    loadTbl = array2table(Lambda, 'VariableNames', ...
        arrayfun(@(i) sprintf('Factor%d',i), 1:nFactors, 'Uni', false));
    loadTbl.Feature = feats(:);
    loadTbl.Communality = communality(:);
    loadTbl = movevars(loadTbl, 'Feature', 'Before', 1);
    writetable(loadTbl, fullfile(statsDir, 'mod07_efa_loadings.csv'));

    % EFA vs PCA comparison (loading correlation)
    fh = figure('Color','w','Position',[80 80 900 400]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile;
    nCompare = min([nFactors, nPC, size(coeff,2)]);  % clamp to available columns
    imagesc(ax1, coeff(:,1:nCompare)'); colorbar(ax1); caxis(ax1,[-1 1]); colormap(ax1, redblue());
    set(ax1,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',35, ...
        'YTick',1:nCompare);
    title(ax1,'PCA Loadings'); grid(ax1,'off'); box(ax1,'on');

    ax2 = nexttile;
    imagesc(ax2, Lambda'); colorbar(ax2); caxis(ax2,[-1 1]); colormap(ax2, redblue());
    set(ax2,'XTick',1:numel(feats),'XTickLabel',feats,'XTickLabelRotation',35, ...
        'YTick',1:nFactors);
    title(ax2,'EFA Loadings (varimax)'); grid(ax2,'off'); box(ax2,'on');
    printpng(fh, fullfile(outDir, 'efa_vs_pca_loadings.png')); close(fh);

    % Factor score trajectory across phases + stats
    % First compute LME p-values per factor
    factorPph = nan(1,nFactors); factorPgr = nan(1,nFactors);
    for fi = 1:nFactors
        fcol = sprintf('F%d', fi);
        Tf = W(:, {'mouse_key','Period','Group', fcol});
        Tf.y = Tf.(fcol);
        Tf.Phase = removecats(categorical(string(Tf.Period)));
        Tf.Grp   = removecats(Tf.Group);
        Tf.mouse = categorical(Tf.mouse_key);
        if height(Tf) < 10, continue; end
        try
            lme = fitlme(Tf, 'y ~ Phase + Grp + (1|mouse)');
            aov = anova(lme);
            factorPph(fi) = safeAnovP(aov,'Phase');
            factorPgr(fi) = safeAnovP(aov,'Grp');
            fprintf('  Factor%d LME: Phase %s, Grp %s\n', fi, fmtP(factorPph(fi)), fmtP(factorPgr(fi)));
        catch, end
    end

    figW = min(1800, max(600, 350*nFactors));
    fh = figure('Color','w','Position',[80 80 figW 460]);
    % Show up to 6 factors in trajectory; if more, split into 2 rows
    nR = ceil(nFactors/4); nCf = ceil(nFactors/nR);
    tiledlayout(nR, nCf, 'TileSpacing','compact','Padding','compact');
    for fi = 1:nFactors
        ax = nexttile; hold(ax,'on');
        fcol = sprintf('F%d',fi);
        muF_grp = struct(); seF_grp = struct();
        for gi = ["Active","Passive"]
            Wg = W(W.Group==gi,:);
            muF = nan(1,numel(periods)); seF = nan(1,numel(periods));
            xoff = 0.05*(gi=="Passive"); c = COL.(char(gi));
            Wgper2 = string(removecats(categorical(Wg.Period)));
            for pi = 1:numel(periods)
                v = Wg.(fcol)(Wgper2==periods(pi));
                v = v(isfinite(v));
                if isempty(v), continue; end
                muF(pi) = mean(v); seF(pi) = std(v)/sqrt(numel(v));
                jit = (rand(numel(v),1)-0.5)*0.08;
                scatter(ax, pi+xoff+jit, v, 14, c, 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility','off');
            end
            muF_grp.(char(gi)) = muF; seF_grp.(char(gi)) = seF;
            errorbar(ax, (1:numel(periods))+xoff, muF, seF, '-o', ...
                'Color',c,'LineWidth',1.8,'MarkerFaceColor',c,'CapSize',5,'DisplayName',char(gi));
        end
        % Per-phase ranksum stars
        for pi = 1:numel(periods)
            vA = W.(fcol)(W.Group=="Active"  & Wper==periods(pi));
            vP = W.(fcol)(W.Group=="Passive" & Wper==periods(pi));
            pRS = safeRanksum(vA, vP);
            addStarsBetweenGroups(ax, pi, muF_grp.Active(pi), muF_grp.Passive(pi), ...
                seF_grp.Active(pi), seF_grp.Passive(pi), pRS);
        end
        set(ax,'XTick',1:numel(periods),'XTickLabel',periods,'XTickLabelRotation',25);
        ylabel(ax, sprintf('Factor %d score',fi));
        title(ax, sprintf('F%d: Phase %s, Grp %s', fi, fmtP(factorPph(fi)), fmtP(factorPgr(fi))), 'FontSize',9);
        grid(ax,'off'); box(ax,'off');
        if fi==nFactors, legend(ax,'Location','northeastoutside'); end
    end
    nAw = nnz(W.Group=="Active"); nPw = nnz(W.Group=="Passive");
    addBottomStatNote(fh, sprintf( ...
        'LME: Factor~Phase+Grp+(1|mouse). Stars: Wilcoxon rank-sum. nA=%d, nP=%d. *p<.05 **p<.01 ***p<.001', nAw, nPw));
    printpng(fh, fullfile(outDir, 'efa_factor_trajectories.png')); close(fh);
end
catch ME7c
    fprintf('  [WARN] Module 7C (EFA) failed: %s\n', ME7c.message);
end

writetable(W, fullfile(statsDir, 'mod07_mouse_period_scores.csv'));
fprintf('  Module 07 done.\n');
end

%% ########################################################################
%  MODULE 08 — Event-Locked Analyses
%  ########################################################################
function module_08_event_locked(csvPath, D, outDir, statsDir)
% Event-locked pupil around lick and reward events from raw frame data

fprintf('  Loading raw CSV for event-locked analysis...\n');
T = readtable(csvPath, 'VariableNamingRule','preserve');
T = ensureStr(T, 'mouse_key');
% Keep mouse_key as string for matching with D.mouse_key

tbCol = pickTimebaseCol(T);
if isempty(tbCol), fprintf('  No time column.\n'); return; end
if ~ismember('Lick_TTL', T.Properties.VariableNames), fprintf('  No Lick_TTL.\n'); return; end
hasPupil = ismember('Diameter_px', T.Properties.VariableNames);
hasReward = ismember('Injector_TTL', T.Properties.VariableNames);

T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;
if hasReward, T.Injector_TTL(isnan(T.Injector_TTL))=0; T.Injector_TTL = T.Injector_TTL > 0.5; end

% Parameters
preWin  = 2;   % seconds before event
postWin = 5;   % seconds after event
dt      = 1/30; % ~30 fps
tAxis   = -preWin:dt:postWin;

[g, km, kd, ks] = findgroups(T.mouse_key, T.day_index, T.session_idx);
nSes = max(g);
fprintf('  Processing %d sessions for event-locked traces...\n', nSes);

% Preallocate result storage
pupilLick = {};  % cell arrays of [1 x nT] traces
pupilRew  = {};

for si = 1:nSes
    idx = (g == si);
    tb  = double(T.(tbCol)(idx));
    ttl = logical(T.Lick_TTL(idx));
    good = isfinite(tb);

    if hasPupil
        pup = double(T.Diameter_px(idx));
        pup = pup(good);
    end
    tb2 = tb(good); ttl2 = ttl(good);

    if numel(tb2) < 100, continue; end

    % Detect lick onsets
    edg = diff([false; ttl2; false]);
    lickOn = find(edg == 1);
    lickTimes = tb2(min(lickOn, numel(tb2)));

    % Detect reward onsets
    rewTimes = [];
    if hasReward
        rew = logical(T.Injector_TTL(idx));
        rew = rew(good);
        edgR = diff([false; rew; false]);
        rewOn = find(edgR == 1);
        rewTimes = tb2(min(rewOn, numel(tb2)));
    end

    mk = string(km(si));
    dayIdx = double(kd(si));
    per = periodOfDay(dayIdx);
    if isundefined(per), continue; end
    grpRows = string(D.mouse_key) == string(km(si));
    if ~any(grpRows), continue; end
    grp = string(D.Group(find(grpRows,1,'first')));

    if hasPupil && numel(lickTimes) >= 3
        traces = extractEventTraces(tb2, pup, lickTimes, tAxis);
        if ~isempty(traces)
            pupilLick(end+1,:) = {mk, dayIdx, char(string(per)), grp, mean(traces,1,'omitnan')}; %#ok<AGROW>
        end
    end

    if hasPupil && numel(rewTimes) >= 1
        traces = extractEventTraces(tb2, pup, rewTimes, tAxis);
        if ~isempty(traces)
            pupilRew(end+1,:) = {mk, dayIdx, char(string(per)), grp, mean(traces,1,'omitnan')}; %#ok<AGROW>
        end
    end
end

COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% --- Plot: Pupil around lick onset by Phase ---
if ~isempty(pupilLick)
    plotEventLockedByPhaseGroup(pupilLick, tAxis, periods, COL, ...
        'Pupil around Lick Onset', 'Pupil (px)', outDir, 'pupil_lick_onset');
end

% --- Plot: Pupil around reward by Phase ---
if ~isempty(pupilRew)
    plotEventLockedByPhaseGroup(pupilRew, tAxis, periods, COL, ...
        'Pupil around Reward Delivery', 'Pupil (px)', outDir, 'pupil_reward');
end

fprintf('  Module 08 done.\n');
end

%% ########################################################################
%  MODULE 09 — Cumulative Curve Fitting
%  ########################################################################
function module_09_cumulative_fit(csvPath, D, outDir, statsDir)
% Fit sigmoid/exponential to cumulative lick count within each session

fprintf('  Loading raw CSV for cumulative fitting...\n');
T = readtable(csvPath, 'VariableNamingRule','preserve');
T = ensureStr(T, 'mouse_key');
% Keep mouse_key as string for matching with D.mouse_key

tbCol = pickTimebaseCol(T);
if isempty(tbCol) || ~ismember('Lick_TTL', T.Properties.VariableNames)
    fprintf('  Missing columns.\n'); return;
end
T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;

[g, km, kd, ks] = findgroups(T.mouse_key, T.day_index, T.session_idx);
nSes = max(g);

fitArr = cell(nSes, 1);
fitIdx = 0;
nSkipFew = 0; nSkipFlat = 0; nFitted = 0; nBadExp = 0; nBadSig = 0;

% Suppress nlinfit warnings (ill-conditioned Jacobian, rank deficiency, iteration limit)
warnIDs = {'stats:nlinfit:IllConditionedJacobian', ...
           'stats:nlinfit:IterationLimitExceeded', ...
           'stats:nlinfit:ModelConstantWRTParam', ...
           'stats:nlinfit:RankDeficient', ...
           'MATLAB:rankDeficientMatrix', ...
           'MATLAB:nearlySingularMatrix'};
prevWarn = cell(numel(warnIDs),1);
for wi = 1:numel(warnIDs)
    prevWarn{wi} = warning('query', warnIDs{wi});
    warning('off', warnIDs{wi});
end
% Blanket suppress all warnings during curve fitting (nlinfit is very chatty)
prevWarnAll = warning('query', 'all');
warning('off', 'all');

for si = 1:nSes
    idx = (g == si);
    tb  = double(T.(tbCol)(idx));
    ttl = logical(T.Lick_TTL(idx));
    good = isfinite(tb);
    tb = tb(good); ttl = ttl(good);
    if numel(tb) < 100, nSkipFew = nSkipFew+1; continue; end

    tb = tb - min(tb);
    edg = diff([false; ttl(:); false]);
    onsets = find(edg == 1);
    if isempty(onsets), nSkipFew = nSkipFew+1; continue; end
    lickTimes = tb(min(onsets, numel(tb)));

    % Quality gate: need >=5 licks to fit a meaningful curve
    nLicks = numel(lickTimes);
    if nLicks < 5, nSkipFew = nSkipFew+1; continue; end

    % Cumulative lick count on a uniform grid
    maxT = max(tb);
    if maxT < 1, nSkipFlat = nSkipFlat+1; continue; end  % session too short
    tGrid = linspace(0, maxT, 200)';
    cumLick = arrayfun(@(t) sum(lickTimes <= t), tGrid);

    Lmax = max(cumLick);
    if Lmax < 5, nSkipFlat = nSkipFlat+1; continue; end  % too few licks for curve

    % Check if curve has any variation (not just a step function at t=0)
    if cumLick(1) >= 0.95*Lmax, nSkipFlat = nSkipFlat+1; continue; end % all licks at start

    opts = statset('MaxIter',1000,'TolFun',1e-8);

    % --- Exponential: y = A*(1 - exp(-k*t)) ---
    % Better initial guess for k: use time to reach 63% of Lmax
    idx63 = find(cumLick >= 0.63*Lmax, 1, 'first');
    if ~isempty(idx63) && tGrid(idx63) > 0
        k0 = 1 / tGrid(idx63);
    else
        k0 = 2 / maxT;
    end
    try
        expModel = @(b,t) b(1) * (1 - exp(-b(2)*t));
        b0 = [Lmax, k0];
        bExp = nlinfit(tGrid, cumLick, expModel, b0, opts);
        yhatExp = expModel(bExp, tGrid);
        R2exp = 1 - sum((cumLick - yhatExp).^2) / sum((cumLick - mean(cumLick)).^2);
        % Sanity: reject if parameters are nonsensical
        if bExp(2) <= 0 || ~isfinite(R2exp), bExp = [NaN NaN]; R2exp = NaN; nBadExp = nBadExp+1; end
    catch
        bExp = [NaN NaN]; R2exp = NaN; nBadExp = nBadExp+1;
    end

    % --- Sigmoid: y = L / (1 + exp(-k*(t - t0))) ---
    % Better initial guess: t0 = time to 50% of Lmax
    idx50 = find(cumLick >= 0.5*Lmax, 1, 'first');
    t0_init = maxT/2;
    if ~isempty(idx50), t0_init = tGrid(idx50); end
    k_sig_init = 4 / maxT;  % rough slope estimate
    try
        sigModel = @(b,t) b(1) ./ (1 + exp(-b(2)*(t - b(3))));
        b0s = [Lmax, k_sig_init, t0_init];
        bSig = nlinfit(tGrid, cumLick, sigModel, b0s, opts);
        yhatSig = sigModel(bSig, tGrid);
        R2sig = 1 - sum((cumLick - yhatSig).^2) / sum((cumLick - mean(cumLick)).^2);
        if bSig(2) <= 0 || ~isfinite(R2sig), bSig = [NaN NaN NaN]; R2sig = NaN; nBadSig = nBadSig+1; end
    catch
        bSig = [NaN NaN NaN]; R2sig = NaN; nBadSig = nBadSig+1;
    end

    nFitted = nFitted + 1;
    fitIdx = fitIdx + 1;
    fitArr{fitIdx} = {char(string(km(si))), double(kd(si)), double(ks(si)), ...
        bExp(1), bExp(2), R2exp, bSig(1), bSig(2), bSig(3), R2sig, Lmax};
end

% Restore warnings
warning(prevWarnAll);
for wi = 1:numel(warnIDs)
    warning(prevWarn{wi}.state, warnIDs{wi});
end

fprintf('  Curve fitting: %d sessions fitted, %d skipped (too few), %d skipped (flat), %d bad exp, %d bad sig.\n', ...
    nFitted, nSkipFew, nSkipFlat, nBadExp, nBadSig);

if fitIdx == 0, fprintf('  No sessions fitted.\n'); return; end

fitTbl = cell2table(vertcat(fitArr{1:fitIdx}), 'VariableNames', ...
    {'mouse_key','day_index','session_idx','exp_A','exp_k','exp_R2', ...
     'sig_L','sig_k','sig_t0','sig_R2','total_licks'});
for v = {'session_idx','exp_A','exp_k','exp_R2','sig_L','sig_k','sig_t0','sig_R2','total_licks'}
    if iscell(fitTbl.(v{1})), fitTbl.(v{1}) = cell2mat(fitTbl.(v{1})); end
end
if iscell(fitTbl.mouse_key), fitTbl.mouse_key = string(fitTbl.mouse_key); end
if iscell(fitTbl.day_index), fitTbl.day_index = cell2mat(fitTbl.day_index); end

% Attach Period & Group using same cohort logic
fitTbl.Period = periodOfDay(fitTbl.day_index);
cohortR = buildNewCohortRoster();
fitTbl = attachCohortGroupAndPair(fitTbl, cohortR);

% --- Fit quality summary ---
nGoodExp = nnz(fitTbl.exp_R2 > 0.5);
nGoodSig = nnz(fitTbl.sig_R2 > 0.5);
fprintf('  Fits with R2 > 0.5:  exponential %d/%d,  sigmoid %d/%d\n', ...
    nGoodExp, height(fitTbl), nGoodSig, height(fitTbl));

% Flag poor fits (R2 < 0.3) by setting parameters to NaN — keeps rows for counting
poorExp = fitTbl.exp_R2 < 0.3 | ~isfinite(fitTbl.exp_R2);
fitTbl.exp_k(poorExp) = NaN; fitTbl.exp_A(poorExp) = NaN;
poorSig = fitTbl.sig_R2 < 0.3 | ~isfinite(fitTbl.sig_R2);
fitTbl.sig_k(poorSig) = NaN; fitTbl.sig_t0(poorSig) = NaN; fitTbl.sig_L(poorSig) = NaN;

COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% Plot: exponential rate (k) by phase/group
fh = figure('Color','w','Position',[80 80 900 450]); hold on;
plotDualBox(fitTbl, 'exp_k', 'Exponential rate k (R^2>0.3 only)', periods, COL);
title('Cumulative Licking: Exponential Rate by Phase');
printpng(fh, fullfile(outDir, 'exp_rate_by_phase.png')); close(fh);

% Plot: sigmoid midpoint (t0) by phase/group
fh = figure('Color','w','Position',[80 80 900 450]); hold on;
plotDualBox(fitTbl, 'sig_t0', 'Sigmoid midpoint t0 (R^2>0.3 only)', periods, COL);
title('Cumulative Licking: Sigmoid Midpoint by Phase');
printpng(fh, fullfile(outDir, 'sig_midpoint_by_phase.png')); close(fh);

writetable(fitTbl, fullfile(statsDir, 'mod09_cumulative_fit_params.csv'));
fprintf('  Module 09 done.\n');
end

%% ########################################################################
%  MODULE 10 — Cross-Modal + Assay Integration
%  ########################################################################
function module_10_crossmodal(D, outDir, statsDir)
COL = groupColors();

% Collect all numeric features + assay columns
behavCols = {'RequirementLast','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
    'bout_n','bout_meanDur_s','rew_freq_per_min','rew_meanDur_s','pupil_mean', ...
    'lick_totalDur_s','bout_totalDur_s','rew_totalDur_s', ...
    'Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min'};
assayCols = {};
if ismember('Immersion_Latency_s', D.Properties.VariableNames), assayCols{end+1}='Immersion_Latency_s'; end
candidates = D.Properties.VariableNames(contains(D.Properties.VariableNames,'TST_Pct') | ...
    contains(D.Properties.VariableNames,'HOT_Pct') | contains(D.Properties.VariableNames,'STRAUB'));
assayCols = [assayCols, candidates];

behavCols = behavCols(ismember(behavCols, D.Properties.VariableNames));
allCols = [behavCols, assayCols];

% Remove columns that have too few values OR near-zero variance
% (TST_Flinching, TST_HindlimbLicking, TST_Jump, HOT_Flinching etc. are mostly NaN or constant)
keepCol = true(1, numel(allCols));
for ci = 1:numel(allCols)
    vals = double(D.(allCols{ci}));
    nFinite = sum(isfinite(vals));
    valsFin = vals(isfinite(vals));
    varVal = 0; if numel(valsFin) > 1, varVal = var(valsFin); end
    if nFinite < 15
        fprintf('  Module 10: dropping %s (only %d finite values).\n', allCols{ci}, nFinite);
        keepCol(ci) = false;
    elseif varVal < 1e-8
        fprintf('  Module 10: dropping %s (near-zero variance).\n', allCols{ci});
        keepCol(ci) = false;
    elseif nFinite < 0.1 * height(D)
        fprintf('  Module 10: dropping %s (%.0f%% missing, below 10%% threshold).\n', ...
            allCols{ci}, 100*(1-nFinite/height(D)));
        keepCol(ci) = false;
    end
end
allCols = allCols(keepCol);

if numel(allCols) < 4
    fprintf('  Skipping Module 10: too few columns (%d).\n', numel(allCols)); return;
end
fprintf('  Module 10: using %d columns for cross-modal analysis.\n', numel(allCols));

% --- 10a. Cross-modal correlation matrix ---
X = nan(height(D), numel(allCols));
for j = 1:numel(allCols), X(:,j) = double(D.(allCols{j})); end
C = corr(X, 'rows','pairwise','type','Spearman');

fh = figure('Color','w','Position',[40 40 max(500, 40*numel(allCols)) max(400, 40*numel(allCols))]);
imagesc(C); colorbar; caxis([-1 1]); colormap(redblue());
set(gca,'XTick',1:numel(allCols),'XTickLabel',allCols,'XTickLabelRotation',40, ...
    'YTick',1:numel(allCols),'YTickLabel',allCols);
title('Cross-Modal Spearman Correlation (All Phases)'); grid off; box on;
printpng(fh, fullfile(outDir, 'crossmodal_correlation.png')); close(fh);

% --- 10a2. Per-phase correlation heatmaps ---
phases = ["Pre","During","Post","Withdrawal","Re-exposure"];
nPh = numel(phases); nC10 = numel(allCols);
figW10 = max(500, 35*nC10); figH10 = max(400, 35*nC10);
for phi = 1:nPh
    phMask = string(D.Period) == phases(phi);
    if sum(phMask) < 5, continue; end
    Xph = nan(sum(phMask), nC10);
    for j = 1:nC10, Xph(:,j) = double(D.(allCols{j})(phMask)); end
    Cph = corr(Xph, 'rows','pairwise','type','Spearman');
    fh = figure('Color','w','Position',[40 40 figW10 figH10]);
    imagesc(Cph); colorbar; caxis([-1 1]); colormap(redblue());
    set(gca,'XTick',1:nC10,'XTickLabel',allCols,'XTickLabelRotation',40, ...
        'YTick',1:nC10,'YTickLabel',allCols,'FontSize',max(5,9-floor(nC10/6)));
    title(sprintf('Spearman Correlation — %s Phase (n=%d)', phases(phi), sum(phMask)));
    grid off; box on;
    printpng(fh, fullfile(outDir, sprintf('correlation_%s.png', lower(char(phases(phi)))))); close(fh);
end
fprintf('  Per-phase correlation heatmaps saved.\n');

% --- 10a3. Per-group correlation heatmaps (Active vs Passive separately) ---
groups10 = ["Active","Passive"];
for gi = 1:numel(groups10)
    grpMask = string(D.Group) == groups10(gi);
    if sum(grpMask) < 5, continue; end
    Xgrp = nan(sum(grpMask), nC10);
    for j = 1:nC10, Xgrp(:,j) = double(D.(allCols{j})(grpMask)); end
    Cgrp = corr(Xgrp, 'rows','pairwise','type','Spearman');
    fh = figure('Color','w','Position',[40 40 figW10 figH10]);
    imagesc(Cgrp); colorbar; caxis([-1 1]); colormap(redblue());
    set(gca,'XTick',1:nC10,'XTickLabel',allCols,'XTickLabelRotation',40, ...
        'YTick',1:nC10,'YTickLabel',allCols,'FontSize',max(5,9-floor(nC10/6)));
    title(sprintf('Spearman Correlation — %s Group (n=%d obs)', groups10(gi), sum(grpMask)));
    grid off; box on;
    printpng(fh, fullfile(outDir, sprintf('correlation_%s.png', lower(char(groups10(gi)))))); close(fh);
end
fprintf('  Per-group correlation heatmaps (Active, Passive) saved.\n');

% --- 10a4. Per-group per-phase correlation heatmaps ---
for gi = 1:numel(groups10)
    grpMask = string(D.Group) == groups10(gi);
    for phi = 1:nPh
        phMask = (string(D.Period) == phases(phi)) & grpMask;
        if sum(phMask) < 5, continue; end
        Xgp = nan(sum(phMask), nC10);
        for j = 1:nC10, Xgp(:,j) = double(D.(allCols{j})(phMask)); end
        Cgp = corr(Xgp, 'rows','pairwise','type','Spearman');
        fh = figure('Color','w','Position',[40 40 figW10 figH10]);
        imagesc(Cgp); colorbar; caxis([-1 1]); colormap(redblue());
        set(gca,'XTick',1:nC10,'XTickLabel',allCols,'XTickLabelRotation',40, ...
            'YTick',1:nC10,'YTickLabel',allCols,'FontSize',max(5,9-floor(nC10/6)));
        title(sprintf('Spearman Correlation — %s, %s Phase (n=%d)', groups10(gi), phases(phi), sum(phMask)));
        grid off; box on;
        printpng(fh, fullfile(outDir, sprintf('correlation_%s_%s.png', lower(char(groups10(gi))), lower(char(phases(phi)))))); close(fh);
    end
end
fprintf('  Per-group per-phase correlation heatmaps saved.\n');

% --- 10b. Partial correlation (controlling for Phase) ---
if numel(allCols) >= 3
    phaseNum = double(D.Period);
    Xp = [X, phaseNum];
    ok = all(isfinite(Xp),2);
    if sum(ok) > 10
        Rpartial = partialcorr(Xp(ok,1:end-1), Xp(ok,end), 'type','Spearman');
        fh = figure('Color','w','Position',[40 40 max(500, 40*numel(allCols)) max(400, 40*numel(allCols))]);
        imagesc(Rpartial); colorbar; caxis([-1 1]); colormap(redblue());
        set(gca,'XTick',1:numel(allCols),'XTickLabel',allCols,'XTickLabelRotation',40, ...
            'YTick',1:numel(allCols),'YTickLabel',allCols);
        title('Partial Correlation (controlling for Phase)'); grid off; box on;
        printpng(fh, fullfile(outDir, 'partial_correlation_ctrl_phase.png')); close(fh);
    end
end

% --- 10c. LME with assay covariates ---
if ~isempty(assayCols) && ismember('RequirementLast', D.Properties.VariableNames)
    for ai = 1:numel(assayCols)
        ac = assayCols{ai};
        T = D(:, {'mouse_key','Period','Group','RequirementLast'});
        T.assay = double(D.(ac));
        T.y = double(T.RequirementLast);
        T = T(isfinite(T.y) & isfinite(T.assay), :);
        T.Phase = removecats(categorical(string(T.Period)));
        T.Grp   = removecats(categorical(string(T.Group)));
        T.mouse = categorical(T.mouse_key);
        if height(T) < 15, continue; end

        try
            lme = fitlme(T, 'y ~ Phase + Grp + assay + (1|mouse)');
            aov = anova(lme);
            fprintf('  LME(PR ~ Phase+Grp+%s): assay p=%.3g\n', ac, safeAnovP(aov,'assay'));
        catch ME
            fprintf('  LME with %s failed: %s\n', ac, ME.message);
        end
    end
end

fprintf('  Module 10 done.\n');
end

%% ########################################################################
%  MODULE 11 — RL Model Scaffold (Actor-Critic TD)
%  ########################################################################
function module_11_rl_model(D, outDir, statsDir)
% Rescorla-Wagner learning model per mouse AND per phase.
% Model: V(t+1) = V(t) + alpha * (R(t) - V(t))
%   alpha = learning rate (how quickly the mouse updates expectations)
%   R(t) = reward on day t (log1p-scaled RequirementLast)
%   V(t) = expected reward (learned value)
%   PE(t) = R(t) - V(t) = prediction error
%
% We fit alpha via maximum likelihood (Gaussian errors on PE).
% We also fit alpha SEPARATELY per phase to see how learning rate changes.

fprintf('  Fitting Rescorla-Wagner model per mouse...\n');
COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

mice = unique(D.mouse_key, 'stable');
if ~ismember('RequirementLast', D.Properties.VariableNames)
    fprintf('  No RequirementLast column. Skipping Module 11.\n'); return;
end

% --- 11A. Global fit: one alpha per mouse (across all days) ---
% IMPORTANT: Skip NaN reward days rather than converting to 0.
% Passive mice have RequirementLast=NaN during "During" phase (no PR task).
% Converting NaN→0 creates artificial reward drops that bias alpha to 0.99.
% Instead, we remove NaN days and fit RW only on days with valid reward data.
% The day indices are preserved for plotting so gaps are visible.
rlRows = {};
dayIdx_per_mouse = {};  % store valid day indices for trajectory plotting
r_per_mouse = {};       % store valid reward sequences
for mi = 1:numel(mice)
    mk = mice(mi);
    Dm = D(D.mouse_key==mk, :);
    Dm = sortrows(Dm, 'day_index');
    grp = string(Dm.Group(1));
    reward = double(Dm.RequirementLast);
    validMask = isfinite(reward) & reward >= 0;
    reward = reward(validMask);
    dayIdxValid = double(Dm.day_index(validMask));
    r = log1p(max(0, reward));
    nT = numel(r);
    if nT < 4, continue; end

    % Fit alpha by minimizing sum of squared prediction errors
    alphaGrid = 0.01:0.01:0.99;
    bestAlpha = 0.1; bestSSE = Inf;
    for ai = 1:numel(alphaGrid)
        [~, pe] = runRW(r, alphaGrid(ai));
        sse = sum(pe.^2);
        if sse < bestSSE, bestSSE = sse; bestAlpha = alphaGrid(ai); end
    end
    [V, pe] = runRW(r, bestAlpha);
    meanPE = mean(abs(pe));
    rlRows(end+1,:) = {char(string(mk)), grp, bestAlpha, bestSSE, meanPE, nT}; %#ok<AGROW>
    dayIdx_per_mouse{end+1} = dayIdxValid;  %#ok<AGROW>
    r_per_mouse{end+1} = r;                 %#ok<AGROW>
end

if isempty(rlRows)
    fprintf('  No mice fitted.\n'); return;
end
rlTbl = cell2table(rlRows, 'VariableNames', {'mouse_key','Group','alpha','SSE','meanAbsPE','N_days'});
for v = {'alpha','SSE','meanAbsPE','N_days'}
    if iscell(rlTbl.(v{1})), rlTbl.(v{1}) = cell2mat(rlTbl.(v{1})); end
end
writetable(rlTbl, fullfile(statsDir, 'mod11_rw_params_global.csv'));

% Plot: alpha by group
fh = figure('Color','w','Position',[80 80 500 430]); hold on;
alphaA = rlTbl.alpha(categorical(rlTbl.Group)=="Active");
alphaP = rlTbl.alpha(categorical(rlTbl.Group)=="Passive");
for gi = ["Active","Passive"]
    vals = rlTbl.alpha(categorical(rlTbl.Group)==gi);
    x = 1 + 0.3*(gi=="Passive");
    c = COL.(char(gi));
    scatter(repmat(x,numel(vals),1), vals, 60, c, 'filled', 'MarkerFaceAlpha',0.7);
    errorbar(x, mean(vals), std(vals)/sqrt(numel(vals)), 'k', 'LineWidth', 1.8, 'CapSize', 10);
end
pAlpha = safeRanksum(alphaA, alphaP);
set(gca,'XTick',[1 1.3],'XTickLabel',{'Active','Passive'},'XLim',[0.6 1.7],'YLim',[0 1]);
ylabel('Learning Rate (\alpha)'); title(sprintf('RW Learning Rate: %s %s', fmtP(pAlpha), starStr(pAlpha)));
grid off; box off;
addBottomStatNote(fh, sprintf(['Rescorla-Wagner: V(t+1)=V(t)+alpha*(R(t)-V(t)). R=log1p(RequirementLast).\n', ...
    'alpha fitted by min SSE over 0.01:0.01:0.99. Ranksum test. nA=%d, nP=%d.'], numel(alphaA), numel(alphaP)));
printpng(fh, fullfile(outDir, 'rw_alpha_by_group.png')); close(fh);

fprintf('  Global alpha: Active=%.2f+/-%.2f, Passive=%.2f+/-%.2f, p=%s\n', ...
    mean(alphaA), std(alphaA), mean(alphaP), std(alphaP), fmtP(pAlpha));

% --- 11B. Per-phase alpha: fit alpha within each phase ---
% Same NaN-skip logic: only use days with valid RequirementLast.
fprintf('  Fitting per-phase alpha...\n');
phaseRows = {};
for mi = 1:numel(mice)
    mk = mice(mi);
    Dm = D(D.mouse_key==mk, :);
    Dm = sortrows(Dm, 'day_index');
    grp = string(Dm.Group(1));
    for pi = 1:numel(periods)
        Dp = Dm(string(Dm.Period)==periods(pi), :);
        reward = double(Dp.RequirementLast);
        validMask = isfinite(reward) & reward >= 0;
        reward = reward(validMask);
        r = log1p(max(0, reward));
        if numel(r) < 2, continue; end
        alphaGrid = 0.01:0.01:0.99;
        bestAlpha = 0.5; bestSSE = Inf;
        for ai = 1:numel(alphaGrid)
            [~, pe] = runRW(r, alphaGrid(ai));
            sse = sum(pe.^2);
            if sse < bestSSE, bestSSE = sse; bestAlpha = alphaGrid(ai); end
        end
        phaseRows(end+1,:) = {char(string(mk)), grp, char(periods(pi)), bestAlpha, bestSSE, numel(r)}; %#ok<AGROW>
    end
end
if ~isempty(phaseRows)
    phaseTbl = cell2table(phaseRows, 'VariableNames', {'mouse_key','Group','Phase','alpha','SSE','N_days'});
    for v = {'alpha','SSE','N_days'}
        if iscell(phaseTbl.(v{1})), phaseTbl.(v{1}) = cell2mat(phaseTbl.(v{1})); end
    end
    writetable(phaseTbl, fullfile(statsDir, 'mod11_rw_params_per_phase.csv'));

    % Plot: alpha trajectory across phases (like effect plot)
    fh = figure('Color','w','Position',[80 80 900 500]); hold on;
    muA = struct(); seA = struct();
    for gi = ["Active","Passive"]
        Tg = phaseTbl(categorical(phaseTbl.Group)==gi, :);
        muV = nan(1,numel(periods)); seV = nan(1,numel(periods));
        c = COL.(char(gi)); xoff = 0.05*(gi=="Passive");
        for pi = 1:numel(periods)
            Tp = Tg(string(Tg.Phase)==periods(pi), :);
            if isempty(Tp), continue; end
            vals = Tp.alpha;
            muV(pi) = mean(vals); seV(pi) = std(vals)/sqrt(numel(vals));
            jit = (rand(numel(vals),1)-0.5)*0.08;
            scatter(pi+xoff+jit, vals, 30, c, 'filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility','off');
        end
        muA.(char(gi)) = muV; seA.(char(gi)) = seV;
        errorbar((1:numel(periods))+xoff, muV, seV, '-o', 'Color', c, ...
            'LineWidth', 1.8, 'MarkerFaceColor', c, 'CapSize', 5, 'DisplayName', char(gi));
    end
    % Ranksum per phase
    for pi = 1:numel(periods)
        vA = phaseTbl.alpha(categorical(phaseTbl.Group)=="Active" & string(phaseTbl.Phase)==periods(pi));
        vP = phaseTbl.alpha(categorical(phaseTbl.Group)=="Passive" & string(phaseTbl.Phase)==periods(pi));
        pRS = safeRanksum(vA, vP);
        addStarsBetweenGroups(gca, pi, muA.Active(pi), muA.Passive(pi), seA.Active(pi), seA.Passive(pi), pRS);
    end
    set(gca,'XTick',1:numel(periods),'XTickLabel',periods,'YLim',[0 1]);
    xlabel('Phase'); ylabel('Learning Rate (\alpha)');
    title('RW Learning Rate Across Phases');
    legend('Location','northeastoutside'); grid off; box off;
    addBottomStatNote(fh, sprintf(['Rescorla-Wagner alpha fitted per mouse per phase (min SSE).\n', ...
        'Stars: Wilcoxon rank-sum. Higher alpha = faster adaptation to reward changes.\n', ...
        'Low alpha in Withdrawal = mouse still expects old reward level (slow extinction).']));
    printpng(fh, fullfile(outDir, 'rw_alpha_per_phase.png')); close(fh);
end

% --- 11C. Value + Prediction Error Trajectory ---
% Uses only valid (non-NaN) reward days. Lines may have gaps where NaN days
% were skipped (e.g., Passive mice during "During" phase).
try
fh = figure('Color','w','Position',[80 80 1000 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% Value trajectory
ax1 = nexttile; hold(ax1,'on');
fittedMice = string(rlTbl.mouse_key);
for mi = 1:numel(fittedMice)
    mk = fittedMice(mi);
    Dm = D(D.mouse_key==mk, :); Dm = sortrows(Dm, 'day_index');
    grp = string(Dm.Group(1));
    reward = double(Dm.RequirementLast);
    validMask = isfinite(reward) & reward >= 0;
    dayIdxV = double(Dm.day_index(validMask));
    r = log1p(max(0, reward(validMask)));
    if numel(r) < 4, continue; end
    bestRow = rlTbl(string(rlTbl.mouse_key)==string(mk), :);
    if isempty(bestRow), continue; end
    [V, ~] = runRW(r, bestRow.alpha);
    c = COL.(char(grp));
    plot(ax1, dayIdxV, V, '-', 'Color', [c 0.5], 'LineWidth', 1.2, 'HandleVisibility','off');
end
xline(ax1, 5.5, 'k:', 'During'); xline(ax1, 10.5, 'k:', 'Post');
xline(ax1, 13.5, 'k:', 'Withdr.'); xline(ax1, 16.5, 'k:', 'Re-exp');
hA = plot(ax1,nan,nan,'-','Color',COL.Active,'LineWidth',2,'DisplayName','Active');
hP = plot(ax1,nan,nan,'-','Color',COL.Passive,'LineWidth',2,'DisplayName','Passive');
legend(ax1,[hA hP],'Location','northeast');
xlabel(ax1,'Day'); ylabel(ax1,'Expected Value V(t)');
title(ax1,'Rescorla-Wagner: Learned Value Trajectory'); grid(ax1,'off'); box(ax1,'off');

% Prediction error trajectory
ax2 = nexttile; hold(ax2,'on');
for mi = 1:numel(fittedMice)
    mk = fittedMice(mi);
    Dm = D(D.mouse_key==mk, :); Dm = sortrows(Dm, 'day_index');
    grp = string(Dm.Group(1));
    reward = double(Dm.RequirementLast);
    validMask = isfinite(reward) & reward >= 0;
    dayIdxV = double(Dm.day_index(validMask));
    r = log1p(max(0, reward(validMask)));
    if numel(r) < 4, continue; end
    bestRow = rlTbl(string(rlTbl.mouse_key)==string(mk), :);
    if isempty(bestRow), continue; end
    [~, pe] = runRW(r, bestRow.alpha);
    c = COL.(char(grp));
    plot(ax2, dayIdxV, pe, '-', 'Color', [c 0.5], 'LineWidth', 1.2, 'HandleVisibility','off');
end
xline(ax2, 5.5, 'k:', 'During'); xline(ax2, 10.5, 'k:', 'Post');
xline(ax2, 13.5, 'k:', 'Withdr.'); xline(ax2, 16.5, 'k:', 'Re-exp');
yline(ax2, 0, 'k--');
xlabel(ax2,'Day'); ylabel(ax2,'Prediction Error PE(t)');
title(ax2,'Rescorla-Wagner: Prediction Error (R-V)'); grid(ax2,'off'); box(ax2,'off');
addBottomStatNote(fh, sprintf(['RW: V(t+1) = V(t) + alpha*(R(t)-V(t)). R = log1p(RequirementLast).\n', ...
    'Top: V(t) = what the mouse expects. Bottom: PE(t) = surprise (positive=better than expected).\n', ...
    'NaN reward days (e.g., Passive During) are skipped — lines show only valid data.\n', ...
    'Large PE at phase transitions = mouse surprised by reward change (e.g., morphine to water).']));
printpng(fh, fullfile(outDir, 'rw_value_pe_trajectory.png')); close(fh);
catch ME_vt
    fprintf('  [WARN] Value/PE trajectory plot failed: %s\n', ME_vt.message);
end

% --- Delete stale TD model plots from previous code version ---
staleFiles = {'rl_value_trajectory.png', 'rl_params_by_group.png'};
for si = 1:numel(staleFiles)
    stalePath = fullfile(outDir, staleFiles{si});
    if isfile(stalePath)
        delete(stalePath);
        fprintf('  Deleted stale TD-model plot: %s\n', staleFiles{si});
    end
end

fprintf('  Module 11 done.\n');
end

%% ########################################################################
%  MODULE 12 — Predictive Modeling (Decoder / Prediction / Cross-Validation)
%  ########################################################################
function module_12_predictive(D, outDir, statsDir)
% Decoders that predict phase / reward type / group from behavioral features.
% Also: predict re-exposure behavior from earlier phases.
% Uses nested leave-one-mouse-out CV + permutation tests.

COL = groupColors();
periods = ["Pre","During","Post","Withdrawal","Re-exposure"];

% --- Feature selection: use all available numeric features ---
allFeats = {'RequirementLast','lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
    'lick_totalDur_s','bout_n','bout_meanDur_s','bout_totalDur_s', ...
    'rew_freq_per_min','rew_meanDur_s','rew_totalDur_s', ...
    'pupil_mean','Requirement_cum','Requirement_speed_per_day','Requirement_speed_per_min'};
allFeats = allFeats(ismember(allFeats, D.Properties.VariableNames));

% Filter out features with too few values
keepF = true(1, numel(allFeats));
for fi = 1:numel(allFeats)
    if sum(isfinite(double(D.(allFeats{fi})))) < 20, keepF(fi) = false; end
end
allFeats = allFeats(keepF);
fprintf('  Module 12: %d features for predictive modeling.\n', numel(allFeats));

if numel(allFeats) < 3
    fprintf('  Skipping Module 12: too few features.\n'); return;
end

% Build feature matrix X (per-mouse-period medians)
mice = unique(D.mouse_key, 'stable');
rows = {};
for i = 1:numel(mice)
    Dm = D(D.mouse_key == mice(i), :);
    gVal = Dm.Group(1);
    for pi = 1:numel(periods)
        Dp = Dm(string(Dm.Period) == periods(pi), :);
        if isempty(Dp), continue; end
        row.mouse = mice(i);
        row.Group = string(gVal);
        row.Period = periods(pi);
        row.RewardType = string(rewardTypeOfDay(double(Dp.day_index(1))));
        for fi = 1:numel(allFeats)
            row.(allFeats{fi}) = median(double(Dp.(allFeats{fi})),'omitnan');
        end
        rows{end+1} = row; %#ok<AGROW>
    end
end
W = struct2table([rows{:}]);
W.mouse = string(W.mouse);
W.Group = categorical(W.Group, {'Passive','Active'});
W.Period = categorical(W.Period, cellstr(periods));

fprintf('  Built predictive table: %d rows (mouse x period).\n', height(W));

% Build numeric matrix
X = nan(height(W), numel(allFeats));
for fi = 1:numel(allFeats)
    vals = W.(allFeats{fi});
    if iscell(vals), X(:,fi) = cellfun(@double, vals); else, X(:,fi) = double(vals); end
end
% Impute NaN with column median
for fi = 1:size(X,2)
    nans = isnan(X(:,fi));
    if any(nans), X(nans,fi) = median(X(~nans,fi),'omitnan'); end
end
% Z-score
Xz = (X - mean(X,1,'omitnan')) ./ max(std(X,0,1,'omitnan'), 1e-10);

% === 12A. Phase Decoder (Predict Phase from features) ===
fprintf('  12A: Phase decoder (LOMO-CV)...\n');
yPhase = double(W.Period);
nPerm = 100;
try
[accPhase, accPhase_perm, confMat_phase] = lomoDecoder(Xz, yPhase, W.mouse, nPerm);
fprintf('    Phase decoder: accuracy = %.1f%% (permutation p = %.4f)\n', accPhase*100, mean(accPhase_perm >= accPhase));

% Confusion matrix — colored, with counts + percentages per cell
fh = figure('Color','w','Position',[80 80 620 520]);
nP5 = numel(periods);
% Normalize by row (true class) to get percentages
rowSums = sum(confMat_phase, 2);
confPct = confMat_phase ./ max(rowSums, 1) * 100;
imagesc(confPct); colorbar;
colormap([ ...
    linspace(1,0.12,64)', linspace(1,0.47,64)', linspace(1,0.71,64)']);  % white→blue
caxis([0 100]);
pLabels = cellstr(periods);
% Add per-class accuracy on Y-axis labels
classAcc = diag(confMat_phase) ./ max(rowSums, 1) * 100;
yLabelsExt = cell(nP5,1);
for k = 1:nP5
    yLabelsExt{k} = sprintf('%s (%.0f%%)', pLabels{k}, classAcc(k));
end
set(gca,'XTick',1:nP5,'XTickLabel',pLabels,'XTickLabelRotation',30, ...
    'YTick',1:nP5,'YTickLabel',yLabelsExt,'FontSize',9);
xlabel('Predicted Phase'); ylabel('True Phase (per-class accuracy)');
title(sprintf('Phase Decoder — Accuracy %.1f%%', accPhase*100));
% Show count + percentage in each cell
for r = 1:nP5
    for c = 1:nP5
        cnt = confMat_phase(r,c);
        pct = confPct(r,c);
        % Use white text on dark cells, black on light
        if pct > 50, tc = 'w'; else, tc = 'k'; end
        text(c, r, sprintf('%d\n(%.0f%%)', cnt, pct), ...
            'HorizontalAlignment','center','FontSize',9, ...
            'FontWeight','bold','Color',tc);
    end
end
grid off; box on;
pPerm = mean(accPhase_perm >= accPhase);
addBottomStatNote(fh, sprintf(['PHASE DECODER: Can we tell which phase a mouse is in from behavior alone?\n', ...
    'Method: Leave-One-Mouse-Out CV with LDA classifier. %d obs, %d features.\n', ...
    'Accuracy=%.1f%% (chance=%.1f%%). Permutation p=%.4f (%d perm).\n', ...
    'Each cell: count (row %%). Diagonal = correct. Color intensity = row-normalized accuracy.'], ...
    height(W), numel(allFeats), accPhase*100, 100/numel(unique(yPhase)), pPerm, nPerm));
printpng(fh, fullfile(outDir, 'confusion_phase_decoder.png')); close(fh);

% Permutation histogram
fh = figure('Color','w','Position',[80 80 500 350]); hold on;
histogram(accPhase_perm*100, 20, 'FaceColor',[0.7 0.7 0.7],'EdgeColor','w');
xline(accPhase*100, 'r-', 'LineWidth', 2);
xlabel('Accuracy (%)'); ylabel('Count');
title(sprintf('Permutation Test — Phase Decoder (p=%.4f)', pPerm));
grid off; box off;
addBottomStatNote(fh, sprintf(['Gray bars=accuracy when labels are randomly shuffled (null distribution).\n', ...
    'Red line=actual accuracy (%.1f%%). If red line is far right of gray bars, behavior truly differs by phase.\n', ...
    'p=%.4f means <%.1f%% of shuffled data beat the real accuracy.'], accPhase*100, pPerm, pPerm*100));
printpng(fh, fullfile(outDir, 'permutation_phase_decoder.png')); close(fh);
catch ME12a
    fprintf('  [WARN] 12A Phase decoder failed: %s\n', ME12a.message);
end

% === 12B. Group Decoder (Predict Active vs Passive from features) ===
fprintf('  12B: Group decoder (LOMO-CV)...\n');
yGroup = double(W.Group);
try
[accGroup, accGroup_perm, confMat_group] = lomoDecoder(Xz, yGroup, W.mouse, nPerm);
fprintf('    Group decoder: accuracy = %.1f%% (permutation p = %.4f)\n', accGroup*100, mean(accGroup_perm >= accGroup));

fh = figure('Color','w','Position',[80 80 500 430]);
nG2 = size(confMat_group,1);
rowSumsG = sum(confMat_group, 2);
confPctG = confMat_group ./ max(rowSumsG, 1) * 100;
imagesc(confPctG); colorbar;
colormap([ ...
    linspace(1,0.12,64)', linspace(1,0.47,64)', linspace(1,0.71,64)']);  % white→blue
caxis([0 100]);
gLabels = {'Passive','Active'};
classAccG = diag(confMat_group) ./ max(rowSumsG, 1) * 100;
yLabelsG = cell(nG2,1);
for k = 1:nG2
    yLabelsG{k} = sprintf('%s (%.0f%%)', gLabels{k}, classAccG(k));
end
set(gca,'XTick',1:nG2,'XTickLabel',gLabels,'YTick',1:nG2,'YTickLabel',yLabelsG,'FontSize',11);
xlabel('Predicted Group'); ylabel('True Group (per-class accuracy)');
title(sprintf('Group Decoder — Accuracy %.1f%%', accGroup*100));
% Show count + pct in each cell, with meaning labels
cellMeaning = {'True Negative\n(Correct Passive)', 'False Positive\n(Passive→Active)'; ...
               'False Negative\n(Active→Passive)', 'True Positive\n(Correct Active)'};
for r = 1:nG2
    for c = 1:nG2
        cnt = confMat_group(r,c);
        pct = confPctG(r,c);
        if pct > 50, tc = 'w'; else, tc = 'k'; end
        text(c, r-0.15, sprintf('%d (%.0f%%)', cnt, pct), ...
            'HorizontalAlignment','center','FontSize',12, ...
            'FontWeight','bold','Color',tc);
        text(c, r+0.2, sprintf(cellMeaning{r,c}), ...
            'HorizontalAlignment','center','FontSize',7, ...
            'Color',tc,'FontAngle','italic');
    end
end
grid off; box on;
pPermG = mean(accGroup_perm >= accGroup);
addBottomStatNote(fh, sprintf(['GROUP DECODER: Can we tell Active vs Passive from behavior alone?\n', ...
    'Method: Leave-One-Mouse-Out CV with LDA. %d obs, %d features.\n', ...
    'Accuracy=%.1f%% (chance=50%%). Permutation p=%.4f.\n', ...
    'Each cell: count (row %%). Diagonal = correct classification.'], ...
    height(W), numel(allFeats), accGroup*100, pPermG));
printpng(fh, fullfile(outDir, 'confusion_group_decoder.png')); close(fh);

fh = figure('Color','w','Position',[80 80 500 350]); hold on;
histogram(accGroup_perm*100, 20, 'FaceColor',[0.7 0.7 0.7],'EdgeColor','w');
xline(accGroup*100, 'r-', 'LineWidth', 2);
xlabel('Accuracy (%)'); ylabel('Count');
title(sprintf('Permutation Test — Group Decoder (p=%.4f)', pPermG));
grid off; box off;
addBottomStatNote(fh, sprintf(['Gray=null distribution (shuffled group labels). Red=actual accuracy (%.1f%%).\n', ...
    'Red line far right means Active vs Passive groups have distinctly different behavior patterns.'], accGroup*100));
printpng(fh, fullfile(outDir, 'permutation_group_decoder.png')); close(fh);
catch ME12b
    fprintf('  [WARN] 12B Group decoder failed: %s\n', ME12b.message);
end

% === 12C. Predict Re-exposure Behavior from Earlier Phases ===
fprintf('  12C: Predict re-exposure from earlier phases...\n');
try
ycolTarget = 'RequirementLast';
if ismember(ycolTarget, allFeats)
    Wre = W(W.Period == "Re-exposure", :);
    if height(Wre) >= 5
        yRe = double(Wre.(ycolTarget));
        yRe(~isfinite(yRe)) = median(yRe,'omitnan');

        % Features: use each mouse's Pre + Post + Withdrawal medians as predictors
        predPeriods = ["Pre","Post","Withdrawal"];
        Xpred = nan(height(Wre), numel(allFeats)*numel(predPeriods));
        predNames = {};
        for pp = 1:numel(predPeriods)
            Wpp = W(W.Period == predPeriods(pp), :);
            for fi = 1:numel(allFeats)
                colIdx = (pp-1)*numel(allFeats) + fi;
                predNames{colIdx} = sprintf('%s_%s', char(predPeriods(pp)), allFeats{fi});
                for mi = 1:height(Wre)
                    mk = Wre.mouse(mi);
                    row = Wpp(Wpp.mouse == mk, :);
                    if ~isempty(row)
                        val = row.(allFeats{fi});
                        if iscell(val), val = val{1}; end
                        Xpred(mi, colIdx) = double(val);
                    end
                end
            end
        end
        % Impute NaN
        for fi = 1:size(Xpred,2)
            nans = isnan(Xpred(:,fi));
            if any(nans) && any(~nans), Xpred(nans,fi) = median(Xpred(~nans,fi),'omitnan'); end
        end
        % Remove near-constant columns
        keepC = std(Xpred,0,1,'omitnan') > 1e-8;
        Xpred = Xpred(:, keepC);
        predNames = predNames(keepC);

        % LOMO regression: predict yRe from earlier features
        yPredAll = nan(numel(yRe),1);
        miceRe = Wre.mouse;
        uMice = unique(miceRe);
        for mi = 1:numel(uMice)
            test = miceRe == uMice(mi);
            train = ~test;
            if sum(train) < 3, continue; end
            Xtrain = Xpred(train,:); yTrain = yRe(train);
            Xtest  = Xpred(test,:);
            % Ridge regression (regularized)
            lambda = 1;
            beta = (Xtrain'*Xtrain + lambda*eye(size(Xtrain,2))) \ (Xtrain'*yTrain);
            yPredAll(test) = Xtest * beta;
        end
        ok = isfinite(yPredAll) & isfinite(yRe);
        if sum(ok) >= 5
            r2 = 1 - sum((yRe(ok)-yPredAll(ok)).^2) / sum((yRe(ok)-mean(yRe(ok))).^2);
            [rho, pRho] = corr(yRe(ok), yPredAll(ok), 'type','Spearman');

            % Permutation test for R2
            nPerm2 = 200;
            r2perm = nan(nPerm2,1);
            for pp = 1:nPerm2
                yShuf = yRe(randperm(numel(yRe)));
                yPredShuf = nan(numel(yShuf),1);
                for mi = 1:numel(uMice)
                    test = miceRe == uMice(mi); train = ~test;
                    if sum(train)<3, continue; end
                    beta = (Xpred(train,:)'*Xpred(train,:) + lambda*eye(size(Xpred,2))) \ (Xpred(train,:)'*yShuf(train));
                    yPredShuf(test) = Xpred(test,:) * beta;
                end
                ok2 = isfinite(yPredShuf) & isfinite(yShuf);
                if sum(ok2)>=3
                    r2perm(pp) = 1 - sum((yShuf(ok2)-yPredShuf(ok2)).^2)/sum((yShuf(ok2)-mean(yShuf(ok2))).^2);
                end
            end
            pR2 = mean(r2perm >= r2, 'omitnan');
            fprintf('    Re-exposure prediction: R2=%.3f, rho=%.3f (p=%.4f), perm p=%.4f\n', r2, rho, pRho, pR2);

            % ── Enhanced Re-exposure Prediction Plots ──
            grp = string(Wre.Group);
            residuals = yRe(ok) - yPredAll(ok);
            mape = mean(abs(residuals) ./ max(abs(yRe(ok)), 1)) * 100;

            % --- PLOT 1: Scatter with mouse labels + residual lines ---
            fh = figure('Color','w','Position',[80 80 700 550]); hold on;
            mn = min([yRe;yPredAll])*0.9; mx = max([yRe;yPredAll])*1.05;
            % Shaded ±5 band around diagonal (good prediction zone)
            bandW = 5;
            fill([mn mx mx mn], [mn-bandW mx-bandW mx+bandW mn+bandW], ...
                [0.9 0.95 0.9], 'EdgeColor','none','FaceAlpha',0.4, ...
                'HandleVisibility','off');
            plot([mn mx],[mn mx],'k--','LineWidth',1.2,'HandleVisibility','off');
            % Residual lines (vertical drop from point to diagonal)
            for i = 1:numel(yRe)
                lc = [0.7 0.7 0.7];
                plot([yRe(i) yRe(i)], [yRe(i) yPredAll(i)], '-', 'Color', lc, ...
                    'LineWidth', 0.8, 'HandleVisibility','off');
            end
            % Scatter points
            scatter(yRe(grp=="Active"), yPredAll(grp=="Active"), 90, COL.Active, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.5,'DisplayName','Active');
            scatter(yRe(grp=="Passive"), yPredAll(grp=="Passive"), 90, COL.Passive, 'filled', ...
                'MarkerEdgeColor','k','LineWidth',0.5,'DisplayName','Passive');
            % Mouse name labels
            for i = 1:height(Wre)
                mkLbl = char(Wre.mouse(i));
                mkLbl = regexprep(mkLbl, '(\d{4})_?(orange)', '$1_org');
                mkLbl = regexprep(mkLbl, '(\d{4})_?(black)', '$1_blk');
                mkLbl = regexprep(mkLbl, '(\d{4})_?(white)', '$1_wht');
                text(yRe(i)+0.5, yPredAll(i)+0.5, mkLbl, 'FontSize', 7, ...
                    'Color', [0.3 0.3 0.3], 'Clipping','on');
            end
            xlim([mn mx]); ylim([mn mx]);
            xlabel('Actual Re-exposure PR Breakpoint');
            ylabel('Predicted Re-exposure PR Breakpoint');
            title(sprintf('Predict Re-exposure from Earlier Phases (R^2=%.2f, \\rho=%.2f)', r2, rho));
            legend('Location','northwest');
            grid off; box off;
            % Add text box with summary stats
            statTxt = sprintf('R^2 = %.3f\n\\rho = %.3f (p=%.4f)\nPerm p = %.4f\nMAPE = %.1f%%', ...
                r2, rho, pRho, pR2, mape);
            text(mx*0.98, mn+2, statTxt, 'FontSize', 8, 'HorizontalAlignment','right', ...
                'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Margin',3);
            addBottomStatNote(fh, sprintf( ...
                ['RE-EXPOSURE PREDICTION: Can earlier behavior predict relapse-like PR breakpoint?\n', ...
                'Green band = ±5 units from perfect prediction. Gray lines = prediction error per mouse.\n', ...
                'R2=%.3f (%.0f%% variance explained). Spearman rho=%.3f (p=%.4f). Perm p=%.4f (%d perm).\n', ...
                'Points near diagonal = well predicted. Labels = mouse ID.'], ...
                r2, r2*100, rho, pRho, pR2, nPerm2));
            printpng(fh, fullfile(outDir, 'predict_reexposure_scatter.png')); close(fh);

            % --- PLOT 2: Side-by-side bar (Actual vs Predicted per mouse) ---
            fh = figure('Color','w','Position',[80 80 900 450]); hold on;
            nMiceRe = height(Wre);
            mouseLabels = cell(nMiceRe,1);
            for i = 1:nMiceRe
                ml = char(Wre.mouse(i));
                ml = regexprep(ml, '(\d{4})_?(orange)', '$1_org');
                ml = regexprep(ml, '(\d{4})_?(black)', '$1_blk');
                ml = regexprep(ml, '(\d{4})_?(white)', '$1_wht');
                mouseLabels{i} = ml;
            end
            % Sort by actual value for cleaner visual
            [~, sortI] = sort(yRe, 'descend');
            barData = [yRe(sortI), yPredAll(sortI)];
            bh = barh(1:nMiceRe, barData, 'grouped');
            bh(1).FaceColor = [0.2 0.6 0.85]; bh(1).DisplayName = 'Actual PR';
            bh(2).FaceColor = [0.95 0.5 0.2]; bh(2).DisplayName = 'Predicted PR';
            % Color-code mouse labels by group
            set(gca, 'YTick', 1:nMiceRe, 'YTickLabel', mouseLabels(sortI), 'FontSize', 8);
            for i = 1:nMiceRe
                gi = sortI(i);
                if string(Wre.Group(gi)) == "Active"
                    clr = COL.Active;
                else
                    clr = COL.Passive;
                end
                % Add error annotation
                errVal = yPredAll(gi) - yRe(gi);
                maxBar = max(barData(i,:));
                text(maxBar + 0.5, i, sprintf('%+.1f', errVal), ...
                    'FontSize', 7, 'Color', [0.4 0.4 0.4]);
            end
            xlabel('PR Breakpoint (RequirementLast)');
            ylabel('Mouse');
            title(sprintf('Re-exposure: Actual vs Predicted (R^2=%.2f)', r2));
            legend('Location','southeast'); grid off; box off;
            addBottomStatNote(fh, sprintf( ...
                ['Each mouse: blue=actual Re-exposure PR, orange=predicted from earlier phases.\n', ...
                'Numbers on right = prediction error (Predicted - Actual). Sorted by actual PR.\n', ...
                'R2=%.3f. Close bars = good prediction. Large gaps = hard-to-predict mice.'], r2));
            printpng(fh, fullfile(outDir, 'predict_reexposure_bars.png')); close(fh);

            % Permutation histogram
            fh = figure('Color','w','Position',[80 80 500 350]); hold on;
            histogram(r2perm, 20, 'FaceColor',[0.7 0.7 0.7],'EdgeColor','w');
            xline(r2, 'r-', 'LineWidth', 2);
            xlabel('R^2'); ylabel('Count');
            title(sprintf('Permutation Test — Re-exposure Prediction (p=%.4f)', pR2));
            grid off; box off;
            addBottomStatNote(fh, sprintf(['Gray=R2 from shuffled outcomes. Red=actual R2 (%.3f).\n', ...
                'Most shuffled R2 are negative (random prediction is worse than guessing the mean).\n', ...
                'p=%.4f: if <0.05, earlier behavior significantly predicts re-exposure relapse.'], r2, pR2));
            printpng(fh, fullfile(outDir, 'permutation_reexposure.png')); close(fh);
        end
    end
end
catch ME12c
    fprintf('  [WARN] 12C Re-exposure prediction failed: %s\n', ME12c.message);
end

% === 12D. Feature Importance (via permutation importance) ===
fprintf('  12D: Feature importance for phase decoder...\n');
try
if exist('accPhase','var') && ~isnan(accPhase)
    baseAcc = accPhase;
    featImp = nan(1, numel(allFeats));
    for fi = 1:numel(allFeats)
        XzShuf = Xz;
        XzShuf(:,fi) = XzShuf(randperm(size(XzShuf,1)), fi);
        [accShuf, ~, ~] = lomoDecoder(XzShuf, yPhase, W.mouse, 0);
        featImp(fi) = baseAcc - accShuf;  % drop in accuracy when shuffled
    end
    [~, sortIdx] = sort(featImp, 'descend');

    fh = figure('Color','w','Position',[80 80 700 max(300, 25*numel(allFeats))]); hold on;
    barh(1:numel(allFeats), featImp(sortIdx), 'FaceColor', [0.3 0.6 0.8]);
    set(gca,'YTick',1:numel(allFeats),'YTickLabel',allFeats(sortIdx),'FontSize',8);
    xlabel('Drop in Accuracy (permutation importance)');
    title('Feature Importance — Phase Decoder');
    grid off; box off;
    addBottomStatNote(fh, sprintf(['FEATURE IMPORTANCE: Which behaviors matter most for distinguishing phases?\n', ...
        'Method: Shuffle one feature at a time, measure accuracy drop. Bigger bar = more important.\n', ...
        'Base accuracy=%.1f%%. Top features drive the phase decoder. Negative = noise (removing helps).\n', ...
        'e.g., lick_freq_per_min at top = lick rate is the most informative behavioral signal for phase.'], baseAcc*100));
    printpng(fh, fullfile(outDir, 'feature_importance_phase.png')); close(fh);

    impTbl = table(allFeats(:), featImp(:), 'VariableNames', {'Feature','PermImportance'});
    impTbl = sortrows(impTbl, 'PermImportance', 'descend');
    writetable(impTbl, fullfile(statsDir, 'mod12_feature_importance.csv'));
end
catch ME12d
    fprintf('  [WARN] 12D Feature importance failed: %s\n', ME12d.message);
end

% Save summary
summTbl = table();
try summTbl.PhaseAcc = accPhase*100; catch, summTbl.PhaseAcc = NaN; end
try summTbl.GroupAcc = accGroup*100; catch, summTbl.GroupAcc = NaN; end
try summTbl.ReexpR2  = r2; catch, summTbl.ReexpR2 = NaN; end
writetable(summTbl, fullfile(statsDir, 'mod12_predictive_summary.csv'));
fprintf('  Module 12 done.\n');
end

%% ########################################################################
%  Module 11 & 12 — Korean Explanation Text File Writer
%  ########################################################################
function write_module_11_12_explanation(baseOut)
% Writes a plain-text explanation of all Module 11 & 12 plots in Korean.
% No external toolbox dependencies — just fprintf to file.

outFile = fullfile(baseOut, 'explain_module_11_12.txt');
fid = fopen(outFile, 'w', 'n', 'UTF-8');
if fid < 0
    fprintf('  [WARN] Cannot open %s for writing.\n', outFile);
    return;
end

w = @(s) fprintf(fid, '%s\n', s);
wn = @() fprintf(fid, '\n');

w('================================================================');
w('  Module 11 & 12 그래프 해설서 (한국어)');
w('================================================================');
w(sprintf('  자동 생성일: %s', datestr(now, 'yyyy-mm-dd HH:MM')));
w('  본 문서는 Module 11(강화학습)과 Module 12(예측 모델링)의');
w('  모든 그래프에 대한 상세 설명을 한국어로 제공합니다.');
wn();

w('================================================================');
w('  MODULE 11: Rescorla-Wagner (RW) 강화학습 모델');
w('================================================================');
wn();
w('모델 공식: V(t+1) = V(t) + alpha * (R(t) - V(t))');
w('  - V(t): 마우스가 t시점에 기대하는 보상의 크기');
w('  - R(t): 실제 받은 보상 = log1p(RequirementLast)');
w('  - alpha: 학습률 (0~1)');
w('  - alpha가 높으면 → 새 보상에 빠르게 적응');
w('  - alpha가 낮으면 → 과거 경험에 더 의존 (느린 업데이트)');
wn();
w('[NaN 처리]');
w('Passive 마우스는 During 시기에 PR task를 하지 않으므로');
w('RequirementLast = NaN. 이 날들은 건너뛰고 유효한 보상');
w('데이터만 사용하여 alpha를 fit합니다.');
wn();

w('----------------------------------------------------------------');
w('  그래프 1: rw_alpha_by_group.png');
w('  RW Learning Rate — 그룹별 비교');
w('----------------------------------------------------------------');
w('  X축: Active / Passive 그룹');
w('  Y축: 학습률 alpha (0~1)');
w('  각 점: 개별 마우스의 학습률 (전체 실험 기간에 걸쳐 하나의 alpha fit)');
w('  검은 에러바: 그룹 평균 ± 표준오차 (SEM)');
w('  제목의 p값: Wilcoxon rank-sum 검정 결과');
wn();
w('[해석]');
w('  - Active alpha 분산이 크면 → 자가투여 경험이 마우스마다');
w('    다른 학습 전략을 유도했음');
w('  - p값 유의 → 두 그룹의 보상 학습 속도가 통계적으로 다름');
w('  - 방법: alpha를 0.01~0.99에서 grid search (min SSE)');
wn();

w('----------------------------------------------------------------');
w('  그래프 2: rw_alpha_per_phase.png');
w('  RW Learning Rate Across Phases (시기별 학습률 변화)');
w('----------------------------------------------------------------');
w('  X축: 실험 시기 (Pre → During → Post → Withdrawal → Re-exposure)');
w('  Y축: 학습률 alpha (0~1)');
w('  빨간 선/점: Active 그룹 (평균 ± SEM), 연한 점 = 개별 마우스');
w('  파란 선/점: Passive 그룹 (평균 ± SEM), 연한 점 = 개별 마우스');
w('  별표: 해당 시기에서 두 그룹 간 Wilcoxon rank-sum 결과');
w('        (n.s. = 유의하지 않음, * p<0.05, ** p<0.01)');
wn();
w('[시기별 해석]');
w('  Pre: 베이스라인. 물만 제공되는 실험 전 기간.');
w('  During: 모르핀 투여. alpha 낮아지면 → 안정적 보상에 적응 완료.');
w('  Post: 모르핀→물 전환. alpha 증가 → 보상 변화에 빠르게 재적응.');
w('  Withdrawal: 금단. alpha 낮으면 → 여전히 이전 보상을 기대 (느린 소거).');
w('  Re-exposure: 재노출. alpha 변화는 재발(relapse) 관련.');
wn();

w('----------------------------------------------------------------');
w('  그래프 3: rw_value_pe_trajectory.png');
w('  Rescorla-Wagner Value & Prediction Error Trajectory');
w('----------------------------------------------------------------');
w('  [상단: V(t) — 학습된 기대 가치]');
w('  X축: 실험일 (Day)');
w('  Y축: V(t) — 마우스가 기대하는 보상 크기');
w('  빨간 선: Active 마우스 (개별), 파란 선: Passive 마우스 (개별)');
w('  점선: 시기 전환점 (During, Post, Withdrawal, Re-exposure)');
wn();
w('  [하단: PE(t) — 예측 오차]');
w('  Y축: PE(t) = R(t) - V(t) = 실제 보상 - 기대 보상');
w('  PE > 0: 기대보다 좋은 보상 (긍정적 놀라움)');
w('  PE < 0: 기대보다 나쁜 보상 (부정적 놀라움)');
w('  PE ≈ 0: 보상이 기대와 일치 (완벽한 예측)');
wn();
w('[해석]');
w('  시기 전환점에서 큰 PE 스파이크가 관찰됨.');
w('  During→Post 전환에서 큰 음의 PE = 모르핀 기대했으나 물을 받음 ("실망").');
w('  Active가 더 큰 PE → 자가투여가 더 강한 보상 기대를 형성했음.');
w('  NaN 보상일(Passive During)은 건너뛰므로, 선에 갭이 있을 수 있음.');
wn();

w('================================================================');
w('  MODULE 12: 예측 모델링 (Predictive Modeling)');
w('================================================================');
wn();
w('핵심 방법: Leave-One-Mouse-Out CV (LOMO-CV)');
w('  한 마우스를 빼고 나머지로 모델 훈련 → 빠진 마우스 예측.');
w('  모든 마우스에 대해 반복하여 정확도 측정.');
w('  해당 마우스의 데이터를 한 번도 본 적 없는 모델이 예측 수행.');
w('사용 특성: RequirementLast, lick 관련, bout 관련, reward 관련,');
w('  pupil, Requirement 누적/속도 등 15개 행동 특성.');
wn();

w('----------------------------------------------------------------');
w('  그래프 4: confusion_phase_decoder.png');
w('  Phase Decoder — 혼동 행렬 (Confusion Matrix)');
w('----------------------------------------------------------------');
w('  행(Y축): 실제 시기 (True Phase)');
w('  열(X축): 모델이 예측한 시기 (Predicted Phase)');
w('  각 셀의 숫자: 횟수 (행 백분율)');
w('    예: 9 (64%) → 실제 Pre인 관찰치 중 64%를 Pre로 정확히 맞춤');
w('  대각선: 정확한 분류 (진한 파랑). 숫자가 클수록 좋음.');
w('  비대각선: 오분류. 어떤 시기를 어떤 시기와 혼동하는지 보여줌.');
w('  Y축 괄호: 해당 시기의 개별 정확도');
w('  색상: 흰색(0%) → 진한 파랑(100%). 행 단위로 정규화.');
wn();
w('[해석]');
w('  우연 수준 = 20% (5개 시기). 정확도 > 20% → 행동이 시기마다 다름.');
w('  특정 시기가 잘 분류 → 그 시기의 행동이 뚜렷이 구별됨.');
w('  두 시기가 자주 혼동 → 행동 패턴이 유사함.');
w('  예: Pre↔During 혼동 → 모르핀 투여가 행동을 크게 안 바꿨을 수 있음.');
w('  예: Withdrawal↔Re-exposure 혼동 → 금단과 재노출의 행동이 유사.');
wn();

w('----------------------------------------------------------------');
w('  그래프 5: permutation_phase_decoder.png');
w('  Phase Decoder — 순열 검정 (Permutation Test)');
w('----------------------------------------------------------------');
w('  X축: 정확도 (%)');
w('  Y축: 빈도 (해당 정확도가 나온 순열 횟수)');
w('  회색 막대: 시기 라벨을 무작위로 섞었을 때의 정확도 분포 (귀무분포)');
w('  빨간 세로선: 실제 데이터의 정확도');
wn();
w('[해석]');
w('  빨간 선이 회색 분포 오른쪽 끝 → 유의미 (행동이 정말로 시기마다 다름)');
w('  p = (무작위 정확도 >= 실제 정확도인 비율)');
w('  p < 0.05 → 통계적으로 유의');
w('  p = 0.0000 → 100번의 순열 중 한 번도 넘지 못함 = 매우 강력한 증거');
wn();

w('----------------------------------------------------------------');
w('  그래프 6: confusion_group_decoder.png');
w('  Group Decoder — 혼동 행렬 (Active vs Passive 분류기)');
w('----------------------------------------------------------------');
w('  행(Y축): 실제 그룹 (True Group)');
w('  열(X축): 예측된 그룹 (Predicted Group)');
w('  각 셀: 횟수 (행 백분율) + 셀 의미');
w('    - 좌상: True Negative (Passive를 Passive로 정확 분류)');
w('    - 우상: False Positive (Passive를 Active로 잘못 분류)');
w('    - 좌하: False Negative (Active를 Passive로 잘못 분류)');
w('    - 우하: True Positive (Active를 Active로 정확 분류)');
wn();
w('[해석]');
w('  우연 수준 = 50%. 정확도 > 80% → 매우 강한 결과.');
w('  높은 정확도 = 자가투여(Active) 경험이 행동을 근본적으로 변화시킴.');
w('  Active→Passive 오분류 = 해당 Active 마우스가 약한 중독 표현형 보임.');
wn();

w('----------------------------------------------------------------');
w('  그래프 7: permutation_group_decoder.png');
w('  Group Decoder — 순열 검정');
w('----------------------------------------------------------------');
w('  회색 막대: 그룹 라벨을 무작위로 섞었을 때의 정확도 분포');
w('  빨간 세로선: 실제 정확도');
w('  회색 분포 중심 ~50% = 정상 (2-class 우연수준)');
w('  빨간 선이 극우단 → Active vs Passive 행동 차이가 매우 확실');
wn();

w('----------------------------------------------------------------');
w('  그래프 8: predict_reexposure_scatter.png');
w('  Re-exposure 행동 예측 (Earlier Phases → Re-exposure)');
w('----------------------------------------------------------------');
w('  X축: 실제 Re-exposure PR breakpoint (관찰값)');
w('  Y축: 예측된 Re-exposure PR breakpoint (모델 예측값)');
w('  빨간 점: Active 마우스, 파란 점: Passive 마우스');
w('  점선 대각선: 완벽한 예측 기준선 (점이 가까울수록 좋음)');
w('  초록 음영: ±5 오차 이내 영역 (예측 양호 구간)');
w('  회색 수직선: 각 마우스의 예측 오차 크기');
w('  텍스트 라벨: 각 마우스 ID');
w('  통계 박스: R², rho, p값, MAPE');
wn();
w('[해석]');
w('  R² > 0 → 이전 행동이 재노출 행동을 예측하는 데 유용함');
w('  R² = 0.27 → 분산의 27%를 설명');
w('  Active 마우스 우상단 = 높은 PR = 강한 재발 경향');
w('  Passive 마우스 좌하단 = 낮은 PR = 약한 동기');
w('  대각선에서 먼 점 = 예측이 어려운 마우스');
w('  핵심 의미: 중독 재발(relapse)을 이전 행동 패턴으로 예측 가능성.');
wn();

w('----------------------------------------------------------------');
w('  그래프 9: predict_reexposure_bars.png');
w('  Re-exposure: Actual vs Predicted (마우스별 막대 비교)');
w('----------------------------------------------------------------');
w('  Y축: 개별 마우스 (실제 PR 기준 내림차순 정렬)');
w('  X축: PR Breakpoint');
w('  파란 막대: 실제 Re-exposure PR');
w('  주황 막대: 모델이 예측한 PR');
w('  오른쪽 숫자: 예측 오차 (Predicted - Actual)');
wn();
w('[해석]');
w('  두 막대가 비슷한 길이 → 해당 마우스를 잘 예측함');
w('  큰 차이 → 해당 마우스의 재노출 행동을 이전 행동으로 예측하기 어려움');
w('  어떤 마우스가 잘 예측되고 안 되는지 한눈에 파악 가능');
wn();

w('----------------------------------------------------------------');
w('  그래프 10: permutation_reexposure.png');
w('  Re-exposure 예측 — 순열 검정');
w('----------------------------------------------------------------');
w('  X축: R² 값');
w('  회색 막대: Re-exposure 값을 무작위로 섞었을 때의 R² 분포');
w('  빨간 세로선: 실제 R²');
wn();
w('[해석]');
w('  회색 분포 대부분 < 0 = 정상 (무작위 예측은 평균보다 나쁨)');
w('  빨간 선이 오른쪽 끝 → 유의미한 예측');
w('  p < 0.05 → 이전 시기의 행동이 재노출 행동을 유의하게 예측');
w('  p >= 0.05 → 샘플 크기 부족으로 검정력 부족 가능');
wn();

w('----------------------------------------------------------------');
w('  그래프 11: feature_importance_phase.png');
w('  Feature Importance — Phase Decoder (특성 중요도)');
w('----------------------------------------------------------------');
w('  Y축: 행동 특성(Feature) 이름');
w('  X축: 정확도 감소량 (Drop in Accuracy)');
w('    해당 특성을 무작위로 섞었을 때 정확도가 얼마나 떨어지는지');
w('  막대가 길수록 → 시기 구분에 더 중요한 특성');
w('  막대 ≈ 0 또는 음수 → 기여하지 않음 (잡음)');
wn();
w('[주요 특성 해석]');
w('  Requirement_cum: 누적 요구량. 시간 경과에 따라 자연히 증가.');
w('  lick_meanDur_s: 평균 핥기 지속시간. 시기마다 변화.');
w('  RequirementLast: PR breakpoint. 보상 동기의 시기별 차이.');
w('  lick_freq_per_min: 핥기 빈도. 구강 행동 빈도가 시기의 지표.');
w('  pupil_mean: 평균 동공 크기. 각성/약물 상태 반영.');
wn();

w('================================================================');
w('  종합 요약');
w('================================================================');
wn();
w('[Module 11 핵심 발견]');
w('  - RW 모델로 각 마우스의 학습률(alpha)을 정량화');
w('  - 시기별 학습률 변화: 금단의 느린 소거, Post의 빠른 재적응');
w('  - 예측 오차(PE)로 시기 전환 시 "놀라움" 크기 정량화');
wn();
w('[Module 12 핵심 발견]');
w('  - Phase Decoder: 행동으로 실험 시기를 우연 이상으로 예측 가능');
w('  - Group Decoder: 행동으로 Active/Passive를 높은 정확도로 구별');
w('    → 자가투여가 행동을 근본적으로 변화시킴');
w('  - Re-exposure 예측: 이전 시기 행동으로 재노출 행동 예측 가능');
w('    → 재발 취약성의 행동적 전조(precursor) 존재 가능성');
w('  - Feature Importance: 핵심 행동 바이오마커 식별');
wn();
w('[임상적 함의]');
w('  중독의 행동적 표현형이 정량적으로 측정 가능하며,');
w('  머신러닝 도구로 중독 위험을 조기 식별하거나');
w('  치료 효과를 모니터링할 수 있는 가능성을 시사.');
w('================================================================');

fclose(fid);
fprintf('  Explanation written: %s\n', outFile);
end

%% ########################################################################
%  LOMO Decoder Helper (Leave-One-Mouse-Out LDA)
%  ########################################################################
function [acc, accPerm, confMat] = lomoDecoder(X, y, mouseID, nPerm)
% Leave-one-mouse-out cross-validation with LDA classifier.
% X: feature matrix (obs x features)
% y: label vector (numeric)
% mouseID: mouse identifier (string/categorical)
% nPerm: number of permutations (0 = skip)

uMice = unique(mouseID);
nTest = numel(y);
yPred = nan(nTest, 1);
uClasses = unique(y(~isnan(y)));
nC = numel(uClasses);

for mi = 1:numel(uMice)
    testIdx = mouseID == uMice(mi);
    trainIdx = ~testIdx;
    if sum(trainIdx) < nC + 1, continue; end
    Xtrain = X(trainIdx,:); yTrain = y(trainIdx);
    Xtest  = X(testIdx,:);

    % LDA: compute class means and shared covariance
    mu = nan(nC, size(X,2));
    for ci = 1:nC
        mu(ci,:) = mean(Xtrain(yTrain==uClasses(ci),:), 1);
    end
    Sw = cov(Xtrain) + 1e-4*eye(size(Xtrain,2));  % regularized

    % Classify by minimum Mahalanobis distance
    testGlobal = find(testIdx(:)');
    for row = 1:size(Xtest,1)
        dists = nan(nC,1);
        for ci = 1:nC
            d = Xtest(row,:) - mu(ci,:);
            dists(ci) = d / Sw * d';
        end
        [~, bestC] = min(dists);
        yPred(testGlobal(row)) = uClasses(bestC);
    end
end

ok = ~isnan(yPred);
acc = mean(yPred(ok) == y(ok));

% Confusion matrix
confMat = zeros(nC, nC);
for r = 1:nC
    for c = 1:nC
        confMat(r,c) = sum(y(ok)==uClasses(r) & yPred(ok)==uClasses(c));
    end
end

% Permutation test
accPerm = nan(max(nPerm,0), 1);
for p = 1:nPerm
    yShuf = y(randperm(numel(y)));
    yPredP = nan(nTest,1);
    for mi = 1:numel(uMice)
        testIdx = mouseID == uMice(mi); trainIdx = ~testIdx;
        if sum(trainIdx) < nC+1, continue; end
        Xtrain = X(trainIdx,:); yTrain = yShuf(trainIdx);
        Xtest  = X(testIdx,:);
        mu = nan(nC, size(X,2));
        for ci = 1:nC
            mask = yTrain==uClasses(ci);
            if any(mask), mu(ci,:) = mean(Xtrain(mask,:),1); else, mu(ci,:) = zeros(1,size(X,2)); end
        end
        Sw = cov(Xtrain) + 1e-4*eye(size(Xtrain,2));
        testGlobal = find(testIdx(:)');
        for row = 1:size(Xtest,1)
            dists = nan(nC,1);
            for ci = 1:nC
                d = Xtest(row,:) - mu(ci,:);
                dists(ci) = d / Sw * d';
            end
            [~, bestC] = min(dists);
            yPredP(testGlobal(row)) = uClasses(bestC);
        end
    end
    okP = ~isnan(yPredP);
    accPerm(p) = mean(yPredP(okP) == yShuf(okP));
end
end

%% ########################################################################
%  Rescorla-Wagner helper: V(t+1) = V(t) + alpha * (R(t) - V(t))
%  ########################################################################
function [V, pe] = runRW(r, alpha)
nT = numel(r);
V = zeros(nT,1);
pe = zeros(nT,1);
V(1) = r(1);  % initialize to first reward
for t = 1:nT-1
    pe(t) = r(t) - V(t);
    V(t+1) = V(t) + alpha * pe(t);
end
pe(nT) = r(nT) - V(nT);
end

%% ########################################################################
%  EVENT-LOCKED HELPERS
%  ########################################################################
function traces = extractEventTraces(tb, signal, eventTimes, tAxis)
% Extract signal snippets around each event time
nEvents = numel(eventTimes);
traces = nan(nEvents, numel(tAxis));

for e = 1:nEvents
    tWin = eventTimes(e) + tAxis;
    if tWin(1) < min(tb) || tWin(end) > max(tb), continue; end

    % Ensure tb is sorted and unique for interp1
    [tbU, ia] = unique(tb, 'stable');
    sigU = signal(ia);
    good = isfinite(tbU) & isfinite(sigU);
    tbU = tbU(good); sigU = sigU(good);
    if numel(tbU) < 2, continue; end
    [tbU, ord] = sort(tbU); sigU = sigU(ord);

    try
        traces(e,:) = interp1(tbU, sigU, tWin, 'linear', NaN);
    catch
        continue;
    end
end

% Remove all-NaN rows
goodRows = any(isfinite(traces), 2);
traces = traces(goodRows, :);
end

function plotEventLockedByPhaseGroup(cellData, tAxis, periods, COL, titleStr, ylab, outDir, prefix)
% cellData: {mouse_key, day_index, period, group, trace}
fh = figure('Color','w','Position',[60 60 1200 550]);
tiledlayout(1, numel(periods), 'TileSpacing','compact','Padding','compact');

for pi = 1:numel(periods)
    ax = nexttile; hold(ax,'on');
    grpTraces = struct();
    for gi = ["Active","Passive"]
        traces = [];
        for r = 1:size(cellData,1)
            if string(cellData{r,3}) == periods(pi) && string(cellData{r,4}) == gi
                traces = [traces; cellData{r,5}]; %#ok<AGROW>
            end
        end
        grpTraces.(char(gi)) = traces;
        if isempty(traces), continue; end
        mu = mean(traces, 1, 'omitnan');
        se = std(traces, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(traces), 1));
        c = COL.(char(gi));
        fill(ax, [tAxis fliplr(tAxis)], [mu-se fliplr(mu+se)], c, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, tAxis, mu, 'Color', c, 'LineWidth', 1.8, 'DisplayName', char(gi));
    end
    xline(ax, 0, 'k--');

    % Ranksum on mean post-event amplitude (0 to +2s) between groups
    postIdx = tAxis >= 0 & tAxis <= 2;
    ampA = []; ampP = [];
    if isfield(grpTraces,'Active') && ~isempty(grpTraces.Active)
        ampA = mean(grpTraces.Active(:, postIdx), 2, 'omitnan');
    end
    if isfield(grpTraces,'Passive') && ~isempty(grpTraces.Passive)
        ampP = mean(grpTraces.Passive(:, postIdx), 2, 'omitnan');
    end
    pRS = safeRanksum(ampA, ampP);
    title(ax, sprintf('%s %s (A:%d, P:%d)', char(periods(pi)), starStr(pRS), numel(ampA), numel(ampP)), 'FontSize', 9);

    if pi==1, ylabel(ax, ylab); end
    xlabel(ax, 'Time (s)');
    grid(ax,'off'); box(ax,'off');
    if pi==numel(periods), legend(ax,'Location','northeastoutside'); end
end
sgtitle(titleStr);
addBottomStatNote(fh, ...
    'Stars in title: Wilcoxon rank-sum on mean post-event amplitude (0-2s). A:n / P:n = sessions. *p<.05 **p<.01 ***p<.001');
printpng(fh, fullfile(outDir, sprintf('%s_by_phase.png', prefix))); close(fh);
end

%% ########################################################################
%  SHARED HELPERS
%  ########################################################################
function feats = pickNumericFeatures(D)
V = D.Properties.VariableNames;
feats = {};
skip = {'mouse_key','day_index','day_name','session_idx','Session_Paradigm','isPassive', ...
        'Cage','Color','Sex','GroupType','PairID','Group','Period','RewardType', ...
        'cage','color','RoleInPair'};
for i = 1:numel(V)
    if any(strcmpi(V{i}, skip)), continue; end
    if isnumeric(D.(V{i})), feats{end+1} = V{i}; end %#ok<AGROW>
end
end

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

function T2 = ensureStr(T2, nm)
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

function p = safeAnovP(aov, termName)
p = NaN;
try
    % Try converting anovaResult to table first (R2024b+)
    tbl = [];
    if ~istable(aov)
        try tbl = dataset2table(aov); catch, end
        if isempty(tbl)
            try
                props = properties(aov);
                if ismember('Term',props) && ismember('pValue',props)
                    terms = string(aov.Term);
                    pvals = double(aov.pValue);
                    idx = find(strcmpi(terms, termName), 1);
                    if isempty(idx), idx = find(contains(lower(terms), lower(string(termName))), 1); end
                    if ~isempty(idx), p = pvals(idx); end
                    return;
                end
            catch, end
        end
        if ~isempty(tbl), aov = tbl; end
    end

    if istable(aov) || (isprop(aov,'Properties'))
        if ismember('Term', aov.Properties.VariableNames)
            terms = string(aov.Term);
        elseif ~isempty(aov.Properties.RowNames)
            terms = string(aov.Properties.RowNames);
        else, return;
        end
        idx = find(strcmpi(terms, termName), 1, 'first');
        if isempty(idx)
            idx = find(contains(lower(terms), lower(string(termName))), 1, 'first');
        end
        if isempty(idx), return; end
        if ismember('pValue', aov.Properties.VariableNames)
            p = double(aov.pValue(idx));
        end
    end
catch, end
end

function s = starStr(p)
if ~isfinite(p), s = ''; return; end
if p < 0.001, s = '***';
elseif p < 0.01, s = '**';
elseif p < 0.05, s = '*';
else, s = 'n.s.';
end
end

function fe = normCoefCols(fe)
% Normalise LME/GLME coefficient table column names across MATLAB versions.
% Ensures columns: Name, Estimate, SE, tStat, DF, pValue.
vn = fe.Properties.VariableNames;
renameMap = { ...
    'Name', {'Name','Term','Row'}; ...
    'Estimate', {'Estimate','Value','Coef','Beta'}; ...
    'SE', {'SE','StdErr','se','StdError','SE_'}; ...
    'tStat', {'tStat','tstat','t','t_Stat','tValue'}; ...
    'DF', {'DF','df','DF1'}; ...
    'pValue', {'pValue','p','Pr___t__','Prob'}};
for ri = 1:size(renameMap,1)
    target = renameMap{ri,1};
    aliases = renameMap{ri,2};
    if ismember(target, vn), continue; end
    for ai = 1:numel(aliases)
        idx = strcmpi(vn, aliases{ai});
        if any(idx)
            fe.Properties.VariableNames{find(idx,1)} = target;
            vn = fe.Properties.VariableNames;
            break;
        end
    end
end
% Add Lower/Upper 95% CI if missing
if ismember('Estimate',fe.Properties.VariableNames) && ismember('SE',fe.Properties.VariableNames)
    if ~ismember('Lower95',fe.Properties.VariableNames)
        fe.Lower95 = fe.Estimate - 1.96 * fe.SE;
        fe.Upper95 = fe.Estimate + 1.96 * fe.SE;
    end
end
end

function aovT = safeAnov2Table(aov)
% Convert anova output (table, dataset, or anovaResult) to a clean table.
try
    if isa(aov, 'dataset'), aovT = dataset2table(aov);
    elseif istable(aov), aovT = aov;
    else
        % MATLAB R2024b+ anovaResult object
        aovT = table();
        aovT.Term = string(aov.Term);
        aovT.FStat = aov.FStat;
        aovT.DF1 = aov.DF1;
        aovT.DF2 = aov.DF2;
        aovT.pValue = aov.pValue;
        return;
    end
    % Normalise: ensure a 'Term' column exists
    vn = aovT.Properties.VariableNames;
    if ~ismember('Term',vn)
        if ismember('Row',vn), aovT.Term = string(aovT.Row);
        else
            rn = aovT.Properties.RowNames;
            if ~isempty(rn), aovT.Term = string(rn); end
        end
    end
catch
    aovT = table();
end
end

function s = fmtP(p)
if ~isfinite(p), s = 'NA'; return; end
if p < 0.001, s = 'p<.001';
elseif p < 0.01, s = sprintf('p=%.3f',p);
else, s = sprintf('p=%.2f',p);
end
end

function addBottomStatNote(fh, noteStr)
% Add a small-font stat annotation at the very bottom of the figure.
% Shrinks the axes slightly to make room, ensuring no overlap.
% Skip axes inside TiledChartLayout (they don't support manual Position).
try
    allAx = findobj(fh, 'Type', 'axes');
    for ai = 1:numel(allAx)
        par = get(allAx(ai), 'Parent');
        if isa(par, 'matlab.graphics.layout.TiledChartLayout'), continue; end
        pos = get(allAx(ai), 'Position');
        if pos(2) < 0.12
            pos(2) = max(pos(2), 0.10);
            pos(4) = min(pos(4), 1-pos(2)-0.05);
            set(allAx(ai), 'Position', pos);
        end
    end
catch, end
annotation(fh, 'textbox', [0.01, 0.002, 0.98, 0.05], ...
    'String', noteStr, 'FontSize', 7, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'FitBoxToText', 'off', 'Margin', 2);
end

function [p, z] = safeRanksum(a, b)
p = NaN; z = NaN;
a = a(isfinite(a)); b = b(isfinite(b));
if numel(a) < 2 || numel(b) < 2, return; end
try
    [p, ~, stats] = ranksum(a, b);
    if isfield(stats,'zval'), z = stats.zval; end
catch, end
end

function addStarsBetweenGroups(ax, xpos, yA, yB, seA, seB, pval)
% Draw bracket + stars between two groups at position xpos
if ~isfinite(pval), return; end
ss = starStr(pval);
if isempty(ss), return; end
yTop = max(yA + seA, yB + seB);
if ~isfinite(yTop), return; end
ylims = ylim(ax);
bump = (ylims(2) - ylims(1)) * 0.04;
yBar = yTop + bump;
text(ax, xpos, yBar + bump*0.5, ss, 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
end

function cmap = redblue()
n = 256;
r = [linspace(0,1,n/2), ones(1,n/2)]';
b = [ones(1,n/2), linspace(1,0,n/2)]';
g = [linspace(0,1,n/2), linspace(1,0,n/2)]';
cmap = [r g b];
end

function plotDualBox(T, ycol, ylab, periods, COL)
off = [-0.12 0.12];
groups = ["Active","Passive"];
muStore = struct(); seStore = struct();
mouseMeansStore = struct('Active',{cell(1,numel(periods))}, 'Passive',{cell(1,numel(periods))});
for gi = 1:numel(groups)
    g = groups(gi); c = COL.(char(g));
    mask = string(T.Group)==g;
    mu = nan(1,numel(periods)); se = nan(1,numel(periods));
    for pi = 1:numel(periods)
        sub = T(mask & string(T.Period)==periods(pi), :);
        if isempty(sub), continue; end
        % Per-mouse mean (one value per mouse, averaged across sessions)
        [G, ~] = findgroups(sub.mouse_key);
        mMeans = splitapply(@(x) mean(x,'omitnan'), double(sub.(ycol)), G);
        mMeans = mMeans(isfinite(mMeans));
        mouseMeansStore.(char(g)){pi} = mMeans;
        if isempty(mMeans), continue; end
        x = pi + off(gi);
        boxchart(repmat(x,numel(mMeans),1), mMeans, 'BoxWidth',0.18, ...
            'MarkerStyle','none','BoxFaceAlpha',0.2,'BoxFaceColor',c);
        scatter(repmat(x,numel(mMeans),1), mMeans, 24, c, 'filled', 'MarkerFaceAlpha',0.65);
        mu(pi) = mean(mMeans); se(pi) = std(mMeans)/sqrt(numel(mMeans));
    end
    muStore.(char(g)) = mu; seStore.(char(g)) = se;
end

% Per-phase ranksum stars (using per-mouse means)
ax = gca;
for pi = 1:numel(periods)
    vA = mouseMeansStore.Active{pi};
    vP = mouseMeansStore.Passive{pi};
    % If Passive has no data in "During" (lick metric) — annotate
    if periods(pi)=="During" && isempty(vP)
        yLim = get(ax, 'YLim');
        text(ax, pi+0.12, yLim(1)+0.05*diff(yLim), 'Passive: no PR', ...
            'FontSize',7, 'Color',[0.5 0.5 0.5], 'FontAngle','italic', ...
            'HorizontalAlignment','center');
        continue;
    end
    pRS = safeRanksum(vA, vP);
    addStarsBetweenGroups(ax, pi, muStore.Active(pi), muStore.Passive(pi), ...
        seStore.Active(pi), seStore.Passive(pi), pRS);
end

% Overall Kruskal-Wallis across phases using per-mouse means
allMM = []; allPer = [];
for pi = 1:numel(periods)
    for gi = 1:numel(groups)
        mm = mouseMeansStore.(char(groups(gi))){pi};
        if ~isempty(mm)
            allMM  = [allMM;  mm(:)]; %#ok<AGROW>
            allPer = [allPer; repmat(periods(pi), numel(mm), 1)]; %#ok<AGROW>
        end
    end
end
pKW = NaN;
if numel(allMM) > 10
    try pKW = kruskalwallis(allMM, allPer, 'off'); catch, end
end

set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
ylabel(ylab);
hA = scatter(nan,nan,30,COL.Active,'filled','DisplayName','Active');
hP = scatter(nan,nan,30,COL.Passive,'filled','DisplayName','Passive');
legend([hA hP],'Location','northeastoutside'); grid off; box off;

% Count mice per group
nAm = 0; nPm = 0;
for pi2 = 1:numel(periods)
    nAm = max(nAm, numel(mouseMeansStore.Active{pi2}));
    nPm = max(nPm, numel(mouseMeansStore.Passive{pi2}));
end
kwStr = '';
if isfinite(pKW), kwStr = sprintf('  KW across phases: %s.', fmtP(pKW)); end
addBottomStatNote(gcf, sprintf( ...
    'Dots=per-mouse means. Stars: Wilcoxon rank-sum. nA=%d, nP=%d.%s *p<.05 **p<.01 ***p<.001', nAm, nPm, kwStr));
end

% ========================= COHORT MAPPING =========================
% Exact same logic as analyze_passive_active_dashboard_dec2 (proven working)

function cohort = buildNewCohortRoster()
% Hard roster:
% 6100: black=Active, orange/red=Passive  (P1)
% 0911: red=Active vs orange=Passive      (P2)
% 0911: white=Active vs black=Passive     (P3)
% 0910: black=Active, orange/red=Passive  (P4)
% 6099: orange=Active vs red=Passive      (P5)
% 6099: black=Active vs white=Passive     (P6; white died day13)

cage  = strings(0,1);
color = strings(0,1);
group = strings(0,1);
pair  = strings(0,1);
role  = strings(0,1);

add("6100","black","Active","P1","Active");
add("6100","orange","Passive","P1","Passive");
add("6100","red","Passive","P1","Passive");

add("0911","red","Active","P2","Active");
add("0911","orange","Passive","P2","Passive");

add("0911","white","Active","P3","Active");
add("0911","black","Passive","P3","Passive");

add("0910","black","Active","P4","Active");
add("0910","orange","Passive","P4","Passive");
add("0910","red","Passive","P4","Passive");

add("6099","orange","Active","P5","Active");
add("6099","red","Passive","P5","Passive");

add("6099","black","Active","P6","Active");
add("6099","white","Passive","P6","Passive");

cohort = table(cage,color,group,pair,role, ...
    'VariableNames', {'cage','color','Group','PairID','RoleInPair'});

    function add(cg, col, gp, pid, rl)
        cage(end+1,1)  = string(cg);
        color(end+1,1) = lower(string(col));
        group(end+1,1) = string(gp);
        pair(end+1,1)  = string(pid);
        role(end+1,1)  = string(rl);
    end
end

function T = attachCohortGroupAndPair(T, cohort)
% Adds Group/PairID/RoleInPair to T using robust parsing of T.mouse_key.
% Rule:
%   1) If mouse_key contains explicit "_a"/"_p" tag, that wins.
%   2) Else use cage+color mapping from cohort table.
%   Any unmapped mouse becomes <missing>.

mk = string(T.mouse_key);
mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));

% Parse cage + color robustly (uses contains, NOT word-boundary regex)
[cg, col] = parseMouseKeyRobust(mk_norm);
T.cage  = cg;
T.color = col;

% Start with explicit A/P tag if present in mouse_key
grp = strings(height(T),1);
grp(contains(mk_norm," a") | endsWith(mk_norm,"a") | ...
    contains(mk_norm," f a") | contains(mk_norm," m a") | contains(mk_norm,"_a")) = "Active";
grp(contains(mk_norm," p") | endsWith(mk_norm,"p") | ...
    contains(mk_norm," f p") | contains(mk_norm," m p") | contains(mk_norm,"_p")) = "Passive";

% Fill remaining by cohort map (cage+color lookup)
needMap = (grp == "");
keyT = T.cage + "|" + T.color;
keyC = cohort.cage + "|" + cohort.color;
[tf, loc] = ismember(keyT, keyC);

grp(needMap & tf) = cohort.Group(loc(needMap & tf));

% Attach PairID/RoleInPair also from mapping
pair = strings(height(T),1);
role = strings(height(T),1);
pair(tf) = cohort.PairID(loc(tf));
role(tf) = cohort.RoleInPair(loc(tf));

% Apply categoricals
T.Group      = categorical(grp, ["Passive","Active"]);
T.PairID     = categorical(pair);
T.RoleInPair = categorical(role, ["Passive","Active"]);
end

function [cage, color] = parseMouseKeyRobust(mk_norm)
% mk_norm is already lowercased with underscores/hyphens -> spaces.
% Uses contains() for color matching (NOT word-boundary regex),
% so "6100red" correctly matches "red".
n = numel(mk_norm);
cage  = strings(n,1);
color = strings(n,1);

colors = ["black","red","orange","white"];

for i = 1:n
    s = mk_norm(i);

    % cage: first 4-digit number
    tok = regexp(s,'(\d{4})','tokens','once');
    if ~isempty(tok), cage(i) = string(tok{1}); end

    % color: substring match (contains), NOT word-boundary regex
    for c = colors
        if contains(s, c)
            color(i) = c;
            break;
        end
    end
end
end
