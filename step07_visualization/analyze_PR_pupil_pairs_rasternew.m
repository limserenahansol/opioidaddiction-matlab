function analyze_PR_pupil_pairs_rasternew()
% analyze_PR_pupil_pairs_raster (REVISED for new cohort + new timeline)
%
% What this script now does (based on your update):
%  1) Uses your NEW cohort mapping (cage/color/sex/group + pair structure).
%  2) Uses your NEW timeline and epoch definitions:
%       - Habituation: Day 1–2 (ignored by default)
%       - Pre (Water PR): Day 3–5
%       - During (Passive training / morphine): Day 6–10
%       - Post (Free morphine PR): Day 11–13
%       - Withdrawal (Water PR): Day 14–16
%       - Re-exposure (Morphine PR): Day 17–18
%  3) Supports “with vs without unreliable switch-days” comparisons:
%       - Default unreliable days: [4, 6, 11, 14] (as you wrote)
%       - Generates outputs for:
%           (A) include all days (Day3–18)
%           (B) exclude unreliable days (Day3–18 minus [4,6,11,14])
%  4) Replaces old hard-coded mouse pairs with your NEW pairs.
%  5) Handles 6099_white death (data only Day1–13) gracefully.
%
% Notes/assumptions (kept minimal and explicit):
%  - Your CSV mouse_key may be formatted like "6100_red" OR "6100red".
%    This script canonicalizes to "6100_red" internally (digits + '_' + color).
%  - If a mouse exists in CSV but not in the cohort map, it is kept but labeled Group="Unknown".
%
% Outputs go under:
%   runDir/figs/PR_pupil_pairs_raster_REVISED/

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

% Normalize key columns
T = ensureString(T,'mouse_key');
T = ensureString(T,'day_name');
T = ensureString(T,'Session_Paradigm');

% Ensure day_index exists
assert(ismember('day_index',T.Properties.VariableNames),'CSV missing day_index');

% TTL numerics
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

%% ---------- NEW: cohort mapping (your provided info) ----------
cohort = build_new_cohort_map();

% Canonicalize mouse_key in T to "####_color"
T.mouse_key = canonicalize_mouse_keys(string(T.mouse_key));

% Attach cohort info into T (Group/Sex/Cage/Color/PairID/AliveLastDay)
T = attach_cohort_to_table(T, cohort);

%% ---------- define timeline / epochs ----------
epochs = define_epochs();             % day ranges per epoch
daysMain = 3:18;                      % ignore day1-2 by default
unreliableDays = [4 6 11 14];         % your list

daysAll_includingSwitch = daysMain;
daysAll_excludingSwitch = setdiff(daysMain, unreliableDays);

%% ---------- per-session metrics, then per-day medians ----------
[S, D] = fast_session_day_metrics_basic(T, runDir);

% Canonicalize keys in S/D too (safety)
S.mouse_key = canonicalize_mouse_keys(string(S.mouse_key));
D.mouse_key = canonicalize_mouse_keys(string(D.mouse_key));

% Attach cohort to S/D
S = attach_cohort_to_table(S, cohort);
D = attach_cohort_to_table(D, cohort);

% Define Group categorical with Unknown included
if ~ismember('Group',D.Properties.VariableNames)
    D.Group = repmat("Unknown",height(D),1);
end
D.Group = categorical(string(D.Group), ["Active","Passive","Unknown"]);

%% ---------- output dir ----------
outDir = fullfile(runDir,'figs','PR_pupil_pairs_raster_REVISED');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Make subfolders for the two day-selection modes
outDirA = fullfile(outDir,'A_includeSwitchDays');
outDirB = fullfile(outDir,'B_excludeUnreliableDays');
if ~exist(outDirA,'dir'), mkdir(outDirA); end
if ~exist(outDirB,'dir'), mkdir(outDirB); end

%% ---------- run BOTH modes: include vs exclude unreliable days ----------
run_mode(D, S, T, cohort, epochs, daysAll_includingSwitch, unreliableDays, true,  outDirA);
run_mode(D, S, T, cohort, epochs, daysAll_excludingSwitch, unreliableDays, false, outDirB);

fprintf('Done. Outputs:\n  %s\n  %s\n', outDirA, outDirB);
end

%% ======================================================================
%                               MODE RUNNER
% ======================================================================
function run_mode(D, S, T, cohort, epochs, daysUse, unreliableDays, includeSwitchDays, outDir)

fprintf('\n=============================\n');
fprintf('MODE: %s\n', tern(includeSwitchDays,'INCLUDE switch days','EXCLUDE unreliable days'));
fprintf('daysUse = [%s]\n', num2str(daysUse));
fprintf('unreliableDays = [%s]\n', num2str(unreliableDays));
fprintf('Output: %s\n', outDir);
fprintf('=============================\n');

% Restrict D/S to daysUse
Duse = D(ismember(D.day_index,daysUse),:);
Suse = S(ismember(S.day_index,daysUse),:);

%% 1) PR spaghetti (RequirementLast)
if ismember('RequirementLast', Duse.Properties.VariableNames)
    % Active only spaghetti across daysUse
    plot_spaghetti_metric(Duse, 'RequirementLast', 'PR / RequirementLast', ...
        sprintf('PR (RequirementLast) - Active only (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        outDir, 'ActiveOnly', daysUse);

    % Active vs Passive panels
    plot_spaghetti_metric_groups(Duse, 'RequirementLast', 'PR / RequirementLast', ...
        sprintf('PR (RequirementLast): Active vs Passive (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        fullfile(outDir,'spaghetti_PR_active_vs_passive.png'), ...
        daysUse, false);
end

%% 2) Pupil + Tail immersion spaghetti
if ismember('pupil_mean', Duse.Properties.VariableNames)
    plot_spaghetti_metric(Duse, 'pupil_mean', 'Pupil diameter (px)', ...
        sprintf('Pupil diameter - all mice (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        outDir, 'All', daysUse);

    plot_spaghetti_metric_groups(Duse, 'pupil_mean', 'Pupil diameter (px)', ...
        sprintf('Pupil: Active vs Passive (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        fullfile(outDir,'spaghetti_pupil_active_vs_passive.png'), ...
        daysUse, true);
end

if ismember('Immersion_Latency_s', Duse.Properties.VariableNames)
    plot_spaghetti_metric(Duse, 'Immersion_Latency_s', 'Tail immersion latency (s)', ...
        sprintf('Tail immersion latency - all mice (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        outDir, 'All', daysUse);

    plot_spaghetti_metric_groups(Duse, 'Immersion_Latency_s', 'Tail immersion latency (s)', ...
        sprintf('Tail immersion: Active vs Passive (%s)', tern(includeSwitchDays,'all days','exclude unreliable')), ...
        fullfile(outDir,'spaghetti_immersion_active_vs_passive.png'), ...
        daysUse, false);

    % Per-day A vs P stats with stars for days where immersion is relevant
    compare_AP_means_per_day_with_stars(Duse, 'Immersion_Latency_s', intersect(daysUse,3:18), outDir);
end

%% 3) Epoch summaries (your requested: pre/during/post/withdrawal/reexposure)
%    - For each metric: epoch mean per mouse, plus group-level visualization
metricsEpoch = {};

% Core metrics
core = {'RequirementLast','pupil_mean','Immersion_Latency_s'};
for i=1:numel(core)
    if ismember(core{i}, D.Properties.VariableNames)
        metricsEpoch{end+1} = core{i}; %#ok<AGROW>
    end
end

% TST/HOT metrics (all columns starting with TST_ or HOT_)
tstCols = D.Properties.VariableNames(startsWith(D.Properties.VariableNames,'TST_'));
hotCols = D.Properties.VariableNames(startsWith(D.Properties.VariableNames,'HOT_'));
% Remove *_File columns
tstCols = tstCols(~endsWith(tstCols,'_File'));
hotCols = hotCols(~endsWith(hotCols,'_File'));
metricsEpoch = [metricsEpoch, tstCols, hotCols];

if ~isempty(metricsEpoch)
    fprintf('Epoch summary for %d metrics...\n', numel(metricsEpoch));
    for i = 1:numel(metricsEpoch)
        mName = metricsEpoch{i};
        if ~ismember(mName, D.Properties.VariableNames), continue; end
        if ~isnumeric(D.(mName)), continue; end
        epoch_summary_plots(D, cohort, epochs, mName, daysUse, outDir);
        epoch_pairwise_active_minus_passive(D, cohort, epochs, mName, daysUse, outDir);
    end
end

%% 4) Licking: spaghetti + cumulative across days (ACTIVE and PASSIVE, as you requested)
% Build per-mouse/day lick metrics from Suse (sessions -> per-day mean)
plot_licking_all_groups_from_S(Suse, daysUse, outDir);

%% 5) Reward & lick rasters + per-day cumulative licking curves (Active only)
% Choose representative days that match your paradigm
rasterDays = intersect([5 10 13 14 16 18], daysUse);

[timeGridMin, meanCumPerDay, meanCumEarly, meanCumLate, sessInfo] = ...
    compute_cumulative_and_events_active(T, D, daysUse);

for dd = rasterDays
    raster_day_plot(sessInfo, dd, 'reward', outDir);
    raster_day_plot(sessInfo, dd, 'lick',   outDir);
end

plot_cumLicks_allDays(timeGridMin, meanCumPerDay, unique(daysUse), outDir);

% “Early vs late” is now epoch-aware (pre+during vs post+withdrawal+reexposure is more meaningful)
plot_cumLicks_epoch_groups(timeGridMin, meanCumPerDay, epochs, outDir);

%% 6) Pair zoom panels (NEW pairs)
pairMetrics = {'pupil_mean','Immersion_Latency_s','RequirementLast', ...
               'TST_Pct_Non_moving','HOT_Pct_Non_moving','HOT_Pct_Licking','HOT_Pct_Flinching'};

pairs = cohort_pairs_as_mousekey_lists(cohort);
pairLabels = cohort_pair_labels(cohort);

for p = 1:numel(pairs)
    mkList = pairs{p};
    for mi = 1:numel(pairMetrics)
        mName = pairMetrics{mi};
        if ~ismember(mName, D.Properties.VariableNames), continue; end
        pair_zoom_plot(D, mkList, pairLabels{p}, mName, outDir, daysUse);
    end
end

end

%% ======================================================================
%                         COHORT / EPOCH DEFINITIONS
% ======================================================================
function cohort = build_new_cohort_map()
% Your cohort definition, canonical keys = "####_color" all lowercase color.

rows = {
    % cage, color, sex, group, pairID, aliveLastDay
    6100,'red',    'F','Passive', 1, 18
    6100,'orange', 'F','Passive', 1, 18
    6100,'black',  'F','Active',  1, 18

    0911,'red',    'F','Active',  2, 18
    0911,'orange', 'F','Passive', 2, 18

    0911,'black',  'F','Passive', 3, 18
    0911,'white',  'F','Active',  3, 18

    0910,'red',    'M','Passive', 4, 18
    0910,'orange', 'M','Passive', 4, 18
    0910,'black',  'M','Active',  4, 18

    6099,'red',    'M','Passive', 5, 18
    6099,'orange', 'M','Active',  5, 18

    6099,'black',  'M','Active',  6, 18
    6099,'white',  'M','Passive', 6, 13   % died after day13
    };

cohort = cell2table(rows, 'VariableNames', ...
    {'Cage','Color','Sex','Group','PairID','AliveLastDay'});

% Make canonical mouse_key
cohort.mouse_key = arrayfun(@(i) sprintf('%04d_%s', cohort.Cage(i), lower(string(cohort.Color(i)))), ...
    (1:height(cohort))', 'UniformOutput', false);
cohort.mouse_key = string(cohort.mouse_key);

end

function epochs = define_epochs()
% Returns struct array with fields: Name, Days

epochs = struct('Name',{},'Days',{});
epochs(1).Name = "Pre_WaterPR";
epochs(1).Days = 3:5;

epochs(2).Name = "During_PassiveTraining_Morphine";
epochs(2).Days = 6:10;

epochs(3).Name = "Post_FreeMorphinePR";
epochs(3).Days = 11:13;

epochs(4).Name = "Withdrawal_WaterPR";
epochs(4).Days = 14:16;

epochs(5).Name = "Reexposure_MorphinePR";
epochs(5).Days = 17:18;
end

function pairs = cohort_pairs_as_mousekey_lists(cohort)
% Each pair returns list of mouse_keys to plot together.
% Pair 1 has 1 active + 2 passive, etc.

pairIDs = unique(cohort.PairID,'stable');
pairs = cell(numel(pairIDs),1);

for i = 1:numel(pairIDs)
    pid = pairIDs(i);
    rows = cohort(cohort.PairID==pid,:);
    pairs{i} = cellstr(rows.mouse_key);
end
end

function labels = cohort_pair_labels(cohort)
pairIDs = unique(cohort.PairID,'stable');
labels = cell(numel(pairIDs),1);
for i=1:numel(pairIDs)
    pid = pairIDs(i);
    rows = cohort(cohort.PairID==pid,:);
    act = rows.mouse_key(rows.Group=="Active");
    pas = rows.mouse_key(rows.Group=="Passive");
    labels{i} = sprintf('Pair %d: Active=%s ; Passive=%s', pid, strjoin(cellstr(act),','), strjoin(cellstr(pas),','));
end
end

%% ======================================================================
%                         ATTACH COHORT INFO
% ======================================================================
function T = attach_cohort_to_table(T, cohort)
% Adds: Cage, Color, Sex, Group, PairID, AliveLastDay
if ~ismember('mouse_key',T.Properties.VariableNames)
    return;
end

mk = canonicalize_mouse_keys(string(T.mouse_key));
T.mouse_key = mk;

% Build lookup
[tf, loc] = ismember(mk, cohort.mouse_key);

T.Cage = nan(height(T),1);
T.Color = repmat("",height(T),1);
T.Sex = repmat("",height(T),1);
T.Group = repmat("Unknown",height(T),1);
T.PairID = nan(height(T),1);
T.AliveLastDay = nan(height(T),1);

T.Cage(tf) = cohort.Cage(loc(tf));
T.Color(tf) = string(cohort.Color(loc(tf)));
T.Sex(tf) = string(cohort.Sex(loc(tf)));
T.Group(tf) = string(cohort.Group(loc(tf)));
T.PairID(tf) = cohort.PairID(loc(tf));
T.AliveLastDay(tf) = cohort.AliveLastDay(loc(tf));

end

function mk2 = canonicalize_mouse_keys(mk)
% Convert various formats to "####_color" if possible.
% Examples:
%   "6100red" -> "6100_red"
%   "6100_red" -> "6100_red"
%   "0911white" -> "0911_white"
% If parsing fails, keeps original string.

mk = string(mk);
mk2 = mk;

for i = 1:numel(mk)
    s = lower(strtrim(mk(i)));
    if strlength(s)==0, continue; end

    % Try pattern: 4 digits then '_' then letters
    tok = regexp(s,'^(\d{4})[_\- ]?([a-z]+)$','tokens','once');
    if ~isempty(tok)
        mk2(i) = sprintf('%s_%s', tok{1}, tok{2});
        continue;
    end

    % Try pattern anywhere inside (last fallback)
    tok = regexp(s,'(\d{4})[_\- ]?([a-z]+)','tokens','once');
    if ~isempty(tok)
        mk2(i) = sprintf('%s_%s', tok{1}, tok{2});
        continue;
    end
end
end

%% ======================================================================
%                           EPOCH SUMMARIES
% ======================================================================
function epoch_summary_plots(D, cohort, epochs, mName, daysUse, outDir)
% For each mouse: compute epoch median across included days.
% Then:
%  - spaghetti across epochs (each mouse)
%  - Active vs Passive panels

% Restrict to daysUse and alive constraint
Duse = D(ismember(D.day_index,daysUse),:);
Duse = apply_alive_filter(Duse);

% Build per-mouse epoch median
E = epoch_table(Duse, epochs, mName);

if isempty(E), return; end

% Plot all mice across epochs
fh = figure('Color','w','Position',[80 80 950 520]); hold on;
mice = unique(E.mouse_key,'stable');
x = 1:numel(epochs);

for i=1:numel(mice)
    ei = E(E.mouse_key==mice(i),:);
    y = nan(1,numel(epochs));
    for k=1:numel(epochs)
        r = ei(ei.Epoch==epochs(k).Name,:);
        if ~isempty(r), y(k) = r.Value; end
    end
    if all(~isfinite(y)), continue; end
    plot(x, y, '-o', 'LineWidth',1, 'MarkerSize',4);
end

% Mean per epoch
G = groupsummary(E, 'Epoch','mean','Value');
yMean = nan(1,numel(epochs));
for k=1:numel(epochs)
    rr = G(G.Epoch==epochs(k).Name,:);
    if ~isempty(rr), yMean(k) = rr.mean_Value; end
end
plot(x, yMean, 'k-', 'LineWidth',3);

set(gca,'XTick',x,'XTickLabel',cellstr(string({epochs.Name})));
xtickangle(25);
ylabel(strrep(mName,'_',' '),'Interpreter','none');
title(sprintf('%s: epoch summary (all mice)', mName),'Interpreter','none');
grid on; box on;

fn = fullfile(outDir, sprintf('epoch_all_%s.png', safeName(mName)));
printpng(fh, fn); close(fh);

% Active vs Passive panels
if ismember('Group',E.Properties.VariableNames)
    fh = figure('Color','w','Position',[80 80 950 880]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    grpNames = ["Active","Passive"];
    for gi=1:2
        nexttile; hold on;
        Eg = E(E.Group==grpNames(gi),:);
        if isempty(Eg)
            text(0.5,0.5,sprintf('No %s data',grpNames(gi)),'Units','normalized'); axis off;
            continue;
        end
        mice = unique(Eg.mouse_key,'stable');
        for i=1:numel(mice)
            ei = Eg(Eg.mouse_key==mice(i),:);
            y = nan(1,numel(epochs));
            for k=1:numel(epochs)
                r = ei(ei.Epoch==epochs(k).Name,:);
                if ~isempty(r), y(k) = r.Value; end
            end
            if all(~isfinite(y)), continue; end
            plot(x, y, '-o', 'LineWidth',1, 'MarkerSize',4);
        end
        Gg = groupsummary(Eg,'Epoch','mean','Value');
        yMean = nan(1,numel(epochs));
        for k=1:numel(epochs)
            rr = Gg(Gg.Epoch==epochs(k).Name,:);
            if ~isempty(rr), yMean(k) = rr.mean_Value; end
        end
        plot(x, yMean, 'k-', 'LineWidth',3);

        set(gca,'XTick',x,'XTickLabel',cellstr(string({epochs.Name})));
        xtickangle(25);
        ylabel(strrep(mName,'_',' '),'Interpreter','none');
        title(sprintf('%s: %s', mName, grpNames(gi)),'Interpreter','none');
        grid on; box on;
    end
    sgtitle(sprintf('Epoch summary: %s (Active vs Passive)', mName),'FontWeight','bold');
    fn = fullfile(outDir, sprintf('epoch_active_vs_passive_%s.png', safeName(mName)));
    printpng(fh, fn); close(fh);
end

end

function epoch_pairwise_active_minus_passive(D, cohort, epochs, mName, daysUse, outDir)
% For each PairID and epoch:
%   diff = ActiveValue - mean(PassiveValues)  (supports 1 or 2 passive)
% Plot per-epoch diffs across pairs + mean +/- SEM.

if ~ismember('PairID',D.Properties.VariableNames) || ~ismember('Group',D.Properties.VariableNames)
    return;
end

Duse = D(ismember(D.day_index,daysUse),:);
Duse = apply_alive_filter(Duse);

E = epoch_table(Duse, epochs, mName);
if isempty(E), return; end

pairIDs = unique(E.PairID(~isnan(E.PairID)),'stable');
if isempty(pairIDs), return; end

x = 1:numel(epochs);
diffMat = nan(numel(pairIDs), numel(epochs));

for pi=1:numel(pairIDs)
    pid = pairIDs(pi);
    Ep = E(E.PairID==pid,:);
    for k=1:numel(epochs)
        Ek = Ep(Ep.Epoch==epochs(k).Name,:);
        if isempty(Ek), continue; end
        a = Ek.Value(Ek.Group=="Active");
        p = Ek.Value(Ek.Group=="Passive");
        if isempty(a) || isempty(p), continue; end
        diffMat(pi,k) = mean(a,'omitnan') - mean(p,'omitnan');
    end
end

% Plot
fh = figure('Color','w','Position',[80 80 980 520]); hold on;
for pi=1:size(diffMat,1)
    y = diffMat(pi,:);
    if all(~isfinite(y)), continue; end
    plot(x, y, '-o', 'LineWidth',1, 'MarkerSize',4);
end

mu = mean(diffMat,1,'omitnan');
se = std(diffMat,0,1,'omitnan') ./ sqrt(sum(isfinite(diffMat),1));
errorbar(x, mu, se, 'k-', 'LineWidth',3, 'CapSize',10);

yline(0,'--k','LineWidth',1);
set(gca,'XTick',x,'XTickLabel',cellstr(string({epochs.Name})));
xtickangle(25);
ylabel(sprintf('%s: Active - Passive', strrep(mName,'_',' ')),'Interpreter','none');
title(sprintf('Pair-wise A-P differences across epochs: %s', mName),'Interpreter','none');
grid on; box on;

fn = fullfile(outDir, sprintf('epoch_pairDiff_AminusP_%s.png', safeName(mName)));
printpng(fh, fn); close(fh);

end

function E = epoch_table(Duse, epochs, mName)
% Returns table with per-mouse per-epoch median Value
% Robust to missing epochs (Dj empty) and missing metadata columns.

if ~ismember(mName, Duse.Properties.VariableNames) || ~isnumeric(Duse.(mName))
    E = table(); 
    return;
end

mice = unique(string(Duse.mouse_key),'stable');
rows = {};

for i=1:numel(mice)
    mk = mice(i);
    Di = Duse(string(Duse.mouse_key)==mk,:);
    if isempty(Di)
        continue;
    end

    % Mouse-level metadata (stable even if an epoch has no rows)
    if ismember('Group',Di.Properties.VariableNames)
        grp = string(Di.Group(1));
    else
        grp = "Unknown";
    end

    if ismember('Sex',Di.Properties.VariableNames)
        sex = string(Di.Sex(1));
    else
        sex = "";
    end

    if ismember('Color',Di.Properties.VariableNames)
        col = string(Di.Color(1));
    else
        col = "";
    end

    if ismember('PairID',Di.Properties.VariableNames)
        pid = double(Di.PairID(1));
    else
        pid = NaN;
    end

    % Epoch loop
    for k=1:numel(epochs)
        dd = epochs(k).Days;
        Dj = Di(ismember(Di.day_index, dd),:);

        if isempty(Dj)
            val = NaN;  % epoch missing for this mouse
        else
            v = Dj.(mName);
            v = v(isfinite(v));
            if isempty(v)
                val = NaN;
            else
                val = median(v,'omitnan');
            end
        end

        rows(end+1,:) = {mk, epochs(k).Name, val, grp, sex, col, pid}; %#ok<AGROW>
    end
end

if isempty(rows)
    E = table();
    return;
end

E = cell2table(rows, 'VariableNames', {'mouse_key','Epoch','Value','Group','Sex','Color','PairID'});
E.Epoch = categorical(E.Epoch, string({epochs.Name}), 'Ordinal', true);
E.Group = categorical(E.Group, ["Active","Passive","Unknown"]);
end


function Duse = apply_alive_filter(Duse)
% Drops rows after AliveLastDay (if available)

if ~ismember('AliveLastDay',Duse.Properties.VariableNames)
    return;
end
alive = Duse.AliveLastDay;
if all(isnan(alive)), return; end

keep = true(height(Duse),1);
has = isfinite(alive);
keep(has) = Duse.day_index(has) <= alive(has);
Duse = Duse(keep,:);
end

%% ======================================================================
%                           PAIR ZOOM
% ======================================================================
function pair_zoom_plot(D, mkList, pairLabel, mName, outDir, daysUse)

Dp = D(ismember(string(D.mouse_key), string(mkList)) & ismember(D.day_index,daysUse), :);
if isempty(Dp) || ~ismember(mName, Dp.Properties.VariableNames), return; end

colors = lines(numel(mkList));
fh = figure('Color','w','Position',[80 80 900 520]); hold on;

for i = 1:numel(mkList)
    mk = string(mkList{i});
    di = Dp(string(Dp.mouse_key)==mk,:);
    if isempty(di), continue; end
    di = di(isfinite(di.(mName)),:);
    if isempty(di), continue; end
    [~,ord] = sort(di.day_index);
    lab = char(mk);
    if ismember('Group',di.Properties.VariableNames)
        g = unique(string(di.Group));
        if ~isempty(g), lab = sprintf('%s (%s)', mk, g(1)); end
    end
    plot(di.day_index(ord), di.(mName)(ord),'-o', ...
        'Color',colors(i,:), 'LineWidth',1.8,'MarkerSize',5, ...
        'DisplayName',lab);
end

xlabel('Day index');
ylabel(strrep(mName,'_',' '),'Interpreter','none');
title({pairLabel; mName},'Interpreter','none');
grid on; box on;
legend('Location','bestoutside');

fn = fullfile(outDir, sprintf('pair_%s_%s.png', safeName(pairLabel), safeName(mName)));
printpng(fh, fn); close(fh);

end

%% ======================================================================
%                   LICKING (ALL GROUPS + CUMULATIVE)
% ======================================================================
function plot_licking_all_groups_from_S(Suse, daysUse, outDir)

if isempty(Suse)
    fprintf('No S data for licking plots.\n'); return;
end

% Build per-mouse/day mean from per-session Suse
lickVars = {'lick_freq_per_min','lick_meanDur_s','lick_medianIEI_s', ...
            'lick_totalDur_s','lick_n','lick_freq_per_10s'};

[g, mk, di] = findgroups(string(Suse.mouse_key), Suse.day_index);
L = table(string(mk), double(di), 'VariableNames',{'mouse_key','day_index'});

% Attach cohort info (if present in Suse)
if ismember('Group',Suse.Properties.VariableNames)
    grp = splitapply(@(x) string(x(1)), string(Suse.Group), g);
    L.Group = categorical(grp, ["Active","Passive","Unknown"]);
else
    L.Group = categorical(repmat("Unknown",height(L),1), ["Active","Passive","Unknown"]);
end

for k = 1:numel(lickVars)
    v = lickVars{k};
    if ismember(v, Suse.Properties.VariableNames)
        L.(v) = splitapply(@(x) mean(x,'omitnan'), Suse.(v), g);
    else
        L.(v) = nan(height(L),1);
    end
end

% Add per-mouse cumulative lick count across days (exclude day1-2 already via daysUse)
L = L(ismember(L.day_index,daysUse),:);
L = add_cumulative_across_days(L);

% Spaghetti across days for BOTH groups + mean lines
lickMetrics = {'lick_n','lick_cumAcrossDays','lick_freq_per_min','lick_freq_per_10s','lick_meanDur_s','lick_medianIEI_s','lick_totalDur_s'};
lickLabels  = {'Total lick count (per day)','Cumulative lick count (across days)','Lick frequency (/min)','Lick frequency (/10 s)', ...
               'Lick mean duration (s)','Lick median IEI (s)','Total lick duration (s)'};

for k=1:numel(lickMetrics)
    ycol = lickMetrics{k};
    if ~ismember(ycol,L.Properties.VariableNames), continue; end

    % all mice
    plot_spaghetti_metric(L, ycol, lickLabels{k}, ...
        sprintf('%s - all mice', lickLabels{k}), outDir, 'All', daysUse);

    % Active vs Passive panels (your request: passive blue, active red)
    plot_spaghetti_metric_groups_colored(L, ycol, lickLabels{k}, ...
        sprintf('%s: Active vs Passive', lickLabels{k}), ...
        fullfile(outDir, sprintf('spaghetti_%s_active_vs_passive.png', safeName(ycol))), ...
        daysUse);
end

end

function L2 = add_cumulative_across_days(L)
L2 = L;
if ~ismember('lick_n',L2.Properties.VariableNames)
    L2.lick_cumAcrossDays = nan(height(L2),1);
    return;
end
L2.lick_cumAcrossDays = nan(height(L2),1);
mice = unique(string(L2.mouse_key),'stable');
for i=1:numel(mice)
    idx = find(string(L2.mouse_key)==mice(i));
    di  = L2.day_index(idx);
    [~,ord] = sort(di);
    lickN = L2.lick_n(idx(ord));
    lickN(~isfinite(lickN)) = 0;
    cumN  = cumsum(lickN);
    L2.lick_cumAcrossDays(idx(ord)) = cumN;
end
end

%% ======================================================================
%                REWARD/LICK EVENTS + CUMULATIVE CURVES
% ======================================================================
function [timeGridMin, meanCumPerDay, meanCumEarly, meanCumLate, sessInfo] = ...
    compute_cumulative_and_events_active(T, D, daysUse)

% Active mice from D cohort mapping
if ismember('Group', D.Properties.VariableNames)
    activeMice = unique(string(D.mouse_key(D.Group=="Active")), 'stable');
else
    activeMice = string.empty(0,1);
end

mkT   = string(T.mouse_key);
dayIx = T.day_index;
T = T( ismember(mkT, activeMice) & ismember(dayIx, daysUse), :);

if isempty(T)
    timeGridMin   = [];
    meanCumPerDay = [];
    meanCumEarly  = [];
    meanCumLate   = [];
    sessInfo      = struct([]);
    return;
end

T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
T.Injector_TTL(isnan(T.Injector_TTL)) = 0;

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
    durMin = (max(t)-min(t))/60;
    maxDurMin = max(maxDurMin, durMin);

    L  = logical(T.Lick_TTL(idx));     L = L(good);
    R  = logical(T.Injector_TTL(idx)); R = R(good);

    sessInfo(k).mouse_key      = mk(k);
    sessInfo(k).day_index      = di(k);
    sessInfo(k).session_idx    = si(k);
    sessInfo(k).lickTimesMin   = (t(L)/60)';
    sessInfo(k).rewardTimesMin = (t(R)/60)';
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

maxDay = max(daysUse);
meanCumPerDay = nan(maxDay, numel(timeGridMin));
sumCum = zeros(maxDay, numel(timeGridMin));
nSess  = zeros(maxDay,1);

for k = 1:numel(sessInfo)
    d = sessInfo(k).day_index;
    if ~ismember(d, daysUse), continue; end
    if isempty(sessInfo(k).lickTimesMin), continue; end
    cvec = arrayfun(@(tt) sum(sessInfo(k).lickTimesMin <= tt), timeGridMin);
    sumCum(d,:) = sumCum(d,:) + cvec;
    nSess(d) = nSess(d) + 1;
end

for d = 1:maxDay
    if nSess(d) > 0
        meanCumPerDay(d,:) = sumCum(d,:) / nSess(d);
    end
end

% Keep legacy outputs (not used by new epoch plot; retained for compatibility)
meanCumEarly = nan(1,numel(timeGridMin));
meanCumLate  = nan(1,numel(timeGridMin));
end

function raster_day_plot(sessInfo, dayIdx, whichType, outDir)
mask = [sessInfo.day_index] == dayIdx;
sess = sessInfo(mask);
if isempty(sess), return; end

fh = figure('Color','w','Position',[80 80 900 520]); hold on;

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
ylabel('Mouse (session index within day)');
title(sprintf('%s raster - Active mice - Day %d', upperFirst(whichType), dayIdx), 'Interpreter','none');
grid on; box on;

set(gca,'YTick',1:numel(sess), 'YTickLabel',labels);
ylim([0.5, numel(sess)+0.5]);

fn = fullfile(outDir, sprintf('raster_%s_day%d.png', whichType, dayIdx));
printpng(fh, fn); close(fh);
end

function plot_cumLicks_allDays(timeGridMin, meanCumPerDay, daysUse, outDir)
if isempty(timeGridMin) || isempty(meanCumPerDay), return; end
fh = figure('Color','w','Position',[80 80 980 560]); hold on;

daysUse = unique(daysUse);
cols = lines(numel(daysUse));

for i=1:numel(daysUse)
    d = daysUse(i);
    if d > size(meanCumPerDay,1), continue; end
    y = meanCumPerDay(d,:);
    if all(~isfinite(y)), continue; end
    plot(timeGridMin, y, 'Color', cols(i,:), 'LineWidth',2, 'DisplayName',sprintf('Day %d',d));
end

xlabel('Session time (min)');
ylabel('Mean cumulative licks (Active mice)');
title('Cumulative licking across session time - Active mice (per day)','Interpreter','none');
grid on; box on;
legend('Location','bestoutside');

fn = fullfile(outDir,'cumLicks_allDays_active.png');
printpng(fh, fn); close(fh);
end

function plot_cumLicks_epoch_groups(timeGridMin, meanCumPerDay, epochs, outDir)
if isempty(timeGridMin) || isempty(meanCumPerDay), return; end

fh = figure('Color','w','Position',[80 80 980 560]); hold on;

for k=1:numel(epochs)
    dd = epochs(k).Days;
    dd = dd(dd<=size(meanCumPerDay,1));
    Y = meanCumPerDay(dd,:);
    y = mean(Y,1,'omitnan');
    if all(~isfinite(y)), continue; end
    plot(timeGridMin, y, 'LineWidth',2.5, 'DisplayName',char(epochs(k).Name));
end

xlabel('Session time (min)');
ylabel('Mean cumulative licks (Active mice)');
title('Cumulative licking by epoch (Active mice)','Interpreter','none');
grid on; box on;
legend('Location','bestoutside');

fn = fullfile(outDir,'cumLicks_byEpoch_active.png');
printpng(fh, fn); close(fh);
end

%% ======================================================================
%                           SPAGHETTI PLOTS
% ======================================================================
function plot_spaghetti_metric(D, ycol, ylab, ttl, outDir, modeStr, restrictDays)
% modeStr: 'All' or 'ActiveOnly'
if ~ismember(ycol, D.Properties.VariableNames), return; end
D = D(ismember(D.day_index, restrictDays), :);

if strcmpi(modeStr,'ActiveOnly') && ismember('Group', D.Properties.VariableNames)
    D = D(D.Group=="Active", :);
end
if isempty(D), return; end

mice = unique(string(D.mouse_key),'stable');

fh = figure('Color','w','Position',[80 80 900 520]); hold on;
for i = 1:numel(mice)
    di = D(string(D.mouse_key)==mice(i),:);
    di = di(isfinite(di.(ycol)),:);
    if isempty(di), continue; end
    [~,ord] = sort(di.day_index);
    plot(di.day_index(ord), di.(ycol)(ord),'-o', 'LineWidth',1.0,'MarkerSize',4);
end

% mean across mice per day
G = groupsummary(D(:,{'day_index',ycol}), 'day_index','mean', ycol);
plot(G.day_index, G.(['mean_' ycol]), 'k-','LineWidth',3);

xlabel('Day index');
ylabel(ylab,'Interpreter','none');
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

fh = figure('Color','w','Position',[80 80 900 920]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

groups = ["Active","Passive"];
for gi = 1:2
    nexttile; hold on;
    if ~ismember('Group', D.Properties.VariableNames)
        text(0.5,0.5,'No Group column','Units','normalized');
        axis off; continue;
    end

    Gd = D(D.Group==groups(gi),:);
    if isempty(Gd)
        text(0.5,0.5,sprintf('No %s mice',groups(gi)),'Units','normalized');
        axis off; continue;
    end

    mice = unique(string(Gd.mouse_key),'stable');
    labX = nan(numel(mice),1);
    labY = nan(numel(mice),1);

    for i = 1:numel(mice)
        di = Gd(string(Gd.mouse_key)==mice(i),:);
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

    Sg = groupsummary(Gd(:,{'day_index',ycol}),'day_index','mean',ycol);
    plot(Sg.day_index, Sg.(['mean_' ycol]), 'k-','LineWidth',3);

    if labelLines
        for i = 1:numel(mice)
            if ~isfinite(labX(i)), continue; end
            text(labX(i)+0.15, labY(i), char(mice(i)), 'FontSize',8, 'Interpreter','none');
        end
    end

    xlabel('Day index');
    ylabel(ylab,'Interpreter','none');
    title(sprintf('%s - %s mice', ylab, groups(gi)),'Interpreter','none');
    grid on; box on;
    xlim([min(restrictDays)-0.5, max(restrictDays)+0.5]);
end

sgtitle(ttl,'FontWeight','bold');
printpng(fh, outFile); close(fh);
end

function plot_spaghetti_metric_groups_colored(D, ycol, ylab, ttl, outFile, restrictDays)
% Passive = blue, Active = red (as you requested)
if ~ismember(ycol, D.Properties.VariableNames), return; end
if ~ismember('Group', D.Properties.VariableNames), return; end

D = D(ismember(D.day_index,restrictDays),:);

fh = figure('Color','w','Position',[80 80 900 920]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

groups = ["Active","Passive"];
fixedColor.Active  = [0.85 0.33 0.10]; % red-ish
fixedColor.Passive = [0.00 0.45 0.74]; % blue-ish

for gi = 1:2
    nexttile; hold on;
    Gd = D(D.Group==groups(gi),:);
    if isempty(Gd)
        text(0.5,0.5,sprintf('No %s mice',groups(gi)),'Units','normalized');
        axis off; continue;
    end

    mice = unique(string(Gd.mouse_key),'stable');
    for i=1:numel(mice)
        di = Gd(string(Gd.mouse_key)==mice(i),:);
        di = di(isfinite(di.(ycol)),:);
        if isempty(di), continue; end
        [~,ord] = sort(di.day_index);
        plot(di.day_index(ord), di.(ycol)(ord), '-o', ...
            'LineWidth',1.0,'MarkerSize',4, 'Color', fixedColor.(char(groups(gi))));
    end

    Sg = groupsummary(Gd(:,{'day_index',ycol}),'day_index','mean',ycol);
    plot(Sg.day_index, Sg.(['mean_' ycol]), 'k-','LineWidth',3);

    xlabel('Day index');
    ylabel(ylab,'Interpreter','none');
    title(sprintf('%s - %s mice', ylab, groups(gi)),'Interpreter','none');
    grid on; box on;
    xlim([min(restrictDays)-0.5, max(restrictDays)+0.5]);
end

sgtitle(ttl,'FontWeight','bold');
printpng(fh, outFile); close(fh);
end

%% ======================================================================
%                       ACTIVE vs PASSIVE PER-DAY STATS
% ======================================================================
function compare_AP_means_per_day_with_stars(D, mName, daysAll, outDir)
fprintf('compare_AP_means_per_day_with_stars: %s\n', mName);

if ~ismember('Group', D.Properties.VariableNames), return; end
if ~ismember(mName, D.Properties.VariableNames), return; end

nDays = numel(daysAll);
meanA = nan(nDays,1); meanP = nan(nDays,1);
semA  = nan(nDays,1); semP  = nan(nDays,1);
nA    = nan(nDays,1); nP    = nan(nDays,1);
p_raw = nan(nDays,1);

for i = 1:nDays
    d  = daysAll(i);
    xa = D.(mName)(D.day_index==d & D.Group=="Active");
    xp = D.(mName)(D.day_index==d & D.Group=="Passive");
    xa = xa(isfinite(xa)); xp = xp(isfinite(xp));

    if numel(xa) >= 2 && numel(xp) >= 2
        nA(i) = numel(xa); nP(i) = numel(xp);
        meanA(i) = mean(xa,'omitnan'); meanP(i) = mean(xp,'omitnan');
        semA(i)  = std(xa,'omitnan')/sqrt(nA(i));
        semP(i)  = std(xp,'omitnan')/sqrt(nP(i));
        p_raw(i) = ranksum(xa,xp);
    end
end

% BH-FDR
p_adj = nan(size(p_raw));
valid = find(isfinite(p_raw) & p_raw>0);
if ~isempty(valid)
    [p_sorted, ord] = sort(p_raw(valid));
    m = numel(p_sorted);
    q = p_sorted .* (m ./ (1:m)');
    q(q>1)=1;
    for k=m-1:-1:1, q(k)=min(q(k),q(k+1)); end
    p_adj(valid(ord)) = q;
end

star = strings(nDays,1);
for i=1:nDays
    if ~isfinite(p_adj(i)), continue; end
    if      p_adj(i) < 0.001, star(i)="***";
    elseif  p_adj(i) < 0.01,  star(i)="**";
    elseif  p_adj(i) < 0.05,  star(i)="*";
    else,                     star(i)="n.s.";
    end
end

resTbl = table(daysAll(:), nA, meanA, semA, nP, meanP, semP, p_raw, p_adj, star, ...
    'VariableNames', {'day_index','N_Active','Mean_Active','SEM_Active', ...
                      'N_Passive','Mean_Passive','SEM_Passive','p_raw','p_adj','Star'});
writetable(resTbl, fullfile(outDir, sprintf('pvals_A_vs_P_%s_ranksum.csv', safeName(mName))));

% Figure
fh = figure('Color','w','Position',[120 120 780 520]); hold on;

ha = errorbar(daysAll, meanA, semA, '-o', 'LineWidth',2, 'MarkerSize',5);
hp = errorbar(daysAll, meanP, semP, '-o', 'LineWidth',2, 'MarkerSize',5);
set(ha,'Color',[0.85 0.33 0.10]); % Active red
set(hp,'Color',[0.00 0.45 0.74]); % Passive blue

xlabel('Day index');
ylabel(strrep(mName,'_',' '),'Interpreter','none');
title(sprintf('%s: Active vs Passive per-day (ranksum + BH-FDR)', mName),'Interpreter','none');
grid on; box on;
legend({'Active','Passive'},'Location','best');

yAll = [meanA-semA; meanA+semA; meanP-semP; meanP+semP];
yAll = yAll(isfinite(yAll));
if isempty(yAll), yMin=0; yMax=1; else, yMin=min(yAll); yMax=max(yAll); end
yRange = max(yMax-yMin, eps);
ylim([yMin-0.05*yRange, yMax+0.35*yRange]);

for i=1:nDays
    if ~isfinite(p_raw(i)) || ~isfinite(p_adj(i)), continue; end
    pStr = tern(p_raw(i)<0.001,'p<0.001',sprintf('p=%.3f',p_raw(i)));
    qStr = tern(p_adj(i)<0.001,'q<0.001',sprintf('q=%.3f',p_adj(i)));
    labelStr = sprintf('%s\n%s, %s', star(i), pStr, qStr);
    yStar = max([meanA(i)+semA(i), meanP(i)+semP(i)]) + 0.08*yRange;
    text(daysAll(i), yStar, labelStr, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',9, ...
        'FontWeight', tern(p_adj(i)<0.05,'bold','normal'));
end

printpng(fh, fullfile(outDir, sprintf('means_A_vs_P_%s.png', safeName(mName))));
close(fh);
end

%% ======================================================================
%                         SESSION/DAY METRICS (FAST)
% ======================================================================
function [S, D] = fast_session_day_metrics_basic(T, runDir)
cacheMat = fullfile(runDir, 'S_D_cache_basic_REVISED.mat');
if exist(cacheMat,'file')
    L = load(cacheMat,'S','D','cache_hash');
    if isfield(L,'cache_hash') && isequal(L.cache_hash, local_hash(T))
        fprintf('Loaded S/D from cache: %s\n', cacheMat);
        S = L.S; D = L.D; return;
    end
end

need = {'mouse_key','day_index','day_name','session_idx','Diameter_px', ...
        'Lick_TTL','Injector_TTL','CamTime_rel_s','PupilTimestamp_s', ...
        'CamTime_s','PlotTime_s_30fps','Session_Paradigm','RequirementLast','Immersion_Latency_s'};
V = T.Properties.VariableNames;

keepExtras = V( contains(V,'TST_','IgnoreCase',true) | contains(V,'HOT_','IgnoreCase',true) );
keepExtras = keepExtras(~endsWith(keepExtras,'_File'));

need = unique([need, keepExtras], 'stable');
need = intersect(need, V, 'stable');
T = T(:, need);

T.mouse_key = categorical(string(T.mouse_key));
T.day_name  = categorical(string(T.day_name));
if ismember('Session_Paradigm',T.Properties.VariableNames)
    T.Session_Paradigm = string(T.Session_Paradigm);
end

T.Lick_TTL = double(T.Lick_TTL)>0.5;
T.Injector_TTL = double(T.Injector_TTL)>0.5;

[g, keys_mouse, keys_day, keys_dayname, keys_sess] = findgroups(T.mouse_key, T.day_index, T.day_name, T.session_idx);
nG = max(g);
fprintf('Computing per-session metrics for %d sessions...\n', nG);

S = table();
S.mouse_key   = string(removecats(keys_mouse));
S.day_index   = double(keys_day);
S.day_name    = removecats(keys_dayname);
S.session_idx = double(keys_sess);

S.RequirementLast = nan(nG,1);
S.SessionMinutes  = nan(nG,1);
S.Session_Paradigm = strings(nG,1);

vars = {'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s','lick_freq_per_10s', ...
        'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s','pupil_mean'};
for v = vars, S.(v{1}) = nan(nG,1); end
for j=1:numel(keepExtras), S.(keepExtras{j}) = nan(nG,1); end

tb_all = pickTimebase_fast(T);
step = max(1, floor(nG/40));

for k=1:nG
    idx = (g==k);
    tb  = tb_all(idx);
    dur = finiteRange_fast(tb);
    S.SessionMinutes(k) = dur/60;

    if ismember('RequirementLast',T.Properties.VariableNames)
        S.RequirementLast(k) = mean(double(T.RequirementLast(idx)),'omitnan');
    end

    if ismember('Session_Paradigm',T.Properties.VariableNames)
        sp = string(T.Session_Paradigm(idx));
        sp = sp(~ismissing(sp));
        if isempty(sp), S.Session_Paradigm(k) = ""; else, S.Session_Paradigm(k) = sp(1); end
    end

    if ismember('Diameter_px',T.Properties.VariableNames)
        S.pupil_mean(k) = mean(double(T.Diameter_px(idx)),'omitnan');
    end

    if ismember('Lick_TTL',T.Properties.VariableNames)
        [n,md,td,iei] = eventMetrics_fast(tb, logical(T.Lick_TTL(idx)));
        S.lick_n(k)=n; S.lick_meanDur_s(k)=md; S.lick_totalDur_s(k)=td; S.lick_medianIEI_s(k)=iei;
        if S.SessionMinutes(k)>0
            S.lick_freq_per_min(k) = n / S.SessionMinutes(k);
            S.lick_freq_per_10s(k) = n / (S.SessionMinutes(k)*6);
        end
    end

    if ismember('Injector_TTL',T.Properties.VariableNames)
        [n,md,td,iri] = eventMetrics_fast(tb, logical(T.Injector_TTL(idx)));
        S.rew_n(k)=n; S.rew_meanDur_s(k)=md; S.rew_totalDur_s(k)=td; S.rew_medianIRI_s(k)=iri;
        if S.SessionMinutes(k)>0
            S.rew_freq_per_min(k)=n/S.SessionMinutes(k);
        end
    end

    if ismember('Immersion_Latency_s',T.Properties.VariableNames)
        % If present per-sample, average it within session
        S.Immersion_Latency_s(k) = mean(double(T.Immersion_Latency_s(idx)),'omitnan');
    end

    for j=1:numel(keepExtras)
        col = keepExtras{j};
        S.(col)(k) = mean(double(T.(col)(idx)),'omitnan');
    end

    if mod(k,step)==0
        fprintf('  %d/%d (%.0f%%)\n', k, nG, 100*k/nG);
    end
end

fprintf('Collapsing to per-day medians...\n');
[g2, mk2, di2, dn2] = findgroups(categorical(S.mouse_key), S.day_index, S.day_name);
D = table(string(removecats(mk2)), double(di2), removecats(dn2), ...
          'VariableNames',{'mouse_key','day_index','day_name'});

baseList = [{'RequirementLast','Immersion_Latency_s','SessionMinutes'}, vars, keepExtras];
baseList = intersect(unique(baseList,'stable'), S.Properties.VariableNames, 'stable');

for v = baseList
    D.(v{1}) = splitapply(@(x) median(x,'omitnan'), S.(v{1}), g2);
end

cache_hash = local_hash(T);
try
    save(cacheMat,'S','D','cache_hash','-v7.3');
    fprintf('Cached S/D to: %s\n', cacheMat);
catch
end
end

%% ======================================================================
%                               SMALL HELPERS
% ======================================================================
function T2 = ensureString(T2, nm)
if ~ismember(nm, T2.Properties.VariableNames)
    T2.(nm) = repmat("",height(T2),1);
    return;
end
if ~isstring(T2.(nm))
    T2.(nm) = string(T2.(nm));
end
end

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
t = double(t(:)); ttl = logical(ttl(:));
good = isfinite(t) & ~isnan(ttl);
t=t(good); ttl=ttl(good);
if numel(t)<2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
dt  = diff(t);
md = median(dt(isfinite(dt)));
if ~isfinite(md), md=1/30; end
t(2:end) = max(t(2:end), t(1:end-1)+md*0.5);
d = diff([false; ttl; false]);
on = find(d==1);
off = find(d==-1)-1;
n = numel(on);
if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
edges = [t; t(end)+md];
segDur = sum(edges(off+1) - edges(on), 2, 'omitnan');
meanDur  = mean(segDur,'omitnan');
totalDur = sum(segDur,'omitnan');
if n>=2, medianIEI = median(diff(t(on)),'omitnan'); else, medianIEI = NaN; end
end

function h = local_hash(T)
% Hash for caching
try
    mk = string(T.mouse_key);
    mk = canonicalize_mouse_keys(mk);
    % crude but stable-ish
    raw = uint64(sum(double(char(join(mk,'|')))));
    raw = raw + uint64(sum(uint64(T.day_index)));
    if ismember('session_idx',T.Properties.VariableNames)
        raw = raw + uint64(sum(uint64(T.session_idx)));
    end
    h = raw;
catch
    h = uint64(now*1e5);
end
end

function s = safeName(nm)
s = regexprep(string(nm),'[^a-zA-Z0-9]+','_');
s = char(s);
end

function printpng(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',180);
catch
    print(fh, fn, '-dpng','-r180');
end
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = upperFirst(s)
s = string(s);
if strlength(s)==0, s=""; return; end
s = lower(s);
s = upper(extractBetween(s,1,1)) + extractAfter(s,1);
end

function v = getOr(T, nm, defaultVal)
if ismember(nm,T.Properties.VariableNames)
    v = T.(nm);
else
    v = repmat(defaultVal,height(T),1);
end
end

function v = getOrNum(T, nm, defaultVal)
if ismember(nm,T.Properties.VariableNames)
    v = T.(nm);
else
    v = repmat(defaultVal,height(T),1);
end
end
