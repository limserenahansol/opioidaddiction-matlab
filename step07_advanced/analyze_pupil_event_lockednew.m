function analyze_pupil_event_lockednew()
% Analyze pupil dynamics aligned to reward / lick events, using explicit cohort mapping.
%
% What this script does:
%  - Loads latest run_* CSV
%  - Attaches explicit cohort info (Group/Sex/Color/Cage/PairID) to each row
%  - Applies death rule: 6099 white only day 1-13
%  - Computes reward-locked and lick-locked pupil traces
%  - Summarizes per epoch (Pre/During/Post/Withdrawal/Reexposure)
%  - Runs 2 modes:
%       A_includeSwitchDays (use all days 3-18)
%       B_excludeSwitchDays (remove unreliable switch days [4 6 11 14])
%
% Output:
%  runDir/figs/pupil_event_locked_newcohort/A_includeSwitchDays
%  runDir/figs/pupil_event_locked_newcohort/B_excludeSwitchDays

%% ---------- parameters ----------
preWin_s   = 2;      % seconds before event
postWin_s  = 2;      % seconds after event
bin_s      = 0.05;    % time bin for interpolation (s)
baselineMode = "zscore";   % "subtract" or "zscore"
% Your epoch definition
epochs = make_epochs();

% Switch days to optionally exclude
unreliableDays = [4 6 11 14];

% Restrict analysis days (ignore habituation day 1-2)
daysAll = 3:18;

%% ---------- locate latest run + load CSV ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder,d(ix).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('Reading: %s\n', csvPath);

T = readtable(csvPath,'VariableNamingRule','preserve');
T = ensureStringCol(T,'mouse_key');
T = ensureStringCol(T,'day_name');

if ~ismember('day_index',T.Properties.VariableNames)
    error('day_index column not found in table.');
end
if ~ismember('session_idx',T.Properties.VariableNames)
    T.session_idx = ones(height(T),1);
end

% Required columns
needCols = {'Diameter_px','Lick_TTL','Injector_TTL'};
for i = 1:numel(needCols)
    if ~ismember(needCols{i},T.Properties.VariableNames)
        error('Missing column %s in table.', needCols{i});
    end
end
T.Lick_TTL(isnan(T.Lick_TTL))         = 0;
T.Injector_TTL(isnan(T.Injector_TTL)) = 0;

%% ---------- build explicit cohort map (NO inference) ----------
cohort = build_cohort_map_new();

% Attach cohort info to T via normalized mouse key
T = attach_cohort_info(T, cohort);

% Enforce death rule (6099_white only day<=13)
T = apply_maxday_rule(T);

% Restrict to analysis days 3-18
T = T(ismember(T.day_index, daysAll), :);

%% ---------- output root ----------
outRoot = fullfile(runDir,'figs','pupil_event_locked_newcohort');
if ~exist(outRoot,'dir'), mkdir(outRoot); end

%% ---------- run two modes ----------
run_mode(T, epochs, unreliableDays, true,  fullfile(outRoot,'A_includeSwitchDays'), preWin_s, postWin_s, bin_s);
run_mode(T, epochs, unreliableDays, false, fullfile(outRoot,'B_excludeSwitchDays'), preWin_s, postWin_s, bin_s);

fprintf('Done. Outputs under:\n  %s\n', outRoot);

end % ===== end main =====


%% =======================================================================
%                               MODE RUNNER
%% =======================================================================

function run_mode(T, epochs, unreliableDays, includeSwitchDays, outDir, preWin_s, postWin_s, bin_s)

if ~exist(outDir,'dir'), mkdir(outDir); end

if includeSwitchDays
    fprintf('\n=============================\nMODE: INCLUDE switch days\n');
else
    fprintf('\n=============================\nMODE: EXCLUDE switch days\n');
end
fprintf('unreliableDays = [%s]\n', num2str(unreliableDays));
fprintf('Output: %s\n', outDir);
fprintf('=============================\n');

Tuse = T;

% Exclude unreliable switch days if requested
if ~includeSwitchDays
    Tuse = Tuse(~ismember(Tuse.day_index, unreliableDays), :);
end

% If after filtering there is no data, stop
if isempty(Tuse)
    warning('No data after filtering; skipping mode output: %s', outDir);
    return;
end

% Compute event-locked traces
fprintf('Computing reward-locked pupil...\n');
[E_rew, meta_rew, tGrid] = extract_pupil_event_locked(Tuse, 'reward', preWin_s, postWin_s, bin_s);

fprintf('Computing lick-locked pupil...\n');
[E_lick, meta_lick, tGrid2] = extract_pupil_event_locked(Tuse, 'lick', preWin_s, postWin_s, bin_s);

% Plot epoch summaries: Active vs Passive per epoch
plot_epoch_summary(E_rew,  meta_rew,  tGrid,  epochs, 'Reward', outDir);
plot_epoch_summary(E_lick, meta_lick, tGrid2, epochs, 'Lick',   outDir);

% Optional: Pair-wise overlays (Active vs Passive within each PairID)
plot_pair_overlays(E_rew,  meta_rew,  tGrid,  epochs, 'Reward', outDir);
plot_pair_overlays(E_lick, meta_lick, tGrid2, epochs, 'Lick',   outDir);

end


%% =======================================================================
%                              COHORT / EPOCHS
%% =======================================================================

function epochs = make_epochs()
% Your definition:
%  Pre        = day 3-5  (water PR)
%  During     = day 6-10 (passive training morphine; active has PR, passive replay)
%  Post       = day 11-13 (all morphine PR freely)
%  Withdrawal = day 14-16 (water, no morphine)
%  Reexposure = day 17-18 (morphine PR)

epochs = struct('Name',{},'Days',{});
epochs(1).Name = "Pre";        epochs(1).Days = 3:5;
epochs(2).Name = "During";     epochs(2).Days = 6:10;
epochs(3).Name = "Post";       epochs(3).Days = 11:13;
epochs(4).Name = "Withdrawal"; epochs(4).Days = 14:16;
epochs(5).Name = "Reexposure"; epochs(5).Days = 17:18;
end


function cohort = build_cohort_map_new()
% Explicit cohort list you provided:
% 6100red     f_s (treat as Passive per your pairing note)
% 6100orange  f_p
% 6100black   f_a
% 0911red     f_a
% 0911orange  f_p
% 0911black   f_p
% 0911white   f_a
% 0910red     m_p
% 0910orange  m_p
% 0910black   m_a
% 6099red     m_p
% 6099orange  m_a
% 6099black   m_a
% (6099white passive died day1-13) -> include for completeness if it exists in CSV

% Define rows: Cage, Color, Sex, Group, PairID, MaxDay(optional)
rows = {
    "6100","red",    "F","Passive", 1, 18
    "6100","orange", "F","Passive", 1, 18
    "6100","black",  "F","Active",  1, 18

    "0911","red",    "F","Active",  2, 18
    "0911","orange", "F","Passive", 2, 18

    "0911","black",  "F","Passive", 3, 18
    "0911","white",  "F","Active",  3, 18

    "0910","black",  "M","Active",  4, 18
    "0910","orange", "M","Passive", 4, 18
    "0910","red",    "M","Passive", 4, 18

    "6099","orange", "M","Active",  5, 18
    "6099","red",    "M","Passive", 5, 18

    "6099","black",  "M","Active",  6, 18
    "6099","white",  "M","Passive", 6, 13  % died day 1-13
    };

cohort = cell2table(rows, 'VariableNames', {'Cage','Color','Sex','Group','PairID','MaxDay'});

% Create canonical mouse_key (with underscore) and normalized key for joining
cohort.mouse_key = cohort.Cage + "_" + cohort.Color;

cohort.mouse_key_norm = normalize_key(cohort.mouse_key);

% Also store alternative format (no underscore) normalization deterministically
cohort.mouse_key_alt = cohort.Cage + cohort.Color;  % e.g., 6100black
cohort.mouse_key_alt_norm = normalize_key(cohort.mouse_key_alt);
end


function T = attach_cohort_info(T, cohort)
% Deterministic join:
% - create normalized key for T.mouse_key
% - match either underscore version or no-underscore version

T.mouse_key = string(T.mouse_key);
T.mouse_key_norm = normalize_key(T.mouse_key);

% Build lookup table where both canonical & alt norms point to same metadata
L1 = cohort(:, {'mouse_key_norm','Cage','Color','Sex','Group','PairID','MaxDay'});
L1.Properties.VariableNames{1} = 'join_key_norm';

L2 = cohort(:, {'mouse_key_alt_norm','Cage','Color','Sex','Group','PairID','MaxDay'});
L2.Properties.VariableNames{1} = 'join_key_norm';

L = [L1; L2];
% Deduplicate by join_key_norm (keep first; both rows identical metadata by construction)
[~, ia] = unique(L.join_key_norm, 'stable');
L = L(ia,:);

% Join
T = outerjoin(T, L, 'LeftKeys','mouse_key_norm', 'RightKeys','join_key_norm', ...
    'MergeKeys', true, 'Type','left');

% Validate: if some rows missing cohort info, keep them but label Unknown
if ~ismember('Group', T.Properties.VariableNames)
    T.Group = repmat("Unknown", height(T), 1);
else
    gg = string(T.Group);
    gg(ismissing(gg)) = "Unknown";
    T.Group = gg;
end
if ~ismember('Sex', T.Properties.VariableNames)
    T.Sex = repmat("", height(T), 1);
else
    ss = string(T.Sex);
    ss(ismissing(ss)) = "";
    T.Sex = ss;
end
if ~ismember('Color', T.Properties.VariableNames)
    T.Color = repmat("", height(T), 1);
else
    cc = string(T.Color);
    cc(ismissing(cc)) = "";
    T.Color = cc;
end
if ~ismember('Cage', T.Properties.VariableNames)
    T.Cage = repmat("", height(T), 1);
else
    cg = string(T.Cage);
    cg(ismissing(cg)) = "";
    T.Cage = cg;
end
if ~ismember('PairID', T.Properties.VariableNames)
    T.PairID = nan(height(T), 1);
else
    pid = double(T.PairID);
    T.PairID = pid;
end
if ~ismember('MaxDay', T.Properties.VariableNames)
    T.MaxDay = nan(height(T), 1);
else
    md = double(T.MaxDay);
    T.MaxDay = md;
end

end


function T = apply_maxday_rule(T)
% Drops rows beyond MaxDay for mice with a finite MaxDay.

if ~ismember('MaxDay', T.Properties.VariableNames)
    return;
end

md = double(T.MaxDay);
keep = true(height(T),1);

finiteMask = isfinite(md) & md > 0;
keep(finiteMask) = T.day_index(finiteMask) <= md(finiteMask);

T = T(keep,:);
end


function k = normalize_key(x)
% Deterministic normalization for joining:
% - lowercase
% - remove all non-alphanumeric
x = string(x);
x = lower(x);
k = regexprep(x, '[^a-z0-9]', '');
end


%% =======================================================================
%                             CORE EXTRACTION
%% =======================================================================

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


function [E, meta, tGrid] = extract_pupil_event_locked(T, whichType, preWin_s, postWin_s, bin_s)
% whichType: 'reward' or 'lick'
% E: nEvents x nTime baseline-subtracted pupil
% meta: struct per event with Group/day_index/mouse_key/PairID/Sex/Color/Cage

switch lower(whichType)
    case 'reward'
        evCol = 'Injector_TTL';
    case 'lick'
        evCol = 'Lick_TTL';
    otherwise
        error('Unknown whichType: %s', whichType);
end

tb_all = pickTimebase_fast(T);
if all(isnan(tb_all))
    error('No valid time base found in table.');
end

tGrid = -preWin_s:bin_s:postWin_s;
nT = numel(tGrid);

E = nan(0,nT);
meta = struct('Group',{},'day_index',{},'mouse_key',{},'PairID',{},'Sex',{},'Color',{},'Cage',{});

baseStart = -2;
baseEnd   = 0;

% group by mouse/day/session
[g, mk, di, si] = findgroups(string(T.mouse_key), T.day_index, T.session_idx); %#ok<ASGLU>
nG = max(g);

for k = 1:nG
    idx = (g==k);

    tb    = double(tb_all(idx));
    diam  = double(T.Diameter_px(idx));
    ttl   = double(T.(evCol)(idx));

    % Session-level metadata (take first non-missing deterministically)
    mk_str = string(mk(k));
    day_k  = double(di(k));

    groupStr = get_first_string(T, idx, 'Group', "Unknown");
    sexStr   = get_first_string(T, idx, 'Sex',   "");
    colStr   = get_first_string(T, idx, 'Color', "");
    cageStr  = get_first_string(T, idx, 'Cage',  "");
    pairID   = get_first_double(T, idx, 'PairID', NaN);

    good = isfinite(tb) & isfinite(diam) & ~isnan(ttl);
    tb   = tb(good);
    diam = diam(good);
    ttl  = ttl(good) > 0.5;

    if numel(tb) < 10
        continue;
    end

    % enforce non-decreasing time deterministically
    dt = diff(tb);
    md = median(dt(isfinite(dt)));
    if isfinite(md) && md > 0
        tb(2:end) = max(tb(2:end), tb(1:end-1) + 0.5*md);
    end

    % event onset indices (rising edges)
    dTTL = diff([false; ttl(:)]);
    onIdx = find(dTTL==1);
    if isempty(onIdx), continue; end

    for j = 1:numel(onIdx)
        t0 = tb(onIdx(j));
        rel = tb - t0;

        inWin = rel>=-preWin_s & rel<=postWin_s;
        if nnz(inWin) < 3
            continue;
        end

        thisT = rel(inWin);
        thisY = diam(inWin);
        % after you compute thisT, thisY (already in the event window)
baseMask = (thisT >= baseStart) & (thisT <= baseEnd);

if nnz(baseMask) >= 5
    mu = mean(thisY(baseMask), 'omitnan');
    sd = std(thisY(baseMask),  'omitnan');
else
    % fallback: use all pre-event samples
    preMask = thisT < 0;
    if nnz(preMask) >= 5
        mu = mean(thisY(preMask), 'omitnan');
        sd = std(thisY(preMask),  'omitnan');
    else
        mu = mean(thisY, 'omitnan');
        sd = std(thisY, 'omitnan');
    end
end

if ~isfinite(sd) || sd < 1e-6
    sd = NaN; % avoid exploding z-score
end

switch lower(string(baselineMode))
    case "zscore"
        thisY = (thisY - mu) ./ sd;
    otherwise % "subtract"
        thisY = thisY - mu;
end

        % Baseline [-2,0] s, else fallback to any pre (<0), else all
        baseMask = thisT >= baseStart & thisT <= baseEnd;
        if nnz(baseMask) >= 3
            baseVal = mean(thisY(baseMask),'omitnan');
        else
            preMask = thisT < 0;
            if nnz(preMask) >= 3
                baseVal = mean(thisY(preMask),'omitnan');
            else
                baseVal = mean(thisY,'omitnan');
            end
        end
        thisY = thisY - baseVal;

        yInterp = interp1(thisT, thisY, tGrid, 'linear', NaN);
        if all(isnan(yInterp)), continue; end

        E(end+1,:) = yInterp; %#ok<AGROW>
        meta(end+1).Group      = groupStr; %#ok<AGROW>
        meta(end).day_index    = day_k;
        meta(end).mouse_key    = mk_str;
        meta(end).PairID       = pairID;
        meta(end).Sex          = sexStr;
        meta(end).Color        = colStr;
        meta(end).Cage         = cageStr;
    end
end

end


function s = get_first_string(T, idx, col, defaultVal)
if ~ismember(col, T.Properties.VariableNames)
    s = string(defaultVal); return;
end
v = string(T.(col)(idx));
v = v(~ismissing(v));
if isempty(v), s = string(defaultVal); else, s = v(1); end
end

function x = get_first_double(T, idx, col, defaultVal)
if ~ismember(col, T.Properties.VariableNames)
    x = defaultVal; return;
end
v = double(T.(col)(idx));
v = v(isfinite(v));
if isempty(v), x = defaultVal; else, x = v(1); end
end


%% =======================================================================
%                                PLOTTING
%% =======================================================================

function plot_epoch_summary(E, meta, tGrid, epochs, labelStr, outDir)
% For each epoch, plot Active vs Passive mean±SEM on the same axes.
% Creates one figure with 5 subplots (one per epoch).

if isempty(E)
    warning('No event-locked data for %s; skipping.', labelStr);
    return;
end

nEv = size(E,1);
G   = strings(nEv,1);
day = nan(nEv,1);

for i = 1:nEv
    G(i)   = string(meta(i).Group);
    day(i) = double(meta(i).day_index);
end

% Colors (fixed)
colActive  = [0.85 0.30 0.30];
colPassive = [0.20 0.60 0.90];

fh = figure('Color','w','Position',[120 80 900 1200]);
tiledlayout(numel(epochs),1,'TileSpacing','compact','Padding','compact');

for ei = 1:numel(epochs)
    nexttile; hold on;

    dd = epochs(ei).Days;
    maskEpoch = ismember(day, dd);

    % Active
    maskA = maskEpoch & strcmpi(G,'Active');
    plot_mean_sem(E(maskA,:), tGrid, colActive);

    % Passive
    maskP = maskEpoch & strcmpi(G,'Passive');
    plot_mean_sem(E(maskP,:), tGrid, colPassive);

    xline(0,'k-','LineWidth',1);
    yline(0,'k:');

    title(sprintf('%s-locked pupil | %s (days %s)', labelStr, epochs(ei).Name, days2str(dd)), ...
        'Interpreter','none');

    ylabel('\Delta pupil (baseline-subtracted)');
    if ei==numel(epochs)
        xlabel('Time from event (s)');
    end
    grid on; box off;

    legend(make_legend(maskA,maskP), 'Location','best', 'Interpreter','none');
end

fn = fullfile(outDir, sprintf('Pupil_%s_locked_epochSummary.png', lower(labelStr)));
exportgraphics(fh, fn, 'Resolution', 300);
close(fh);

end


function L = make_legend(maskA, maskP)
na = nnz(maskA);
np = nnz(maskP);
L = {sprintf('Active (N events=%d)',na), sprintf('Passive (N events=%d)',np)};
end


function plot_mean_sem(D, tGrid, c)
% D: nEvents x nTime
if isempty(D)
    % still keep placeholder so legend counts make sense; do nothing
    return;
end
mTrace = mean(D,1,'omitnan');
sTrace = std(D,0,1,'omitnan') ./ sqrt(size(D,1));

x = tGrid(:)';
m = mTrace(:)';
s = sTrace(:)';

xx = [x fliplr(x)];
yy = [m+s fliplr(m-s)];
patch(xx, yy, c, 'FaceAlpha',0.18, 'EdgeColor','none');
plot(x, m, 'Color', c, 'LineWidth', 2);
end


function plot_pair_overlays(E, meta, tGrid, epochs, labelStr, outDir)
% For each epoch, and each PairID, overlay Active vs Passive mean traces.
% Saves one figure per epoch.

if isempty(E)
    return;
end

% Extract meta arrays
nEv = size(E,1);
pairID = nan(nEv,1);
G      = strings(nEv,1);
day    = nan(nEv,1);

for i=1:nEv
    pairID(i) = double(meta(i).PairID);
    G(i)      = string(meta(i).Group);
    day(i)    = double(meta(i).day_index);
end

uPairs = unique(pairID(isfinite(pairID)));
if isempty(uPairs)
    return;
end

colActive  = [0.85 0.30 0.30];
colPassive = [0.20 0.60 0.90];

for ei = 1:numel(epochs)
    dd = epochs(ei).Days;
    maskEpoch = ismember(day, dd);

    fh = figure('Color','w','Position',[120 120 1000 700]);
    tiledlayout(ceil(numel(uPairs)/2),2,'TileSpacing','compact','Padding','compact');

    for pi = 1:numel(uPairs)
        pid = uPairs(pi);

        nexttile; hold on;

        maskPair = maskEpoch & (pairID==pid);

        maskA = maskPair & strcmpi(G,'Active');
        maskP = maskPair & strcmpi(G,'Passive');

        plot_mean_only(E(maskA,:), tGrid, colActive);
        plot_mean_only(E(maskP,:), tGrid, colPassive);

        xline(0,'k-','LineWidth',1);
        yline(0,'k:');
        title(sprintf('%s | %s | PairID %d', labelStr, epochs(ei).Name, pid), 'Interpreter','none');
        grid on; box off;

        na = nnz(maskA); np = nnz(maskP);
        legend({sprintf('Active (N=%d)',na), sprintf('Passive (N=%d)',np)}, ...
            'Location','best', 'Interpreter','none');
    end

    fn = fullfile(outDir, sprintf('Pupil_%s_locked_pairOverlays_%s.png', ...
        lower(labelStr), safeName(char(epochs(ei).Name))));
    exportgraphics(fh, fn, 'Resolution', 250);
    close(fh);
end

end


function plot_mean_only(D, tGrid, c)
if isempty(D), return; end
mTrace = mean(D,1,'omitnan');
plot(tGrid, mTrace, 'Color', c, 'LineWidth', 2);
end


function s = days2str(v)
v = v(:)';
if isempty(v), s = ""; return; end
s = sprintf('%d-', v);
s = s(1:end-1);
end


function s = safeName(nm)
s = regexprep(string(nm),'[^a-zA-Z0-9]+','_');
s = char(s);
end


%% =======================================================================
%                               SMALL UTILS
%% =======================================================================

function T = ensureStringCol(T, nm)
if ~ismember(nm, T.Properties.VariableNames)
    T.(nm) = repmat("",height(T),1);
else
    if ~isstring(T.(nm))
        T.(nm) = string(T.(nm));
    end
end
end
