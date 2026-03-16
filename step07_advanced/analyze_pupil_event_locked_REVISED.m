function analyze_pupil_event_locked_REVISED()
% Analyze pupil dynamics aligned to reward / lick events with:
% - explicit cohort pairing map (NO guessing)
% - epochs: Pre, During, Post, Withdrawal, Reexposure
% - option to include/exclude unreliable switch days (4,6,11,14)
% - SEM shown on all plots (mean ± SEM)
% - NO discarding "unmatched events" (we only report count differences)
%
% Input:
%   latest run_* folder under:
%     K:\addiction_concate_Dec_2025\longitudinal_outputs
%   requires ALL_mice_longitudinal.csv
%
% Output:
%   runDir\figs\pupil_event_locked_REVISED\...

%% ---------------- user parameters ----------------
preWin_s   = 2;      % seconds before event
postWin_s  = 2;      % seconds after event
bin_s      = 0.2;    % interpolation bin size
useZscore  = false;  % false = baseline-subtract, true = z-score using baseline mean/std

% Epoch definition (your timeline)
epochs = struct('Name',{},'Days',{});
epochs(end+1) = struct('Name',"Pre",        'Days',3:5);
epochs(end+1) = struct('Name',"During",     'Days',6:10);
epochs(end+1) = struct('Name',"Post",       'Days',11:13);
epochs(end+1) = struct('Name',"Withdrawal", 'Days',14:16);
epochs(end+1) = struct('Name',"Reexposure", 'Days',17:18);

unreliableDays = [4 6 11 14];

% Cohort map (explicit)
% mouse_key must match what is in ALL_mice_longitudinal.csv
% If your CSV uses something else, adjust make_mouse_key() below.
M = build_cohort_map();  % table with mouse_key, Cage, Color, Sex, Group, PairID

%% ---------------- locate latest run and load ----------------
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

assert(ismember('day_index',T.Properties.VariableNames), 'day_index column not found.');
if ~ismember('session_idx',T.Properties.VariableNames)
    T.session_idx = ones(height(T),1);
end

needCols = {'Diameter_px','Lick_TTL','Injector_TTL'};
for i=1:numel(needCols)
    assert(ismember(needCols{i},T.Properties.VariableNames), 'Missing column %s.', needCols{i});
end
T.Lick_TTL(isnan(T.Lick_TTL))         = 0;
T.Injector_TTL(isnan(T.Injector_TTL)) = 0;

%% ---------------- attach cohort info (Group/PairID/Sex) ----------------
T.mouse_key = string(T.mouse_key);

% Join on mouse_key (explicit, no guessing)
T = outerjoin(T, M(:,{'mouse_key','Group','PairID','Sex','Cage','Color'}), ...
    'Keys','mouse_key','MergeKeys',true,'Type','left');

% Hard fail if any mouse_key missing from map (you said: no guessing)
if any(ismissing(T.Group))
    missingKeys = unique(T.mouse_key(ismissing(T.Group)));
    error("Some mouse_key values are missing from the cohort map. Add them to build_cohort_map():\n%s", ...
        join(missingKeys,newline));
end

T.Group  = string(T.Group);
T.PairID = double(T.PairID);
T.Sex    = string(T.Sex);

%% ---------------- output dir ----------------
outRoot = fullfile(runDir,'figs','pupil_event_locked_REVISED');
if ~exist(outRoot,'dir'), mkdir(outRoot); end

%% ---------------- run both modes: include vs exclude unreliable switch days ----------------
modes = { ...
    struct('Name',"A_includeSwitchDays",  'dropDays',[]), ...
    struct('Name',"B_excludeSwitchDays",  'dropDays',unreliableDays) ...
};

eventTypes = ["reward","lick"];

for mi=1:numel(modes)
    modeName = modes{mi}.Name;
    dropDays = modes{mi}.dropDays;

    outDir = fullfile(outRoot, modeName);
    if ~exist(outDir,'dir'), mkdir(outDir); end

    fprintf('\n=============================\n');
    fprintf('MODE: %s\n', modeName);
    if isempty(dropDays)
        fprintf('  dropDays = (none)\n');
    else
        fprintf('  dropDays = [%s]\n', num2str(dropDays));
    end
    fprintf('Output: %s\n', outDir);
    fprintf('=============================\n');

    % apply day drop
    Tuse = T;
    if ~isempty(dropDays)
        Tuse = Tuse(~ismember(Tuse.day_index, dropDays), :);
    end

    for et = eventTypes
        whichType = et;
        fprintf('\nComputing %s-locked pupil...\n', whichType);

        [E, meta, tGrid, counts] = extract_event_locked_by_mouse(Tuse, whichType, preWin_s, postWin_s, bin_s, useZscore);

        % Always save count diagnostics (THIS answers "why N differs")
        save_counts_report(counts, outDir, whichType);

        % Plot pooled across mice per epoch (mean ± SEM across mice)
        plot_epoch_group_mean_sem(E, meta, tGrid, epochs, whichType, outDir, useZscore);

        % Plot per PairID per epoch (mean ± SEM across events within each mouse, then mean ± SEM across mice)
        plot_pairs_epoch(E, meta, tGrid, epochs, whichType, outDir, useZscore);
    end
end

fprintf('\nDone. Output under:\n  %s\n', outRoot);

end % main


%% =====================================================================
function M = build_cohort_map()
% Explicit cohort definition (NO inference).
% mouse_key MUST match your CSV's mouse_key.
%
% Convention assumed here: "CAGE_COLOR" (e.g., "6100_red").
% If your CSV differs, edit make_mouse_key() below.

rows = {};

% 6100: orange+red passive, black active
rows(end+1,:) = {"6100_red",    "6100","red",    "f","Passive", 1};
rows(end+1,:) = {"6100_orange", "6100","orange", "f","Passive", 1};
rows(end+1,:) = {"6100_black",  "6100","black",  "f","Active",  1};

% 0911: red active vs orange passive
rows(end+1,:) = {"0911_red",    "0911","red",    "f","Active",  2};
rows(end+1,:) = {"0911_orange", "0911","orange", "f","Passive", 2};

% 0911: white active vs black passive
rows(end+1,:) = {"0911_white",  "0911","white",  "f","Active",  3};
rows(end+1,:) = {"0911_black",  "0911","black",  "f","Passive", 3};

% 0910: black active vs orange+red passive
rows(end+1,:) = {"0910_black",  "0910","black",  "m","Active",  4};
rows(end+1,:) = {"0910_orange", "0910","orange", "m","Passive", 4};
rows(end+1,:) = {"0910_red",    "0910","red",    "m","Passive", 4};

% 6099: orange active vs red passive
rows(end+1,:) = {"6099_orange", "6099","orange", "m","Active",  5};
rows(end+1,:) = {"6099_red",    "6099","red",    "m","Passive", 5};

% 6099: black active vs white passive (white died day1-13; keep as Passive but later we will just have fewer days)
rows(end+1,:) = {"6099_black",  "6099","black",  "m","Active",  6};
rows(end+1,:) = {"6099_white",  "6099","white",  "m","Passive", 6};

M = cell2table(rows, 'VariableNames', {'mouse_key','Cage','Color','Sex','Group','PairID'});
M.mouse_key = string(M.mouse_key);
M.Cage      = string(M.Cage);
M.Color     = string(M.Color);
M.Sex       = string(M.Sex);
M.Group     = string(M.Group);
M.PairID    = double(M.PairID);
end


%% =====================================================================
function T = ensureStringCol(T, nm)
if ~ismember(nm, T.Properties.VariableNames)
    T.(nm) = repmat("",height(T),1);
else
    if ~isstring(T.(nm)), T.(nm) = string(T.(nm)); end
end
end


%% =====================================================================
function tb = pickTimebase_fast(T)
cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
tb = nan(height(T),1);
for i=1:numel(cands)
    if ismember(cands{i},T.Properties.VariableNames)
        v = double(T.(cands{i}));
        if any(isfinite(v)), tb = v; return; end
    end
end
end


%% =====================================================================
function [E, meta, tGrid, counts] = extract_event_locked_by_mouse(T, whichType, preWin_s, postWin_s, bin_s, useZscore)
% Extract event-locked traces but STORE them at mouse level:
% - For each mouse/day/session: detect TTL rising edges
% - For each event: extract window [-pre,+post] and baseline-normalize (baseline within [-2,0])
% - Return:
%   E(i).mouse_key
%   E(i).Group
%   E(i).PairID
%   E(i).day_index
%   E(i).eventTraces  (nEvents x nT)
%
% counts: diagnostic table for TTL counts per mouse/day/session

switch lower(whichType)
    case 'reward'
        evCol = 'Injector_TTL';
    case 'lick'
        evCol = 'Lick_TTL';
    otherwise
        error('Unknown whichType: %s', whichType);
end

tb_all = pickTimebase_fast(T);
assert(~all(isnan(tb_all)), 'No valid time base found in table.');

tGrid = -preWin_s:bin_s:postWin_s;
nT    = numel(tGrid);

% Baseline window you requested: [-2,0] within the extracted [-2,+2] window
baseStart = -2;
baseEnd   = 0;

% group by mouse/day/session
[g, mk, di, si] = findgroups(string(T.mouse_key), T.day_index, T.session_idx); %#ok<ASGLU>
nG = max(g);

E = struct('mouse_key',{},'Group',{},'PairID',{},'day_index',{},'eventTraces',{});
meta = struct('mouse_key',{},'Group',{},'PairID',{},'day_index',{});

% counts table (diagnostics)
counts = table('Size',[0 8], ...
    'VariableTypes',{'string','string','double','double','double','double','double','double'}, ...
    'VariableNames',{'mouse_key','Group','PairID','day_index','session_idx','nRows','nTTL_onsets','nFinitePupil'});

for k=1:nG
    idx = (g==k);

    mk_str = string(mk(k));
    day_k  = double(di(k));
    sess_k = double(si(k));

    tb    = tb_all(idx);
    diam  = double(T.Diameter_px(idx));
    ttl   = double(T.(evCol)(idx)) > 0.5;

    grp   = string(T.Group(idx));   grp = grp(~ismissing(grp));
    pid   = double(T.PairID(idx));  pid = pid(~isnan(pid));

    assert(~isempty(grp) && ~isempty(pid), 'Missing Group/PairID after mapping (should not happen).');
    grp = grp(1);
    pid = pid(1);

    good = isfinite(tb) & isfinite(diam);
    nFinite = nnz(good);

    counts(end+1,:) = {mk_str, grp, pid, day_k, sess_k, nnz(idx), nnz(find(diff([false; ttl(:)])==1)), nFinite}; %#ok<AGROW>

    tb   = tb(good);
    diam = diam(good);
    ttl  = ttl(good);

    if numel(tb) < 10, continue; end

    % enforce monotonic time
    dt = diff(tb);
    md = median(dt(isfinite(dt)));
    if ~isempty(md) && md>0
        tb(2:end) = max(tb(2:end), tb(1:end-1)+0.5*md);
    end

    dTTL = diff([false; ttl(:)]);
    onIdx = find(dTTL==1);
    if isempty(onIdx), continue; end

    traces = nan(0,nT);

    for j=1:numel(onIdx)
        t0  = tb(onIdx(j));
        rel = tb - t0;

        inWin = rel>=-preWin_s & rel<=postWin_s;
        if nnz(inWin) < 3, continue; end

        thisT = rel(inWin);
        thisY = diam(inWin);

        % baseline in [-2,0]
        baseMask = thisT>=baseStart & thisT<=baseEnd;

        if nnz(baseMask) >= 3
            mu = mean(thisY(baseMask),'omitnan');
            sd = std(thisY(baseMask),0,'omitnan');
        else
            % fallback to any pre-event samples
            preMask = thisT < 0;
            mu = mean(thisY(preMask),'omitnan');
            sd = std(thisY(preMask),0,'omitnan');
        end

        if ~isfinite(mu), continue; end

        if useZscore
            if ~isfinite(sd) || sd==0
                % if sd invalid, fallback to baseline subtract (still no discard)
                thisY = thisY - mu;
            else
                thisY = (thisY - mu) ./ sd;
            end
        else
            thisY = thisY - mu;
        end

        yInterp = interp1(thisT, thisY, tGrid, 'linear', NaN);
        if all(isnan(yInterp)), continue; end

        traces(end+1,:) = yInterp; %#ok<AGROW>
    end

    if isempty(traces), continue; end

    E(end+1).mouse_key   = mk_str; %#ok<AGROW>
    E(end).Group         = grp;
    E(end).PairID        = pid;
    E(end).day_index     = day_k;
    E(end).eventTraces   = traces;

    meta(end+1).mouse_key = mk_str; %#ok<AGROW>
    meta(end).Group       = grp;
    meta(end).PairID      = pid;
    meta(end).day_index   = day_k;
end
end


%% =====================================================================
function save_counts_report(counts, outDir, whichType)
% Save diagnostics to understand WHY counts differ between active/passive.
% This directly answers your "why N differs" question.

fn = fullfile(outDir, sprintf('%s_event_count_diagnostics.csv', lower(whichType)));
writetable(counts, fn);

% also print a quick summary
fprintf('Saved %s event count diagnostics:\n  %s\n', whichType, fn);

% helpful console warning: show large mismatches by PairID/day
try
    U = counts(:,{'PairID','day_index','Group','nTTL_onsets'});
    % pivot-ish
    pairs = unique(U.PairID);
    for p=pairs(:)'
        Up = U(U.PairID==p,:);
        days = unique(Up.day_index);
        for d=days(:)'
            Ud = Up(Up.day_index==d,:);
            nA = sum(Ud.nTTL_onsets(Ud.Group=="Active"));
            nP = sum(Ud.nTTL_onsets(Ud.Group=="Passive"));
            if nA==0 && nP==0, continue; end
            if nA~=nP
                fprintf('COUNT MISMATCH: %s %s PairID=%d day=%d  ActiveTTL=%d  PassiveTTL=%d\n', ...
                    lower(whichType), 'onsets', p, d, nA, nP);
            end
        end
    end
catch
end
end


%% =====================================================================
function plot_epoch_group_mean_sem(E, meta, tGrid, epochs, whichType, outDir, useZscore)
% mean ± SEM across MICE (not events) for Active and Passive within each epoch
% For each mouse: average its event traces -> one mean trace per mouse
% Then across mice: mean ± SEM

colActive  = [0.85 0.30 0.30];
colPassive = [0.20 0.60 0.90];

for e=1:numel(epochs)
    daysE = epochs(e).Days;

    % build per-mouse mean traces for this epoch
    mice = struct('Group',{},'mouse_key',{},'trace',{});

    for i=1:numel(E)
        if ~ismember(E(i).day_index, daysE), continue; end
        tr = mean(E(i).eventTraces,1,'omitnan'); % average across events within mouse/day/session group
        if all(~isfinite(tr)), continue; end

        mice(end+1).Group = string(E(i).Group); %#ok<AGROW>
        mice(end).mouse_key = string(E(i).mouse_key);
        mice(end).trace = tr;
    end

    if isempty(mice), continue; end

    fh = figure('Color','w','Position',[120 120 1100 450]); hold on;

    hA = []; hP = [];
    nA = 0; nP = 0;

    % Active
    A = vertcat_if_any(mice, "Active");
    if ~isempty(A)
        mA = mean(A,1,'omitnan');
        sA = std(A,0,1,'omitnan') ./ sqrt(size(A,1));
        patch([tGrid fliplr(tGrid)], [mA+sA fliplr(mA-sA)], colActive, 'FaceAlpha',0.18,'EdgeColor','none');
        hA = plot(tGrid, mA, 'Color', colActive, 'LineWidth',2.5);
        nA = size(A,1);
    end

    % Passive
    P = vertcat_if_any(mice, "Passive");
    if ~isempty(P)
        mP = mean(P,1,'omitnan');
        sP = std(P,0,1,'omitnan') ./ sqrt(size(P,1));
        patch([tGrid fliplr(tGrid)], [mP+sP fliplr(mP-sP)], colPassive, 'FaceAlpha',0.18,'EdgeColor','none');
        hP = plot(tGrid, mP, 'Color', colPassive, 'LineWidth',2.5);
        nP = size(P,1);
    end

    xline(0,'k-','LineWidth',1);
    yline(0,'k:');

    xlabel('Time from event (s)');
    if useZscore
        ylabel('Pupil (z-score; baseline -2..0s)');
    else
        ylabel('\Delta pupil (baseline-subtracted; -2..0s)');
    end

    ttl = sprintf('%s-locked pupil | %s (mean \\pm SEM across mice)', cap1(whichType), epochs(e).Name);
    title(ttl,'FontWeight','bold');

    % Legend: explicit handle-based, unambiguous color mapping
    legH = [];
    legT = {};
    if ~isempty(hA), legH(end+1)=hA; legT{end+1}=sprintf('Active (N mice=%d)', nA); end %#ok<AGROW>
    if ~isempty(hP), legH(end+1)=hP; legT{end+1}=sprintf('Passive (N mice=%d)', nP); end %#ok<AGROW>
    if ~isempty(legH)
        legend(legH, legT, 'Location','best');
    end

    grid on; box off;

    fn = fullfile(outDir, sprintf('%s_locked_%s_meanSEM_acrossMice.png', lower(whichType), lower(epochs(e).Name)));
    exportgraphics(fh, fn, 'Resolution',300);
    close(fh);
end
end


%% =====================================================================
function plot_pairs_epoch(E, meta, tGrid, epochs, whichType, outDir, useZscore)
% For each epoch, make a multi-panel figure showing PairID-specific A vs P.
% We do NOT attempt to match events. We compute mean ± SEM across mice
% within that pair/group, using per-mouse event-mean traces.

colActive  = [0.85 0.30 0.30];
colPassive = [0.20 0.60 0.90];

pairIDs = unique([E.PairID]);
pairIDs = pairIDs(:)';

for e=1:numel(epochs)
    daysE = epochs(e).Days;

    fh = figure('Color','w','Position',[50 50 1400 900]);
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    for pi=1:numel(pairIDs)
        pid = pairIDs(pi);

        nexttile; hold on;

        % collect per-mouse mean traces for this pair & epoch
        A_tr = [];
        P_tr = [];

        rawCountA = 0; rawCountP = 0;

        for i=1:numel(E)
            if E(i).PairID ~= pid, continue; end
            if ~ismember(E(i).day_index, daysE), continue; end

            trMouse = mean(E(i).eventTraces,1,'omitnan');  % per mouse/day/session mean over events
            if all(~isfinite(trMouse)), continue; end

            if E(i).Group == "Active"
                A_tr(end+1,:) = trMouse; %#ok<AGROW>
                rawCountA = rawCountA + size(E(i).eventTraces,1);
            else
                P_tr(end+1,:) = trMouse; %#ok<AGROW>
                rawCountP = rawCountP + size(E(i).eventTraces,1);
            end
        end

        % plot A
        hA = [];
        if ~isempty(A_tr)
            mA = mean(A_tr,1,'omitnan');
            sA = std(A_tr,0,1,'omitnan') ./ sqrt(size(A_tr,1));
            patch([tGrid fliplr(tGrid)], [mA+sA fliplr(mA-sA)], colActive, 'FaceAlpha',0.18,'EdgeColor','none');
            hA = plot(tGrid, mA, 'Color', colActive, 'LineWidth',2.5);
        end

        % plot P
        hP = [];
        if ~isempty(P_tr)
            mP = mean(P_tr,1,'omitnan');
            sP = std(P_tr,0,1,'omitnan') ./ sqrt(size(P_tr,1));
            patch([tGrid fliplr(tGrid)], [mP+sP fliplr(mP-sP)], colPassive, 'FaceAlpha',0.18,'EdgeColor','none');
            hP = plot(tGrid, mP, 'Color', colPassive, 'LineWidth',2.5);
        end

        xline(0,'k-','LineWidth',1);
        yline(0,'k:');
        grid on; box off;

        if isempty(A_tr) || isempty(P_tr)
            title(sprintf('%s | %s | PairID %d (missing A/P)', cap1(whichType), epochs(e).Name, pid), 'FontWeight','bold');
        else
            title(sprintf('%s | %s | PairID %d', cap1(whichType), epochs(e).Name, pid), 'FontWeight','bold');
        end

        % add diagnostic text (NO discard)
        txt = sprintf('mouse means: A=%d, P=%d\nraw TTL onsets: A=%d, P=%d', ...
            size(A_tr,1), size(P_tr,1), rawCountA, rawCountP);
        text(0.02,0.05,txt,'Units','normalized','FontSize',9);

        % explicit legend with handles
        legH = [];
        legT = {};
        if ~isempty(hA), legH(end+1)=hA; legT{end+1}=sprintf('Active (N mice=%d)', size(A_tr,1)); end %#ok<AGROW>
        if ~isempty(hP), legH(end+1)=hP; legT{end+1}=sprintf('Passive (N mice=%d)', size(P_tr,1)); end %#ok<AGROW>
        if ~isempty(legH)
            legend(legH, legT, 'Location','best');
        end

        xlabel('Time (s)');
        if useZscore
            ylabel('Pupil z');
        else
            ylabel('\Delta pupil');
        end
    end

    fn = fullfile(outDir, sprintf('%s_locked_%s_pairs_meanSEM.png', lower(whichType), lower(epochs(e).Name)));
    exportgraphics(fh, fn, 'Resolution',300);
    close(fh);
end
end


%% =====================================================================
function X = vertcat_if_any(mice, groupName)
% Collect traces into matrix for a group
tr = [];
for i=1:numel(mice)
    if mice(i).Group ~= groupName, continue; end
    tr(end+1,:) = mice(i).trace; %#ok<AGROW>
end
X = tr;
end


%% =====================================================================
function s = cap1(x)
x = char(x);
s = [upper(x(1)) x(2:end)];
end
