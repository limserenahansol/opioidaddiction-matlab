function analyze_pupil_event_locked()
% Analyze pupil dynamics aligned to reward / lick events.
% - Uses latest run_* under I:\addiction_concate_Nov\longitudinal_outputs
% - Needs S_D_cache_basic.mat from analyze_PR_pupil_pairs_raster (for Active/Passive)
% - Outputs figures under runDir/figs/pupil_event_locked
%
% Windows (can edit):
preWin_s   = 2;     % seconds before event
postWin_s  = 2;     % seconds after event
bin_s      = 0.2;   % time bin for interpolation (s)
earlyDays  = 3:5;   % early epoch
lateDays   = 6:10;  % late epoch

%% --- locate latest run + load CSV ---
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

% ensure session_idx exists (else treat as single session per day)
if ~ismember('session_idx',T.Properties.VariableNames)
    T.session_idx = ones(height(T),1);
end

% TTL and pupil columns
needCols = {'Diameter_px','Lick_TTL','Injector_TTL'};
for i = 1:numel(needCols)
    if ~ismember(needCols{i},T.Properties.VariableNames)
        error('Missing column %s in table.', needCols{i});
    end
end
T.Lick_TTL(isnan(T.Lick_TTL))         = 0;
T.Injector_TTL(isnan(T.Injector_TTL)) = 0;

%% --- load S/D cache and classify Active vs Passive mice ---
cacheMat = fullfile(runDir, 'S_D_cache_basic.mat');
assert(exist(cacheMat,'file')>0, ...
    'S_D_cache_basic.mat not found in %s. Run analyze_PR_pupil_pairs_raster first.', runDir);

L = load(cacheMat,'S','D');
S = L.S;
D = L.D;

D.HadPassive = classifyHadPassive(S, D);
D.Group = categorical(D.HadPassive, [false true], {'Active','Passive'});

% map mouse_key -> Group
map = unique(D(:,{'mouse_key','Group'}));
map.mouse_key = string(map.mouse_key);
map.Group     = string(map.Group);

% attach Group to T by mouse_key
T.mouse_key = string(T.mouse_key);
T = outerjoin(T, map, 'Keys','mouse_key', ...
    'MergeKeys', true, 'Type','left');

% any missing Group -> treat as Active
if ~ismember('Group',T.Properties.VariableNames)
    T.Group = repmat("Active",height(T),1);
else
    G = string(T.Group);
    G(ismissing(G)) = "Active";
    T.Group = G;
end

%% --- output dir ---
outDir = fullfile(runDir,'figs','pupil_event_locked');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% --- reward-locked pupil ---
fprintf('Computing reward-locked pupil dynamics...\n');
[E_rew, meta_rew, tGrid] = extract_pupil_event_locked(T, ...
    'reward', preWin_s, postWin_s, bin_s);

plot_pupil_event_locked(E_rew, meta_rew, tGrid, ...
    earlyDays, lateDays, 'Reward', outDir);

%% --- lick-locked pupil ---
fprintf('Computing lick-locked pupil dynamics...\n');
[E_lick, meta_lick, tGrid2] = extract_pupil_event_locked(T, ...
    'lick', preWin_s, postWin_s, bin_s);

% (time grid should be identical, but just in case:)
if ~isequal(tGrid, tGrid2)
    warning('Reward-locked and lick-locked time grids differ; using lick grid for labels.');
end

plot_pupil_event_locked(E_lick, meta_lick, tGrid2, ...
    earlyDays, lateDays, 'Lick', outDir);

fprintf('Done. Event-locked pupil figures are under:\n  %s\n', outDir);
end  % ===== end main function =====


%% ================= helper functions (local) =================

function T = ensureStringCol(T, nm)
    if ~ismember(nm, T.Properties.VariableNames)
        T.(nm) = repmat("",height(T),1);
    else
        if ~isstring(T.(nm))
            T.(nm) = string(T.(nm));
        end
    end
end

%% --- Active/Passive classifier (same logic as other script) ---
function HadPassive = classifyHadPassive(S, D)
    mice = unique(S.mouse_key,'stable');
    % mouse_key in S is categorical; convert to string for comparison
    miceStr = string(removecats(mice));
    DmkStr  = string(removecats(D.mouse_key));

    HadPassive = false(height(D),1);
    havePar = ismember('Session_Paradigm', S.Properties.VariableNames);

    for i = 1:numel(miceStr)
        rS  = string(S.mouse_key)==miceStr(i);
        ip  = S.isPassive(rS);
        di  = S.day_index(rS);

        par = strings(nnz(rS),1);
        if havePar
            p = S.Session_Paradigm(rS);
            if iscategorical(p), p = string(p); end
            par = string(p);
        end

        flag = false;
        if any(ip==1,'all')
            flag = true;
        elseif havePar && any(contains(lower(par),'passive'),'all')
            flag = true;
        elseif any(ismember(double(di),6:10) & (ip>=1 | isnan(ip)),'all')
            flag = true;
        end

        Dmask = DmkStr==miceStr(i);
        HadPassive(Dmask) = flag;
    end
end

%% --- timebase picker ---
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

%% --- core extractor: event-locked pupil traces ---
function [E, meta, tGrid] = extract_pupil_event_locked(T, whichType, preWin_s, postWin_s, bin_s)
% whichType: 'reward' or 'lick'
%
% Outputs:
%   E    : nEvents x nTime matrix of baseline-subtracted pupil
%   meta : struct array with .Group, .day_index, .mouse_key
%   tGrid: time vector (s)

    if ~ismember('Diameter_px',T.Properties.VariableNames)
        error('Diameter_px column not found.');
    end
    switch lower(whichType)
        case 'reward'
            evCol = 'Injector_TTL';
        case 'lick'
            evCol = 'Lick_TTL';
        otherwise
            error('Unknown whichType: %s', whichType);
    end
    if ~ismember(evCol,T.Properties.VariableNames)
        error('Event column %s not found.', evCol);
    end

    tb_all = pickTimebase_fast(T);
    if all(isnan(tb_all))
        error('No valid time base found in table.');
    end

    % group by mouse/day/session
    [g, mk, di, si] = findgroups(string(T.mouse_key), T.day_index, T.session_idx); %#ok<ASGLU>
    nG = max(g);

    tGrid = -preWin_s:bin_s:postWin_s;
    nT = numel(tGrid);
    E = nan(0,nT);
    meta = struct('Group',{}, 'day_index',{}, 'mouse_key',{});

    % Fixed baseline window: [-2, 0] s
    baseStart = -2;
    baseEnd   = 0;

    for k = 1:nG
        idx = (g==k);

        tb    = tb_all(idx);
        diam  = double(T.Diameter_px(idx));
        ttl   = double(T.(evCol)(idx));
        group = string(T.Group(idx));
        group = group(~ismissing(group));
        if isempty(group), groupStr = "Active"; else, groupStr = group(1); end

        mk_str = string(mk(k));
        day_k  = di(k);

        good = isfinite(tb) & isfinite(diam) & ~isnan(ttl);
        tb   = tb(good);
        diam = diam(good);
        ttl  = ttl(good) > 0.5;

        if numel(tb) < 10
            continue;
        end

        % enforce monotonic time
        dt = diff(tb);
        md = median(dt(isfinite(dt)));
        if ~isempty(md) && md>0
            tb(2:end) = max(tb(2:end), tb(1:end-1)+0.5*md);
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

            % baseline using [-2, 0] s
            baseMask = thisT >= baseStart & thisT <= baseEnd;

            if nnz(baseMask) >= 3
                baseVal = nanmean(thisY(baseMask));
            else
                % Fallback: use whatever we have in the pre-event window
                preMask = thisT < 0;
                if nnz(preMask) >= 3
                    baseVal = nanmean(thisY(preMask));
                else
                    baseVal = nanmean(thisY);  % last resort
                end
            end

            thisY = thisY - baseVal;

            % interpolate onto uniform grid
            yInterp = interp1(thisT, thisY, tGrid, 'linear', NaN);
            if all(isnan(yInterp))
                continue;
            end

            E(end+1,:) = yInterp; %#ok<AGROW>
            meta(end+1).Group      = groupStr; %#ok<AGROW>
            meta(end).day_index    = day_k;
            meta(end).mouse_key    = mk_str;
        end
    end
end

%% --- plotting function: combined + separate Active/Passive ---
%% --- plotting function: combined + separate Active/Passive ---
function plot_pupil_event_locked(E, meta, tGrid, earlyDays, lateDays, labelStr, outDir)
% Plot event-locked pupil traces:
%   - rows of E = individual events
%   - meta(k).Group = 'Active' / 'Passive'
%   - meta(k).day_index = day index
%   - earlyDays, lateDays = vectors of day indices
%   - labelStr = 'Reward' or 'Lick' (for title / filename)

    if isempty(E)
        warning('plot_pupil_event_locked:Empty', ...
            'No event-locked pupil data for %s; skipping plot.', labelStr);
        return;
    end

    nEv = size(E,1);

    % Extract metadata into arrays
    G   = strings(nEv,1);
    day = nan(nEv,1);
    for i = 1:nEv
        if isfield(meta(i),'Group')
            G(i) = string(meta(i).Group);
        else
            G(i) = "Active";
        end
        if isfield(meta(i),'day_index')
            day(i) = double(meta(i).day_index);
        else
            day(i) = NaN;
        end
    end

    groups = {'Active','Passive'};
    epochs = {'Early','Late'};

    % Colors similar to your bar plots
    colActive  = [0.85 0.30 0.30];
    colPassive = [0.20 0.60 0.90];

    % Line styles: solid = early, dashed = late
    lsEpoch = {'-','--'};   % lsEpoch{1} = Early, lsEpoch{2} = Late

    %% ========== Combined plot: Active + Passive ==========

    figure('Position',[200 150 900 650],'Color','w'); hold on;

    legendEntriesCombined = {};

    for gi = 1:numel(groups)
        for ei = 1:numel(epochs)
            gName = groups{gi};
            eName = epochs{ei};

            % epoch mask
            switch eName
                case 'Early'
                    mEpoch = ismember(day, earlyDays);
                    daysVec = earlyDays;
                case 'Late'
                    mEpoch = ismember(day, lateDays);
                    daysVec = lateDays;
                otherwise
                    mEpoch = true(size(day));
                    daysVec = unique(day(~isnan(day)));
            end

            mGroup = strcmpi(G, gName);
            mask   = mGroup & mEpoch;

            if ~any(mask)
                continue;
            end

            D = E(mask,:);              % events x time
            mTrace = mean(D,1,'omitnan');
            sTrace = std(D,0,1,'omitnan') ./ sqrt(size(D,1));

            % pick color
            if strcmpi(gName,'Active')
                c = colActive;
            else
                c = colPassive;
            end

            % plot mean ± SEM (shaded)
            x = tGrid(:)';
            mTrace = mTrace(:)';
            sTrace = sTrace(:)';

            % shaded error
            xx = [x fliplr(x)];
            yy = [mTrace+sTrace fliplr(mTrace-sTrace)];
            hPatch = patch(xx, yy, c, ...
                'FaceAlpha',0.15, 'EdgeColor','none');
            uistack(hPatch,'bottom');

            % main line
            plot(x, mTrace, 'Color', c, ...
                'LineWidth',2, 'LineStyle', lsEpoch{ei});

            legendEntriesCombined{end+1} = sprintf('%s %s (days %s)', ...
                gName, eName, vec2str(daysVec)); %#ok<AGROW>
        end
    end

    yline(0,'k:');
    xline(0,'k-','LineWidth',1);

    xlabel('Time from event (s)');
    ylabel('\Delta pupil diameter (a.u., baseline-subtracted)');
    ttl = sprintf('%s-locked pupil (event-aligned, baseline-corrected)', labelStr);
    title(ttl, 'FontWeight','bold');

    if ~isempty(legendEntriesCombined)
        legend(legendEntriesCombined,'Location','best','Interpreter','none');
    end

    grid on; box off;

    % Save combined
    if nargin >= 7 && ~isempty(outDir)
        if ~exist(outDir,'dir'), mkdir(outDir); end
        fn = sprintf('Pupil_%s_locked_combined.png', lower(labelStr));
        outPNG = fullfile(outDir, fn);
        exportgraphics(gcf, outPNG, 'Resolution',300);
        fprintf('Saved combined %s-locked pupil plot to:\n  %s\n', labelStr, outPNG);
    end

    %% ========== Separate plots: Active-only and Passive-only ==========

    for gi = 1:numel(groups)
        gName = groups{gi};
        groupMask = strcmpi(G, gName);
        if ~any(groupMask)
            continue;
        end

        % color for this group
        if strcmpi(gName,'Active')
            c = colActive;
        else
            c = colPassive;
        end

        figure('Position',[200 150 900 650],'Color','w'); hold on;

        legendEntriesGroup = {};

        for ei = 1:numel(epochs)
            eName = epochs{ei};

            switch eName
                case 'Early'
                    mEpoch = ismember(day, earlyDays);
                    daysVec = earlyDays;
                case 'Late'
                    mEpoch = ismember(day, lateDays);
                    daysVec = lateDays;
                otherwise
                    mEpoch = true(size(day));
                    daysVec = unique(day(~isnan(day)));
            end

            mask = groupMask & mEpoch;
            if ~any(mask)
                continue;
            end

            D = E(mask,:);  % events x time
            mTrace = mean(D,1,'omitnan');
            sTrace = std(D,0,1,'omitnan') ./ sqrt(size(D,1));

            x = tGrid(:)';
            mTrace = mTrace(:)';
            sTrace = sTrace(:)';

            % shaded error
            xx = [x fliplr(x)];
            yy = [mTrace+sTrace fliplr(mTrace-sTrace)];
            hPatch = patch(xx, yy, c, ...
                'FaceAlpha',0.15, 'EdgeColor','none');
            uistack(hPatch,'bottom');

            % main line, style encodes epoch
            plot(x, mTrace, 'Color', c, ...
                 'LineWidth',2, 'LineStyle', lsEpoch{ei});

            legendEntriesGroup{end+1} = sprintf('%s %s (days %s)', ...
                gName, eName, vec2str(daysVec)); %#ok<AGROW>
        end

        yline(0,'k:');
        xline(0,'k-','LineWidth',1);

        xlabel('Time from event (s)');
        ylabel('\Delta pupil diameter (a.u., baseline-subtracted)');
        ttl2 = sprintf('%s-locked pupil (%s only)', labelStr, gName);
        title(ttl2, 'FontWeight','bold');

        if ~isempty(legendEntriesGroup)
            legend(legendEntriesGroup,'Location','best','Interpreter','none');
        end
        grid on; box off;

        % Save per-group
        if nargin >= 7 && ~isempty(outDir)
            if ~exist(outDir,'dir'), mkdir(outDir); end
            fn = sprintf('Pupil_%s_locked_%s.png', lower(labelStr), gName);
            outPNG = fullfile(outDir, fn);
            exportgraphics(gcf, outPNG, 'Resolution',300);
            fprintf('Saved %s-only %s-locked pupil plot to:\n  %s\n', ...
                gName, labelStr, outPNG);
        end
    end
end

function s = vec2str(v)
    v = v(:)';          % row
    s = sprintf('%d-', v);
    s = s(1:end-1);     % remove trailing dash
end


