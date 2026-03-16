function run_lick_MEGA_followup2()
% FOLLOW-UP 2: visualize MEGA CSVs with short labels (mouse ID only)

%% locate latest run
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up=1:5
        p = fullfile(cand,'longitudinal_outputs'); if exist(p,'dir'), rootTry=p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*')); assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]); runDir = fullfile(D(ix).folder, D(ix).name);
outDir = fullfile(runDir,'figs','lick_MEGA_followup2'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% load MEGA outputs
S_sess  = safeRead(fullfile(runDir,'figs','lick_MEGA','per_session_features.csv'));
S_trial = safeRead(fullfile(runDir,'figs','lick_MEGA','per_trial_features.csv'));
S_spec  = safeRead(fullfile(runDir,'figs','lick_MEGA','C_rhythm_spectral_metrics.csv'));
S_ac    = safeRead(fullfile(runDir,'figs','lick_MEGA','C_rhythm_autocorr_metrics.csv'));
Pstart  = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','psth_start_group.csv'));
Prew    = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','psth_reward_group.csv'));
Drpr    = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','reward_ramp_pause.csv'));

%% sanitize key columns
S_sess  = ensureCols(S_sess);
S_trial = ensureCols(S_trial);
S_spec  = ensureCols(S_spec);
S_ac    = ensureCols(S_ac);

%% session rates (mouse labels only)
boxAndJitterByGroup(S_sess,'licks_per_min','Licks / min', outDir,'A_session_licks_per_min.png');
boxAndJitterByBin  (S_sess,'licks_per_min','Licks / min', outDir,'A_session_licks_per_min_by_bin.png');
spaghettiPerMouse  (S_sess,'licks_per_min', outDir,'S_session_spaghetti_licks_per_min.png');

%% rhythm summaries
if ~isempty(S_spec)
    quickViolin(S_spec,'peak_hz','Peak Hz',outDir,'C_spec_peak_hz.png');
    quickViolin(S_spec,'band_power_3_12','Band power 3–12 Hz',outDir,'C_spec_band312.png');
    quickViolin(S_spec,'spectral_entropy','Spectral entropy',outDir,'C_spec_entropy.png');
end
if ~isempty(S_ac)
    quickViolin(S_ac,'ac_peak_lag','AC peak lag (s)',outDir,'C_ac_lag.png');
    quickViolin(S_ac,'ac_peak_height','AC peak height',outDir,'C_ac_height.png');
end

%% PSTH (group) — robust to “wide” CSVs, no reliance on a 't' column
%% PSTH (group) — start-only, split panels + difference
if ~isempty(Pstart)
    plotPSTHtableSplitDiff(Pstart, ...
        'Start-aligned lick rate (0–5 s)', ...
        fullfile(outDir,'PSTH_start_group_SPLIT_and_DIFFERENCE.png'));
end
% (skip reward-aligned per your request)
if ~isempty(Drpr), writetable(Drpr, fullfile(outDir,'reward_ramp_pause_copy.csv')); end

fprintf('\nFOLLOWUP2 complete. Outputs in:\n  %s\n', outDir);
end

%% ================= helpers ============================
function plotPSTHtableSplitDiff(P, ttl, saveTo)
    if isempty(P), return; end

    % groups
    if ismember('Group',P.Properties.VariableNames)
        G = categorical(string(P.Group));
        G = setcats(G, {'ActiveOnly','HadPassive'});
        groups = categories(G);
        groups = groups(~cellfun(@isempty,groups));
    else
        G = categorical(repmat("All",height(P),1));
        groups = {"All"};
    end
    assert(numel(groups)>=1,'No groups in PSTH table.');

    % pull vectors for each group (wide or vector-in-cell)
    T = cell(1,numel(groups)); MU = T; SE = T;
    for i = 1:numel(groups)
        r = (G == groups{i});
        t  = rowFromWide(P,r,'t_');
        mu = rowFromWide(P,r,'mean_rate_');
        se = rowFromWide(P,r,'sem_rate_');
        if isempty(t) && ismember('t',P.Properties.VariableNames)
            ri = find(r,1,'first');
            t  = parseVec(P{ri,'t'});
            mu = parseVec(P{ri,'mean_rate'});
            if ismember('sem_rate',P.Properties.VariableNames)
                se = parseVec(P{ri,'sem_rate'});
            else
                se = zeros(size(mu));
            end
        end
        T{i}=t(:); MU{i}=mu(:);
        if isempty(se), se=zeros(size(mu)); end
        SE{i}=se(:);
    end

    C = [0 0.45 0.74; 0.85 0.33 0.10; 0.3 0.3 0.3]; % colors
    ymax = 1.05*max(cellfun(@(m,s) max(m+s), MU, SE));

    fig = figure('Color','w','Position',[100 100 900 650]);
    tl  = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % panel 1: first group
    nexttile; hold on
    i = 1;
    fill([T{i}; flipud(T{i})], [MU{i}-SE{i}; flipud(MU{i}+SE{i})], C(i,:), ...
        'FaceAlpha',0.18,'EdgeColor','none');
    plot(T{i}, MU{i}, 'LineWidth',2.2, 'Color', C(i,:));
    xline(0,'k:'); ylim([0 ymax]);
    ylabel('Lick rate (Hz)'); title(sprintf('%s — %s', ttl, groups{i}));
    grid on; box on

    % panel 2: second group (if present)
    if numel(groups) >= 2
        j = 2;
        nexttile; hold on
        fill([T{j}; flipud(T{j})], [MU{j}-SE{j}; flipud(MU{j}+SE{j})], C(j,:), ...
            'FaceAlpha',0.18,'EdgeColor','none');
        plot(T{j}, MU{j}, 'LineWidth',2.2, 'Color', C(j,:));
        xline(0,'k:'); ylim([0 ymax]);
        ylabel('Lick rate (Hz)'); title(groups{j});
        grid on; box on

        % panel 3: difference (HadPassive − ActiveOnly)
        % align time vectors if needed
        tC = T{1};
        if numel(T{1})~=numel(T{2}) || any(abs(T{1}-T{2})>1e-12)
            tmin = max(min(T{1}), min(T{2}));
            tmax = min(max(T{1}), max(T{2}));
            tC   = linspace(tmin,tmax,min(numel(T{1}),numel(T{2})))';
            mu1  = interp1(T{1}, MU{1}, tC, 'linear','extrap');
            se1  = interp1(T{1}, SE{1}, tC, 'linear','extrap');
            mu2  = interp1(T{2}, MU{2}, tC, 'linear','extrap');
            se2  = interp1(T{2}, SE{2}, tC, 'linear','extrap');
        else
            mu1 = MU{1}; se1 = SE{1}; mu2 = MU{2}; se2 = SE{2};
        end
        muD = mu2 - mu1; seD = sqrt(se1.^2 + se2.^2);

        nexttile; hold on
        fill([tC; flipud(tC)], [muD-seD; flipud(muD+seD)], C(3,:), ...
            'FaceAlpha',0.15,'EdgeColor','none');
        plot(tC, muD, 'Color', C(3,:), 'LineWidth',2.2);
        yline(0,'k:'); xline(0,'k:');
        xlabel('Time from trial start (s)'); ylabel('\Delta rate (P − A)');
        title('Difference'); grid on; box on
    else
        xlabel('Time from trial start (s)');
    end

    savepng(fig, saveTo); close(fig);
end

function T = safeRead(p)
    if exist(p,'file')>0, T = readtable(p,'VariableNamingRule','preserve'); else, T = table(); end
end

function S = ensureCols(S)
    if isempty(S), return; end
    if ismember('mouse_key',S.Properties.VariableNames) && ~isstring(S.mouse_key)
        S.mouse_key = string(S.mouse_key);
    end
    if ismember('GroupMouse',S.Properties.VariableNames)
        if ~iscategorical(S.GroupMouse), S.GroupMouse = categorical(string(S.GroupMouse)); end
        S.GroupMouse = setcats(S.GroupMouse, {'ActiveOnly','HadPassive'});
    end

    % --- NEW: robust bin handling with 15/16/25 bucket ---
    binsWanted = {'D3-5','D6-8','D9-11','D12-14','D15-16-25','<undef>'};

    if ismember('day_index',S.Properties.VariableNames)
        % Recompute from day_index if available
        S.bin5 = dayBin5_15_16_25(S.day_index);
    elseif ismember('bin5',S.Properties.VariableNames)
        % Keep existing bin5, but rename and set categories
        if ~iscategorical(S.bin5), S.bin5 = categorical(string(S.bin5)); end
        if any(strcmp('D15-16', categories(S.bin5)))
            S.bin5 = renamecats(S.bin5,'D15-16','D15-16-25');
        end
        S.bin5 = setcats(S.bin5, binsWanted);
    else
        % No bin info at all -> fill as <undef>
        S.bin5 = categorical(repmat("<undef>",height(S),1), binsWanted);
    end
end
function B = dayBin5_15_16_25(day_index)
% Map days to bins: 3–5, 6–8, 9–11, 12–14, and {15,16,25}
    lab = strings(size(day_index)); lab(:) = "<undef>";
    d = double(day_index);
    lab(d>=3  & d<=5 )  = "D3-5";
    lab(d>=6  & d<=8 )  = "D6-8";
    lab(d>=9  & d<=11 ) = "D9-11";
    lab(d>=12 & d<=14 ) = "D12-14";
    lab(d==15 | d==16 | d==25) = "D15-16-25";
    B = categorical(lab, ["D3-5","D6-8","D9-11","D12-14","D15-16-25","<undef>"]);
end


function s = mouselabel(mk)
% "7597_black" -> "7597 black"
    s = regexprep(string(mk), '[_\-]+', ' ');
end

function boxAndJitterByGroup(S, col, ylab, outDir, fname)
    if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember('GroupMouse',S.Properties.VariableNames), return, end
    cats = categories(S.GroupMouse); cats=cats(~cellfun(@isempty,cats)); if isempty(cats), return, end
    fig = figure('Color','w','Position',[90 90 740 520]); hold on
    for i=1:numel(cats)
        r = (S.GroupMouse == cats{i});   % <-- direct compare to category name (no categorical(...))
        x = i + 0.35*(rand(sum(r),1)-0.5); y = double(S.(col)(r));
        scatter(x, y, 18, 'filled','MarkerFaceAlpha',0.45);
        if ismember('mouse_key',S.Properties.VariableNames)
            lab = mouselabel(S.mouse_key(r));
            text(x+0.02, y, lab, 'FontSize',6, 'Interpreter','none');
        end
    end
    xlim([0.5 numel(cats)+0.5]); set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
    ylabel(ylab); title([ylab ' by Group']); grid on; box on
    savepng(fig, fullfile(outDir, fname)); close(fig);
end

function boxAndJitterByBin(S, col, ylab, outDir, fname)
    if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember('bin5',S.Properties.VariableNames), return, end
    binsAll = {'D3-5','D6-8','D9-11','D12-14','D15-16-25'};
    bins = binsAll(ismember(binsAll, categories(S.bin5)));   % show only present bins
    fig = figure('Color','w','Position',[90 90 900 520]); hold on
    for i=1:numel(bins)
        r = (S.bin5 == bins{i});
        x = i + 0.35*(rand(sum(r),1)-0.5); y = double(S.(col)(r));
        scatter(x, y, 18, 'filled','MarkerFaceAlpha',0.45);
        if ismember('mouse_key',S.Properties.VariableNames)
            lab = regexprep(string(S.mouse_key(r)),'[_\-]+',' ');
            text(x+0.02, y, lab, 'FontSize',6, 'Interpreter','none');
        end
    end
    xlim([0.5 numel(bins)+0.5]); set(gca,'XTick',1:numel(bins),'XTickLabel',bins);
    ylabel(ylab); title([ylab ' by Day-bin']); grid on; box on
    savepng(fig, fullfile(outDir, fname)); close(fig);
end


function spaghettiPerMouse(S, col, outDir, fname)
    if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember('mouse_key',S.Properties.VariableNames), return, end
    fig = figure('Color','w','Position',[90 90 900 520]); hold on
    ms = unique(S.mouse_key,'stable');
    for i=1:numel(ms)
        r = S.mouse_key==ms(i);
        x = double(S.day_index(r)); y = double(S.(col)(r));
        [x,ord] = sort(x); y=y(ord);
        plot(x,y,'-o','LineWidth',1,'MarkerSize',3);
        if ~isempty(y)
            [~,mxi] = max(y);
            text(x(mxi)+0.05, y(mxi), mouselabel(ms(i)), 'FontSize',7, 'Interpreter','none');
        end
    end
    xlabel('Day'); ylabel(strrep(col,'_','\_')); title(['Spaghetti per mouse: ' strrep(col,'_','\_')]); grid on; box on
    savepng(fig, fullfile(outDir, fname)); close(fig);
end

function quickViolin(T, col, ylab, outDir, fname)
    if isempty(T) || ~ismember(col,T.Properties.VariableNames) || ~ismember('GroupMouse',T.Properties.VariableNames), return, end
    cats = categories(T.GroupMouse); cats = cats(~cellfun(@isempty,cats)); if isempty(cats), return, end
    fig = figure('Color','w','Position',[100 100 700 420]); hold on
    for i=1:numel(cats)
        y = double(T.(col)(T.GroupMouse==cats{i}));
        if isempty(y), continue, end
        [yk,xk] = kde1(y,200); yk = yk/max(yk)*0.35;  % density
        patch([i-yk, fliplr(i+yk)], [xk, fliplr(xk)], 0.92*[1 1 1], 'EdgeColor','none'); % violin
        scatter(i + 0.25*(rand(numel(y),1)-0.5), y, 8, 'filled','MarkerFaceAlpha',0.35);
        plot([i-0.18 i+0.18], [median(y) median(y)], 'k-', 'LineWidth',1.8);
    end
    set(gca,'XTick',1:numel(cats),'XTickLabel',cats); ylabel(ylab); title(ylab); grid on; box on
    savepng(fig, fullfile(outDir,fname)); close(fig);
end

function plotPSTHtableWide(P, ttl, saveTo)
% Robust to wide CSVs. Never assumes a variable named 't'.
    if isempty(P), return, end
    if ismember('Group',P.Properties.VariableNames)
        G = categorical(string(P.Group)); G = setcats(G, {'ActiveOnly','HadPassive'});
        groups = categories(G); groups = groups(~cellfun(@isempty,groups));
    else
        G = categorical(repmat("All",height(P),1)); groups = {"All"};
    end

    fig = figure('Color','w','Position',[100 100 820 480]); hold on
    cols = lines(max(3,numel(groups)));
    for i=1:numel(groups)
        r = (G == groups{i});
        if ~any(r), continue, end
        t  = rowFromWide(P, r, 't_');
        mu = rowFromWide(P, r, 'mean_rate_');
        se = rowFromWide(P, r, 'sem_rate_');
        if isempty(t) && ismember('t',P.Properties.VariableNames)
            ri = find(r,1,'first');
            t  = parseVec(P{ri,'t'});  mu = parseVec(P{ri,'mean_rate'});  se = parseVec(P{ri,'sem_rate'});
        end
        if isempty(t) || isempty(mu), continue, end
        if isempty(se), se = zeros(size(mu)); end
        fill([t(:); flipud(t(:))]', [mu(:)-se(:); flipud(mu(:)+se(:))]', cols(i,:), 'FaceAlpha',0.15,'EdgeColor','none');
        plot(t, mu, 'LineWidth',2,'Color',cols(i,:));
    end
    xlabel('Time (s)'); ylabel('Lick rate (Hz)'); title(ttl); grid on; box on
    legend(groups,'Location','best');
    savepng(fig, saveTo); close(fig);
end

function v = parseVec(x)
    if isnumeric(x), v = x(:)'; return, end
    s = char(string(x)); s = strrep(s,'[',''); s = strrep(s,']',''); s = strrep(s,';',',');
    v = str2num(s); v = v(:)'; %#ok<ST2NM>
end

function v = rowFromWide(T, mask, prefix)
    cols = T.Properties.VariableNames;
    k = startsWith(cols, prefix); cols = cols(k);
    if isempty(cols), v = []; return, end
    suf = regexp(cols, ['^' regexptranslate('escape',prefix) '(\d+)$'], 'tokens','once');
    suf = cellfun(@(c) str2double(c{1}), suf); [~,ord] = sort(suf); cols = cols(ord);
    r = find(mask,1,'first'); v = zeros(1,numel(cols));
    for i=1:numel(cols), v(i) = double(T.(cols{i})(r)); end
end

function [y,x] = kde1(v, n)
    v=v(isfinite(v)); if numel(v)<3, x=linspace(min(v)-eps,max(v)+eps,20); y=ones(size(x)); return, end
    s = std(v); h = 1.06*s* numel(v)^(-1/5);
    x = linspace(min(v)-3*h, max(v)+3*h, n); y=zeros(size(x));
    for i=1:numel(v), y = y + exp(-0.5*((x-v(i))/h).^2); end
    y = y/(numel(v)*h*sqrt(2*pi));
end

function savepng(fh, fn)
    set(fh,'PaperPositionMode','auto');
    try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end
