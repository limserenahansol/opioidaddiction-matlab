function analyze_passive_active_dashboard1()
% One-stop cohort visualization + nonparametric stats for addiction_concate longitudinal outputs.

%% ---------- locate latest run + load ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*')); assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix]   = max([d.datenum]);
runDir   = fullfile(d(ix).folder,d(ix).name);
csvPath  = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s',csvPath);
fprintf('Reading: %s\n', csvPath);
T = readtable(csvPath,'VariableNamingRule','preserve');

% Normalise types
T = ensureString(T,'mouse_key');
T = ensureString(T,'day_name');
T = ensureString(T,'Session_Paradigm');
if ~ismember('isPassive',T.Properties.VariableNames)
    T.isPassive = nan(height(T),1);
elseif ~isnumeric(T.isPassive)
    T.isPassive = double(T.isPassive);
end

% TTL numerics (keep raw; threshold later)
if ismember('Lick_TTL',T.Properties.VariableNames),     T.Lick_TTL(isnan(T.Lick_TTL)) = 0; end
if ismember('Injector_TTL',T.Properties.VariableNames), T.Injector_TTL(isnan(T.Injector_TTL)) = 0; end

%% ---------- exclude only 6872 black + 8606 forange ----------
if ismember('mouse_key',T.Properties.VariableNames)
    mk      = string(T.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    drop    = (contains(mk_norm,"6872") & contains(mk_norm,"black")) | ...
              (contains(mk_norm,"8606") & contains(mk_norm,"forange"));
    if any(drop)
        fprintf('Excluding %d rows from: %s\n', nnz(drop), strjoin(unique(cellstr(mk(drop)))', ', '));
        T(drop,:) = [];
    end
end

%% ---------- per-session metrics, then per-day medians ----------
[S, D] = fast_session_day_metrics(T, runDir);   % caching + speedups + carries extras
% --- PUPIL-ONLY exclusion: 7597 black ---
if ismember('mouse_key', D.Properties.VariableNames) && ismember('pupil_mean', D.Properties.VariableNames)
    mk      = string(D.mouse_key);
    mk_norm = lower(regexprep(strtrim(mk),'[_\-]+',' '));
    mask7597 = contains(mk_norm,"7597") & contains(mk_norm,"black");
    if any(mask7597)
        fprintf('Pupil-only exclude: setting pupil_mean=NaN for %d rows from 7597 black.\n', nnz(mask7597));
        D.pupil_mean(mask7597) = NaN;   % affects only pupil; all other metrics stay intact
    end
end
%% ---------- cohort label ----------
D.HadPassive = classifyHadPassive(S, D);

%% ---------- epochs / groups ----------
D.Epoch = epochOfDay(double(D.day_index));
D.Group = categorical(D.HadPassive, [false true], {'Active','Passive'});

%% ---------- metric list (dynamic) ----------
baseMetrics = { ...
    'RequirementLast','lick_freq_per_min','rew_freq_per_min', ...
    'lick_meanDur_s','rew_meanDur_s', ...
    'lick_totalDur_s','rew_totalDur_s', ...
    'lick_medianIEI_s','rew_medianIRI_s','pupil_mean'};
baseLabels  = { ...
    'Requirement','Lick freq (/min)','Reward freq (/min)', ...
    'Lick mean dur (s)','Reward mean dur (s)', ...
    'Lick total dur (s)','Reward total dur (s)', ...
    'Lick median IEI (s)','Reward median IRI (s)','Pupil mean (px)'};

% Optional extras: Immersion; prefer *Nonmoving* for TST/HOT if present.
extraM = {}; extraL = {};
immCol = colLike(D,'Immersion_Latency_s');
if ~isempty(immCol), extraM{end+1}=immCol; extraL{end+1}='Immersion latency (s)'; end
tstCols = D.Properties.VariableNames(contains(D.Properties.VariableNames,'TST_Pct'));
hotCols = D.Properties.VariableNames(contains(D.Properties.VariableNames,'HOT_Pct'));
tstPick = pickPctColumn(tstCols, 'Nonmoving'); if ~isempty(tstPick), extraM{end+1}=tstPick; extraL{end+1}='TST nonmoving (%)'; end
hotPick = pickPctColumn(hotCols, 'Nonmoving'); if ~isempty(hotPick), extraM{end+1}=hotPick; extraL{end+1}='HOT nonmoving (%)'; end

metrics = [baseMetrics, extraM];
labels  = [baseLabels,  extraL];
keepIdx = ismember(metrics, D.Properties.VariableNames);
metrics = metrics(keepIdx); labels = labels(keepIdx);

%% ---------- output dir ----------
outDir = fullfile(runDir,'figs','dashboard');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---------- main plots + stats ----------
allStats = table();
for k = 1:numel(metrics)
    ycol = metrics{k}; ylab = labels{k};
    W = perMouseEpochTable(D, ycol);     % one row per mouse×epoch (median across days)
    if isempty(W), continue; end

    % Dual-panel epoch lines
    fh = figure('Color','w','Position',[80 60 900 900]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    dualPanelEpochLines(W, ylab);
    printpng(fh, fullfile(outDir, sprintf('dual_epoch_%s.png', safeName(ycol)))); close(fh);

    % Change-from-Pre strip/box (per group) with p-value annotations
    fh = figure('Color','w','Position',[100 80 800 500]);
   % NEW (choose 'holm' | 'fdr' | 'none')
changeFromPreStrip(W, ylab, 'none');  % annotate Passive vs Active per-epoch (Holm)
% Repeated-measures one-way ANOVA figure (per group, Δ from Pre)
rmOneWayWithinGroupFigure(W, ylab, ycol, outDir, true);  % true = delta-from-Pre

    printpng(fh, fullfile(outDir, sprintf('delta_from_pre_%s.png', safeName(ycol)))); close(fh);

ST  = runStats(W, ycol);
STd = runDeltaStats(W, ycol);
ST.metric  = repmat(string(ycol), height(ST),1);
ST.label   = repmat(string(ylab), height(ST),1);
STd.label  = repmat(string(ylab), height(STd),1);

STrm = runRMAnovaWithinGroup(W, ycol, true);   % true = delta-from-Pre
STrm.label  = repmat(string(ylab), height(STrm),1);

allStats = [allStats; ST; STd; STrm];  %#ok<AGROW>      %#ok<AGROW>
end

% Global multiple-comparison corrections across ALL tests for ALL metrics
if ~isempty(allStats) && ismember('p', allStats.Properties.VariableNames)
    allStats.p_FDR  = fdrBH(allStats.p);
    allStats.p_Holm = holmBonferroni(allStats.p);
    writetable(allStats, fullfile(outDir,'stats_summary.csv'));
    disp(allStats(:,{'metric','test','effect','level','p','p_FDR','p_Holm','N'}));
end

% Correlation + PCA examples (per mouse-epoch, core behavior metrics)
core = {'lick_freq_per_min','rew_freq_per_min','lick_meanDur_s','lick_medianIEI_s','pupil_mean'};
core = core(ismember(core, D.Properties.VariableNames));
if numel(core) >= 3
    Wcore = perMouseEpochTable(D, core);   % one row/mouse×epoch, columns=core metrics
    corrPCA(Wcore, outDir);
end

fprintf('Done. Plots + stats in:\n  %s\n', outDir);
end

%% ================= helpers =================
function T2 = ensureString(T2, nm)
    if ~ismember(nm, T2.Properties.VariableNames), T2.(nm) = repmat("",height(T2),1); return; end
    if ~isstring(T2.(nm)), T2.(nm) = string(T2.(nm)); end
end
function pick = pickPctColumn(cols, preferKey)
    pick = '';
    if isempty(cols), return; end
    pctCols = cols(contains(cols,'Pct'));
    if isempty(pctCols), return; end
    pref = pctCols(contains(lower(pctCols), lower(preferKey)));
    if ~isempty(pref), pick = pref{1}; else, pick = pctCols{1}; end
end
function nm = colLike(T, key)
    nm = '';
    if istable(T)
        V = T.Properties.VariableNames;
        hit = find(strcmpi(V,key) | contains(lower(V), lower(key)), 1, 'first');
        if ~isempty(hit), nm = V{hit}; end
    end
end
function s = safeName(nm), s = regexprep(nm,'[^a-zA-Z0-9]+','_'); end
function printpng(fh, fn)
    set(fh,'PaperPositionMode','auto');
    try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
end

function E = epochOfDay(d)
    e = strings(size(d));
    e(d>=1  & d<=5 ) = "Pre";
    e(d>=6  & d<=8 ) = "During";
    e(d>=9  & d<=11) = "Post";
    e(d>=12 & d<=14) = "Withdrawal";
    e(d>=15        ) = "Re-exposure";
    E = categorical(e, ["Pre","During","Post","Withdrawal","Re-exposure"], 'Ordinal',true);
end

function HadPassive = classifyHadPassive(S, D)
    mice = unique(S.mouse_key,'stable');
    HadPassive = false(height(D),1);
    havePar = ismember('Session_Paradigm', S.Properties.VariableNames);
    for i = 1:numel(mice)
        rS  = S.mouse_key==mice(i);
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
        elseif any(ismember(double(di),6:8) & (ip>=1 | isnan(ip)),'all')
            flag = true;
        end
        Dmask = D.mouse_key==mice(i);
        HadPassive(Dmask) = flag;
    end
end

%% ---------- fast per-session → per-day with extras + caching ----------
function [S, D] = fast_session_day_metrics(T, runDir)
    cacheMat = fullfile(runDir, 'S_D_cache.mat');
    if exist(cacheMat,'file')
        L = load(cacheMat,'S','D','cache_hash');
        if isfield(L,'cache_hash') && isequal(L.cache_hash, local_hash(T))
            fprintf('Loaded S/D from cache: %s\n', cacheMat);
            S = L.S; D = L.D; return;
        end
    end

    tic;
    % keep only columns we actually use + EXTRAS
    need = {'mouse_key','day_index','day_name','session_idx','Diameter_px', ...
            'Lick_TTL','Injector_TTL','CamTime_rel_s','PupilTimestamp_s', ...
            'CamTime_s','PlotTime_s_30fps','Session_Paradigm','isPassive','RequirementLast'};
    V = T.Properties.VariableNames;
    keepExtras = V( contains(V,'Immersion_Latency','IgnoreCase',true) | ...
                    contains(V,'TST_','IgnoreCase',true)            | ...
                    contains(V,'HOT_','IgnoreCase',true) );
    need = unique([need, keepExtras], 'stable');
    need = intersect(need, V, 'stable');
    T = T(:, need);

    % compact types
    if ~isstring(T.mouse_key),        T.mouse_key = string(T.mouse_key); end
    if ~isstring(T.day_name),         T.day_name  = string(T.day_name);  end
    if ismember('Session_Paradigm', T.Properties.VariableNames) && ~isstring(T.Session_Paradigm)
        T.Session_Paradigm = string(T.Session_Paradigm);
    end
    T.mouse_key  = categorical(T.mouse_key);
    T.day_name   = categorical(T.day_name);
    if ismember('Session_Paradigm', T.Properties.VariableNames)
        T.Session_Paradigm = categorical(T.Session_Paradigm);
    end
    if ismember('isPassive', T.Properties.VariableNames) && ~isnumeric(T.isPassive)
        T.isPassive = double(T.isPassive);
    end

    % TTLs to logical
    if ismember('Lick_TTL', T.Properties.VariableNames)
        T.Lick_TTL(isnan(T.Lick_TTL)) = 0; T.Lick_TTL = T.Lick_TTL > 0.5;
    end
    if ismember('Injector_TTL', T.Properties.VariableNames)
        T.Injector_TTL(isnan(T.Injector_TTL)) = 0; T.Injector_TTL = T.Injector_TTL > 0.5;
    end

    % group by session
    [g, keys_mouse, keys_day, keys_dayname, keys_sess] = findgroups( ...
        T.mouse_key, T.day_index, T.day_name, T.session_idx);
    nG = max(g);
    fprintf('Computing per-session metrics for %d sessions...\n', nG);

    % preallocate S
    S = table();
    S.mouse_key        = removecats(keys_mouse);
    S.day_index        = double(keys_day);
    S.day_name         = removecats(keys_dayname);
    S.session_idx      = double(keys_sess);
    S.RequirementLast  = nan(nG,1);
    S.isPassive        = nan(nG,1);
    S.SessionMinutes   = nan(nG,1);
    S.Session_Paradigm = strings(nG,1);
    vars = {'lick_n','lick_freq_per_min','lick_meanDur_s','lick_totalDur_s','lick_medianIEI_s', ...
            'rew_n','rew_freq_per_min','rew_meanDur_s','rew_totalDur_s','rew_medianIRI_s','pupil_mean'};
    for v = vars, S.(v{1}) = nan(nG,1); end
    for j=1:numel(keepExtras), S.(keepExtras{j}) = nan(nG,1); end

    tb_all = pickTimebase_fast(T);
    lastTick = tic; step = max(1, floor(nG/50));
    for k = 1:nG
        idx = (g == k);
        tb  = tb_all(idx);
        dur = finiteRange_fast(tb);
        S.SessionMinutes(k) = dur/60;

        if ismember('RequirementLast',T.Properties.VariableNames)
            S.RequirementLast(k) = mean(double(T.RequirementLast(idx)),'omitnan');
        end
        if ismember('isPassive',T.Properties.VariableNames)
            ip = double(T.isPassive(idx)); ip = ip(isfinite(ip));
            if ~isempty(ip), S.isPassive(k) = mode(round(ip)); end
        end
        if ismember('Session_Paradigm', T.Properties.VariableNames)
            sp = T.Session_Paradigm(idx);
            if iscategorical(sp), sp = string(mode(sp(~ismissing(sp))));
            else, sp = string(sp); sp = sp(~ismissing(sp)); if isempty(sp), sp=""; else, sp=sp(1); end
            end
            S.Session_Paradigm(k) = sp;
        end

        if ismember('Diameter_px',T.Properties.VariableNames)
            S.pupil_mean(k) = mean(double(T.Diameter_px(idx)),'omitnan');
        end
        if ismember('Lick_TTL',T.Properties.VariableNames)
            [n,md,td,iei] = eventMetrics_fast(tb, logical(T.Lick_TTL(idx)));
            S.lick_n(k)=n; S.lick_meanDur_s(k)=md; S.lick_totalDur_s(k)=td; S.lick_medianIEI_s(k)=iei;
            if S.SessionMinutes(k)>0, S.lick_freq_per_min(k)=n/S.SessionMinutes(k); end
        end
        if ismember('Injector_TTL',T.Properties.VariableNames)
            [n,md,td,iri] = eventMetrics_fast(tb, logical(T.Injector_TTL(idx)));
            S.rew_n(k)=n; S.rew_meanDur_s(k)=md; S.rew_totalDur_s(k)=td; S.rew_medianIRI_s(k)=iri;
            if S.SessionMinutes(k)>0, S.rew_freq_per_min(k)=n/S.SessionMinutes(k); end
        end
        for j=1:numel(keepExtras)
            col = keepExtras{j};
            if ismember(col, T.Properties.VariableNames)
                S.(col)(k) = mean(double(T.(col)(idx)),'omitnan');
            end
        end
        if mod(k, step)==0 && toc(lastTick)>0.5
            fprintf('  %d/%d (%.0f%%)\n', k, nG, 100*k/nG); lastTick = tic;
        end
    end
    fprintf('Per-session metrics done in %.1f s.\n', toc);

    % per-day collapse (median across sessions)
    fprintf('Collapsing to per-day medians...\n');
    [g2, mk2, di2, dn2] = findgroups(S.mouse_key, S.day_index, S.day_name);
    D = table(removecats(mk2), double(di2), removecats(dn2), ...
              'VariableNames',{'mouse_key','day_index','day_name'});
    baseList  = [{'RequirementLast','isPassive','SessionMinutes'}, vars];
    extraList = intersect(keepExtras, S.Properties.VariableNames, 'stable');
    list = unique([baseList, extraList], 'stable');
    for v = list
        D.(v{1}) = splitapply(@(x) median(x,'omitnan'), S.(v{1}), g2);
    end

    % cache
    cache_hash = local_hash(T(:, {'mouse_key','day_index','session_idx'}));
    try, save(cacheMat,'S','D','cache_hash','-v7.3'); fprintf('Cached S/D to: %s\n', cacheMat); catch, end
end

% --- tiny fast helpers
function tb = pickTimebase_fast(T)
    cands = {'CamTime_rel_s','PupilTimestamp_s','CamTime_s','PlotTime_s_30fps'};
    tb = nan(height(T),1);
    for i=1:numel(cands)
        if ismember(cands{i}, T.Properties.VariableNames)
            v = double(T.(cands{i})); if any(isfinite(v)), tb = v; return; end
        end
    end
end
function r = finiteRange_fast(x)
    x = double(x(:)); x = x(isfinite(x));
    if isempty(x), r = 0; else, r = max(x)-min(x); end
end
function [n, meanDur, totalDur, medianIEI] = eventMetrics_fast(t, ttl)
    t   = double(t(:)); ttl = logical(ttl(:));
    good = isfinite(t) & ~isnan(ttl); t=t(good); ttl=ttl(good);
    if numel(t)<2, n=0; meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
    dt  = diff(t); md = median(dt(isfinite(dt))); if ~isfinite(md), md=1/30; end
    t(2:end) = max(t(2:end), t(1:end-1)+md*0.5);
    d   = diff([false; ttl; false]); on = find(d==1); off = find(d==-1)-1;
    n   = numel(on);
    if n==0, meanDur=NaN; totalDur=0; medianIEI=NaN; return; end
    edges  = [t; t(end)+md];
    segDur = sum(edges(off+1) - edges(on), 2, 'omitnan');
    meanDur  = mean(segDur,'omitnan'); totalDur = sum(segDur,'omitnan');
    if n>=2, medianIEI = median(diff(t(on)),'omitnan'); else, medianIEI = NaN; end
end
function rmOneWayWithinGroupFigure(W, ylab, ycol, outDir, useDelta)
% Draws two panels (Passive / Active) showing within-group changes across time
% and annotates: omnibus RM-ANOVA p and Tukey-Kramer post-hoc sig pairs.

    if nargin<5, useDelta = true; end
    E = ["During","Post","Withdrawal","Re-exposure"];
    groups = ["Passive","Active"];
    safeE  = matlab.lang.makeValidName(cellstr(E));   % {'During','Post','Withdrawal','Re_exposure'}
    xMap   = containers.Map(cellstr(E), num2cell(1:numel(E)));

    fh = figure('Color','w','Position',[90 90 980 420]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    for gi = 1:numel(groups)
        G = groups(gi);
        nexttile; hold on; title(sprintf('%s — %s', ylab, G));

        % ----- build long and wide tables -----
        if useDelta
            base = groupsummary(W(W.Group==G & W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
            base.Properties.VariableNames{'median_value'}='base';
            Tg = outerjoin(W(W.Group==G & ismember(W.Epoch,E),:), base(:,{'mouse_key','base'}), ...
                           'Keys','mouse_key','MergeKeys',true);
            Tg.value = Tg.value - Tg.base;
        else
            Tg = W(W.Group==G & ismember(W.Epoch,E),:);
        end

        % Boxcharts (all available mice per epoch)
        for ei = 1:numel(E)
            Ge = Tg(Tg.Epoch==E(ei) & isfinite(Tg.value),:);
            if isempty(Ge), continue; end
            b = boxchart(repmat(ei,height(Ge),1), Ge.value, 'BoxWidth',0.25, ...
                         'MarkerStyle','none','BoxFaceAlpha',0.22);
            if G=="Passive", b.BoxFaceColor=[0 0.45 0.74]; else, b.BoxFaceColor=[0.85 0.33 0.10]; end
            scatter(repmat(ei,height(Ge),1), Ge.value, 26, 'k', 'filled', 'MarkerFaceAlpha',0.75);
        end
        set(gca,'XTick',1:numel(E),'XTickLabel',E); grid on; box on;
        yline(0,'k:'); ylabel(sprintf('%s – Pre', ylab));

        % ----- fit RM model on complete cases only (needed for ranova) -----
        Wg   = unstack(Tg(:,{'mouse_key','Epoch','value'}), 'value', 'Epoch');
        have = safeE(ismember(safeE, Wg.Properties.VariableNames));
        ok   = all(isfinite(Wg{:,have}),2);
        Wc   = Wg(ok, [{'mouse_key'}, have]);

        if height(Wc) < 3
            text(0.02,0.95,'N<3 complete mice','Units','normalized','FontWeight','bold');
            continue
        end

        % map the wide-column names back to original epoch labels as a COLUMN vector
orig = strings(numel(have),1);
for i = 1:numel(have)
    orig(i) = E(strcmp(have{i}, safeE));
end
within = table(categorical(orig,'Ordinal',true), 'VariableNames', {'Epoch'});

        rm  = fitrm(Wc, sprintf('%s ~ 1', strjoin(have,',')), 'WithinDesign', within);
        RA  = ranova(rm,'WithinModel','Epoch');

        % robust term label extraction
        termLabels = ranovaTermLabels(RA);
        idx = find(strcmpi(termLabels,'Epoch') | contains(lower(termLabels),'epoch'), 1, 'first');
        pOmni = NaN;
        if ~isempty(idx) && ismember('pValue', RA.Properties.VariableNames)
            pOmni = RA.pValue(idx);
        end
        % show omnibus p in the panel title
        t = title(sprintf('%s — %s (RM-ANOVA p=%.3g)', ylab, G, pOmni)); set(t,'FontWeight','bold');

        % Tukey-Kramer post-hoc; draw only significant pairs (alpha=0.05)
        try
            C = multcompare(rm,'Epoch','ComparisonType','tukey-kramer');
        catch
            C = multcompare(rm,'Epoch','ComparisonType','bonferroni');
        end
        vnames = string(C.Properties.VariableNames);
        L1 = 'LevelA'; L2 = 'LevelB';
        if ~ismember(L1,vnames)
            i1 = find(endsWith(vnames,'_1'),1,'first'); i2 = find(endsWith(vnames,'_2'),1,'first');
            if ~isempty(i1), L1 = char(vnames(i1)); end
            if ~isempty(i2), L2 = char(vnames(i2)); end
        end
        if ismember('pValue',vnames), pv = C.pValue;
        elseif ismember('pValueAdj',vnames), pv = C.pValueAdj;
        else, pv = NaN(height(C),1);
        end

        yl = ylim; baseY = yl(2); pad = 0.07*diff(yl);
        used = zeros(1,numel(E));   % height increments by span
        for r = 1:height(C)
            a = string(C.(L1)(r)); b = string(C.(L2)(r));
            if ~isKey(xMap, char(a)) || ~isKey(xMap, char(b)), continue; end
            if ~isfinite(pv(r)) || pv(r) >= 0.05, continue; end
            x1 = xMap(char(a)); x2 = xMap(char(b));
            span = abs(x2-x1);  if span==0, span=1; end
            used(span) = used(span) + 1;
            y = baseY + pad*(used(span)-0.5 + 0.8*(span-1));
            plot([x1 x1 x2 x2],[y-0.4*pad y y y-0.4*pad],'k-','LineWidth',0.9);
            text(mean([x1 x2]), y+0.05*pad, starStr(pv(r)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9,'FontWeight','bold');
            yl(2) = max(yl(2), y+pad);
        end
        ylim(yl);
    end

    fn = fullfile(outDir, sprintf('rm_oneway_%s.png', safeName(ycol)));
    set(fh,'PaperPositionMode','auto');
    try, exportgraphics(fh, fn, 'Resolution',180); catch, print(fh, fn, '-dpng','-r180'); end
    close(fh);
end

function h = local_hash(Tkeys)
    try
        raw = [uint32(double(grp2idx(categorical(Tkeys.mouse_key)))), ...
               uint32(double(Tkeys.day_index)), uint32(double(Tkeys.session_idx))];
        h = uint64(sum(uint64(raw(:)).*1664525 + 1013904223));
    catch, h = now; end
end
function STrm = runRMAnovaWithinGroup(W, ycol, useDelta)
% Repeated-measures one-way ANOVA across time within each group.
% Within factor = Epoch (During, Post, Withdrawal, Re-exposure).
% If useDelta==true, subtract per-mouse Pre first (recommended).

    if nargin<3, useDelta = true; end
    E = ["During","Post","Withdrawal","Re-exposure"];
    rows = {};

    for G = ["Passive","Active"]

        % ----- build per-mouse wide table -----
        if useDelta
            base = groupsummary(W(W.Group==G & W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
            base.Properties.VariableNames{'median_value'} = 'base';
            Tg = outerjoin(W(W.Group==G & ismember(W.Epoch,E),:), base(:,{'mouse_key','base'}), ...
                           'Keys','mouse_key','MergeKeys',true);
            Tg.value = Tg.value - Tg.base;
        else
            Tg = W(W.Group==G & ismember(W.Epoch,E),:);
        end

        Wg    = unstack(Tg(:,{'mouse_key','Epoch','value'}), 'value', 'Epoch');  % rows=mice, cols=epochs
        safeE = matlab.lang.makeValidName(cellstr(E));                           % {'During','Post','Withdrawal','Re_exposure'}
        have  = safeE(ismember(safeE, Wg.Properties.VariableNames));
        if numel(have) < 3
            rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Epoch', NaN, height(Wg)}; %#ok<AGROW>
            continue
        end
        mask = all(isfinite(Wg{:,have}),2);
        Wc   = Wg(mask, [{'mouse_key'}, have]);
        if height(Wc) < 3
            rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Epoch', NaN, height(Wc)}; %#ok<AGROW>
            continue
        end

        % within-design table with original epoch order/labels
        orig = strings(numel(have),1);
        for i=1:numel(have)
            orig(i) = E(strcmp(have{i}, safeE));
        end
        within = table(categorical(orig,'Ordinal',true), 'VariableNames', {'Epoch'});

        % ----- fit repeated-measures and test omnibus effect -----
        measStr = strjoin(have, ',');
        rm  = fitrm(Wc, sprintf('%s ~ 1', measStr), 'WithinDesign', within);
        RA  = ranova(rm, 'WithinModel', 'Epoch');

        % Get the row corresponding to 'Epoch' regardless of MATLAB version
        termLabels = ranovaTermLabels(RA);  % helper below
        idx = find(strcmpi(termLabels,'Epoch') | contains(lower(termLabels),'epoch'), 1, 'first');
        if isempty(idx)
            % fallback: first non-error row
            idx = find(~contains(lower(termLabels),'error'), 1, 'first');
        end
        if isempty(idx) || ~ismember('pValue', RA.Properties.VariableNames)
            pOmni = NaN;
        else
            pOmni = RA.pValue(idx);
        end
        rows(end+1,:) = {'RM-ANOVA','Group='+string(G),'Epoch', pOmni, height(Wc)}; %#ok<AGROW>

        % ----- post-hoc multiple comparisons across epochs -----
        try
            C = multcompare(rm, 'Epoch', 'ComparisonType', 'tukey-kramer');
        catch
            C = multcompare(rm, 'Epoch', 'ComparisonType', 'bonferroni');
        end

        % level column names vary by release; pick them robustly
        L1 = 'LevelA'; L2 = 'LevelB';
        vnames = string(C.Properties.VariableNames);
        if ~ismember(L1, vnames)
            i1 = find(endsWith(vnames,'_1'), 1, 'first');
            i2 = find(endsWith(vnames,'_2'), 1, 'first');
            if ~isempty(i1) && ~isempty(i2)
                L1 = char(vnames(i1)); L2 = char(vnames(i2));
            end
        end
        % p-value column can also differ
        if ismember('pValue', vnames), pv = C.pValue;
        elseif ismember('pValueAdj', vnames), pv = C.pValueAdj;
        else, pv = NaN(height(C),1);
        end

        for r = 1:height(C)
            lvl = string(C.(L1)(r)) + "|" + string(C.(L2)(r));
            rows(end+1,:) = {'RM-ANOVA posthoc','Group='+string(G), lvl, pv(r), height(Wc)}; %#ok<AGROW>
        end
    end

    STrm = cell2table(rows, 'VariableNames', {'test','effect','level','p','N'});
    STrm.metric = repmat(string(ycol), height(STrm), 1);
end

% Helper to read labels from ranova result across MATLAB versions
function labs = ranovaTermLabels(tbl)
    if ismember('Term', tbl.Properties.VariableNames)
        labs = string(tbl.Term);
    elseif ~isempty(tbl.Properties.RowNames)
        labs = string(tbl.Properties.RowNames);
    else
        % last-resort: try to find a “term-like” variable
        v = tbl.Properties.VariableNames;
        cand = v(contains(lower(v),'term') | contains(lower(v),'effect') | contains(lower(v),'within'));
        if ~isempty(cand)
            labs = string(tbl.(cand{1}));
        else
            labs = strings(height(tbl),1);
        end
    end
end

%% ---------- plotting & stats ----------
function W = perMouseEpochTable(D, ycols)
    if ischar(ycols) || isstring(ycols), ycols = cellstr(string(ycols)); end
    E = epochOfDay(double(D.day_index));
    D2 = D; D2.Epoch = E;

    if numel(ycols)==1
        y = ycols{1};
        G = groupsummary(D2(:,{'mouse_key','Group','Epoch',y}), {'mouse_key','Group','Epoch'}, 'median', y);
        W = renamevars(G, "median_"+y, "value");
    else % wide for correlation/PCA
        key = unique(D2(:,{'mouse_key','Group','Epoch'}),'rows','stable');
        W = key;
        for j=1:numel(ycols)
            y = ycols{j};
            G = groupsummary(D2(:,{'mouse_key','Group','Epoch',y}), {'mouse_key','Group','Epoch'}, 'median', y);
            W = outerjoin(W, renamevars(G, "median_"+y, y), 'Keys',{'mouse_key','Group','Epoch'}, 'MergeKeys',true);
        end
    end
end

function COL = palette()
    COL.passive = [0 0.45 0.74];       % blue
    COL.active  = [0.85 0.33 0.10];    % orange
    COL.grey    = 0.7*[1 1 1];         % light gray for individuals
    COL.epochOrder = {'Pre','During','Post','Withdrawal','Re-exposure'};
    COL.epochMark  = {'o','^','s','d','p'};
end

function dualPanelEpochLines(W, ylab)
    COL = palette(); epochs = COL.epochOrder;
    nexttile; hold on; title(sprintf('HadPassive — %s', ylab));
    drawPanel(W(W.Group=="Passive",:),  epochs, COL, ylab);
    nexttile; hold on; title(sprintf('Active-only (NoPassive) — %s', ylab));
    drawPanel(W(W.Group=="Active",:),   epochs, COL, ylab);
end

function drawPanel(Wg, epochs, COL, ylab)
    if isempty(Wg)
        text(0.02,0.95,'N=0 mice','Units','normalized');
        set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs); xlim([0.7 numel(epochs)+0.3]); grid on; box on; ylabel(ylab); return;
    end
    mice = unique(Wg.mouse_key,'stable');
    for i=1:numel(mice)
        s = Wg(Wg.mouse_key==mice(i),:);
        [x,ord] = sort(categorical(s.Epoch,epochs,'Ordinal',true)); y = s.value(ord);
        plot(double(x), y, '-o', 'Color',COL.grey, 'MarkerSize',3, 'LineWidth',0.9, ...
             'MarkerFaceColor',COL.grey*0.9, 'MarkerEdgeColor','none');
    end
    M = groupsummary(Wg, 'Epoch', 'mean', 'value');
    E = groupsummary(Wg, 'Epoch', @(x) std(x,'omitnan')./sqrt(sum(isfinite(x))), 'value');
    x = double(categorical(M.Epoch,epochs,'Ordinal',true)); y = M.mean_value; e = E.fun1_value;
    fill([x; flipud(x)], [y-e; flipud(y+e)], [0 0 0], 'FaceAlpha',0.12, 'EdgeColor','none');
    plot(x, y, 'k-', 'LineWidth',2.0); errorbar(x, y, e, 'k','LineWidth',1.1,'CapSize',8);
    set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs); xlim([0.7 numel(epochs)+0.3]); grid on; box on; ylabel(ylab);
    text(0.02,0.95,sprintf('N=%d mice', numel(mice)),'Units','normalized','FontWeight','bold');
end

function changeFromPreStrip(W, ylab, adjustMethod)
% adjustMethod: 'holm' (default), 'fdr', or 'none'
    if nargin<3 || isempty(adjustMethod), adjustMethod = 'holm'; end

    COL    = palette();
    epochs = COL.epochOrder;
    W = W(~isnan(W.value),:);

    % per-mouse Pre baseline and deltas
    base = groupsummary(W(W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
    base.Properties.VariableNames{'median_value'} = 'base';
    W = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key', 'MergeKeys',true);
    W.delta = W.value - W.base;

    hold on
    groups = {'Passive','Active'};
    off    = [-0.12 0.12];
    hG(1) = scatter(nan,nan,36,'o','filled','MarkerFaceColor',COL.passive,'MarkerEdgeColor','k','DisplayName','Passive');
    hG(2) = scatter(nan,nan,36,'o','filled','MarkerFaceColor',COL.active, 'MarkerEdgeColor','k','DisplayName','Active');

    topY = -inf(1,numel(epochs));
    for g=1:numel(groups)
        G = W(W.Group==groups{g},:);
        thisCol = COL.passive; if g==2, thisCol = COL.active; end
        for e=1:numel(epochs)
            Ge = G(G.Epoch==epochs{e} & isfinite(G.delta),:);
            if isempty(Ge), continue; end
            x = e + off(g);
            b = boxchart(repmat(x,height(Ge),1), Ge.delta, 'BoxWidth',0.18, 'MarkerStyle','none', 'BoxFaceAlpha',0.18);
            if isprop(b,'BoxFaceColor'),     b.BoxFaceColor     = thisCol; end
            if isprop(b,'BoxEdgeColor'),     b.BoxEdgeColor     = thisCol; end
            if isprop(b,'WhiskerLineColor'), b.WhiskerLineColor = 0.3*[1 1 1]; end
            scatter(repmat(x,height(Ge),1), Ge.delta, 24, 'filled', 'MarkerFaceColor',thisCol, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.75);
            topY(e) = max(topY(e), max(Ge.delta));
        end
    end

    set(gca,'XTick',1:numel(epochs),'XTickLabel',epochs);
    xlim([0.5 numel(epochs)+0.5]); yline(0,'k:'); grid on; box on;
    ylabel(sprintf('%s – Pre', ylab)); legend(hG,'Location','best');
    title('Change from Pre (per mouse)');

    % ---------- Passive vs Active per epoch: rank-sum, with optional correction ----------
    pRaw  = nan(1,numel(epochs)); valid = false(1,numel(epochs));
    for e = 1:numel(epochs)
        if strcmp(epochs{e},'Pre'), continue; end
        Pa = W(W.Group=="Passive" & W.Epoch==epochs{e} & isfinite(W.delta),:);
        Ac = W(W.Group=="Active"  & W.Epoch==epochs{e} & isfinite(W.delta),:);
        x = Pa.delta; y = Ac.delta;
        if numel(x)>=2 && numel(y)>=2
            pRaw(e) = ranksum(x,y);   % <-- RANK-SUM as requested
            valid(e)= true;
        end
    end
    pAdj = pRaw;
    switch lower(string(adjustMethod))
        case "holm", if any(valid), pAdj(valid) = holmBonferroni(pRaw(valid)); end
        case "fdr",  if any(valid), pAdj(valid) = fdrBH(pRaw(valid));          end
        otherwise    % 'none' -> leave raw p's
    end

    yl = ylim; rngY = diff(yl); pad = 0.06*rngY;
    for e = 1:numel(epochs)
        if ~valid(e), continue; end
        x1 = e + off(1);  x2 = e + off(2);
        y  = max(topY(e), yl(2)-2*pad) + pad;
        plot([x1 x1 x2 x2],[y-0.5*pad y y y-0.5*pad],'k-','LineWidth',1);
        tag = upper(char(adjustMethod)); if tag=="NONE", tag="raw"; end
        txt = sprintf('%s (%s %.3g)', starStr(pAdj(e)), tag, pAdj(e));
        text(mean([x1 x2]), y+0.2*pad, txt, 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontSize',9, 'FontWeight','bold');
        yl(2) = max(yl(2), y+0.8*pad);
    end
    ylim(yl);
end


% ---------- Nonparametric “simple” stats on raw per-epoch medians ----------
function ST = runStats(W, ycol)
% Within-group paired (Pre vs Post/Withdrawal/Re-exposure): signrank
% Between-group unpaired at each epoch: ranksum
    T = W(:,{'mouse_key','Group','Epoch','value'});
    T.mouse_key = categorical(T.mouse_key);
    T.Group     = categorical(T.Group);
    T.Epoch     = categorical(T.Epoch, {'Pre','During','Post','Withdrawal','Re-exposure'}, 'Ordinal',true);

    rows = {};

    % Within-group paired: Pre vs {Post, Withdrawal, Re-exposure}
    G = ["Passive","Active"];
    Ecmp = ["Post","Withdrawal","Re-exposure"];
    for g = 1:numel(G)
        for e = 1:numel(Ecmp)
            A = groupsummary(T(T.Group==G(g) & T.Epoch=="Pre",:),  'mouse_key','median','value');
            B = groupsummary(T(T.Group==G(g) & T.Epoch==Ecmp(e),:),'mouse_key','median','value');
            M = outerjoin(A(:,{'mouse_key','median_value'}), B(:,{'mouse_key','median_value'}), ...
                          'Keys','mouse_key','MergeKeys',true);
            xa = M.median_value_left; xb = M.median_value_right;
            good = isfinite(xa) & isfinite(xb);
            p = NaN; N = nnz(good);
            if N >= 3, p = signrank(xa(good), xb(good)); end
            rows(end+1,:) = {'WithinNP','Group='+string(G(g)),'Pre|'+Ecmp(e), p, N}; %#ok<AGROW>
        end
    end

    % Between-group at each epoch
    Eall = ["Pre","During","Post","Withdrawal","Re-exposure"];
    for e = 1:numel(Eall)
        Pa = T(T.Group=="Passive" & T.Epoch==Eall(e) & isfinite(T.value),:);
        Ac = T(T.Group=="Active"  & T.Epoch==Eall(e) & isfinite(T.value),:);
        x = Pa.value; y = Ac.value;
        p = NaN; N = numel(x)+numel(y);
        if numel(x)>=2 && numel(y)>=2, p = ranksum(x,y); end
        rows(end+1,:) = {'BetweenNP','Epoch='+string(Eall(e)),'Passive vs Active', p, N}; %#ok<AGROW>
    end

    ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
    ST.metric = repmat(string(ycol), height(ST),1);
end

% ---------- Nonparametric stats on delta-from-Pre (ALL pairs) ----------
function ST = runDeltaStats(W, ycol)
    % Build delta table (per mouse)
    base = groupsummary(W(W.Epoch=="Pre",:), 'mouse_key', 'median', 'value');
    base.Properties.VariableNames{'median_value'} = 'base';
    T = outerjoin(W, base(:,{'mouse_key','base'}), 'Keys','mouse_key', 'MergeKeys',true);
    T.delta = T.value - T.base;

    E = ["During","Post","Withdrawal","Re-exposure"];
    T = T(ismember(string(T.Epoch), E) & isfinite(T.delta),:);

    rows = {};
    % within-group, paired signrank across epoch pairs
    G = ["Passive","Active"];
    for g = 1:numel(G)
        Tg = T(T.Group==G(g),:);
        for i = 1:numel(E)-1
            for j = i+1:numel(E)
                A = outerjoin( ...
                    Tg(Tg.Epoch==E(i), {'mouse_key','delta'}), ...
                    Tg(Tg.Epoch==E(j), {'mouse_key','delta'}), ...
                    'Keys','mouse_key','MergeKeys',true);
                xa = A.delta_left; xb = A.delta_right;
                good = isfinite(xa) & isfinite(xb); N = nnz(good); p = NaN;
                if N >= 3, p = signrank(xa(good), xb(good)); end
                rows(end+1,:) = {'DeltaWithinNP', 'Group='+string(G(g)), E(i)+"|"+E(j), p, N}; %#ok<AGROW>
            end
        end
    end
    % between groups at each epoch (ranksum)
    for e = 1:numel(E)
        Pa = T(T.Group=="Passive" & T.Epoch==E(e),:);
        Ac = T(T.Group=="Active"  & T.Epoch==E(e),:);
        x = Pa.delta; y = Ac.delta; x=x(isfinite(x)); y=y(isfinite(y)); p = NaN; N = numel(x)+numel(y);
        if numel(x)>=2 && numel(y)>=2, p = ranksum(x, y); end
        rows(end+1,:) = {'DeltaBetweenNP', 'Epoch='+string(E(e)), 'Passive vs Active', p, N}; %#ok<AGROW>
    end
    % all cell-pairs (unpaired ranksum)
    cells = unique(T(:,{'Group','Epoch'}),'rows','stable');
    for i = 1:height(cells)-1
        for j = i+1:height(cells)
            A = T(T.Group==cells.Group(i) & T.Epoch==cells.Epoch(i),:);
            B = T(T.Group==cells.Group(j) & T.Epoch==cells.Epoch(j),:);
            x = A.delta; y = B.delta; x=x(isfinite(x)); y=y(isfinite(y)); p = NaN; N = numel(x)+numel(y);
            if numel(x)>=2 && numel(y)>=2, p = ranksum(x, y); end
            lab = string(cells.Group(i))+"|"+string(cells.Epoch(i))+" vs "+ ...
                  string(cells.Group(j))+"|"+string(cells.Epoch(j));
            rows(end+1,:) = {'DeltaAllPairsNP', 'CellPair', lab, p, N}; %#ok<AGROW>
        end
    end
    ST = cell2table(rows, 'VariableNames',{'test','effect','level','p','N'});
    ST.metric = repmat(string(ycol), height(ST),1);
end

function s = starStr(p)
    if ~isfinite(p), s = 'n/a'; return; end
    if p < 1e-4, s = '****';
    elseif p < 1e-3, s = '***';
    elseif p < 1e-2, s = '**';
    elseif p < 0.05, s = '*';
    else, s = 'n.s.';
    end
end

function q = fdrBH(p)
    p = double(p); m = numel(p); [ps,idx] = sort(p); q = nan(size(p));
    ranks = (1:m)'; adj = ps.*m./ranks;
    for i=m-1:-1:1, adj(i)=min(adj(i),adj(i+1)); end
    q(idx) = min(adj,1);
end

function p_holm = holmBonferroni(p)
    p = double(p(:)); [ps,idx] = sort(p); m = numel(p);
    adj = (m - (1:m)' + 1) .* ps;             % (m-k+1)*p_(k)
    for k = 2:m, adj(k) = max(adj(k), adj(k-1)); end
    p_holm = zeros(size(p)); p_holm(idx) = min(adj,1);
end

% ---------- correlation & PCA ----------
function corrPCA(Wwide, outDir)
    COL = palette();
    dat = Wwide; key = dat(:,1:3); X = dat{:,4:end};
    vars = dat.Properties.VariableNames(4:end);

    % correlation heatmap
    C = corr(X,'rows','pairwise');
    fh = figure('Color','w','Position',[60 60 650 550]);
    imagesc(C); axis image; colorbar; caxis([-1 1]);
    set(gca,'XTick',1:numel(vars),'XTickLabel',vars,'XTickLabelRotation',45);
    set(gca,'YTick',1:numel(vars),'YTickLabel',vars);
    title('Metric correlation (per mouse×epoch)'); grid on;
    printpng(fh, fullfile(outDir,'corr_heatmap.png')); close(fh);

    % PCA scatter: color = group, marker = epoch
    Xz = X; for j=1:size(Xz,2), Xz(:,j) = (Xz(:,j)-nanmean(Xz(:,j)))./nanstd(Xz(:,j)); end
    Xz(~isfinite(Xz))=0;
    [~,score,~,~,expl] = pca(Xz);

    fh = figure('Color','w','Position',[80 80 740 540]); hold on;
    grp = string(key.Group); epo = string(key.Epoch);
    colors = zeros(numel(grp),3);
    colors(grp=="Passive",:) = repmat(COL.passive, sum(grp=="Passive"), 1);
    colors(grp=="Active", :) = repmat(COL.active,  sum(grp=="Active"),  1);
    for e = 1:numel(COL.epochOrder)
        mask = (epo==COL.epochOrder{e});
        if any(mask)
            scatter(score(mask,1), score(mask,2), 46, colors(mask,:), 'filled', ...
                    'Marker', COL.epochMark{e}, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.85);
        end
    end
    hG(1) = scatter(nan,nan,60,'o','filled','MarkerFaceColor',COL.passive,'MarkerEdgeColor','k','DisplayName','Passive');
    hG(2) = scatter(nan,nan,60,'o','filled','MarkerFaceColor',COL.active, 'MarkerEdgeColor','k','DisplayName','Active');
    hE = gobjects(numel(COL.epochOrder),1);
    for e=1:numel(COL.epochOrder)
        hE(e) = plot(nan,nan,'k','LineStyle','none','Marker',COL.epochMark{e}, 'MarkerSize',8,'DisplayName',COL.epochOrder{e});
    end
    lg1 = legend(hG, 'Location','northeastoutside'); title(lg1,'Group');
    lg2 = legend(hE, 'Location','southeastoutside');  title(lg2,'Epoch');
    uistack(lg1,'top'); uistack(lg2,'top');

    xlabel(sprintf('PC1 (%.1f%%)', expl(1))); ylabel(sprintf('PC2 (%.1f%%)', expl(2)));
    grid on; box on; title('PCA of mouse×epoch metric vectors');
    printpng(fh, fullfile(outDir,'pca_scatter.png')); close(fh);
end
