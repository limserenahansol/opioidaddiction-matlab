function run_pr_licking_full_analysis
% run_pr_licking_full_analysis — full-data PR lick rate vs trial, by phase and group, + stats
%
% 1) Within- and between-trial metrics vs PR (Trial); Active vs Passive; per-phase figures.
% 2) For each, plots: Hz, licks/min, licks/10s, and a fourth view (within: licks/s over full
%    camera trial span; between: raw lick count in the ITI). Sanitize + IQR censor spikes.
%
% Defaults to ALL longitudinal data: latest run_* + figs/lick_MEGA/per_trial_features.csv
% Set USE_DEMO = true to run the small demo only.

%% ========== config ==========
USE_DEMO = false;              % true = demo only | false = all data (latest run_* lick MEGA + enrich between-trial)
% Step 7 EFA and lick MEGA must come from the SAME run_* (e.g. run_009):
%   <run_*>/figs/lick_MEGA/per_trial_features.csv
%   <run_*>/ALL_mice_longitudinal.csv
% Step 7 AI table is under <run_*>/figs/.../AI_mouse_summary.csv (used by run_step9_one_click).
OVERRIDE_RUN = 'K:\addiction_concate_Dec_2025\longitudinal_outputs\run_009';            % e.g. 'K:\longitudinal_outputs\run_009' — set if pr_step08_locate_data should not pick latest
SKIP_PLOTS = false;          % true = no PNG/FIG/boxplots (fast, tables only); false = full figures (default for analysis)
MAX_TRIAL_PLOT = 45;           % x-axis cap for PR trial index
WRITE_ENRICHED = true;         % save per_trial_enriched_step08.csv
SKIP_ENRICH_IF_CACHED = true;  % reuse enriched CSV in output if present
% QC: lick_rate = lick_count / time_to_reward; near-zero time explodes values and slopes.
% Recompute within-trial rate with a floor on TTR; cap licks/s for stats/plots (not a biology claim).
PR_LICK_MIN_TTR_S   = 0.02;    % s, minimum denominator for within-trial lick / s
PR_LICK_MAX_HZ      = 50;     % licks/s cap, within-trial
PR_LICK_MAX_HZ_BTWN = 22;     % stricter cap for between-trial (typical licking << this)
% Per (phase, PR-trial) IQR fence: censor y > Q3 + k*IQR (removes 1–2 subject spikes; needs min N)
PR_LICK_IQR_K       = 1.5;
PR_LICK_IQR_MIN_N   = 6;      % skip IQR censoring for that slice if fewer finite points

here = fileparts(mfilename('fullpath'));
if USE_DEMO
    outDir = fullfile(here, 'output_demo');
else
    outDir = fullfile(here, 'output_full');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if USE_DEMO
    % Toy CSV: G1 / G2 groups, no lick MEGA
    TT2 = readtable(fullfile(here, 'demo_data', 'demo_pr_trial_table.csv'));
    fprintf('USE_DEMO: %s\n', fullfile(here, 'demo_data', 'demo_pr_trial_table.csv'));
else
    p = pr_step08_locate_data(OVERRIDE_RUN);
    if exist(p.trial_csv, 'file') ~= 2
        error(['Missing %s\nSet USE_DEMO = true or run run_lick_mega_pipeline_new first, ', ...
            'or set OVERRIDE_RUN.'], p.trial_csv);
    end
    fprintf('Step08: runDir = %s\n', p.runDir);
    TT = readtable(p.trial_csv);
    if ~ismember('lick_rate', TT.Properties.VariableNames)
        error('per_trial_features.csv must contain lick_rate');
    end
    if ~ismember('Trial', TT.Properties.VariableNames)
        error('per_trial_features.csv must contain Trial');
    end
    T = readtable(p.long_csv);
    assert(ismember('Lick_TTL', T.Properties.VariableNames), 'ALL_mice_longitudinal must have Lick_TTL for between-trial');
    tb = pr_step08_pickTimebase(T);
    T.Lick_TTL(isnan(T.Lick_TTL)) = 0;
    T.Lick_TTL = T.Lick_TTL > 0.5;
    fEnr = fullfile(outDir, 'per_trial_enriched_step08.csv');
    if SKIP_ENRICH_IF_CACHED && exist(fEnr, 'file') == 2
        fprintf(['Loading cached enriched table: %s\n', ...
            '  (delete to recompute between-trial or add new columns like BetweenTrialNlicks)\n'], fEnr);
        TT2 = readtable(fEnr);
    else
        fprintf('Enriching between-trial lick rates (this may take a few minutes)...\n');
        TT2 = pr_step08_enrich_between_lick(T, tb, TT);
        if WRITE_ENRICHED
            writetable(TT2, fEnr);
            fprintf('Wrote: %s\n', fEnr);
        end
    end
    if ~ismember('BetweenTrialNlicks', TT2.Properties.VariableNames)
        fprintf('Adding BetweenTrialNlicks / refreshing between-trial features...\n');
        TT2 = pr_step08_enrich_between_lick(T, tb, TT);
        if WRITE_ENRICHED
            writetable(TT2, fEnr);
        end
    end
    fprintf('Enriching trial time span (full trial width on camera)...\n');
    TT2 = pr_step08_enrich_trial_span(T, tb, TT2);
end

U = pr_step08_normalize_trial_table(TT2);
if ~ismember('WithinTrialLickRate', U.Properties.VariableNames)
    error('normalize_trial_table must set WithinTrialLickRate');
end
if ~USE_DEMO && ismember('time_to_reward', TT2.Properties.VariableNames) && ismember('lick_count', TT2.Properties.VariableNames)
    ttr = max(double(TT2.time_to_reward), PR_LICK_MIN_TTR_S);
    U.WithinTrialLickRate = double(TT2.lick_count) ./ ttr;
    fprintf('Within-trial lick rate recomputed as lick_count / max(time_to_reward, %.3g s).\n', PR_LICK_MIN_TTR_S);
end
U = pr_step08_sanitize_lick_rate_table(U, PR_LICK_MAX_HZ, PR_LICK_MAX_HZ_BTWN);

% Phase & group
if USE_DEMO
    U.Group = pr_step08_cell_group(U);
    U.phase = repmat({'All'}, height(U), 1);
else
    if ~ismember('Period', U.Properties.VariableNames)
        U.Period = repmat("Pre", height(U), 1);
    end
    U.Group = pr_step08_cell_group(U);
    U.phase = cellstr(string(U.Period));
    ux = ismember(U.phase, {'<undef>', '<excluded>'});
    U.phase(ux) = {'(drop)'};
end

% Per (phase, PR-trial) high-end outliers (drives one-mouse spikes in group means)
U = pr_step08_censor_trialwise_iqr(U, 'WithinTrialLickRate', PR_LICK_IQR_K, PR_LICK_IQR_MIN_N);
if ismember('BetweenTrialLickRate', U.Properties.VariableNames)
    U = pr_step08_censor_trialwise_iqr(U, 'BetweenTrialLickRate', PR_LICK_IQR_K, PR_LICK_IQR_MIN_N);
end

% Derived rate views: /min, /10 s, and (within) licks / full camera trial width
U = pr_step08_derive_lick_scales(U, PR_LICK_MIN_TTR_S, PR_LICK_MAX_HZ, PR_LICK_MAX_HZ_BTWN);
if ismember('WithinTrialLickPerWhole', U.Properties.VariableNames)
    U = pr_step08_censor_trialwise_iqr(U, 'WithinTrialLickPerWhole', PR_LICK_IQR_K, PR_LICK_IQR_MIN_N);
end
if ismember('BetweenTrialNlicks', U.Properties.VariableNames)
    U = pr_step08_censor_trialwise_iqr(U, 'BetweenTrialNlicks', PR_LICK_IQR_K, PR_LICK_IQR_MIN_N);
end
hasB = ismember('BetweenTrialLickRate', U.Properties.VariableNames) && any(isfinite(U.BetweenTrialLickRate));

% Step 9 (run_step9_one_click) needs this file — write before the heavy plot loop (recoverable if figs error).
pr_step08_write_mousemean_pr_table(U, outDir);

if ~SKIP_PLOTS
    if USE_DEMO
        phases = {'All'};
    else
        phases = {'All', 'Pre', 'During', 'Post', 'Withdrawal', 'Reexposure'};
    end
    % {column in U, file tag piece, y-label, do slopes (only for Hz base columns)}
    rateBlocks = { ...
        {'WithinTrialLickRate',     'W_Hz',     'Within-trial lick rate (Hz)', true}, ...
        {'WithinTrialLickPerMin',   'W_lpm',    'Within-trial lick rate (licks / min)', false}, ...
        {'WithinTrialLickPer10s',   'W_p10s',   'Within-trial lick rate (licks / 10 s)', false}, ...
        {'WithinTrialLickPerWhole', 'W_whole',  'Within-trial lick density (1/s, full trial window)', false} ...
        };
    if hasB
        rateBlocks = [rateBlocks, { ...
            {'BetweenTrialLickRate',   'B_Hz',     'Between-trial lick rate (Hz)', true}, ...
            {'BetweenTrialLickPerMin',  'B_lpm',    'Between-trial lick rate (licks / min)', false}, ...
            {'BetweenTrialLickPer10s',  'B_p10s',   'Between-trial lick rate (licks / 10 s)', false}, ...
            {'BetweenTrialNlicks',      'B_nLicks', 'Licks in between-trial interval (count)', false} ...
            }];
    end
    allStats = table();
    allSlopes = table();
    for ph = 1:numel(phases)
        ptag = phases{ph};
        if strcmp(ptag, 'All')
            sub = U(~strcmp(U.phase, '(drop)'), :);
        else
            sub = U(strcmp(U.phase, ptag), :);
        end
        if height(sub) < 4
            fprintf('Skipping phase %s (n=%d rows too small)\n', ptag, height(sub));
            continue
        end
        for rb = 1:numel(rateBlocks)
            spec = rateBlocks{rb};
            rknm = spec{1};
            ttag = spec{2};
            ylab = spec{3};
            doSlopes = spec{4};
            if ~ismember(rknm, sub.Properties.VariableNames)
                continue
            end
            y = sub.(rknm);
            if ~any(isfinite(y))
                continue
            end
            S = sub; S.y = y; %#ok<SAGROW>
            tag = sprintf('%s_%s', ptag, ttag);
            pr_step08_figure_trialcurves(S, 'y', ptag, ylab, outDir, tag, MAX_TRIAL_PLOT);
            st = pr_step08_trialwise_stats(S, 'y', ptag, sprintf('%s_%s', strrep(rknm, 'TrialLickRate', 'TR'), ttag));
            if height(st) > 0
                if isempty(allStats) || width(allStats) < 1
                    allStats = st;
                else
                    allStats = [allStats; st];
                end
            end
            if doSlopes
                sl = pr_step08_slopes(S, 'y', ptag, rknm);
                if height(sl) > 0
                    if isempty(allSlopes) || width(allSlopes) < 1
                        allSlopes = sl;
                    else
                        allSlopes = [allSlopes; sl];
                    end
                end
            end
        end
    end
else
    fprintf('SKIP_PLOTS: skipping figures, trialwise stats, and slope outputs.\n');
    allStats = table();
    allSlopes = table();
end

if height(allStats) > 0
    writetable(allStats, fullfile(outDir, 'PR_lick_stats_trialwise.csv'));
end
if height(allSlopes) > 0
    writetable(allSlopes, fullfile(outDir, 'PR_lick_slopes_by_mouse.csv'));
    pr_step08_figure_slope_bars(allSlopes, outDir);
end
fprintf('Step 08 full analysis done. Output: %s\n', outDir);
end

% pr_step08_write_mousemean_pr_table is in pr_step08_write_mousemean_pr_table.m

% ---------------- timebase (same logic as lick MEGA) ----------------
function tb = pr_step08_pickTimebase(T)
cands = {'CamTime_rel_s', 'PupilTimestamp_s', 'CamTime_s', 'PlotTime_s_30fps', 'Frame'};
tb = nan(height(T), 1);
for i = 1:numel(cands)
    c = cands{i};
    if ismember(c, T.Properties.VariableNames)
        v = double(T.(c));
        if any(isfinite(v))
            tb = v;
            if strcmpi(c, 'Frame')
                tb = v / 30;
            end
            return
        end
    end
end
end

function U = pr_step08_sanitize_lick_rate_table(U, maxHz, maxHzBetween)
% Set implausible lick rates to NaN so means/SEM/slopes are not driven by 1e12 junk.
if nargin < 3, maxHzBetween = maxHz; end
colMax = { {'WithinTrialLickRate', maxHz}, {'BetweenTrialLickRate', maxHzBetween} };
for c = 1:numel(colMax)
    nm = colMax{c}{1};
    lim = colMax{c}{2};
    if ~ismember(nm, U.Properties.VariableNames)
        continue
    end
    y = U.(nm);
    y = double(y(:));
    bad = isfinite(y) & (y < 0 | y > lim);
    n = sum(bad);
    if n > 0
        y(bad) = NaN;
        U.(nm) = y;
        fprintf('Sanitized %s: set %d row(s) to NaN (|rate| out of [0, %.0f] 1/s).\n', nm, n, lim);
    end
end
end

function U = pr_step08_censor_trialwise_iqr(U, colname, kIqr, minN)
y = U.(colname);
y = double(y(:));
nCen = 0;
if ~ismember('phase', U.Properties.VariableNames)
    return
end
ph = U.phase;
tr = double(U.Trial);
phU = unique(ph, 'stable');
trU = unique(tr);
for ip = 1:numel(phU)
    for it = 1:numel(trU)
        r = strcmp(ph, phU{ip}) & (tr == trU(it)) & isfinite(y);
        ix = find(r);
        vals = y(ix);
        vals = vals(isfinite(vals));
        nv = numel(vals);
        if nv < minN
            continue
        end
        q1 = prctile(vals, 25);
        q3 = prctile(vals, 75);
        iqr = q3 - q1;
        if iqr < 1e-9
            iqr = max(1e-6, 0.1 * max(abs(median(vals, 'omitnan')), 1e-6));
        end
        cap = q3 + kIqr * iqr;
        for j = 1:numel(ix)
            jj = ix(j);
            if isfinite(y(jj)) && y(jj) > cap
                y(jj) = NaN;
                nCen = nCen + 1;
            end
        end
    end
end
U.(colname) = y;
if nCen > 0
    fprintf('IQR censor %s: set %d point(s) to NaN (Q3 + %.1f*IQR by phase and PR trial).\n', colname, nCen, kIqr);
end
end

function U = pr_step08_derive_lick_scales(U, minTtr, capHz, capBtw)
% licks/min, licks/10s, within-trial density on full camera span, ITI lick count for plots.
W = U.WithinTrialLickRate;
W = double(W(:));
U.WithinTrialLickPerMin   = 60  * W;
U.WithinTrialLickPer10s   = 10  * W;
if ismember('lick_count', U.Properties.VariableNames) && ismember('trial_span_s', U.Properties.VariableNames)
    lc  = double(U.lick_count(:));
    sp  = max(double(U.trial_span_s(:)), minTtr);
    pWh = lc ./ sp;
    bad = isfinite(pWh) & pWh > capHz;
    pWh(bad) = NaN;
    pWh(~isfinite(W)) = NaN;
    U.WithinTrialLickPerWhole = pWh;
else
    u = nan(size(W));
    u(:) = W;
    u(~isfinite(W)) = NaN;
    U.WithinTrialLickPerWhole = u;
end
if ismember('BetweenTrialLickRate', U.Properties.VariableNames)
    B = double(U.BetweenTrialLickRate(:));
    U.BetweenTrialLickPerMin  = 60 * B;
    U.BetweenTrialLickPer10s  = 10 * B;
else
    U.BetweenTrialLickPerMin  = nan(height(U), 1);
    U.BetweenTrialLickPer10s  = nan(height(U), 1);
end
if ismember('BetweenTrialNlicks', U.Properties.VariableNames)
    nL = double(U.BetweenTrialNlicks(:));
    if ismember('BetweenTrialLickRate', U.Properties.VariableNames)
        bR = double(U.BetweenTrialLickRate(:));
        nL(~isfinite(bR)) = NaN;
    end
    U.BetweenTrialNlicks = nL;
end
end

function U = pr_step08_normalize_trial_table(TT)
U = TT;
U.mouse_key = cellstr(string(TT.mouse_key));
U.Trial = double(U.Trial);
if ismember('lick_count', TT.Properties.VariableNames)
    U.lick_count = double(TT.lick_count);
end
if ismember('trial_span_s', TT.Properties.VariableNames)
    U.trial_span_s = double(TT.trial_span_s);
end
if ismember('BetweenTrialNlicks', TT.Properties.VariableNames)
    U.BetweenTrialNlicks = double(TT.BetweenTrialNlicks);
end
if ismember('lick_rate', TT.Properties.VariableNames)
    U.WithinTrialLickRate = double(TT.lick_rate);
end
if ismember('WithinTrialLickRate', TT.Properties.VariableNames)
    U.WithinTrialLickRate = double(TT.WithinTrialLickRate);
end
if ismember('Period', U.Properties.VariableNames)
    U.Period = string(TT.Period);
else
    U.Period = repmat("Pre", height(U), 1);
end
end

function g = pr_step08_cell_group(U)
if ismember('Group', U.Properties.VariableNames)
    v = U.Group;
    if iscell(v)
        g = cellfun(@(x) char(x), v, 'uni', 0);
    else
        g = cellstr(string(v));
    end
else
    g = repmat({'Unknown'}, height(U), 1);
end
end

function pr_step08_figure_trialcurves(S, ycol, ptag, ylabTitle, outDir, tag, maxT)
% Individuals (thin) + group mean ± SEM
S.mouse_key = cellstr(string(S.mouse_key));
S.Group = cellstr(string(S.Group));
G = S.(ycol);
tr = double(S.Trial);
gstr = S.Group;
[ug, ~, ~] = unique(gstr, 'stable');
cmap = [0.25, 0.45, 0.85; 0.90, 0.35, 0.30; 0.35, 0.65, 0.35; 0.6, 0.3, 0.7; 0.2, 0.55, 0.5];
f = figure('Color', 'w', 'Position', [80, 80, 960, 540]);
hold on
mice = unique(S.mouse_key, 'stable');
for mi = 1:numel(mice)
    r = strcmp(S.mouse_key, mice{mi});
    ttr = tr(r);
    yy = G(r);
    tmu = unique(ttr);
    mu = nan(size(tmu));
    for u = 1:numel(tmu)
        f2 = ttr == tmu(u) & isfinite(yy);
        mu(u) = mean(yy(f2), 'omitnan');
    end
    [tS, o] = sort(tmu);
    mS = mu(o);
    keep = tS <= maxT;
    tS = tS(keep);
    mS = mS(keep);
    if isempty(tS)
        continue
    end
    plot(tS, mS, '-', 'Color', [0.80, 0.80, 0.80], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for g = 1:numel(ug)
    col = cmap(mod(g - 1, size(cmap, 1)) + 1, :);
    tlist = unique(tr(strcmp(gstr, ug{g})));
    tlist = tlist(tlist <= maxT);
    tlist = tlist(:);
    mmu = nan(size(tlist));
    sse = nan(size(tlist));
    for ti = 1:numel(tlist)
        rr = strcmp(gstr, ug{g}) & (tr == tlist(ti)) & isfinite(G);
        mmu(ti) = mean(G(rr), 'omitnan');
        n = sum(rr);
        sse(ti) = std(G(rr), 0, 'omitnan') / max(1, sqrt(n));
    end
    if isempty(tlist)
        continue
    end
    ok = isfinite(tlist) & isfinite(mmu) & isfinite(sse);
    tlist = tlist(ok);
    mmu = mmu(ok);
    sse = sse(ok);
    if isempty(tlist)
        continue
    end
    % Name-Value avoids ambiguous line-spec parsing; column vectors for errorbar
    errorbar(tlist, mmu, sse, 'LineStyle', '-', 'Color', col, 'LineWidth', 2, ...
        'CapSize', 0, 'DisplayName', char(ug{g}));
end
xlabel('PR trial index (breakpoint step)');
ylabel(ylabTitle, 'Interpreter', 'none');
title(sprintf('%s | %s', ptag, strrep(ylabTitle, '_', ' ')));
grid on
box off
legend('Location', 'best');
fn = fullfile(outDir, sprintf('PR_lick_trialCurve_%s.png', tag));
try
    saveas(f, fn);
    saveas(f, strrep(fn, '.png', '.fig'));
catch
    saveas(f, fn);
end
close(f);
end

function st = pr_step08_trialwise_stats(S, ycol, ptag, rlabel)
% Per PR trial: two main groups (Active/Passive) or first two label levels
G = S.(ycol);
tr = double(S.Trial);
gL = cellstr(string(S.Group));
iA = strcmpi(gL, 'active');
iP = strcmpi(gL, 'passive');
if sum(iA) < 1 || sum(iP) < 1
    ugS = unique(gL, 'stable');
    if numel(ugS) >= 2
        iA = strcmp(gL, ugS{1});
        iP = strcmp(gL, ugS{2});
    else
        iA = false(size(gL));
        iP = false(size(gL));
    end
end
tlist = unique(tr);
st = table();
for t = 1:numel(tlist)
    r = tr == tlist(t) & isfinite(G);
    yA = G(r & iA);
    yP = G(r & iP);
    yA = yA(isfinite(yA));
    yP = yP(isfinite(yP));
    p_t = NaN;
    p_m = NaN;
    if numel(yA) >= 2 && numel(yP) >= 2
        [~, p_t] = ttest2(yA, yP, 'Vartype', 'unequal');
    end
    if numel(yA) >= 1 && numel(yP) >= 1
        p_m = ranksum(yA, yP);
    end
    row = table({char(ptag)}, {char(rlabel)}, tlist(t), mean(yA, 'omitnan'), pr_step08_sem(yA), numel(yA), ...
        mean(yP, 'omitnan'), pr_step08_sem(yP), numel(yP), p_t, p_m, ...
        'VariableNames', {'phase', 'rate', 'PR_trial', 'mean_g1', 'sem_g1', 'n_g1', ...
        'mean_g2', 'sem_g2', 'n_g2', 'p_ttest2', 'p_ranksum'});
    if height(st) < 1
        st = row;
    else
        st = [st; row];
    end
end
end

function s = pr_step08_sem(x)
s = std(x, 0, 'omitnan') / max(1, sqrt(numel(x(isfinite(x)))));
end

function sl = pr_step08_slopes(S, ycol, ptag, rlabel)
% One slope per mouse: polyfit(Trial, y, 1) on that mouse's mean y per trial
G = S.(ycol);
tr = double(S.Trial);
mk = cellstr(string(S.mouse_key));
gm = cellstr(string(S.Group));
mice = unique(mk, 'stable');
sl = table();
for k = 1:numel(mice)
    r = strcmp(mk, mice{k});
    tlist = unique(tr(r));
    mmu = nan(size(tlist));
    for t = 1:numel(tlist)
        f = r & (tr == tlist(t)) & isfinite(G);
        mmu(t) = mean(G(f), 'omitnan');
    end
    v = isfinite(mmu) & (tlist(:) <= 80);
    if sum(v) < 3
        p1 = [NaN, NaN];
    else
        p1 = polyfit(tlist(v), mmu(v), 1);
    end
    ix = find(r, 1, 'first');
    gg = gm{ix};
    row = table(mice(k), {gg}, {ptag}, {char(rlabel)}, p1(1), p1(2), ...
        'VariableNames', {'mouse_key', 'Group', 'phase', 'rate', 'slope', 'intercept'});
    if height(sl) < 1
        sl = row;
    else
        sl = [sl; row];
    end
end
end

function pr_step08_figure_slope_bars(SL, outDir)
if isempty(SL) || height(SL) < 1
    return
end
rates = unique(SL.rate, 'stable');
phases = unique(SL.phase, 'stable');
for ri = 1:numel(rates)
    for pi = 1:numel(phases)
        s = SL(strcmpi(SL.rate, rates{ri}) & strcmp(SL.phase, phases{pi}), :);
        s.Group = cellstr(string(s.Group));
        if height(s) < 4
            continue
        end
        gU = unique(s.Group, 'stable');
        a = s.slope(strcmpi(s.Group, 'Active') | cellfun(@(x) contains(lower(x), 'act'), s.Group));
        b = s.slope(strcmpi(s.Group, 'Passive') | cellfun(@(x) contains(lower(x), 'pass'), s.Group));
        a = a(isfinite(a));
        b = b(isfinite(b));
        if (numel(a) < 1 || numel(b) < 1) && numel(gU) >= 2
            a = s.slope(strcmp(s.Group, gU{1}));
            b = s.slope(strcmp(s.Group, gU{2}));
            a = a(isfinite(a));
            b = b(isfinite(b));
            l1 = gU{1};
            l2 = gU{2};
        else
            l1 = 'g1';
            l2 = 'g2';
        end
        if numel(a) < 1 || numel(b) < 1
            continue
        end
        f = figure('Color', 'w', 'Position', [100, 100, 500, 480]);
        hold on
        boxplot([a; b], [ones(size(a)); 2*ones(size(b))], 'Labels', {{l1, l2}});
        ylabel('Slope of lick rate vs PR trial (Hz / step)');
        title(sprintf('Slopes %s / %s', rates{ri}, phases{pi}));
        grid on
        try
            p = ranksum(a, b);
            text(0.1, 0.95, sprintf('ranksum p=%.4g (n1=%d n2=%d)', p, numel(a), numel(b)), 'Units', 'normalized');
        end
        fn = fullfile(outDir, sprintf('slope_boxplot_%s_%s.png', strrep(rates{ri}, '_', ''), phases{pi}));
        try
            saveas(f, fn);
        end
        close(f);
    end
end
end
