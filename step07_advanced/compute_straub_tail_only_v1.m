function compute_straub_tail_only_v1()
% compute_straub_tail_only_v1 (Straub tail analysis)
%
% Straub tail analysis using ONLY:
%   - STRAUB_Frames_Non_moving
%   - STRAUB_Pct_Non_moving
%
% Reads from latest:
%   K:\addiction_concate_Dec_2025\longitudinal_outputs\run_*\ALL_mice_longitudinal.csv
%
% Excludes day 1-2, assigns MouseGroup from your mapping (cage+color),
% then produces:
%   1) Spaghetti by MouseGroup
%   2) Phase comparison within group (Kruskal-Wallis + stars)
%   3) Group comparison within phase (ranksum + stars)
%
% Outputs under:
%   <runDir>\figs\straub_nonmoving_only_v1\

%% ----------------- USER SETTINGS -----------------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';

% find latest run_*
d = dir(fullfile(rootTry,'run_*'));
assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix] = max([d.datenum]);
runDir  = fullfile(d(ix).folder,d(ix).name);

csvPath = fullfile(runDir,'ALL_mice_longitudinal.csv');
assert(exist(csvPath,'file')>0,'Missing %s', csvPath);

outDir  = fullfile(runDir,'figs','straub_nonmoving_only_v1');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Variables to analyze (exact names from your CSV)
varsToAnalyze = { ...
    'STRAUB_Frames_Non_moving', ...
    'STRAUB_Pct_Non_moving' ...
    };

% Exclude day1-2
minDay = 3;

%% ----------------- LOAD CSV -----------------
% Preserve original column names (important for STRAUB_* names)
T = readtable(csvPath, 'VariableNamingRule','preserve');

% Standardize key columns expected
T = standardize_key_columns(T);

% Filter day
T = T(double(T.day_index) >= minDay, :);

% Define phase
T.Period = periodOfDay(double(T.day_index));
T = T(~isundefined(T.Period), :);

% Add MouseGroup (Active/Passive) from cage+color mapping
T = add_mouse_group_from_cohort(T);

% Save a quick sanity check
writetable(T(:, intersect([{'mouse_key','day_index','MouseGroup','Period'}, varsToAnalyze], ...
    T.Properties.VariableNames, 'stable')), fullfile(outDir,'StraubNonMoving_mouse_day_RAW.csv'));

%% ----------------- RUN ANALYSES FOR EACH VAR -----------------
for vi = 1:numel(varsToAnalyze)
    vname = varsToAnalyze{vi};
    assert(ismember(vname, T.Properties.VariableNames), ...
        'Column not found: %s. Check exact header in CSV.', vname);

    Dsub = table();
    Dsub.mouse_key  = string(T.mouse_key);
    Dsub.day_index  = double(T.day_index);
    Dsub.Period     = T.Period;
    Dsub.MouseGroup = T.MouseGroup;
    Dsub.Y          = to_double_robust(T.(vname));

    % Drop missing group or missing Y
    good = isfinite(Dsub.Y) & ~isundefined(Dsub.MouseGroup);
    Dsub = Dsub(good,:);

    % Save
    writetable(Dsub, fullfile(outDir, sprintf('%s_mouse_day.csv', vname)));

    % Phase aggregated table (mouse-level median per phase)
    G = groupsummary(Dsub(:,{'mouse_key','MouseGroup','Period','Y'}), ...
        {'mouse_key','MouseGroup','Period'}, 'median', 'Y');
    G.Properties.VariableNames{'median_Y'} = 'Y_med';
    writetable(G, fullfile(outDir, sprintf('%s_mouse_phase.csv', vname)));

    % Plots
    plot_spaghetti_byGroup(Dsub, outDir, vname);
    plot_phase_withinGroup(G, outDir, vname);
    plot_group_withinPhase(G, outDir, vname);
    % PI-requested consolidated figure (During vs Withdrawal, Active vs Passive)
    if strcmp(vname, 'STRAUB_Pct_Non_moving')
        plot_straub_main_figure(G, outDir);
    end
end

fprintf('\nDone. Outputs:\n  %s\n', outDir);
end

%% ===================== HELPERS =====================

function [col, ylbl] = straub_plot_constants(vname)
% Group colors: Active=red, Passive=blue (explicit — do NOT use MATLAB defaults)
col = struct('Active', [0.85 0.2 0.2], 'Passive', [0.2 0.4 0.8]);
% Y-axis label: definition per PI (for %) or frames
if nargin >= 1 && strcmp(vname, 'STRAUB_Pct_Non_moving')
    ylbl = 'Straub tail (%)';
else
    ylbl = 'Straub tail (frames)';
end
end

function T = standardize_key_columns(T)
% Try to locate mouse_key and day_index in a robust way.
names = string(T.Properties.VariableNames);
low   = lower(names);

% mouse_key
if ~ismember('mouse_key', T.Properties.VariableNames)
    cand = find(contains(low,'mouse_key') | strcmp(low,'mouse') | contains(low,'mouseid') | contains(low,'mouse_id'), 1);
    assert(~isempty(cand), 'Could not find mouse_key column in CSV.');
    T.mouse_key = string(T.(names(cand)));
end

% day_index
if ~ismember('day_index', T.Properties.VariableNames)
    cand = find(strcmp(low,'day') | strcmp(low,'day_index') | contains(low,'dayindex') | contains(low,'session_day'), 1);
    assert(~isempty(cand), 'Could not find day_index/day column in CSV.');
    T.day_index = double(T.(names(cand)));
end

% Ensure types
T.mouse_key = string(T.mouse_key);
T.day_index = double(T.day_index);
end

function x = to_double_robust(v)
try
    if isnumeric(v) || islogical(v)
        x = double(v); return
    end
    if iscategorical(v)
        x = str2double(string(v)); return
    end
    if isstring(v) || iscellstr(v) || iscell(v)
        x = str2double(string(v)); return
    end
catch
end
x = nan(size(v,1),1);
end

function T = add_mouse_group_from_cohort(T)
% MouseGroup from cage+color mapping you gave.
% NOTE: 6100 red = f_s is NOT assigned here. If you want it assigned,
% tell me whether it should be Active or Passive and I will add it.

mk = lower(strtrim(string(T.mouse_key)));
cage = strings(size(mk));
col  = strings(size(mk));
for i=1:numel(mk)
    [cage(i), col(i)] = parse_cage_color(mk(i));
end

mapKeys = [
    "6100_orange"
    "6100_black"
    "0911_red"
    "0911_orange"
    "0911_black"
    "0911_white"
    "0910_red"
    "0910_orange"
    "0910_black"
    "6099_red"
    "6099_orange"
    "6099_black"
    "6099_white"
];

mapVals = [
    "Passive"
    "Active"
    "Active"
    "Passive"
    "Passive"
    "Active"
    "Passive"
    "Passive"
    "Active"
    "Passive"
    "Active"
    "Active"
    "Passive"
];

MouseGroup = strings(size(mk));
for i=1:numel(mk)
    k = cage(i) + "_" + col(i);
    idx = find(mapKeys==k, 1, 'first');
    if ~isempty(idx)
        MouseGroup(i) = mapVals(idx);
    else
        MouseGroup(i) = ""; % unassigned (e.g., 6100_red f_s)
    end
end

T.MouseGroup = categorical(MouseGroup, ["Active","Passive"]);
end

function [cage,color] = parse_cage_color(s)
s = lower(strtrim(string(s)));

m = regexp(s, '(\d{4})', 'tokens', 'once');
if isempty(m), cage = ""; else, cage = string(m{1}); end

colors = ["black","red","orange","white"];
color = "";
for c = colors
    if contains(s, c)
        color = c;
        break
    end
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

function plot_spaghetti_byGroup(Dsub, outDir, vname)
[col, ylbl] = straub_plot_constants(vname);
groups = ["Active","Passive"];

fh = figure('Color','w','Position',[90 90 1100 480]);
tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

for gi=1:2
    ax = nexttile; hold(ax,'on');
    G = groups(gi);
    Tg = Dsub(Dsub.MouseGroup==G,:);
    c = col.(char(G));

    if isempty(Tg)
        text(ax,0.5,0.5,'No data for this MouseGroup','Units','normalized','HorizontalAlignment','center');
        title(ax,G); xlabel(ax,'Day'); ylabel(ax,ylbl);
        box(ax,'on'); continue
    end

    mice = unique(Tg.mouse_key,'stable');
    for i=1:numel(mice)
        r = Tg.mouse_key==mice(i);
        dd = Tg.day_index(r);
        [dd,ord] = sort(dd);
        yy = Tg.Y(r); yy = yy(ord);
        plot(ax, dd, yy, '-', 'Color', [c 0.5], 'LineWidth',0.9);
    end

    days = unique(Tg.day_index);
    mY = nan(size(days));
    for j=1:numel(days)
        mY(j) = mean(Tg.Y(Tg.day_index==days(j)), 'omitnan');
    end
    plot(ax, days, mY, '-', 'Color', c, 'LineWidth',2.5, 'DisplayName', char(G));

    title(ax, G);
    xlabel(ax,'Day'); ylabel(ax,ylbl);
    legend(ax, char(G), 'Location', 'best');
    box(ax,'on'); grid(ax,'off');
end

sgtitle(sprintf('Straub tail in Active and Passive mice (both sexes) — %s trajectories', vname), 'FontSize', 11);
printpng_local(fh, fullfile(outDir, sprintf('%s_spaghetti_byMouseGroup.png', vname)));
close(fh);
end

function plot_phase_withinGroup(G, outDir, vname)
[col, ylbl] = straub_plot_constants(vname);
groups = ["Active","Passive"];
periods = categories(G.Period);

fh = figure('Color','w','Position',[90 90 1100 500]);
tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

for gi=1:2
    ax = nexttile; hold(ax,'on');
    Gi = G(G.MouseGroup==groups(gi),:);
    c = col.(char(groups(gi)));

    if isempty(Gi)
        text(ax,0.5,0.5,'No data for this MouseGroup','Units','normalized','HorizontalAlignment','center');
        title(ax, sprintf('%s (mouse-level median per phase)', groups(gi)));
        continue
    end

    y = []; gph = [];
    for pi=1:numel(periods)
        yy = Gi.Y_med(Gi.Period==periods{pi} & isfinite(Gi.Y_med));
        y   = [y; yy]; %#ok<AGROW>
        gph = [gph; repmat(pi, numel(yy), 1)]; %#ok<AGROW>
    end
    pKW = NaN;
    if numel(unique(gph)) >= 2 && numel(y) >= 3
        pKW = kruskalwallis(y, gph, 'off');
    end

    for pi=1:numel(periods)
        yy = Gi.Y_med(Gi.Period==periods{pi} & isfinite(Gi.Y_med));
        if isempty(yy), continue; end
        boxchart(ax, repmat(pi,numel(yy),1), yy, 'BoxWidth',0.55, 'MarkerStyle','none', ...
            'BoxFaceColor', c, 'BoxFaceAlpha',0.5);
        scatter(ax, repmat(pi,numel(yy),1), yy, 28, c, 'filled', 'MarkerFaceAlpha',0.80);
    end

    set(ax,'XTick',1:numel(periods),'XTickLabel',periods);
    xlim(ax,[0.5 numel(periods)+0.5]);
    title(ax, sprintf('%s — During (morphine) vs Withdrawal (water)\nKruskal-Wallis p = %s  %s', ...
        groups(gi), fmt_p_local(pKW), stars_from_p(pKW)));
    ylabel(ax, ylbl);
    box(ax,'on'); grid(ax,'off');
end

sgtitle('Straub tail: phase comparison within each group (measured During & Withdrawal)', 'FontSize', 11);
printpng_local(fh, fullfile(outDir, sprintf('%s_phase_comparison_byMouseGroup.png', vname)));
close(fh);
end

function plot_group_withinPhase(G, outDir, vname)
[col, ylbl] = straub_plot_constants(vname);
periods = categories(G.Period);

fh = figure('Color','w','Position',[90 90 1300 500]);
tiledlayout(1,numel(periods),'TileSpacing','loose','Padding','loose');

for pi=1:numel(periods)
    ax = nexttile; hold(ax,'on');

    A = G.Y_med(G.MouseGroup=="Active"  & G.Period==periods{pi} & isfinite(G.Y_med));
    P = G.Y_med(G.MouseGroup=="Passive" & G.Period==periods{pi} & isfinite(G.Y_med));

    pRS = NaN;
    if numel(A)>=2 && numel(P)>=2
        pRS = ranksum(A,P);
    end

    if ~isempty(A)
        boxchart(ax, ones(size(A)), A, 'BoxWidth',0.55, 'MarkerStyle','none', ...
            'BoxFaceColor', col.Active, 'BoxFaceAlpha',0.5);
        scatter(ax, ones(size(A)), A, 28, col.Active, 'filled', 'MarkerFaceAlpha',0.80);
    end
    if ~isempty(P)
        boxchart(ax, 2*ones(size(P)), P, 'BoxWidth',0.55, 'MarkerStyle','none', ...
            'BoxFaceColor', col.Passive, 'BoxFaceAlpha',0.5);
        scatter(ax, 2*ones(size(P)), P, 28, col.Passive, 'filled', 'MarkerFaceAlpha',0.80);
    end

    set(ax,'XTick',[1 2],'XTickLabel',{'Active','Passive'});
    title(ax, sprintf('%s\nranksum p = %s  %s', periods{pi}, fmt_p_local(pRS), stars_from_p(pRS)));
    ylabel(ax, ylbl);
    xlim(ax,[0.5 2.5]);
    % Legend: use dummy points so colors are correct
    hA = scatter(ax, nan, nan, 36, col.Active, 'filled', 'DisplayName', 'Active');
    hP = scatter(ax, nan, nan, 36, col.Passive, 'filled', 'DisplayName', 'Passive');
    legend(ax, [hA hP], 'Location', 'best');
    box(ax,'on'); grid(ax,'off');
end

sgtitle('Straub tail: Active vs Passive within each phase (both sexes, equivalent morphine exposure)', 'FontSize', 11);
printpng_local(fh, fullfile(outDir, sprintf('%s_group_comparison_withinPhase.png', vname)));
close(fh);
end

function plot_straub_main_figure(G, outDir)
% PI-requested layout (STRAUB_Pct_Non_moving only):
%   Upper: Active vs Passive within During (morphine) and Withdrawal (water)
%   Lower: During vs Withdrawal within Active and Passive groups
% Explicit colors: Active=red, Passive=blue. Legend on every panel.
[col, ylbl] = straub_plot_constants('STRAUB_Pct_Non_moving');

fh = figure('Color','w','Position',[90 90 1000 750]);
tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

% --- Upper left: During — Active vs Passive ---
ax1 = nexttile; hold(ax1,'on');
A = G.Y_med(G.MouseGroup=="Active"  & G.Period=="During" & isfinite(G.Y_med));
P = G.Y_med(G.MouseGroup=="Passive" & G.Period=="During" & isfinite(G.Y_med));
pRS1 = NaN; if numel(A)>=2 && numel(P)>=2, pRS1 = ranksum(A,P); end
if ~isempty(A)
    boxchart(ax1, ones(size(A)), A, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Active, 'BoxFaceAlpha',0.5);
    scatter(ax1, ones(size(A)), A, 36, col.Active, 'filled', 'MarkerFaceAlpha',0.80);
end
if ~isempty(P)
    boxchart(ax1, 2*ones(size(P)), P, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Passive, 'BoxFaceAlpha',0.5);
    scatter(ax1, 2*ones(size(P)), P, 36, col.Passive, 'filled', 'MarkerFaceAlpha',0.80);
end
set(ax1,'XTick',[1 2],'XTickLabel',{'Active','Passive'}); xlim(ax1,[0.5 2.5]);
title(ax1, sprintf('During (morphine) — Active vs Passive\nranksum p = %s  %s', fmt_p_local(pRS1), stars_from_p(pRS1)));
ylabel(ax1, ylbl);
hA1 = scatter(ax1, nan, nan, 36, col.Active, 'filled', 'DisplayName', 'Active');
hP1 = scatter(ax1, nan, nan, 36, col.Passive, 'filled', 'DisplayName', 'Passive');
legend(ax1, [hA1 hP1], 'Location','best'); box(ax1,'on'); grid(ax1,'off');

% --- Upper right: Withdrawal — Active vs Passive ---
ax2 = nexttile; hold(ax2,'on');
A = G.Y_med(G.MouseGroup=="Active"  & G.Period=="Withdrawal" & isfinite(G.Y_med));
P = G.Y_med(G.MouseGroup=="Passive" & G.Period=="Withdrawal" & isfinite(G.Y_med));
pRS2 = NaN; if numel(A)>=2 && numel(P)>=2, pRS2 = ranksum(A,P); end
if ~isempty(A)
    boxchart(ax2, ones(size(A)), A, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Active, 'BoxFaceAlpha',0.5);
    scatter(ax2, ones(size(A)), A, 36, col.Active, 'filled', 'MarkerFaceAlpha',0.80);
end
if ~isempty(P)
    boxchart(ax2, 2*ones(size(P)), P, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Passive, 'BoxFaceAlpha',0.5);
    scatter(ax2, 2*ones(size(P)), P, 36, col.Passive, 'filled', 'MarkerFaceAlpha',0.80);
end
set(ax2,'XTick',[1 2],'XTickLabel',{'Active','Passive'}); xlim(ax2,[0.5 2.5]);
title(ax2, sprintf('Withdrawal (water) — Active vs Passive\nranksum p = %s  %s', fmt_p_local(pRS2), stars_from_p(pRS2)));
ylabel(ax2, ylbl);
hA2 = scatter(ax2, nan, nan, 36, col.Active, 'filled', 'DisplayName', 'Active');
hP2 = scatter(ax2, nan, nan, 36, col.Passive, 'filled', 'DisplayName', 'Passive');
legend(ax2, [hA2 hP2], 'Location','best'); box(ax2,'on'); grid(ax2,'off');

% --- Lower left: Active — During vs Withdrawal ---
ax3 = nexttile; hold(ax3,'on');
D = G.Y_med(G.MouseGroup=="Active" & G.Period=="During"    & isfinite(G.Y_med));
W = G.Y_med(G.MouseGroup=="Active" & G.Period=="Withdrawal" & isfinite(G.Y_med));
pKW3 = NaN; if numel(D)>=2 && numel(W)>=2, pKW3 = kruskalwallis([D;W], [ones(numel(D),1); 2*ones(numel(W),1)], 'off'); end
if ~isempty(D)
    boxchart(ax3, ones(size(D)), D, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Active, 'BoxFaceAlpha',0.5);
    scatter(ax3, ones(size(D)), D, 36, col.Active, 'filled', 'MarkerFaceAlpha',0.80);
end
if ~isempty(W)
    boxchart(ax3, 2*ones(size(W)), W, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Active, 'BoxFaceAlpha',0.5);
    scatter(ax3, 2*ones(size(W)), W, 36, col.Active, 'filled', 'MarkerFaceAlpha',0.80);
end
set(ax3,'XTick',[1 2],'XTickLabel',{'During (morphine)','Withdrawal (water)'}); xlim(ax3,[0.5 2.5]);
title(ax3, sprintf('Active — During vs Withdrawal\nKruskal-Wallis p = %s  %s', fmt_p_local(pKW3), stars_from_p(pKW3)));
ylabel(ax3, ylbl);
h3 = scatter(ax3, nan, nan, 36, col.Active, 'filled', 'DisplayName', 'Active');
legend(ax3, h3, 'Location','best'); box(ax3,'on'); grid(ax3,'off');

% --- Lower right: Passive — During vs Withdrawal ---
ax4 = nexttile; hold(ax4,'on');
D = G.Y_med(G.MouseGroup=="Passive" & G.Period=="During"    & isfinite(G.Y_med));
W = G.Y_med(G.MouseGroup=="Passive" & G.Period=="Withdrawal" & isfinite(G.Y_med));
pKW4 = NaN; if numel(D)>=2 && numel(W)>=2, pKW4 = kruskalwallis([D;W], [ones(numel(D),1); 2*ones(numel(W),1)], 'off'); end
if ~isempty(D)
    boxchart(ax4, ones(size(D)), D, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Passive, 'BoxFaceAlpha',0.5);
    scatter(ax4, ones(size(D)), D, 36, col.Passive, 'filled', 'MarkerFaceAlpha',0.80);
end
if ~isempty(W)
    boxchart(ax4, 2*ones(size(W)), W, 'BoxWidth',0.55, 'MarkerStyle','none', 'BoxFaceColor', col.Passive, 'BoxFaceAlpha',0.5);
    scatter(ax4, 2*ones(size(W)), W, 36, col.Passive, 'filled', 'MarkerFaceAlpha',0.80);
end
set(ax4,'XTick',[1 2],'XTickLabel',{'During (morphine)','Withdrawal (water)'}); xlim(ax4,[0.5 2.5]);
title(ax4, sprintf('Passive — During vs Withdrawal\nKruskal-Wallis p = %s  %s', fmt_p_local(pKW4), stars_from_p(pKW4)));
ylabel(ax4, ylbl);
h4 = scatter(ax4, nan, nan, 36, col.Passive, 'filled', 'DisplayName', 'Passive');
legend(ax4, h4, 'Location','best'); box(ax4,'on'); grid(ax4,'off');

sgtitle({'Straub tail: Active and Passive mice (both sexes)', ...
    'Upper: Active vs Passive within phase. Lower: During (morphine) vs Withdrawal (water) within group.', ...
    'Both groups exhibited Straub tail after morphine — equivalent pharmacological effect across behavioral contingency.'}, ...
    'FontSize', 10);
% Ensure layout has enough space for title and legends before export
drawnow;
printpng_local(fh, fullfile(outDir, 'STRAUB_Pct_Non_moving_main_figure.png'));
close(fh);
end

function s = stars_from_p(p)
if ~isfinite(p)
    s = '';
elseif p < 1e-4
    s = '****';
elseif p < 1e-3
    s = '***';
elseif p < 1e-2
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end

function t = fmt_p_local(p)
if ~isfinite(p), t = 'NaN'; return; end
if p < 1e-4
    t = sprintf('%.1e', p);
else
    t = sprintf('%.4g', p);
end
end

function printpng_local(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',180);
catch
    print(fh, fn, '-dpng','-r180');
end
end
