function run_lick_MEGA_followup2_new()
% FOLLOW-UP 2 (REVISED): visualize MEGA CSVs with:
% - new cohort mapping (cage/color/sex/active-passive + pair ID)
% - period bins: pre/during/post/withdrawal/reexposure
% - optional "clean" version excluding transition days (D4, D6, D11, D14, D17)

%% ================= USER OPTIONS =========================
% Transition days you said are less reliable:
TRANSITION_DAYS = [4, 6, 11, 14, 17];

% Period definitions you specified:
% D1-2 habituation (ignored => <undef>)
% Pre: D3-5
% During: D6-10
% Post: D11-13
% Withdrawal: D14-16
% Reexposure: D17-18

%% ================= locate latest run ====================
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
if ~exist(rootTry,'dir')
    here = pwd; cand = here;
    for up=1:5
        p = fullfile(cand,'longitudinal_outputs');
        if exist(p,'dir'), rootTry=p; break; end
        cand = fileparts(cand);
    end
end
D = dir(fullfile(rootTry,'run_*'));
assert(~isempty(D),'No run_* under %s', rootTry);
[~,ix] = max([D.datenum]);
runDir = fullfile(D(ix).folder, D(ix).name);
outDir = fullfile(runDir,'figs','lick_MEGA_followup2');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ================= cohort mapping =======================
% ASSUMPTION: your "f_s" for 6100red means "female, passive-like"
% (because you described 6100 red as passive in the pairing).
% If "s" means something else, change ActPass for that row below.

cohort = buildCohortTable();  % table with mouse_key, Cage, Color, Sex, ActPass, PairID

%% ================= load MEGA outputs ====================
S_sess  = safeRead(fullfile(runDir,'figs','lick_MEGA','per_session_features.csv'));
S_trial = safeRead(fullfile(runDir,'figs','lick_MEGA','per_trial_features.csv'));
S_spec  = safeRead(fullfile(runDir,'figs','lick_MEGA','C_rhythm_spectral_metrics.csv'));
S_ac    = safeRead(fullfile(runDir,'figs','lick_MEGA','C_rhythm_autocorr_metrics.csv'));
Pstart  = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','psth_start_group.csv'));
Prew    = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','psth_reward_group.csv')); %#ok<NASGU> (kept for completeness)
Drpr    = safeRead(fullfile(runDir,'figs','lick_MEGA','D_psth','reward_ramp_pause.csv'));

%% ================= sanitize + attach cohort + bins ======
S_sess  = ensureCols(S_sess, cohort, TRANSITION_DAYS);
S_trial = ensureCols(S_trial, cohort, TRANSITION_DAYS);
S_spec  = ensureCols(S_spec, cohort, TRANSITION_DAYS);
S_ac    = ensureCols(S_ac, cohort, TRANSITION_DAYS);

%% ================= session rates =========================
boxAndJitterByGroup(S_sess,'licks_per_min','Licks / min', outDir,'A_session_licks_per_min.png');

% Period plots (include all days in period)
boxAndJitterByPeriod(S_sess,'licks_per_min','Licks / min', outDir, ...
    'A_session_licks_per_min_by_PERIOD_including_transition.png', 'period');

% Period plots (excluding transition days: D4, D6, D11, D14, D17)
boxAndJitterByPeriod(S_sess,'licks_per_min','Licks / min', outDir, ...
    'A_session_licks_per_min_by_PERIOD_EXCLUDING_transition.png', 'period_clean');

spaghettiPerMousePeriod(S_sess,'licks_per_min', outDir, ...
    'S_session_spaghetti_licks_per_min.png', TRANSITION_DAYS);

%% ================= rhythm summaries ======================
if ~isempty(S_spec)
    quickViolin(S_spec,'peak_hz','Peak Hz',outDir,'C_spec_peak_hz.png');
    quickViolin(S_spec,'band_power_3_12','Band power 3–12 Hz',outDir,'C_spec_band312.png');
    quickViolin(S_spec,'spectral_entropy','Spectral entropy',outDir,'C_spec_entropy.png');
end
if ~isempty(S_ac)
    quickViolin(S_ac,'ac_peak_lag','AC peak lag (s)',outDir,'C_ac_lag.png');
    quickViolin(S_ac,'ac_peak_height','AC peak height',outDir,'C_ac_height.png');
end

%% ================= PSTH (group) ==========================
% Robust to “wide” CSVs, no reliance on a 't' column
% Start-aligned only, split panels + difference (Passive − Active)
if ~isempty(Pstart)
    plotPSTHtableSplitDiff(Pstart, ...
        'Start-aligned lick rate (0–5 s)', ...
        fullfile(outDir,'PSTH_start_group_SPLIT_and_DIFFERENCE.png'));
end

% (skip reward-aligned per your request; keep copy of ramp/pause)
if ~isempty(Drpr)
    writetable(Drpr, fullfile(outDir,'reward_ramp_pause_copy.csv'));
end

fprintf('\nFOLLOWUP2 (REVISED) complete. Outputs in:\n  %s\n', outDir);
end

%% ================= cohort builder =========================
function cohort = buildCohortTable()
% mouse_key format assumed in your pipeline: "CAGE_color" (e.g., "6100_red")
% If your CSV uses "6100red" (no underscore), we normalize in attachCohort().

rows = {
% mouse_key_raw   Sex  ActPass  PairID
'6100_red',      'F', 'Passive', 1   % user wrote f_s; treated as passive (see note above)
'6100_orange',   'F', 'Passive', 1
'6100_black',    'F', 'Active',  1

'0911_red',      'F', 'Active',  2
'0911_orange',   'F', 'Passive', 2

'0911_black',    'F', 'Passive', 3
'0911_white',    'F', 'Active',  3

'0910_red',      'M', 'Passive', 4
'0910_orange',   'M', 'Passive', 4
'0910_black',    'M', 'Active',  4

'6099_red',      'M', 'Passive', 5
'6099_orange',   'M', 'Active',  5

'6099_black',    'M', 'Active',  6
'6099_white',    'M', 'Passive', 6  % user: died after day13 (script will just plot what exists)
};

mouse_key = string(rows(:,1));
Sex      = categorical(string(rows(:,2)), {'F','M'});
ActPass  = categorical(string(rows(:,3)), {'Active','Passive'});
PairID   = double(cell2mat(rows(:,4)));

% parse Cage + Color from mouse_key
Cage  = strings(size(mouse_key));
Color = strings(size(mouse_key));
for i=1:numel(mouse_key)
    parts = split(mouse_key(i), "_");
    if numel(parts)>=2
        Cage(i)  = parts(1);
        Color(i) = parts(2);
    else
        % fallback (should not happen here)
        Cage(i)  = mouse_key(i);
        Color(i) = "";
    end
end

cohort = table(mouse_key, Cage, Color, Sex, ActPass, PairID);
end

%% ================= helpers ================================
function T = safeRead(p)
if exist(p,'file')>0
    T = readtable(p,'VariableNamingRule','preserve');
else
    T = table();
end
end

function S = ensureCols(S, cohort, transitionDays)
if isempty(S), return; end

% --- normalize mouse_key to string if present
if ismember('mouse_key',S.Properties.VariableNames) && ~isstring(S.mouse_key)
    S.mouse_key = string(S.mouse_key);
end

% --- attach cohort metadata where possible
S = attachCohort(S, cohort);

% --- normalize group label columns (support older labels)
% GroupMouse preferred: Active/Passive
if ismember('GroupMouse', S.Properties.VariableNames)
    S.GroupMouse = normalizeGroupLabels(S.GroupMouse);
elseif ismember('ActPass', S.Properties.VariableNames)
    S.GroupMouse = categorical(string(S.ActPass), {'Active','Passive'});
end

% --- period bins
% period: includes all days in each period
% period_clean: same but transition days set to "<exclude>" so we can drop them
if ismember('day_index', S.Properties.VariableNames)
    d = double(S.day_index);
    S.period       = dayToPeriod(d, false, transitionDays);
    S.period_clean = dayToPeriod(d, true,  transitionDays);
else
    S.period       = categorical(repmat("<undef>",height(S),1), periodCats(false));
    S.period_clean = categorical(repmat("<undef>",height(S),1), periodCats(true));
end
end

function cats = periodCats(withExclude)
base = ["Pre_D3-5","During_D6-10","Post_D11-13","Withdrawal_D14-16","Reexposure_D17-18","<undef>"];
if withExclude
    cats = [base(1:5), "<exclude>", base(6)];
else
    cats = base;
end
cats = cellstr(cats);
end

function P = dayToPeriod(day_index, excludeTransitions, transitionDays)
lab = strings(size(day_index)); lab(:) = "<undef>";
d = double(day_index);

lab(d>=3  & d<=5 )  = "Pre_D3-5";
lab(d>=6  & d<=10 ) = "During_D6-10";
lab(d>=11 & d<=13 ) = "Post_D11-13";
lab(d>=14 & d<=16 ) = "Withdrawal_D14-16";
lab(d>=17 & d<=18 ) = "Reexposure_D17-18";

if excludeTransitions
    lab(ismember(d, transitionDays)) = "<exclude>";
    P = categorical(lab, periodCats(true));
else
    P = categorical(lab, periodCats(false));
end
end

function G = normalizeGroupLabels(x)
% Accepts string/categorical/numeric; maps old labels to Active/Passive if present
sx = string(x);

% common legacy mapping:
sx(strcmpi(sx,"ActiveOnly"))  = "Active";
sx(strcmpi(sx,"HadPassive"))  = "Passive";
sx(strcmpi(sx,"PassiveOnly")) = "Passive";

G = categorical(sx, {'Active','Passive'});
end

function S = attachCohort(S, cohort)
if isempty(S) || isempty(cohort), return; end
if ~ismember('mouse_key', S.Properties.VariableNames), return; end

% Normalize keys:
% - allow "6100red" or "6100_red" etc.
mk = normalizeMouseKey(S.mouse_key);
ck = normalizeMouseKey(cohort.mouse_key);

% match
[tf,loc] = ismember(mk, ck);

% add columns (fill missing with <undef>)
S.Cage   = strings(height(S),1);
S.Color  = strings(height(S),1);
S.Sex    = categorical(repmat("U",height(S),1), {'F','M','U'});
S.ActPass= categorical(repmat("U",height(S),1), {'Active','Passive','U'});
S.PairID = nan(height(S),1);

S.Cage(tf)    = cohort.Cage(loc(tf));
S.Color(tf)   = cohort.Color(loc(tf));
S.Sex(tf)     = addcats(cohort.Sex(loc(tf)), "U"); % ensure U exists
S.Sex(tf)     = cohort.Sex(loc(tf));
S.ActPass(tf) = addcats(cohort.ActPass(loc(tf)), "U");
S.ActPass(tf) = cohort.ActPass(loc(tf));
S.PairID(tf)  = cohort.PairID(loc(tf));

% If we found a match and GroupMouse exists, prefer ActPass as ground-truth
if any(tf)
    S.GroupMouse = categorical(string(S.ActPass), {'Active','Passive','U'});
    S.GroupMouse = removecats(S.GroupMouse, 'U'); % keep only Active/Passive if possible
end
end

function mk = normalizeMouseKey(mkIn)
s = string(mkIn);
s = strtrim(lower(s));
s = regexprep(s,'[\-]+','_');      % dashes -> underscore
s = regexprep(s,'\s+','_');        % spaces -> underscore

% If someone used "6100red" without underscore, try to insert after first 4 digits
% Pattern: 4 digits + letters
s2 = s;
for i=1:numel(s)
    tok = regexp(s(i), '^(\d{4})([a-z]+)$', 'tokens','once');
    if ~isempty(tok)
        s2(i) = tok{1} + "_" + tok{2};
    end
end
mk = s2;
end

function s = mouselabel(mk)
% "6100_red" -> "6100 red"
s = regexprep(string(mk), '[_\-]+', ' ');
end

%% ================= plotting ==============================
function boxAndJitterByGroup(S, col, ylab, outDir, fname)
if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember('GroupMouse',S.Properties.VariableNames)
    return
end
cats = categories(S.GroupMouse);
cats = cats(~cellfun(@isempty,cats));
if isempty(cats), return, end

fig = figure('Color','w','Position',[90 90 740 520]); hold on
for i=1:numel(cats)
    r = (S.GroupMouse == cats{i});
    x = i + 0.35*(rand(sum(r),1)-0.5);
    y = double(S.(col)(r));
    scatter(x, y, 18, 'filled','MarkerFaceAlpha',0.45);
    if ismember('mouse_key',S.Properties.VariableNames)
        lab = mouselabel(S.mouse_key(r));
        text(x+0.02, y, lab, 'FontSize',6, 'Interpreter','none');
    end
end
xlim([0.5 numel(cats)+0.5]);
set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
ylabel(ylab); title([ylab ' by Group']); grid on; box on
savepng(fig, fullfile(outDir, fname)); close(fig);
end

function boxAndJitterByPeriod(S, col, ylab, outDir, fname, periodVar)
if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember(periodVar,S.Properties.VariableNames)
    return
end

P = S.(periodVar);

% drop <undef> and (for clean) drop <exclude>
keep = (P ~= "<undef>");
if any(categories(P) == "<exclude")
    keep = keep & (P ~= "<exclude>");
end
S2 = S(keep,:); P2 = S2.(periodVar);

binsOrder = {'Pre_D3-5','During_D6-10','Post_D11-13','Withdrawal_D14-16','Reexposure_D17-18'};
bins = binsOrder(ismember(binsOrder, categories(P2)));
if isempty(bins), return, end

fig = figure('Color','w','Position',[90 90 980 520]); hold on
for i=1:numel(bins)
    r = (P2 == bins{i});
    x = i + 0.35*(rand(sum(r),1)-0.5);
    y = double(S2.(col)(r));
    scatter(x, y, 18, 'filled','MarkerFaceAlpha',0.45);
    if ismember('mouse_key',S2.Properties.VariableNames)
        lab = mouselabel(S2.mouse_key(r));
        text(x+0.02, y, lab, 'FontSize',6, 'Interpreter','none');
    end
end
xlim([0.5 numel(bins)+0.5]);
set(gca,'XTick',1:numel(bins),'XTickLabel',bins);
ylabel(ylab); title([ylab ' by Period (' periodVar ')']); grid on; box on
savepng(fig, fullfile(outDir, fname)); close(fig);
end

function spaghettiPerMousePeriod(S, col, outDir, fname, transitionDays)
if isempty(S) || ~ismember(col,S.Properties.VariableNames) || ~ismember('mouse_key',S.Properties.VariableNames)
    return
end

fig = figure('Color','w','Position',[90 90 980 540]); hold on
ms = unique(S.mouse_key,'stable');

for i=1:numel(ms)
    r = (S.mouse_key==ms(i));
    x = double(S.day_index(r));
    y = double(S.(col)(r));
    [x,ord] = sort(x); y=y(ord);
    plot(x,y,'-o','LineWidth',1,'MarkerSize',3);

    % label at max
    if ~isempty(y)
        [~,mxi] = max(y);
        text(x(mxi)+0.06, y(mxi), mouselabel(ms(i)), 'FontSize',7, 'Interpreter','none');
    end
end

% period boundary markers (visual aid)
xline(3,'k:'); xline(6,'k:'); xline(11,'k:'); xline(14,'k:'); xline(17,'k:');

% transition day markers (visual aid)
for d = transitionDays(:)'
    xline(d,'k--','LineWidth',0.6);
end

xlabel('Day'); ylabel(strrep(col,'_','\_'));
title(['Spaghetti per mouse: ' strrep(col,'_','\_')]);
grid on; box on
savepng(fig, fullfile(outDir, fname)); close(fig);
end

function quickViolin(T, col, ylab, outDir, fname)
if isempty(T) || ~ismember(col,T.Properties.VariableNames), return, end
if ~ismember('GroupMouse',T.Properties.VariableNames), return, end

cats = categories(T.GroupMouse);
cats = cats(~cellfun(@isempty,cats));
if isempty(cats), return, end

fig = figure('Color','w','Position',[100 100 700 420]); hold on
for i=1:numel(cats)
    y = double(T.(col)(T.GroupMouse==cats{i}));
    if isempty(y), continue, end
    [yk,xk] = kde1(y,200);
    yk = yk/max(yk)*0.35;
    patch([i-yk, fliplr(i+yk)], [xk, fliplr(xk)], 0.92*[1 1 1], 'EdgeColor','none');
    scatter(i + 0.25*(rand(numel(y),1)-0.5), y, 8, 'filled','MarkerFaceAlpha',0.35);
    plot([i-0.18 i+0.18], [median(y) median(y)], 'k-', 'LineWidth',1.8);
end
set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
ylabel(ylab); title(ylab); grid on; box on
savepng(fig, fullfile(outDir,fname)); close(fig);
end

%% ================= PSTH helpers ==========================
function plotPSTHtableSplitDiff(P, ttl, saveTo)
if isempty(P), return; end

% groups
if ismember('Group',P.Properties.VariableNames)
    Graw = normalizeGroupLabels(P.Group);
    G = Graw;
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

C = [0 0.45 0.74; 0.85 0.33 0.10; 0.3 0.3 0.3];
ymax = 1.05*max(cellfun(@(m,s) max(m+s), MU, SE));

fig = figure('Color','w','Position',[100 100 900 650]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% panel 1
nexttile; hold on
i = 1;
fill([T{i}; flipud(T{i})], [MU{i}-SE{i}; flipud(MU{i}+SE{i})], C(i,:), ...
    'FaceAlpha',0.18,'EdgeColor','none');
plot(T{i}, MU{i}, 'LineWidth',2.2, 'Color', C(i,:));
xline(0,'k:'); ylim([0 ymax]);
ylabel('Lick rate (Hz)'); title(sprintf('%s — %s', ttl, groups{i}));
grid on; box on

% panel 2 + difference
if numel(groups) >= 2
    j = 2;
    nexttile; hold on
    fill([T{j}; flipud(T{j})], [MU{j}-SE{j}; flipud(MU{j}+SE{j})], C(j,:), ...
        'FaceAlpha',0.18,'EdgeColor','none');
    plot(T{j}, MU{j}, 'LineWidth',2.2, 'Color', C(j,:));
    xline(0,'k:'); ylim([0 ymax]);
    ylabel('Lick rate (Hz)'); title(groups{j});
    grid on; box on

    % difference (Passive − Active) if labels are Active/Passive
    % otherwise (group2 - group1)
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

    % choose order for Passive - Active if possible
    if any(strcmp(groups,'Active')) && any(strcmp(groups,'Passive'))
        muA = MU{strcmp(groups,'Active')};
        seA = SE{strcmp(groups,'Active')};
        muP = MU{strcmp(groups,'Passive')};
        seP = SE{strcmp(groups,'Passive')};
        % align if needed
        tA = T{strcmp(groups,'Active')};
        tP = T{strcmp(groups,'Passive')};
        if numel(tA)~=numel(tP) || any(abs(tA-tP)>1e-12)
            tmin = max(min(tA), min(tP));
            tmax = min(max(tA), max(tP));
            tC   = linspace(tmin,tmax,min(numel(tA),numel(tP)))';
            muA  = interp1(tA, muA, tC, 'linear','extrap');
            seA  = interp1(tA, seA, tC, 'linear','extrap');
            muP  = interp1(tP, muP, tC, 'linear','extrap');
            seP  = interp1(tP, seP, tC, 'linear','extrap');
        else
            tC = tA;
        end
        muD = muP - muA;
        seD = sqrt(seA.^2 + seP.^2);
        diffTitle = 'Difference (Passive − Active)';
    else
        muD = mu2 - mu1;
        seD = sqrt(se1.^2 + se2.^2);
        diffTitle = 'Difference (Group2 − Group1)';
    end

    nexttile; hold on
    fill([tC; flipud(tC)], [muD-seD; flipud(muD+seD)], C(3,:), ...
        'FaceAlpha',0.15,'EdgeColor','none');
    plot(tC, muD, 'Color', C(3,:), 'LineWidth',2.2);
    yline(0,'k:'); xline(0,'k:');
    xlabel('Time from trial start (s)'); ylabel('\Delta rate (Hz)');
    title(diffTitle); grid on; box on
else
    xlabel('Time from trial start (s)');
end

savepng(fig, saveTo); close(fig);
end

function v = parseVec(x)
if isnumeric(x), v = x(:)'; return, end
s = char(string(x));
s = strrep(s,'[',''); s = strrep(s,']',''); s = strrep(s,';',',');
v = str2num(s); %#ok<ST2NM>
v = v(:)';
end

function v = rowFromWide(T, mask, prefix)
cols = T.Properties.VariableNames;
k = startsWith(cols, prefix); cols = cols(k);
if isempty(cols), v = []; return, end
suf = regexp(cols, ['^' regexptranslate('escape',prefix) '(\d+)$'], 'tokens','once');
suf = cellfun(@(c) str2double(c{1}), suf);
[~,ord] = sort(suf); cols = cols(ord);
r = find(mask,1,'first');
v = zeros(1,numel(cols));
for i=1:numel(cols)
    v(i) = double(T.(cols{i})(r));
end
end

function [y,x] = kde1(v, n)
v=v(isfinite(v));
if numel(v)<3
    x=linspace(min(v)-eps,max(v)+eps,20);
    y=ones(size(x));
    return
end
s = std(v);
h = 1.06*s* numel(v)^(-1/5);
x = linspace(min(v)-3*h, max(v)+3*h, n);
y=zeros(size(x));
for i=1:numel(v)
    y = y + exp(-0.5*((x-v(i))/h).^2);
end
y = y/(numel(v)*h*sqrt(2*pi));
end

function savepng(fh, fn)
set(fh,'PaperPositionMode','auto');
try
    exportgraphics(fh, fn, 'Resolution',180);
catch
    print(fh, fn, '-dpng','-r180');
end
end
