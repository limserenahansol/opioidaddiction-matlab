function analyze_addiction_score_efa_mousefit()
% Mouse-level EFA/PCA projection (separate pipeline)
% Produces outputs under figs/addiction_score_efa_mousefit/

%% ---------- locate latest run + load ----------
rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
d = dir(fullfile(rootTry,'run_*')); assert(~isempty(d),'No run_* under %s',rootTry);
[~,ix]   = max([d.datenum]);
runDir   = fullfile(d(ix).folder,d(ix).name);
cacheMat = fullfile(runDir, 'S_D_cache.mat');
assert(exist(cacheMat,'file')>0, 'Missing %s', cacheMat);
L = load(cacheMat,'D');
D = L.D;

% Exclude day 1-2
D = D(D.day_index >= 3, :);

% Group assignment (explicit mouse map)
mouseInfo = build_mouse_map();
map = containers.Map(cellstr(mouseInfo.mouse_key), cellstr(string(mouseInfo.MouseGroup)));
mg = strings(height(D),1);
for r=1:height(D)
    k = char(string(D.mouse_key(r)));
    if isKey(map,k)
        mg(r) = string(map(k));
    else
        mg(r) = "Unknown";
    end
end
D.Group = categorical(mg, {'Active','Passive','Unknown'});

outDir = fullfile(runDir,'figs','addiction_score_efa_mousefit');
if ~exist(outDir,'dir'), mkdir(outDir); end
noteFile = fullfile(outDir,'mousefit_note.txt');
fidNote = fopen(noteFile,'w');
fprintf(fidNote,'Mousefit EFA/PCA notes\n');

%% ---------- features ----------
exclude = ["day_index","session_idx","SessionMinutes","isPassive"];
exclude = [exclude, "mouse_key","day_name","Group","Period","Session_Paradigm"];
vn = string(D.Properties.VariableNames);
numMask = varfun(@isnumeric, D, 'OutputFormat','uniform');
useVars = setdiff(vn(numMask), exclude, 'stable');

%% ---------- mouse-level matrix ----------
mk = unique(D.mouse_key,'stable');
Xmouse = nan(numel(mk), numel(useVars));
for j=1:numel(useVars)
    v = useVars(j);
    for i=1:numel(mk)
        r = (D.mouse_key==mk(i));
        Xmouse(i,j) = median(double(D.(v)(r)),'omitnan');
    end
end

missFeat = mean(~isfinite(Xmouse),1);
keepFeat = missFeat <= 0.70;
Xmouse = Xmouse(:, keepFeat);
useVars = useVars(keepFeat);
if size(Xmouse,2) < 4
    v = var(Xmouse,0,1,'omitnan');
    [~,ord] = sort(v,'descend','MissingPlacement','last');
    ord = ord(1:min(6, numel(ord)));
    Xmouse = Xmouse(:,ord);
    useVars = useVars(ord);
end

missMouse = mean(~isfinite(Xmouse),2);
keepMouse = missMouse <= 0.70;
Xmouse = Xmouse(keepMouse,:);
mk = mk(keepMouse);

for j=1:size(Xmouse,2)
    col = Xmouse(:,j);
    med = median(col,'omitnan');
    col(~isfinite(col)) = med;
    Xmouse(:,j) = col;
end

[Xz, muFeat, sdFeat] = zscore_with_store(Xmouse);

% Drop near-constant features
v = var(Xz,0,1);
keep = v > 1e-6;
Xz = Xz(:,keep);
muFeat = muFeat(keep);
sdFeat = sdFeat(keep);
useVars = useVars(keep);

%% ---------- try EFA, fallback PCA ----------
targetExpl = 0.85;
nf = 0;
Lambda = []; Psi = []; Fmouse = [];

maxFactors = min([12, size(Xz,2)-1, size(Xz,1)-1]);
for k=1:maxFactors
    [ok, Lm, Ps, msg] = tryFactoranMouse(Xz, k);
    if ~ok
        fprintf(fidNote,'  nf=%d failed: %s\n', k, msg);
        continue
    end
    commonExpl = sum(1 - Ps) / numel(Ps);
    fprintf(fidNote,'  nf=%d OK -> common variance %.1f%%\n', k, 100*commonExpl);
    Lambda = Lm; Psi = Ps; nf = k;
    if commonExpl >= targetExpl, break; end
end

usePCA = false;
if nf == 0
    usePCA = true;
    fprintf(fidNote,'EFA failed; using PCA fallback.\n');
    [coeff, score, ~, ~, expl] = pca(Xz);
    cumExpl = cumsum(expl)/100;
    nf = find(cumExpl>=targetExpl,1,'first');
    if isempty(nf), nf = min(3, size(Xz,2)-1); end
    Fmouse = score(:,1:nf);
else
    Fmouse = score_factors_regression(Xz, Lambda, Psi);
end

[FmouseZ, muF, sdF] = zscore_with_store(Fmouse);
AI_mouse_raw = mean(FmouseZ, 2, 'omitnan');
AI_mouse = rescale01(AI_mouse_raw);
mouseTbl = table(mk, AI_mouse, AI_mouse_raw, 'VariableNames', {'mouse_key','AI_mouse','AI_mouse_z'});

%% ---------- project to mouse×day ----------
Dsub = D(:, intersect([{'mouse_key','day_index','day_name','Group'}, cellstr(useVars)], ...
    D.Properties.VariableNames, 'stable'));
Xday = Dsub{:, useVars};
for j=1:size(Xday,2)
    col = Xday(:,j);
    col(~isfinite(col)) = muFeat(j);
    Xday(:,j) = col;
end
XdayZ = (Xday - muFeat) ./ sdFeat;
XdayZ(~isfinite(XdayZ)) = 0;

if usePCA
    Fday = XdayZ * coeff(:,1:nf);
else
    Fday = score_factors_regression(XdayZ, Lambda, Psi);
end
FdayZ = (Fday - muF) ./ sdF;
FdayZ(~isfinite(FdayZ)) = 0;
AI_day_raw = mean(FdayZ, 2, 'omitnan');
AI_day = rescale01(AI_day_raw);

out = Dsub(:, {'mouse_key','day_index','day_name','Group'});
out.AddictionScore = AI_day;
out.AddictionScore_z = AI_day_raw;
for k=1:nf
    out.("Factor"+k) = Fday(:,k);
end

writetable(out, fullfile(outDir,'addiction_score_mousefit_by_mouse_day.csv'));
writetable(mouseTbl, fullfile(outDir,'addiction_score_mousefit_by_mouse.csv'));

fprintf(fidNote,'Success: nf=%d, nMice=%d, nFeat=%d, method=%s\n', ...
    nf, size(Xmouse,1), numel(useVars), ternary(usePCA,'PCA','EFA'));
fclose(fidNote);

writeMousefitWorkflow(outDir, useVars, nf, usePCA);

% Basic plots
plotAddictionSpaghetti(out, outDir);
out.Period = periodOfDay(double(out.day_index));
plotAddictionByPhase(out, outDir);

fprintf('Mousefit outputs saved to:\n  %s\n', outDir);
end

%% ---- helpers ----
function [Z, mu, sd] = zscore_with_store(X)
X = double(X);
mu = median(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(sd==0 | ~isfinite(sd)) = 1;
Z = (X - mu) ./ sd;
Z(~isfinite(Z)) = 0;
end

function y = rescale01(x)
x = double(x);
good = isfinite(x);
if ~any(good)
    y = 0.5*ones(size(x)); return
end
lo = min(x(good));
hi = max(x(good));
if hi==lo
    y = 0.5*ones(size(x));
else
    y = (x - lo) / (hi - lo);
end
y(~isfinite(y)) = 0.5;
end

function F = score_factors_regression(Xz, Lambda, Psi)
Xz = double(Xz);
Psi = double(Psi(:));
k = size(Lambda,2);
invPsi = diag(1 ./ max(Psi, 1e-6));
A = (Lambda' * invPsi * Lambda + eye(k));
W = (A \ (Lambda' * invPsi)); % k x p
F = (W * Xz')'; % n x k
F(~isfinite(F)) = 0;
end

function [ok, Lambda, Psi, msg] = tryFactoranMouse(Xz, nf)
ok = false; Lambda=[]; Psi=[]; msg='';
try
    [Lambda, Psi] = factoran(Xz, nf, 'rotate','varimax', 'scores','regression');
    ok = true; return
catch ME
    msg = ME.message;
end
try
    Xj = Xz + 1e-6*randn(size(Xz));
    [Lambda, Psi] = factoran(Xj, nf, 'rotate','varimax', 'scores','regression');
    ok = true; return
catch ME
    msg = [msg ' | jitter: ' ME.message];
end
try
    R = corrcoef(Xz);
    R = (R + R')/2;
    R = R + 1e-3*eye(size(R,1));
    [Lambda, Psi] = factoran(R, nf, 'nobs', size(Xz,1), 'rotate','varimax', 'scores','regression');
    ok = true; return
catch ME
    msg = [msg ' | R+ridge: ' ME.message];
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

function plotAddictionSpaghetti(out, outDir)
fig = figure('Color','w','Position',[80 80 1000 520]); hold on
groups = ["Passive","Active"];
cols = struct('Passive',[0.2 0.4 0.9],'Active',[0.9 0.3 0.2]);
for gi=1:numel(groups)
    g = groups(gi);
    mk = unique(out.mouse_key(string(out.Group)==g),'stable');
    for i=1:numel(mk)
        r = (out.mouse_key==mk(i)) & (string(out.Group)==g);
        [d,ord] = sort(out.day_index(r));
        y = double(out.AddictionScore(r)); y=y(ord);
        plot(d,y,'-o','LineWidth',0.8,'MarkerSize',3,'Color',cols.(g));
    end
end
xlabel('Day'); ylabel('Addiction score (0-1)');
title('Mousefit Addiction score over days');
grid off; box off
exportgraphics(fig, fullfile(outDir,'mousefit_addiction_spaghetti.png'), 'Resolution', 200);
close(fig);
end

function plotAddictionByPhase(out, outDir)
fig = figure('Color','w','Position',[80 80 900 520]); hold on
periods = {'Pre','During','Post','Withdrawal','Re-exposure'};
groups = ["Passive","Active"];
off = [-0.18 0.18];
cols = struct('Passive',[0.2 0.4 0.9],'Active',[0.9 0.3 0.2]);
for gi=1:numel(groups)
    g = groups(gi);
    for pi=1:numel(periods)
        r = string(out.Group)==g & out.Period==periods{pi};
        if ~any(r), continue; end
        x = pi + off(gi);
        b = boxchart(repmat(x,sum(r),1), out.AddictionScore(r), 'BoxWidth',0.25, 'MarkerStyle','none','BoxFaceAlpha',0.25);
        b.BoxFaceColor = cols.(g);
    end
end
set(gca,'XTick',1:numel(periods),'XTickLabel',periods);
xlabel('Period'); ylabel('Addiction score (0-1)');
title('Mousefit Addiction score by phase');
grid off; box off
exportgraphics(fig, fullfile(outDir,'mousefit_addiction_by_phase.png'), 'Resolution', 200);
close(fig);
end

function y = ternary(cond,a,b)
if cond, y=a; else, y=b; end
end

function writeMousefitWorkflow(outDir, useVars, nf, usePCA)
fid = fopen(fullfile(outDir,'mousefit_workflow.txt'),'w');
fprintf(fid,'Mousefit EFA/PCA workflow (v3-style)\\n');
fprintf(fid,'-----------------------------------\\n');
fprintf(fid,'Inputs: S_D_cache.mat (table D), day>=3 only\\n');
fprintf(fid,'Mouse grouping: explicit mapping via build_mouse_map()\\n\\n');

fprintf(fid,'1) Mouse-level aggregation\\n');
fprintf(fid,'   - One row per mouse, median across days per feature\\n');
fprintf(fid,'   - Drop features with >70%% missing (fallback: top-variance features)\\n');
fprintf(fid,'   - Drop mice with >70%% missing\\n');
fprintf(fid,'   - Impute remaining missing with feature medians\\n');
fprintf(fid,'   - Z-score (median/std) and drop near-constant features\\n\\n');

fprintf(fid,'2) Factor model selection\\n');
fprintf(fid,'   - Target common variance >= 85%%\\n');
fprintf(fid,'   - EFA varimax attempted; fallback to PCA if EFA fails\\n');
fprintf(fid,'   - Selected nf = %d (method=%s)\\n\\n', nf, ternary(usePCA,'PCA','EFA'));

fprintf(fid,'3) Mouse-level AI\\n');
fprintf(fid,'   - Compute factor scores\\n');
fprintf(fid,'   - AI_mouse_raw = mean(z-scored factors)\\n');
fprintf(fid,'   - AI_mouse = rescale to 0-1\\n\\n');

fprintf(fid,'4) Project to mouse×day\\n');
fprintf(fid,'   - Impute missing day features with mouse-level medians\\n');
fprintf(fid,'   - Standardize with mouse-level mu/sd\\n');
fprintf(fid,'   - Project to factor space (EFA loadings or PCA coeff)\\n');
fprintf(fid,'   - AI_day_raw = mean(z-scored factor scores)\\n');
fprintf(fid,'   - AI_day = rescale to 0-1\\n\\n');

fprintf(fid,'Features used (%d):\\n', numel(useVars));
for i=1:numel(useVars)
    fprintf(fid,'  - %s\\n', useVars(i));
end
fclose(fid);
end

function mouseInfo = build_mouse_map()
keys  = ["6100_red","6100_orange","6100_black", ...
         "0911_red","0911_orange","0911_black","0911_white", ...
         "0910_red","0910_orange","0910_black", ...
         "6099_red","6099_orange","6099_black","6099_white"];

Sex = ["F","F","F",  "F","F","F","F",  "M","M","M",  "M","M","M","M"];
Group = ["Passive","Passive","Active", ...
         "Active","Passive","Passive","Active", ...
         "Passive","Passive","Active", ...
         "Passive","Active","Active","Passive"];
PairID = [1,1,1, 2,2,3,3, 4,4,4, 5,5,6,6]';

mouseInfo = table(string(keys(:)), categorical(Sex(:),{'F','M'}), ...
                  categorical(Group(:),{'Active','Passive'}), PairID, ...
                  'VariableNames',{'mouse_key','Sex','MouseGroup','PairID'});
end
