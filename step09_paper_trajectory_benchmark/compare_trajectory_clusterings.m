function S = compare_trajectory_clusterings(pathEfa, pathPr, varargin)
% compare_trajectory_clusterings — overlap between two mouse-level clusterings
% Usage:
%   S = compare_trajectory_clusterings(efaCsv, prCsv, 'outDir', outD)
% Each CSV should contain mouse_key and a cluster column named one of:
%   'cluster', 'trajectory_cluster', 'efa_cluster', 'Cluster', 'clust'
%
% Outputs: adjusted Rand index, NMI, confusion matrix CSV/PNG if outDir given

p = inputParser;
addParameter(p, 'outDir', '', @ischar);
addParameter(p, 'efaLabelVar', '', @ischar);
addParameter(p, 'prLabelVar', '', @ischar);
parse(p, varargin{:});
outD = p.Results.outDir;
efaCol = p.Results.efaLabelVar;
prCol = p.Results.prLabelVar;

E = readtable(pathEfa);
P = readtable(pathPr);
% Match Step-7 key rules (slashes/spaces->underscore, color lowercased) or innerjoin gives n=0
E.mouse_key = cellstr(step09_normalize_mouse_key(E.mouse_key));
P.mouse_key = cellstr(step09_normalize_mouse_key(P.mouse_key));
vE = pickClustVar(E, efaCol);
vP = pickClustVar(P, prCol);
E.cl = relabel_categorical_to_int(vE);
P.cl = relabel_categorical_to_int(vP);
E2 = E(:, {'mouse_key', 'cl'});
E2.Properties.VariableNames{'cl'} = 'cl_efa';
P2 = P(:, {'mouse_key', 'cl'});
P2.Properties.VariableNames{'cl'} = 'cl_pr';
[~, ia] = unique(E2.mouse_key, 'stable');
E2 = E2(ia, :);
[~, ib] = unique(P2.mouse_key, 'stable');
P2 = P2(ib, :);
J = innerjoin(E2, P2, 'Keys', 'mouse_key');
if isempty(J) || height(J) < 1
    pr_warn_key_mismatch(E2, P2);
    error('step09:KeyMismatch', ...
        ['No overlapping mouse_key after join (n=0). ', ...
         'Check that EFA and PR tables use the same run/cohort, and see printed key samples above. ']);
end

c1 = J.cl_efa;
c2 = J.cl_pr;
n = numel(c1);
C = confusionmat(c1, c2);
[ari, ~] = pr_adjusted_rand_index(c1, c2);
nmi = pr_normalized_mutual_info(c1, c2);
S = struct();
S.n_matched_mice = n;
S.adj_Rand = ari;
S.NMI = nmi;
S.confusion = C;
fprintf('Matched mice: n=%d  |  ARI=%.4f  NMI=%.4f\n', n, ari, nmi);

if ~isempty(outD) && exist(outD, 'dir') ~= 7
    mkdir(outD);
end
if ~isempty(outD)
    writetable(array2table(C, 'VariableNames', ...
        arrayfun(@(x) sprintf('pr_%d', x), 1:size(C,2), 'uni', 0), ...
        'RowNames', arrayfun(@(x) sprintf('efa_%d', x), 1:size(C,1), 'uni', 0)), ...
        fullfile(outD, 'confusion_cluster_overlap.csv'), 'WriteRowNames', true);
    f = figure('Color', 'w', 'Position', [100, 100, 520, 420]);
    imagesc(C);
    colormap(flipud(gray));
    colorbar
    xlabel('PR-trajectory cluster (or B)');
    ylabel('EFA cluster (or A)');
    title(sprintf('Confusion  n=%d  ARI=%.3f NMI=%.3f', n, ari, nmi));
    try
        saveas(f, fullfile(outD, 'confusion_cluster_overlap.png'));
        close(f);
    end
    mtab = table(S.n_matched_mice, S.adj_Rand, S.NMI, 'VariableNames', {'n_mice', 'ARI', 'NMI'});
    writetable(mtab, fullfile(outD, 'compare_metrics.csv'));
end
end

function pr_warn_key_mismatch(E2, P2)
uE = unique(cellstr(string(E2.mouse_key)), 'stable');
uP = unique(cellstr(string(P2.mouse_key)), 'stable');
fprintf(2, '\n--- mouse_key mismatch (join n=0) ---\n');
ne = min(8, numel(uE));
np = min(8, numel(uP));
fprintf(2, 'EFA file: n=%d unique keys. First 8: %s\n', numel(uE), strjoin(uE(1:ne), ' | '));
fprintf(2, 'PR  file: n=%d unique keys. First 8: %s\n', numel(uP), strjoin(uP(1:np), ' | '));
nBoth = numel(intersect(uE, uP));
fprintf(2, 'Intersection size (after normalize): %d\n', nBoth);
fprintf(2, '-------------------------------------\n\n');
end

function v = pickClustVar(T, userCol)
cands = {'trajectory_cluster', 'cluster', 'efa_cluster', 'Cluster', 'clust'};
if ~isempty(userCol) && ismember(userCol, T.Properties.VariableNames)
    v = T.(userCol);
    return
end
for i = 1:numel(cands)
    if ismember(cands{i}, T.Properties.VariableNames)
        v = T.(cands{i});
        return
    end
end
error('No cluster column found. Expected one of: %s', strjoin(cands, ', '));
end

function x = relabel_categorical_to_int(v)
if iscell(v)
    [u, ~, x] = unique(v, 'stable');
    return
end
if iscategorical(v)
    x = double(v);
    return
end
x = double(v);
end

function [ARI, s] = pr_adjusted_rand_index(c1, c2)
% Adjusted Rand (Hubert & Arabie)
c1 = c1(:);
c2 = c2(:);
N = numel(c1);
k1 = max(c1);
k2 = max(c2);
M = zeros(k1, k2);
for i = 1:N
    M(c1(i), c2(i)) = M(c1(i), c2(i)) + 1;
end
a = sum(M, 2);
b = sum(M, 1);
b = b(:)'; % row
nc2 = N * (N - 1) / 2;
if nc2 < 1
    ARI = NaN; s = [];
    return
end
sumM2 = 0;
for g = 1:k1
    for h = 1:k2
        sumM2 = sumM2 + pr_comb2(M(g, h));
    end
end
sumA2 = 0;
for g = 1:k1
    sumA2 = sumA2 + pr_comb2(a(g));
end
sumB2 = 0;
for h = 1:k2
    sumB2 = sumB2 + pr_comb2(b(h));
end
expected = (sumA2 * sumB2) / max(nc2, 1);
maxIndex = 0.5 * (sumA2 + sumB2);
num = sumM2 - expected;
den = maxIndex - expected;
if abs(den) < eps
    ARI = 0;
else
    ARI = num / den;
end
s = struct('M', M);
end

function nmi = pr_normalized_mutual_info(c1, c2)
c1 = c1(:);
c2 = c2(:);
N = numel(c1);
k1 = max(c1);
k2 = max(c2);
P = zeros(k1, k2);
for i = 1:N
    P(c1(i), c2(i)) = P(c1(i), c2(i)) + 1;
end
P = P / N;
P1 = sum(P, 2);
P2 = sum(P, 1);
H1 = -sum(P1(P1 > 0) .* log(P1(P1 > 0) + 1e-20));
H2 = -sum(P2(P2 > 0) .* log(P2(P2 > 0) + 1e-20));
mi = 0;
for a = 1:k1
    for b2 = 1:k2
        p = P(a, b2);
        if p > 0 && (P1(a) * P2(b2) > 0)
            mi = mi + p * log(p / (P1(a) * P2(b2)) + 1e-20);
        end
    end
end
den = 0.5 * (H1 + H2);
if den < 1e-20
    nmi = 0;
else
    nmi = mi / den;
end
end

function y = pr_comb2(n)
% binomial n choose 2, 0 if n<2
n = double(n(1));
if n < 2
    y = 0;
else
    y = n * (n - 1) / 2;
end
end
