function f = pr_step08_rebuild_mousemean_from_enriched(enrichedPath, outDir)
% pr_step08_rebuild_mousemean_from_enriched — build MouseMean_LickRate_by_PRLevel.csv
%   from an existing per_trial_enriched_step08.csv (same hard caps as Step 8).
% Use when Step 8 wrote enriched data but not MouseMean, or to refresh quickly.
%
%   f = pr_step08_rebuild_mousemean_from_enriched()   % default: .../step08_pr_trajectory_downstream/output_full/...
%   f = pr_step08_rebuild_mousemean_from_enriched('D:\data\per_trial_enriched_step08.csv', 'D:\data\output_full')
%
% The strings 'C:\path\to\...' in docs are placeholders — use your real paths or call with
% no arguments (after addpath to step08) so this script’s output_full is used.
%
% outputDir: folder to receive MouseMean_LickRate_by_PRLevel.csv (not the file name).

if nargin < 1 || isempty(enrichedPath)
    here = fileparts(mfilename('fullpath'));
    enrichedPath = fullfile(here, 'output_full', 'per_trial_enriched_step08.csv');
end
if exist(enrichedPath, 'file') ~= 2
    error([ 'per_trial enriched table not found:\n  %s\n\n', ...
        'If you used example text like C:\\path\\to\\..., replace it with a real .csv path, or run ', ...
        'with no arguments:\n  f = pr_step08_rebuild_mousemean_from_enriched\n', ...
        'after: addpath(''...\\step08_pr_trajectory_downstream'')\n' ], enrichedPath);
end
if nargin < 2 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'output_full');
end
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

PR_LICK_MIN_TTR_S = 0.02;
PR_LICK_MAX_HZ = 50;
PR_LICK_MAX_HZ_BTWN = 22;

TT = readtable(enrichedPath);
if ismember('time_to_reward', TT.Properties.VariableNames) && ismember('lick_count', TT.Properties.VariableNames)
    ttr = max(double(TT.time_to_reward), PR_LICK_MIN_TTR_S);
    w = double(TT.lick_count) ./ ttr;
elseif ismember('lick_rate', TT.Properties.VariableNames) || ismember('WithinTrialLickRate', TT.Properties.VariableNames)
    if ismember('lick_rate', TT.Properties.VariableNames)
        w = double(TT.lick_rate);
    else
        w = double(TT.WithinTrialLickRate);
    end
else
    error('enriched table needs lick_count+time_to_reward, or lick_rate / WithinTrialLickRate');
end
w = w(:);
w(isfinite(w) & (w < 0 | w > PR_LICK_MAX_HZ)) = NaN;
if ismember('BetweenTrialLickRate', TT.Properties.VariableNames)
    b = double(TT.BetweenTrialLickRate(:));
    b(isfinite(b) & (b < 0 | b > PR_LICK_MAX_HZ_BTWN)) = NaN;
else
    b = nan(height(TT), 1);
end
U = table;
U.mouse_key = cellstr(string(TT.mouse_key));
U.Trial = double(TT.Trial);
U.WithinTrialLickRate = w;
U.BetweenTrialLickRate = b;
if ismember('Group', TT.Properties.VariableNames)
    g = TT.Group;
    if iscell(g)
        U.Group = g;
    else
        U.Group = cellstr(string(g));
    end
else
    U.Group = repmat({'All'}, height(U), 1);
end
if ismember('Period', TT.Properties.VariableNames)
    U.phase = cellstr(string(TT.Period));
    ux = ismember(U.phase, {'<undef>', '<excluded>'});
    U.phase(ux) = {'(drop)'};
else
    U.phase = repmat({'All'}, height(U), 1);
end
pr_step08_write_mousemean_pr_table(U, outDir);
f = fullfile(outDir, 'MouseMean_LickRate_by_PRLevel.csv');
end
