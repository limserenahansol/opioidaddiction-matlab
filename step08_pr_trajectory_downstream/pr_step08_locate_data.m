function p = pr_step08_locate_data(overrideRunDir)
% p = pr_step08_locate_data — paths to latest longitudinal run for Step 8
% p.runDir, p.trial_csv (lick_MEGA), p.long_csv (ALL_mice_longitudinal)

if nargin < 1, overrideRunDir = []; end
if ~isempty(overrideRunDir) && ischar(overrideRunDir) && exist(overrideRunDir, 'dir')
    runDir = char(overrideRunDir);
else
    rootTry = 'K:\addiction_concate_Dec_2025\longitudinal_outputs';
    if ~exist(rootTry, 'dir')
        here = fileparts(mfilename('fullpath'));
        cand = here;
        for up = 1:8
            pth = fullfile(cand, 'longitudinal_outputs');
            if exist(pth, 'dir')
                rootTry = pth;
                break
            end
            cand = fileparts(cand);
        end
    end
    D = dir(fullfile(rootTry, 'run_*'));
    if isempty(D)
        error('No run_* under %s', rootTry);
    end
    [~, ix] = max([D.datenum]);
    runDir = fullfile(D(ix).folder, D(ix).name);
end
p = struct;
p.runDir = runDir;
p.lick_mega = fullfile(runDir, 'figs', 'lick_MEGA');
p.trial_csv = fullfile(p.lick_mega, 'per_trial_features.csv');
p.long_csv = fullfile(runDir, 'ALL_mice_longitudinal.csv');
end
