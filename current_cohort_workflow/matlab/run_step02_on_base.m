function outRoot = run_step02_on_base(base, step02File)
% run_step02_on_base  Run the ORIGINAL Step 02 script on any BASE folder without editing it.
%
%   outRoot = run_step02_on_base('D:\cohort\step02_base')
%   outRoot = run_step02_on_base(base, '<repo>\step02_build_longitudinal\Longitudinal_final_....m')
%
% The original Step 02 (step02_build_longitudinal/Longitudinal_final_trialrequire_HOTTST_passive_final_handle_
% nomatchingstrabu.m) is a script whose first lines are "clear; clc;" and a hard-coded BASE = '...'. This wrapper
% reads that file, replaces ONLY the BASE assignment line with the folder given here, writes the result to a temporary
% script and runs it. Every other line is used unchanged, so the output (BASE\longitudinal_outputs\run_###\
% ALL_mice_longitudinal.csv, session_inventory.csv, issues.txt, summary.txt, per_mouse\) is what Step 02 itself makes.
%
% BASE layout Step 02 expects:  BASE\dayN\<cage>\<mouse>\concat_out_<stamp>\combined_pupil_digital.csv
%                               plus the session's *.jsonl in the mouse folder (or the concat_out, cage or day folder)
%                               plus, when they exist, ManualScoringResults*tst*/*hot*/*straub*.mat and the
%                               tail-immersion file in the mouse folder.
% Use prepare_step02_base.py to build that layout from the pupil/TST package output and the task logs.

if nargin < 2 || isempty(step02File)
    here = fileparts(mfilename('fullpath'));
    step02File = fullfile(here, '..', '..', 'step02_build_longitudinal', ...
        'Longitudinal_final_trialrequire_HOTTST_passive_final_handle_nomatchingstrabu.m');
end
assert(isfile(step02File), 'Step 02 script not found: %s', step02File);
assert(isfolder(base), 'BASE folder not found: %s', base);

src = fileread(step02File);
pat = '(?m)^\s*BASE\s*=\s*''[^'']*''\s*;';
nHit = numel(regexp(src, pat, 'match'));
assert(nHit == 1, 'expected exactly one "BASE = ''...'';" line in %s, found %d', step02File, nHit);
baseQ = strrep(char(base), '''', '''''');
src = regexprep(src, pat, sprintf('BASE = ''%s'';  %% set by run_step02_on_base', strrep(baseQ, '\', '\\')), 'once');

tmpDir = tempname;   % (no Java needed; newer MATLAB releases may run without a Java runtime)
mkdir(tmpDir);
tmp = fullfile(tmpDir, 'step02_tmp_run.m');
fid = fopen(tmp, 'w', 'n', 'UTF-8'); fwrite(fid, src, 'char'); fclose(fid);
fprintf('Step 02 (original logic) on BASE = %s\n', base);
before = dir(fullfile(base, 'longitudinal_outputs', 'run_*'));
run_isolated(tmp);   % the script starts with "clear": run it in its own function workspace
after = dir(fullfile(base, 'longitudinal_outputs', 'run_*'));
new = setdiff({after.name}, {before.name});
if isempty(new), outRoot = ''; else, outRoot = fullfile(base, 'longitudinal_outputs', new{end}); end
fprintf('Step 02 output: %s\n', outRoot);
try, rmdir(tmpDir, 's'); catch, end
end

function run_isolated(scriptFile)
% "clear" inside the Step 02 script only clears this function's workspace
run(scriptFile);
end
