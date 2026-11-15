function [RUN_ROOT, info] = pupil_align_ttl_sync(DIGITAL_FILE, PUPIL_FILE, OUT_PARENT, FRAMES_CSV, RUN_NAME)
% pupil_align_ttl_sync  The package's pupil_align_ttl.m (original Step 01a alignment) with ONE change: which camera
% pulse belongs to which video frame is taken from sync_sessions.py (<session>_frames.csv), not from the END.
%   [RUN_ROOT, info] = pupil_align_ttl_sync(digitalCsv, pupilCsv, outParent, framesCsv)
%   [RUN_ROOT, info] = pupil_align_ttl_sync(..., runName)   % output folder name (default concat_out_<stamp>_syncaligned)
%
% Why (current cohort, 2026-11-15): when a video has more or fewer frames than camera TTL edges, the original matches
% them at the END (front rows NaN, or the first camera edges dropped). On the current rig the video starts with the
% first camera pulse (the recording stopped early, or extra frames were written at the end), so the END match shifts
% licks and rewards against the pupil by the difference (several seconds). A motion-vs-lick cross-correlation of the
% video confirmed lag 0 for the start match in the sessions checked. sync_sessions.py decides the match from the video
% file creation time (or a manual video_offsets.csv) and writes it to <session>_frames.csv.
%
% Everything else is the original code, unchanged: reading digital.csv and the pupil CSV, the TTL edge strategy, the
% lick / injector state read at each camera edge, the columns of combined_pupil_digital.csv, lick_onsets.csv,
% injector_onsets.csv, the figure and summary.txt (plus a few added lines that say how the rows were matched).
% The output folder is concat_out_<stamp>_syncaligned under OUT_PARENT. The package's own folder is not touched.
% Use it only for sessions the user approved (pupil_tst_final\alignment_fixes.csv, tools\fix_concat_alignment.py).

% -------- Load digital CSV --------
try
    Tdig = readtable(DIGITAL_FILE, 'VariableNamingRule','preserve');
catch
    Tdig = readtable(DIGITAL_FILE);
end

timeCol   = pickVar(Tdig, {'Time [s]','Time[s]','Time_s','Time'});
cameraCol = pickVar(Tdig, {'camera','Camera','Channel 6','Channel_6'});
lickCol   = pickVar(Tdig, {'lick','Lick','Channel 7'});

% injector is OPTIONAL; try to find, else placeholder
injCol = ''; hasInjector = false;
try
    % CURRENT-COHORT EDIT 2026-11-14: the current rig's Saleae export names Channel 4 'reward' (added last, so old files are unchanged)
    injCol = pickVar(Tdig, {'injector','Injector','injection','Inj','Channel 4','Channel_4','reward','Reward'});
    hasInjector = true;
catch
end

D = table;
D.Time_s = numify(Tdig.(timeCol));
D.camera = logical(numify(Tdig.(cameraCol)));
D.lick   = logical(numify(Tdig.(lickCol)));
if hasInjector
    D.injector = logical(numify(Tdig.(injCol)));
else
    D.injector = nan(height(Tdig),1);   % for export consistency
end

% -------- Raw TTL rising edges (may be doubled) --------
cam = D.camera(:) > 0.5;
idxRiseRaw = find(cam & [false; ~cam(1:end-1)]);
if isempty(idxRiseRaw)
    idxRiseRaw = find(cam); % fallback if edges couldn't be formed
end
tRiseRaw = D.Time_s(idxRiseRaw);
isiRaw   = diff(tRiseRaw);
ttlHzRaw = 1/median(isiRaw, 'omitnan');

% -------- Load pupil results --------
try
    Tpup = readtable(PUPIL_FILE, 'VariableNamingRule','preserve');
catch
    Tpup = readtable(PUPIL_FILE);
end

tsCol    = pickVar(Tpup, {'Timestamp_seconds','Timestamp','time_s','time'});
diamCol  = pickVar(Tpup, {'SmoothedDiameter_px','Diameter_pixels','Diameter'});
cxCol    = pickVar(Tpup, {'CenterX_px','CenterX'});
cyCol    = pickVar(Tpup, {'CenterY_px','CenterY'});
scoreCol = pickVar(Tpup, {'UsedScoreFallback','UsedScore','Score'});

p_ts    = numify(Tpup.(tsCol));
p_diam  = numify(Tpup.(diamCol));
p_cx    = numify(Tpup.(cxCol));
p_cy    = numify(Tpup.(cyCol));
p_score = numify(Tpup.(scoreCol));

nPupil  = numel(p_ts);
nTTLraw = numel(idxRiseRaw);
ratioRaw = nTTLraw / max(1, nPupil);

% -------- Decide strategy using ONLY TTL & pupil --------
if isempty(isiRaw)
    estFramePeriod = NaN;
else
    q75 = prctile(isiRaw, 75);
    longISI = isiRaw(isiRaw >= q75);                   % keep longer spacings
    estFramePeriod = median(longISI, 'omitnan');       % robust real frame period
end

tolRel = 0.08;   % 8% tolerance for "match"
strategy = 'raw_direct_match';
idxFrames = idxRiseRaw;

if abs(ratioRaw - 1) <= tolRel
    strategy = 'raw_direct_match';
elseif abs(ratioRaw - 2) <= 0.25
    if ~isnan(estFramePeriod) && estFramePeriod > 0
        minISI = 0.7 * estFramePeriod;                % refractory window
    else
        minISI = 0.7 * (2*median(isiRaw,'omitnan'));  % fallback
    end
    keep = false(size(idxRiseRaw));
    last_t = -inf;
    for k = 1:numel(idxRiseRaw)
        if tRiseRaw(k) - last_t >= minISI
            keep(k) = true; last_t = tRiseRaw(k);
        end
    end
    idxFrames = idxRiseRaw(keep);
    if numel(idxFrames) > 1.4*nPupil && abs(numel(idxFrames)/nPupil - 2) < 0.3
        idxFrames = idxFrames(1:2:end);
        strategy = 'debounce_then_halve';
    else
        strategy = 'debounce_time_guard';
    end
else
    strategy = 'raw_mismatch_end_align';
end

nFrames = numel(idxFrames);

% what the original pupil_align_ttl does with these counts (reported only)
if nPupil == nFrames
    origMode = 'start_1to1'; origOffset = 0;
elseif nPupil > nFrames
    origMode = 'end_align_camera_to_pupil'; origOffset = nPupil - nFrames;
else
    origMode = 'truncate_camera_to_last_pupil'; origOffset = nFrames - nPupil;
end

% -------- Build mapping from sync_sessions.py's <session>_frames.csv (replaces the END alignment) --------
N = nPupil;
CamTime_s     = nan(N,1);
Camera_TTL    = nan(N,1);
Lick_TTL      = nan(N,1);
Injector_TTL  = nan(N,1);

% (all checks come BEFORE the output folder is made, so a failed check leaves nothing behind)
Tfr = readtable(FRAMES_CSV, 'VariableNamingRule','preserve');
frm = numify(Tfr.frame);           % 0-based video frame of each session camera pulse (NaN: pulse without video)
tfr = numify(Tfr.t_saleae_s);      % time of that pulse on the Saleae clock (= Time [s] of digital.csv)
use = ~isnan(frm);
if ~any(use)
    error('%s: no camera pulse has a video frame (the frame column is empty)', FRAMES_CSV);
end
if any(frm(use) < 0 | frm(use) >= N | frm(use) ~= round(frm(use)))
    error('%s: frame numbers outside 0..%d (the video has %d frames)', FRAMES_CSV, N-1, N);
end
rowOf = frm(use) + 1;              % pupil row (1-based) = video frame + 1
tOf   = tfr(use);
if numel(unique(rowOf)) ~= numel(rowOf) || any(diff(rowOf) <= 0)
    error('%s: frame numbers are not strictly increasing', FRAMES_CSV);
end
tSel = D.Time_s(idxFrames);        % the camera edges chosen by the strategy above (unchanged)
if numel(tSel) > 1
    pos = interp1(tSel, (1:numel(tSel))', tOf, 'nearest', 'extrap');
else
    pos = ones(size(tOf));
end
edgeDt = tSel(pos) - tOf;
TOL_S = 0.0005;                    % 0.5 ms: both times come from the same digital.csv (frames.csv keeps 6 decimals)
if any(abs(edgeDt) > TOL_S) || any(diff(pos) <= 0)
    error('%s: %d of %d frame times have no camera edge in %s within %.1f ms (or the order differs)', ...
        FRAMES_CSV, sum(abs(edgeDt) > TOL_S), numel(edgeDt), DIGITAL_FILE, 1000*TOL_S);
end
mapRows = rowOf;
selIdx  = idxFrames(pos);
alignmentMode = 'sync_frames_csv';
offset = pos(1) - 1;               % camera edges before video frame 1 (0: the video starts with the first edge)
nNoPulse    = N - numel(mapRows);  % pupil rows without a camera edge (digital fields NaN)
nEdgeUnused = nFrames - numel(mapRows);

% -------- Prepare output folder + logging --------
outDir = OUT_PARENT;
if ~exist(outDir,'dir'), mkdir(outDir); end
stamp   = datestr(now,'yyyymmdd_HHMMSS');
if nargin >= 5 && ~isempty(RUN_NAME)
    RUN_ROOT = fullfile(outDir, char(RUN_NAME));   % fix_concat_alignment.py: a temporary name until its checks pass
else
    RUN_ROOT = fullfile(outDir, ['concat_out_' stamp '_syncaligned']);   % SYNC COPY: name says how rows were matched
end
if exist(RUN_ROOT,'dir'), error('output folder already exists: %s', RUN_ROOT); end
mkdir(RUN_ROOT);

diary(fullfile(RUN_ROOT, 'log.txt'));
fprintf('DIGITAL: %s\n', DIGITAL_FILE);
fprintf('  TTL rising edges (raw)=%d, ttlHzRaw≈%.3f, ratioRaw (TTL/pupil)=%.3f\n', nTTLraw, ttlHzRaw, ratioRaw);
fprintf('PUPIL file: %s | rows=%d\n', PUPIL_FILE, nPupil);
fprintf('Estimated frame period from TTL (s): %.6f\n', estFramePeriod);
fprintf('Chosen strategy: %s | TTL frames used=%d\n', strategy, nFrames);
fprintf('Alignment: video frame k <-> camera pulse of frame k in %s (sync_sessions.py)\n', FRAMES_CSV);
fprintf('  rows with a camera edge: %d (rows %d-%d) | rows without: %d | camera edges not used: %d | max |edge - frames.csv| = %.2g s\n', ...
    numel(mapRows), mapRows(1), mapRows(end), nNoPulse, nEdgeUnused, max(abs(edgeDt)));
fprintf('  (the original pupil_align_ttl would use: %s, offset %d)\n', origMode, origOffset);

CamTime_s(mapRows)   = D.Time_s(selIdx);
Camera_TTL(mapRows)  = double(D.camera(selIdx));
Lick_TTL(mapRows)    = double(D.lick(selIdx));
if hasInjector
    Injector_TTL(mapRows) = double(D.injector(selIdx));
end

% Relative time based on first non-NaN camera time
CamTime_rel_s = CamTime_s;
firstValid = find(~isnan(CamTime_s), 1, 'first');
if ~isempty(firstValid)
    CamTime_rel_s = CamTime_s - CamTime_s(firstValid);
end

% -------- Choose plot time base (always 30 fps as requested) --------
FPS_PLOT = 30;
PlotTime_s = ((1:N)' - 1) / FPS_PLOT;    % 0, 1/30, 2/30, ...

% -------- Assemble final table (one row per pupil frame) --------
Frame             = (1:N).';
PupilTimestamp_s  = p_ts(:);
Diameter_px       = p_diam(:);
CenterX_px        = p_cx(:);
CenterY_px        = p_cy(:);
UsedScoreFallback = p_score(:);

Tout = table(Frame, PupilTimestamp_s, Diameter_px, CenterX_px, CenterY_px, UsedScoreFallback, ...
             CamTime_s, CamTime_rel_s, Camera_TTL, Lick_TTL, Injector_TTL);
Tout.PlotTime_s_30fps = PlotTime_s;      % for plotting only

% -------- Save combined CSV --------
outCSV = fullfile(RUN_ROOT, 'combined_pupil_digital.csv');
writetable(Tout, outCSV);
fprintf('Saved combined table: %s\n', outCSV);

% -------- Save lick onsets --------
lick_bin = Lick_TTL; lick_bin(isnan(lick_bin)) = 0; lick_bin = lick_bin > 0.5;
lick_prev = [false; lick_bin(1:end-1)];
lick_on_idx = find(lick_bin & ~lick_prev);
LickOn = table(lick_on_idx, ...
               PupilTimestamp_s(lick_on_idx), ...
               CamTime_s(lick_on_idx), ...
               CamTime_rel_s(lick_on_idx), ...
               'VariableNames', {'Frame','PupilTimestamp_s','CamTime_s','CamTime_rel_s'});
outLick = fullfile(RUN_ROOT, 'lick_onsets.csv');
writetable(LickOn, outLick);
fprintf('Saved lick onsets: %s (n=%d)\n', outLick, height(LickOn));

% -------- Save injector onsets (if present) --------
if hasInjector
    inj_bin = Injector_TTL; inj_bin(isnan(inj_bin)) = 0; inj_bin = inj_bin > 0.5;
    inj_prev = [false; inj_bin(1:end-1)];
    inj_on_idx = find(inj_bin & ~inj_prev);
    InjOn = table(inj_on_idx, ...
                  PupilTimestamp_s(inj_on_idx), ...
                  CamTime_s(inj_on_idx), ...
                  CamTime_rel_s(inj_on_idx), ...
                  'VariableNames', {'Frame','PupilTimestamp_s','CamTime_s','CamTime_rel_s'});
    outInj = fullfile(RUN_ROOT, 'injector_onsets.csv');
    writetable(InjOn, outInj);
    fprintf('Saved injector onsets: %s (n=%d)\n', outInj, height(InjOn));
else
    inj_on_idx = [];
    outInj = 'N/A';
end

% -------- Plot (3-panel: pupil, stacked TTL tracks, onset raster) --------
x = PlotTime_s;          % fixed 30 fps time for figures
y = Diameter_px;

lick_plot = Lick_TTL;     lick_plot(isnan(lick_plot)) = 0;  lick_plot = lick_plot > 0.5;
inj_plot  = Injector_TTL; inj_plot(isnan(inj_plot))  = 0;  inj_plot  = inj_plot  > 0.5;

lick_track = double(lick_plot) * 1;       % Lick at y=1
inj_track  = double(inj_plot)  * 1 + 1;   % Injector at y=2

fig = figure('Color','w','Position',[100 100 1300 700]);
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(x, y, '-', 'LineWidth', 1.1);
grid on; box off;
xlabel('Time (s)'); ylabel('Pupil diameter (px)');
title('Pupil diameter');
xlim([x(1) x(end)]);

nexttile;
stairs(x, lick_track, '-', 'LineWidth', 1.0); hold on;
if hasInjector
    stairs(x, inj_track,  '--', 'LineWidth', 1.0);
end
grid on; box off;
xlabel('Time (s)'); ylabel('TTL tracks');
yticks([0 1 2]); yticklabels({'0','Lick','Injector'});
ylim([-0.2 2.5]); xlim([x(1) x(end)]);
if hasInjector, legend({'Lick','Injector'}, 'Location','best'); else, legend({'Lick'}, 'Location','best'); end
title('Stacked TTL tracks (Lick @1, Injector @2)');

nexttile; hold on;
if ~isempty(lick_on_idx)
    stem(x(lick_on_idx), ones(size(lick_on_idx))*1, 'filled', 'LineWidth', 0.9, 'MarkerSize', 3);
end
if ~isempty(inj_on_idx)
    stem(x(inj_on_idx),  ones(size(inj_on_idx))*2, 'filled', 'LineWidth', 0.9, 'MarkerSize', 3);
end
grid on; box off;
xlabel('Time (s)'); ylabel('Onsets');
yticks([1 2]); yticklabels({'Lick','Injector'});
ylim([0.5 2.5]); xlim([x(1) x(end)]);
title('Onset raster (event timing)');

outPNG = fullfile(RUN_ROOT, 'plot_pupil_lick_injector_STACKED_RASTER.png');
outFIG = fullfile(RUN_ROOT, 'plot_pupil_lick_injector_STACKED_RASTER.fig');
try
    exportgraphics(fig, outPNG, 'Resolution', 300);
catch
    saveas(fig, outPNG);
end
savefig(fig, outFIG);
close(fig);
fprintf('Saved figure: %s and %s\n', outPNG, outFIG);

% -------- Write summary.txt --------
summaryFile = fullfile(RUN_ROOT, 'summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Concatenation summary (no video)\n');
fprintf(fid, 'Timestamp: %s\n', datestr(now, 31));
fprintf(fid, 'Digital file: %s\n', DIGITAL_FILE);
fprintf(fid, '  TTL rising edges (raw)=%d, ttlHzRaw≈%.3f, TTL/pupil ratio=%.3f\n', nTTLraw, ttlHzRaw, ratioRaw);
fprintf(fid, 'Estimated frame period from TTL (s): %.6f\n', estFramePeriod);
fprintf(fid, 'Chosen strategy: %s\n', strategy);
fprintf(fid, 'Pupil file: %s\n', PUPIL_FILE);
fprintf(fid, 'nPupil rows: %d\n', nPupil);
fprintf(fid, 'nCamera frames used: %d\n', nFrames);
fprintf(fid, 'Alignment mode: %s\n', alignmentMode);
fprintf(fid, 'Offset: %d\n', offset);
if ~isempty(firstValid)
    fprintf(fid, 'First valid camera time (s): %.6f\n', CamTime_s(firstValid));
else
    fprintf(fid, 'First valid camera time (s): NaN\n');
end
fprintf(fid, 'Has injector column: %d\n', hasInjector);
fprintf(fid, 'Output CSV: %s\n', outCSV);
fprintf(fid, 'Lick onsets CSV: %s\n', outLick);
fprintf(fid, 'Injector onsets CSV: %s\n', outInj);
fprintf(fid, 'Figure PNG: %s\n', outPNG);
fprintf(fid, 'Figure FIG: %s\n', outFIG);
fprintf(fid, 'Alignment source: %s (sync_sessions.py; video frame k <-> camera pulse of frame k)\n', FRAMES_CSV);
fprintf(fid, 'Rows with a camera edge: %d (rows %d-%d); rows without: %d; camera edges not used: %d\n', ...
    numel(mapRows), mapRows(1), mapRows(end), nNoPulse, nEdgeUnused);
fprintf(fid, 'Max |camera edge time - frames.csv time| (s): %.3g\n', max(abs(edgeDt)));
fprintf(fid, 'The original pupil_align_ttl would use: %s, offset %d\n', origMode, origOffset);
fclose(fid);
fprintf('Wrote summary: %s\n', summaryFile);

% also print durations for sanity
ttl_dt = diff(CamTime_rel_s(~isnan(CamTime_rel_s)));
if ~isempty(ttl_dt)
    fprintf('TTL-derived median rate ≈ %.3f Hz | duration ≈ %.1f s\n', ...
        1/median(ttl_dt), CamTime_rel_s(find(~isnan(CamTime_rel_s),1,'last')));
end
fprintf('Plot (fixed) 30 fps duration ≈ %.1f s\n', PlotTime_s(end));

diary off;

info = struct('nTTLraw', nTTLraw, 'ttlHzRaw', ttlHzRaw, 'ratioRaw', ratioRaw, 'strategy', strategy, ...
    'nPupil', nPupil, 'nFrames', nFrames, 'alignmentMode', alignmentMode, 'offset', offset, ...
    'nLickOnsets', numel(lick_on_idx), 'nInjectorOnsets', numel(inj_on_idx), ...
    'nRowsNoPulse', nNoPulse, 'nEdgesUnused', nEdgeUnused, 'maxEdgeDt', max(abs(edgeDt)), ...
    'origMode', origMode, 'origOffset', origOffset);
end

% -------- Helpers (unchanged) --------
function vname = pickVar(T, candidates)
    VN = string(T.Properties.VariableNames);
    VNl = lower(VN);
    for c = 1:numel(candidates)
        patt = lower(string(candidates{c}));
        idx = find(VNl == patt, 1);
        if ~isempty(idx), vname = char(VN(idx)); return; end
        patt2 = replace(patt, ["[","]","_"], "");
        idx = find(contains(VNl, patt2), 1);
        if ~isempty(idx), vname = char(VN(idx)); return; end
    end
    error('Could not find any of: %s\nAvailable: %s', ...
          strjoin(string(candidates), ', '), strjoin(string(VN), ', '));
end

function x = numify(v)
    if islogical(v), x = double(v); return; end
    if iscell(v),    x = str2double(string(v)); return; end
    x = double(v);
end
