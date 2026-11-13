function R = run_pupil_one_click(varargin)
% run_pupil_one_click  ONE call: pupil diameter for every mouse x day, aligned to the Saleae TTLs (2026-11-12).
%
%   run_pupil_one_click                              % everything with the defaults below
%   run_pupil_one_click('DryRun', true)              % only list what it found / would do (nothing written)
%   run_pupil_one_click('Days', [19 20], 'Mice', {'0910_black'})
%   run_pupil_one_click('VideoRoot',<videos>, 'SessionRoot',<sessions>, 'OutRoot',<outputs>)   % else one_click_paths.json
%
% Folder layout it expects (the one used for the Dec-2025 cohort):
%   VideoRoot\dayN\<cage>\<color>\*.avi (+ roi.mat)      cage = 4 digits, color = black / orange / red / white ...
%   SessionRoot\dayN\<cage>\<color>\digital.csv          (Saleae export; Ch6 camera, Ch7 lick, Ch4 injector)
%   If digital.csv is not there, it is searched anywhere under VideoRoot\dayN\ as ...\<cage>\<color>\digital.csv
%   (the Dec-2025 TTL folders were named ttl, tt, day4Jan2, Dec29day1, ...).
%
% For each session:
%   1) camera TTL rising edges from digital.csv
%   2) pupil video = the .avi whose frame count is closest to the TTL count (flagged if off by > 2 frames)
%   3) ROI = roi.mat in the video folder; if missing: ROIMode 'ask' (default) collects ALL missing ROIs at the start
%      (one rectangle per session, previous day's ROI shown, Enter keeps it), 'previous' reuses the same mouse's last ROI
%   4) pupil_track_session with the general U-Net  ->  OutRoot\dayN\<cage>\<color>\pupil_outputs_general\
%   5) pupil_align_ttl (the original alignment)    ->  OutRoot\dayN\<cage>\<color>\concat_out_<stamp>\
%   Sessions with both outputs already present are skipped ('Redo', true to redo).
% At the end: OutRoot\one_click_report_<stamp>.csv  (one row per session, QC numbers, FLAGS column) + log.
%
% Paths: pass them, or put them in one_click_paths.json next to this file (copy one_click_paths.example.json; the .json is
% git-ignored so no local path is committed).
% IMPORTANT: Step 02 counts EVERY concat_out* folder in a session folder as a separate session. So keep OutRoot a separate
% folder. Use OutRoot = SessionRoot only for NEW sessions that have no concat_out yet - the runner refuses to add a second
% concat_out to a session folder that already has one.

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('VideoRoot',   '');           % '' = from one_click_paths.json (local, not in git)
p.addParameter('SessionRoot', '');
p.addParameter('OutRoot',     '');
p.addParameter('Model',       fullfile(here, 'model', 'general_pupil_unet.mat'));
p.addParameter('Days',        []);           % [] = every dayN folder
p.addParameter('Mice',        {});           % {} = all, else {'0910_black', ...}
p.addParameter('ROIMode',     'ask');        % 'ask' | 'previous'
p.addParameter('DryRun',      false);
p.addParameter('Redo',        false);
p.addParameter('MaxFrames',   Inf);          % for quick tests
p.parse(varargin{:});
o = p.Results;
o = fill_paths(o, here);
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
if ~o.DryRun && ~exist(o.OutRoot, 'dir'), mkdir(o.OutRoot); end
if ~o.DryRun, diary(fullfile(o.OutRoot, ['one_click_log_' stamp '.txt'])); end
cleanup = onCleanup(@() diary('off'));
fprintf('=== run_pupil_one_click %s ===\nVideoRoot %s\nSessionRoot %s\nOutRoot %s\nModel %s\n', stamp, o.VideoRoot, o.SessionRoot, o.OutRoot, o.Model);

% ---------------- 1. discover sessions ----------------
L = discover(o);
fprintf('found %d sessions (video folders with .avi)\n', height(L));
if isempty(L), R = L; return; end

% ---------------- 2. TTL + video choice ----------------
for i = 1:height(L)
    L.digital(i) = string(find_digital(o, L.day(i), L.cage(i), L.color(i)));
    if strlength(L.digital(i)) > 0, L.ttl_edges(i) = count_camera_edges(L.digital(i)); end
    vids = dir(fullfile(L.video_dir(i), '*.avi')); nf = zeros(numel(vids),1); fr = zeros(numel(vids),1);
    for v = 1:numel(vids), [nf(v), fr(v)] = video_frames(fullfile(vids(v).folder, vids(v).name)); end
    L.videos(i) = strjoin(string({vids.name}) + "(" + string(nf') + " fr, " + string(fr') + " fps)", '; ');
    if ~isempty(vids)
        % pupil camera = frame count closest to the camera TTLs (243/244 Dec-2025 sessions); tie -> the 40-fps-header file
        if L.ttl_edges(i) > 0, cost = abs(nf - L.ttl_edges(i)); else, cost = -nf; end
        cost = cost + 0.5*(abs(fr - 40) > 0.5);
        [~, b] = min(cost);
        L.video(i) = string(fullfile(vids(b).folder, vids(b).name)); L.video_frames(i) = nf(b);
        L.frame_minus_ttl(i) = nf(b) - L.ttl_edges(i);
        L.camera_ambiguous(i) = numel(vids) > 1 && sum(abs(nf - nf(b)) <= 2) > 1;
    end
    rf = fullfile(L.video_dir(i), 'roi.mat');
    if exist(rf, 'file'), r = load(rf); L.roi(i) = string(num2str(double(r.roi(:)'))); L.roi_source(i) = "roi.mat"; end
    fprintf('  %-6s %s_%-7s TTL %6d | video %-34s %6d | roi %s\n', L.day(i), L.cage(i), L.color(i), L.ttl_edges(i), ...
        shortname(L.video(i)), L.video_frames(i), L.roi(i));
end

% ---------------- 3. ROIs that are missing ----------------
need = strlength(L.roi) == 0 & strlength(L.video) > 0;
if any(need)
    fprintf('%d sessions have no roi.mat (ROIMode = %s)\n', nnz(need), o.ROIMode);
    L = fill_rois(L, need, o);
end

done_ = false(height(L),1);
for i = 1:height(L)
    od = fullfile(o.OutRoot, L.day(i), L.cage(i), L.color(i));
    done_(i) = exist(fullfile(od, 'pupil_outputs_general', 'qc.json'), 'file') && ~isempty(dir(fullfile(od, 'concat_out_*')));
end
L.already_done = done_;
if o.DryRun
    R = L; fprintf('\nDryRun: %d to run, %d already done, %d without TTL, %d without ROI\n', ...
        nnz(~done_), nnz(done_), nnz(L.ttl_edges <= 0), nnz(strlength(L.roi) == 0)); return;
end

% ---------------- 4-5. track + align ----------------
M = load(o.Model); fprintf('model loaded (%s, input %d)\n', o.Model, M.INPUT);
for i = 1:height(L)
    tag = sprintf('%s %s_%s', L.day(i), L.cage(i), L.color(i));
    if L.already_done(i) && ~o.Redo, L.status(i) = "skipped (done before)"; continue; end
    if strlength(L.video(i)) == 0 || strlength(L.roi(i)) == 0, L.status(i) = "no video or ROI"; continue; end
    od = fullfile(o.OutRoot, L.day(i), L.cage(i), L.color(i));
    fprintf('\n[%d/%d] %s\n', i, height(L), tag);
    try
        Q = pupil_track_session(L.video(i), str2num(L.roi(i)), M, fullfile(od, 'pupil_outputs_general'), 'MaxFrames', o.MaxFrames); %#ok<ST2NM>
        L.n_frames(i) = Q.nFrames; L.frac_no_blob(i) = Q.fracNoBlob; L.frac_rejected(i) = Q.fracRejected;
        L.median_diam_px(i) = Q.medianDiameter; L.track_s(i) = Q.seconds;
        if strlength(L.digital(i)) > 0
            if strcmpi(strip(string(o.OutRoot),'right','\'), strip(string(o.SessionRoot),'right','\')) && ...
                    ~isempty(dir(fullfile(od, 'concat_out_*')))
                error('refusing to add a 2nd concat_out to %s (Step 02 would count it as another session)', od);
            end
            [rr, info] = pupil_align_ttl(char(L.digital(i)), char(Q.cleanCsv), char(od));
            L.concat_out(i) = string(rr); L.align_mode(i) = string(info.alignmentMode); L.align_offset(i) = info.offset;
            L.n_lick_onsets(i) = info.nLickOnsets;
        end
        L.status(i) = "ok";
    catch ME
        L.status(i) = "FAILED: " + string(ME.message); fprintf(2, '  FAILED: %s\n', ME.message);
    end
end

% ---------------- 6. report ----------------
L.flags = strings(height(L),1);
key = L.cage + "_" + L.color;
for i = 1:height(L)
    f = strings(0);
    if L.ttl_edges(i) <= 0, f(end+1) = "no TTL"; end %#ok<AGROW>
    if abs(L.frame_minus_ttl(i)) > 2 && L.ttl_edges(i) > 0, f(end+1) = sprintf("video-TTL %+d frames", L.frame_minus_ttl(i)); end %#ok<AGROW>
    if L.camera_ambiguous(i), f(end+1) = "2 videos with the same frame count - check which is the pupil camera"; end %#ok<AGROW>
    if L.roi_source(i) ~= "roi.mat" && strlength(L.roi_source(i)) > 0, f(end+1) = "ROI " + L.roi_source(i); end %#ok<AGROW>
    if L.frac_no_blob(i) > 0.05, f(end+1) = sprintf("no pupil in %.0f%% frames", 100*L.frac_no_blob(i)); end %#ok<AGROW>
    if L.frac_rejected(i) > 0.10, f(end+1) = sprintf("%.0f%% frames rejected", 100*L.frac_rejected(i)); end %#ok<AGROW>
    mm = median(L.median_diam_px(key == key(i) & L.median_diam_px > 0), 'omitnan');
    if L.median_diam_px(i) > 0 && abs(L.median_diam_px(i)/mm - 1) > 0.35, f(end+1) = sprintf("diameter %.0f%% off this mouse's median", 100*(L.median_diam_px(i)/mm-1)); end %#ok<AGROW>
    if startsWith(L.status(i), "FAILED"), f(end+1) = "FAILED"; end %#ok<AGROW>
    L.flags(i) = strjoin(f, '; ');
end
rep = fullfile(o.OutRoot, ['one_click_report_' stamp '.csv']);
writetable(L, rep);
fprintf('\nreport: %s\n  ok %d | skipped %d | failed %d | flagged %d\n', rep, nnz(L.status == "ok"), ...
    nnz(startsWith(L.status, "skipped")), nnz(startsWith(L.status, "FAILED")), nnz(strlength(L.flags) > 0));
R = L;
end

% ======================================================================
function o = fill_paths(o, here)
% empty roots come from one_click_paths.json next to this file (copy one_click_paths.example.json and edit it)
f = fullfile(here, 'one_click_paths.json'); keys = {'VideoRoot','video_root'; 'SessionRoot','session_root'; 'OutRoot','out_root'};
need = cellfun(@(k) isempty(o.(k)), keys(:,1));
if ~any(need), return; end
assert(exist(f, 'file') == 2, 'Pass VideoRoot/SessionRoot/OutRoot, or create %s from one_click_paths.example.json', f);
c = jsondecode(fileread(f));
for i = find(need)'
    v = c.(keys{i,2}); assert(~startsWith(v, '<'), 'Set %s in %s', keys{i,2}, f); o.(keys{i,1}) = v;
end
end

function L = discover(o)
rows = {};
dd = dir(fullfile(o.VideoRoot, 'day*'));
for d = dd(:)'
    tok = regexp(d.name, '^day(\d+)$', 'tokens', 'once');
    if ~d.isdir || isempty(tok), continue; end
    di = str2double(tok{1});
    if ~isempty(o.Days) && ~ismember(di, o.Days), continue; end
    cc = dir(fullfile(d.folder, d.name));
    for c = cc(:)'
        if ~c.isdir || isempty(regexp(c.name, '^\d{4}$', 'once')), continue; end
        kk = dir(fullfile(c.folder, c.name));
        for k = kk(:)'
            if ~k.isdir || startsWith(k.name, '.') || isempty(dir(fullfile(k.folder, k.name, '*.avi'))), continue; end
            if ~isempty(o.Mice) && ~ismember([c.name '_' k.name], o.Mice), continue; end
            rows(end+1,:) = {string(d.name), di, string(c.name), string(k.name), string(fullfile(k.folder, k.name))}; %#ok<AGROW>
        end
    end
end
if isempty(rows), L = table(); return; end
L = cell2table(rows, 'VariableNames', {'day','day_index','cage','color','video_dir'});
L = sortrows(L, {'day_index','cage','color'});
n = height(L); z = zeros(n,1); s = strings(n,1);
L.digital = s; L.ttl_edges = z; L.videos = s; L.video = s; L.video_frames = z; L.frame_minus_ttl = z;
L.roi = s; L.roi_source = s; L.status = s; L.n_frames = z; L.frac_no_blob = z; L.frac_rejected = z;
L.median_diam_px = z; L.track_s = z; L.concat_out = s; L.align_mode = s; L.align_offset = z; L.n_lick_onsets = z;
L.camera_ambiguous = false(n,1);
end

function f = find_digital(o, day, cage, color)
f = fullfile(o.SessionRoot, day, cage, color, 'digital.csv');
if exist(f, 'file'), return; end
hits = dir(fullfile(o.VideoRoot, day, '**', cage, color, 'digital.csv'));
if isempty(hits), f = ''; return; end
[~, b] = max([hits.datenum]); f = fullfile(hits(b).folder, hits(b).name);
end

function n = count_camera_edges(f)
opts = detectImportOptions(f, 'VariableNamingRule', 'preserve');
T = readtable(f, opts);
vn = T.Properties.VariableNames;
c = find(contains(lower(vn), 'channel 6') | strcmpi(vn, 'camera'), 1);
if isempty(c), n = -1; return; end
cam = T{:, c} > 0.5;
n = nnz(cam & [false; ~cam(1:end-1)]);
end

function [n, fps] = video_frames(f)
% Duration x FrameRate: same as NumFrames within 1 frame on these AVIs, but < 1 s instead of ~10 s per 3 GB file
try
    v = VideoReader(f); n = round(v.Duration * v.FrameRate); fps = v.FrameRate;
catch
    n = -1; fps = NaN;
end
end

function s = shortname(f)
[~, s, e] = fileparts(char(f)); s = [s e];
end

function L = fill_rois(L, need, o)
idx = find(need);
for i = idx'
    key = L.cage(i) + "_" + L.color(i);
    prev = find(L.cage == L.cage(i) & L.color == L.color(i) & strlength(L.roi) > 0 & L.day_index < L.day_index(i), 1, 'last');
    sugg = []; if ~isempty(prev), sugg = str2num(L.roi(prev)); end %#ok<ST2NM>
    if strcmpi(o.ROIMode, 'previous') || o.DryRun
        if ~isempty(sugg) && ~o.DryRun
            L.roi(i) = string(num2str(sugg)); L.roi_source(i) = "reused from " + L.day(prev);
        end
        continue
    end
    v = VideoReader(L.video(i)); v.CurrentTime = v.Duration/2; F = readFrame(v);
    fig = figure('Name', sprintf('ROI %s %s', L.day(i), key), 'Color', 'w');
    imshow(F, 'InitialMagnification', 'fit');
    pos = []; src = "";
    try
        if ~isempty(sugg)
            title(sprintf(['%s %s: the box is the previous day''s ROI. Move/resize it over the eye and DOUBLE-CLICK it.\n' ...
                'Close the window to keep it as it is.'], L.day(i), key), 'Interpreter', 'none');
            r = drawrectangle('Position', sugg, 'StripeColor', 'm');
        else
            title(sprintf('%s %s: drag a box around the eye, then DOUBLE-CLICK it.', L.day(i), key), 'Interpreter', 'none');
            r = drawrectangle('StripeColor', 'm');
        end
        wait(r); pos = round(r.Position); src = "drawn " + string(datetime('today'));
    catch
        % window closed: keep the previous day's ROI (if there is one)
    end
    if isvalid(fig), close(fig); end
    if isempty(pos) || any(pos(3:4) < 5)
        if isempty(sugg), fprintf(2, '  no ROI for %s %s - session skipped\n', L.day(i), key); continue; end
        pos = sugg; src = "kept from " + L.day(prev);
    end
    roi = pos; save(fullfile(L.video_dir(i), 'roi.mat'), 'roi');
    L.roi(i) = string(num2str(roi)); L.roi_source(i) = "roi.mat (" + src + ")";
end
end

