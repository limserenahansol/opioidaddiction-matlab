function eye_track_session(video, crop, modelFile, outDir, varargin)
% eye_track_session  Whole-eye size in every frame of one session (9670f black, 2026-11-17).
% crop = [x0 y0 w h] in 0-based frame pixels (analysis/wholeeye/crops.json) - the same crop as in training.
% Video frame k = camera pulse k (frames are read in order). Output in outDir:
%   wholeeye_frames.csv  frame, found, area_px, major_px, minor_px, height_px, width_px, cx, cy (full frame, 0-based),
%                        prob (mean eye probability inside the outline)
%   qc/f<k>_img.png, qc/f<k>_mask.png  24 frames spread over the session, for the example images
%   wholeeye_qc.json     frames read, found %, median area
p = inputParser; p.addParameter('Batch', 256); p.addParameter('Threads', 6); p.addParameter('NQC', 24);
p.parse(varargin{:}); o = p.Results;
if o.Threads > 0, maxNumCompThreads(o.Threads); end
M = load(modelFile); S = M.INPUT; net = M.net;
v = VideoReader(video);
x0 = crop(1); y0 = crop(2); w = crop(3); h = crop(4);
rows = y0 + (1:h); cols = x0 + (1:w);           % 0-based crop -> 1-based MATLAB indices
nEst = max(1, round(v.Duration * 60));          % IC Capture stamps frames 1/60 s apart
qcK = unique(round(linspace(1, max(1, nEst - 30), o.NQC)));
if ~exist(fullfile(outDir, 'qc'), 'dir'), mkdir(fullfile(outDir, 'qc')); end
R = zeros(0, 10);
buf = zeros(S, S, 1, o.Batch, 'uint8'); bk = zeros(o.Batch, 1); nb = 0; k = 0;
qcCrop = containers.Map('KeyType', 'double', 'ValueType', 'any');
t0 = tic;
while hasFrame(v)
    F = readFrame(v); k = k + 1;
    C = F(rows, cols, :); if size(C, 3) > 1, C = rgb2gray(C); end
    nb = nb + 1; buf(:, :, 1, nb) = imresize(eye_prep(C), [S S], 'bilinear'); bk(nb) = k;
    if any(qcK == k), qcCrop(k) = C; end
    if nb == o.Batch, R = [R; flush(buf(:, :, :, 1:nb), bk(1:nb))]; nb = 0; end %#ok<AGROW>
end
if nb > 0, R = [R; flush(buf(:, :, :, 1:nb), bk(1:nb))]; end
T = array2table(R, 'VariableNames', {'frame', 'found', 'area_px', 'major_px', 'minor_px', 'height_px', ...
    'width_px', 'cx', 'cy', 'prob'});
writetable(T, fullfile(outDir, 'wholeeye_frames.csv'));
Q = struct('video', video, 'crop', crop, 'model', modelFile, 'frames_read', k, ...
    'found_pct', 100 * mean(T.found), 'median_area_px', median(T.area_px(T.found > 0)), ...
    'minutes', toc(t0) / 60, 'when', char(datetime('now')));
fid = fopen(fullfile(outDir, 'wholeeye_qc.json'), 'w'); fprintf(fid, '%s', jsonencode(Q, 'PrettyPrint', true)); fclose(fid);
fprintf('[eye] %s: %d frames, found %.1f %%, median area %.0f px, %.1f min\n', video, k, Q.found_pct, ...
    Q.median_area_px, Q.minutes);

    function out = flush(B, ks)
        P = predict(net, B, 'MiniBatchSize', 64);
        out = zeros(numel(ks), 10);
        for j = 1:numel(ks)
            pm = imresize(P(:, :, 2, j), [h w], 'bilinear');
            m = pm > 0.5;
            m = imopen(m, strel('disk', 1)); m = imfill(m, 'holes');
            if any(m(:)), m = bwareafilt(m, 1); end
            r = [ks(j), 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
            if nnz(m) >= 30
                st = regionprops(m, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'BoundingBox', 'Centroid');
                r = [ks(j), 1, st.Area, st.MajorAxisLength, st.MinorAxisLength, st.BoundingBox(4), ...
                    st.BoundingBox(3), x0 + st.Centroid(1) - 1, y0 + st.Centroid(2) - 1, mean(pm(m))];
            end
            out(j, :) = r;
            if isKey(qcCrop, ks(j))
                imwrite(qcCrop(ks(j)), fullfile(outDir, 'qc', sprintf('f%d_img.png', ks(j))));
                imwrite(m, fullfile(outDir, 'qc', sprintf('f%d_mask.png', ks(j))));
            end
        end
    end
end
