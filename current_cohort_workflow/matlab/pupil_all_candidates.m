function pupil_all_candidates(videoFile, roiFile, modelFile, pkgPupilDir, ftDir, outCsv, varargin)
% pupil_all_candidates  Every pupil candidate of every video frame, for the best-effort tracking of hard eyes (2026-11-16).
%
% The package's tracker keeps one blob per frame (hybrid_detections.csv: the best candidate only). For eyes where it often
% picks the wrong blob (narrowed / inflamed / cloudy eyes), the choice has to be redone with the whole candidate list. This
% function writes that list. It changes nothing: it only reads the video, the session's roi.mat and the model.
%
% Per frame, exactly as the package does it:
%   pupil_mask_frame(frames, roi, model, 'IrisCorrection', true)   (package code; U-Net mask + probability map)
%   pupil_hybrid_candidates(crop, mask)                             (pupil_finetune code; U-Net blobs + dark-region chains,
%                                                                    the candidates shown on the labelling sheets)
% Per candidate it adds the bounding-box width / height, major / minor axis, and the U-Net probability inside it.
%
% Output CSV (one row per candidate; a frame without any candidate gets one row with cand = 0):
%   frame (0-based = row of results2.csv), cand, src (1 = U-Net blob, 2 = dark region), area, diam (equivalent diameter),
%   width, height, major, minor, orient, cx, cy (ROI-crop pixels, MATLAB 1-based), circ, sol, ecc, con, dark, dglint,
%   dctr, unet_ov, maxp, meanp
%
% usage (MATLAB -batch):
%   pupil_all_candidates(video, roiMat, modelMat, '<package>\pupil', '<current_cohort_workflow>\pupil_finetune', outCsv)
%   pupil_all_candidates(..., 'MaxFrames', 600, 'Batch', 256)
p = inputParser;
p.addParameter('MaxFrames', Inf);
p.addParameter('Batch', 256);
p.parse(varargin{:});
o = p.Results;
addpath(pkgPupilDir); addpath(ftDir);
M = load(modelFile);
r = load(roiFile); roi = double(r.roi(:)');
v = VideoReader(videoFile);
nTot = min(round(v.Duration * v.FrameRate), o.MaxFrames);
tmp = [outCsv '.part'];
fid = fopen(tmp, 'w');
fprintf(fid, ['frame,cand,src,area,diam,width,height,major,minor,orient,cx,cy,circ,sol,ecc,con,dark,dglint,dctr,' ...
    'unet_ov,maxp,meanp\n']);
k = 0; t0 = tic;
while hasFrame(v) && k < nTot
    F = {};
    while hasFrame(v) && numel(F) < o.Batch && k + numel(F) < nTot
        F{end + 1} = readFrame(v); %#ok<AGROW>
    end
    [masks, ~, crops, ~, ~, ~, prob] = pupil_mask_frame(F, roi, M, 'IrisCorrection', true);
    for j = 1:numel(F)
        fr = k + j - 1;
        cand = pupil_hybrid_candidates(crops{j}, masks{j});
        if isempty(cand)
            fprintf(fid, '%d,0,0,,,,,,,,,,,,,,,,,,,\n', fr);
            continue
        end
        P = prob{j};
        for c = 1:numel(cand)
            m = cand(c).m;
            R = regionprops(m, 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Area');
            [~, b] = max([R.Area]); R = R(b);
            pv = P(m);
            src = 1 + double(~strcmp(cand(c).src, 'unet'));
            fprintf(fid, '%d,%d,%d,%d,%.3f,%.2f,%.2f,%.3f,%.3f,%.2f,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f,%.4f,%.3f,%.4f,%.4f\n', ...
                fr, c, src, cand(c).area, cand(c).diam, R.BoundingBox(3), R.BoundingBox(4), R.MajorAxisLength, ...
                R.MinorAxisLength, R.Orientation, cand(c).cx, cand(c).cy, cand(c).circ, cand(c).sol, cand(c).ecc, ...
                cand(c).con, cand(c).dark, cand(c).dglint, cand(c).dctr, cand(c).unet_ov, max(pv), mean(pv));
        end
    end
    k = k + numel(F);
    if mod(k, 2048) < o.Batch
        fprintf('%s: %d / %d frames (%.1f min)\n', outCsv, k, nTot, toc(t0) / 60);
    end
end
fclose(fid);
movefile(tmp, outCsv, 'f');
fprintf('done %s: %d frames in %.1f min\n', outCsv, k, toc(t0) / 60);
end
