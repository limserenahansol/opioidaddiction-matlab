 % =======================================================================
% Pupil extraction (ROI-only) with POINT SEEDS + WATERSHED candidate
% - You click: center + 4 edge points (top/right/bottom/left) on 2 frames
%   (~10% & ~90%). During Step-1 it WILL ask for a 3rd seed on the
%   "hardest so far" frame, even if Step-1 looks good.
% - Same ROI is used for all seeds. Clicks are constrained to the ROI view.
% - If still struggling after 3rd seed, it asks one more guide:
%   click 1 pixel INSIDE pupil + 1 pixel just OUTSIDE (same ROI view).
% - Seeds build REF (center/radii/polarity/threshold band)
% - Step-1: generate pseudo-labels from dark/bright/watershed candidates
%   with area/shape/border/contrast/temporal gates + diversity filter
% - Step-2: train U-Net on Step-1 pairs (or fallback) with class weights
% =======================================================================
clear; clc; close all; rng default;
% ===================== PATHS (EDIT) =====================

% BASE = 'E:\Aug_new_addiction\addiction_aug20';
% VIDEO = fullfile(BASE,'side','Fr1_6872_suc_d2_side.avi');
BASE      = 'H:\addiction_dec29\day18\6099\';
VIDEO     = fullfile(BASE,'black','6099BBBB0039 26-01-22 11-01-38.avi');
   % <-- change per video
ROI_FILE  = fullfile(BASE,'black','roi.mat');
% ===================== RUN FOLDERS ==================
% ===================== RUN FOLDERS ======================
STAMP = datestr(now,'yyyymmdd_HHMMSS');
RUN_ROOT = fullfile(BASE,'training_data1',['run_' STAMP]);
DIRS.img    = fullfile(RUN_ROOT,'images');
DIRS.lab    = fullfile(RUN_ROOT,'labels');
DIRS.qc     = fullfile(RUN_ROOT,'qc');
DIRS.audit  = fullfile(RUN_ROOT,'audit_step1');
DIRS.model  = fullfile(BASE,'trainpupil');
DIRS.seedimg= fullfile(RUN_ROOT,'seeds_images'); % normalized seed frames
DIRS.seedpts= fullfile(RUN_ROOT,'seeds_points'); % point-click json/mat
cellfun(@(p) ~exist(p,'dir') && mkdir(p), ...
  {DIRS.img,DIRS.lab,DIRS.qc,DIRS.audit,DIRS.model,DIRS.seedimg,DIRS.seedpts});
fprintf('=== RUN: %s ===\n', RUN_ROOT);

% ===================== DEBUG SWITCHES =====================
DEBUG.ENABLE        = true;
DEBUG.SAVE_ALL_CANDS= true;
DEBUG.PRINT_EVERY   = 25;
DEBUG.OUT           = fullfile(DIRS.qc,'debug'); if ~exist(DEBUG.OUT,'dir'), mkdir(DEBUG.OUT); end

% ===================== STEP-1 CONFIG (defaults, relaxed) ====================
S1 = struct();
S1.SAVE_EVERY_N_FRAMES   = 1;
S1.STOP_AFTER_PAIRS      = 800;

% NEW: spatial/size constraints tied to seeds
S1.SEARCH_RADIUS_PX      = [];        % [] -> auto from REF (≈1.6*max(rx,ry))
S1.MIN_AREA_FRAC_OF_SEED = 1.16;
S1.MAX_AREA_REL_TO_SEED  = 6;

% preprocessing / normalization
S1.USE_CLAHE      = true;
S1.USE_STRETCHLIM = false;
S1.ADAPTIVE_CLIP  = 0.01;

% polarity percentiles (bases)
S1.PERCENTILE_DARK = 42;
S1.BRIGHT_PCT      = 92;

% morphology
S1.OPEN_RADIUS     = 1;
S1.CLOSE_RADIUS    = 0;
S1.SHRINK_MAX_PX   = 1;
S1.RING_WIDTH_PX   = 6;

% shape gates
S1.MIN_CIRCULARITY = 0.10;
S1.MIN_SOLIDITY    = 0.88;
S1.MAX_ECCENTRICITY= 0.99;

% border + temporal gates
S1.BORDER_MARGIN_PX   = 2;S1.MIN_AREA_FRAC_OF_SEED = 0.20;      % allow much smaller than seed ellipse
S1.MAX_AREA_REL_TO_SEED  = 1.62;
S1.JUMP_MAX_PX        = 95;
S1.AREA_CHANGE_MAX_FRAC = 0.85;

% ROI handling
S1.AUTO_ROI_PAD_PX  = 120;
S1.MANUAL_ROI       = true;
S1.MANUAL_ROI_MARGIN= 40;

% ring contrast (base)
S1.DELTA_I_MIN      = 1.30;

% glint & polarity
S1.SUPPRESS_GLINT   = true;
S1.GLINT_PCT_MAX    = 0.50;
S1.USE_BRIGHT_POL   = false;

% area gates (auto-scaled by ROI)
S1.AREA_MIN_FRAC      = 0.0015;
S1.AREA_MAX_FRAC      = 0.80;
S1.MIN_PUPIL_AREA_ABS = 120;

% center tolerance (fraction of min(roi dims)); updated after REF
S1.CENTER_TOL_FRAC  = 0.70;

% diversity
S1.MIN_MEDIAN_ABS_DIFF = 2;
S1.DIVERSITY_MIN_DA    = 0.02;

% QC cadence
S1.QC_EVERY_SAVED   = 2;
S1.QC_GALLERY_MAX   = 40;

% ===== When / how to ask for the 3rd seed (DURING Step-1) =====
S1.THIRD = struct();
S1.THIRD.ENABLE          = true;
S1.THIRD.ALWAYS_ASK      = true;   % << ALWAYS prompt once mid-pass
S1.THIRD.FORCE_AFTER_SEEN= [];     % [] => auto ~15% of STOP_AFTER
S1.THIRD.MIN_PROCESSED   = 120;    % (still used for auto-bad-case trigger)
S1.THIRD.MIN_FAILS       = 50;
S1.THIRD.MIN_FAIL_RATE   = 0.65;

% ===== EXTRA GUIDE (after 3rd seed, still struggling) =====
S1.GUIDE = struct();
S1.GUIDE.ENABLE                 = true;
S1.GUIDE.MIN_AFTER_THIRD_FRAMES = 100;
S1.GUIDE.MIN_AFTER_THIRD_FAILS  = 35;
S1.GUIDE.MIN_AFTER_THIRD_RATE   = 0.60;
S1.GUIDE.DISK_RADIUS_PX         = 6;
S1.GUIDE.BAND_FRAC              = 0.40;

S1.JUMP_MAX_PX           = 40;    % was 40
S1.CENTER_TOL_FRAC       = 0.62;  % was 0.60
S1.MIN_AREA_FRAC_OF_SEED = 0.41;  % was 0.45 (allow foreshortened pupil)
S1.MAX_AREA_REL_TO_SEED  = 5;   % was 1.9 (allow dilation/foreshortening)
S1.MIN_CIRCULARITY       = 0.50;  % was 0.50 (oblique ellipse)
S1.MAX_ECCENTRICITY      = 0.96;  % was 0.95

% --- size & search tied to seeds ---
S1.MIN_AREA_FRAC_OF_SEED = 0.26;   % keep
S1.MAX_AREA_REL_TO_SEED  = 2;   % was 2.0 → tighter

% --- global ROI-relative area caps ---
S1.AREA_MAX_FRAC      = 0.8;      % was 0.80

% --- morphology / ring ---

S1.CLOSE_RADIUS  = 0;              % keep it off
S1.SHRINK_MAX_PX = 5;              % was 1 → allow 1–3 px shrink trial
S1.RING_WIDTH_PX = 5.5;              % was 6 → narrower contrast band

% --- shape gates ---
S1.MIN_CIRCULARITY = 0.52;         % was 0.50
S1.MAX_ECCENTRICITY= 0.94;         % was 0.98
S1.MIN_SOLIDITY    = 0.86;         % was 0.88

% --- border / jump ---
S1.BORDER_MARGIN_PX = 2;           % was 2
S1.JUMP_MAX_PX      = 50;          % keep

% --- contrast / glint ---
S1.DELTA_I_MIN   = 5;            % higher base (will be raised again by REF)
S1.GLINT_PCT_MAX = 0.30;           % was 0.50 → reject bright contamination



% ===================== STEP-2 CONFIG =====================
S2 = struct();
S2.encoderDepth     = 5;
S2.numFirstFilters  = 32;
S2.miniBatch        = 4;
S2.maxEpochs        = 120;
S2.baseLR           = 3e-4;
S2.weightClipMax    = 5;
S2.minTrainFrac     = 0.005;

% ===================== ROI LOAD / CREATE =================
vr = VideoReader(VIDEO);
numFramesApprox = max(1, floor(vr.Duration * vr.FrameRate));

if exist(ROI_FILE,'file')
  R = load(ROI_FILE); roi = R.roi;
  fprintf('ROI loaded: [x=%d y=%d w=%d h=%d]\n', roi);
else
  midFrameNum = max(1,floor(numFramesApprox/2));
  vr.CurrentTime = (midFrameNum-1)/vr.FrameRate; F0 = readFrame(vr);
  [H,W,~] = size(F0);
  figure('Name','Draw ROI around the eye (double-click to confirm)','Color','w');
  imshow(F0,'InitialMagnification','fit'); title('Draw ROI around the eye, double-click to confirm.');
  try
    if exist('drawrectangle','file')
      r = drawrectangle('StripeColor','m'); pos = round(wait(r));
    else
      h = imrect; pos = round(wait(h));
    end
  catch
    if exist('r','var') && isvalid(r), pos = round(r.Position);
    else, pos = round(getPosition(h));
    end
  end
  close(gcf);
  if isempty(pos), error('ROI drawing cancelled.'); end
  m = max(0,S1.MANUAL_ROI_MARGIN);
  x1 = max(1, pos(1)-m); y1 = max(1, pos(2)-m);
  x2 = min(W, pos(1)+pos(3)-1 + m);
  y2 = min(H, pos(2)+pos(4)-1 + m);
  roi = [x1, y1, x2-x1+1, y2-y1+1];
  save(ROI_FILE,'roi');
  fprintf('ROI saved: [x=%d y=%d w=%d h=%d]\n', roi);
end
assert(exist('roi','var')==1 && numel(roi)==4);

lockedCenter = [roi(3)/2, roi(4)/2]; % temp; replaced by REF.center later
roiArea = roi(3)*roi(4);
S1.MIN_PUPIL_AREA = max(S1.MIN_PUPIL_AREA_ABS, round(S1.AREA_MIN_FRAC * roiArea));
S1.MAX_PUPIL_AREA = round(S1.AREA_MAX_FRAC * roiArea);
S1.CENTER_TOL_PX  = round(S1.CENTER_TOL_FRAC * min(roi(3),roi(4)));
fprintf('Area gates: [%d .. %d] px; center tol=%d px\n',S1.MIN_PUPIL_AREA,S1.MAX_PUPIL_AREA,S1.CENTER_TOL_PX);

% =========== POINT SEEDS (2 fixed; 3rd requested later DURING Step-1) ===========
PTS.enable   = true;
PTS.frames   = [];  % [] -> start with 2 frames (10% and 90%)
PTS.order    = {'center','top','right','bottom','left'};

REF = []; lastCenter=[]; lastArea=NaN; startFrameNum=0;

% Difficulty trackers
DIFF.frames = []; DIFF.bestScore = []; DIFF.difficulty = [];
DIFF2.frames = []; DIFF2.bestScore = []; DIFF2.difficulty = [];

% 3rd seed & guide states
THIRD = struct('frame',NaN,'clicked',false);
GUIDE = struct('frame',NaN,'clicked',false);
STATS = struct('processedSinceThird',0,'failsSinceThird',0);

% ---- choose first 2 seed frames (10% and 90%) ----
if PTS.enable
  if isempty(PTS.frames)
    base2 = unique(max(2, round([0.1 0.9]*numFramesApprox)));
    if numel(base2) < 2, base2 = [2, max(3,numFramesApprox-2)]; end
    PTS.frames = base2(:).';
  else
    PTS.frames = unique(PTS.frames);
    if numel(PTS.frames) > 2, PTS.frames = PTS.frames(1:2); end
  end

  % ---- collect points on the 2 selected seed frames ONLY ----
  vr_pts = VideoReader(VIDEO);
  seeds = struct('frame',[],'center',[],'edges',[],'rx',[],'ry',[], ...
                 'meanIn',[],'meanOut',[],'contrast',[],'perGrad',[]);
  for i = 1:numel(PTS.frames)
    f = PTS.frames(i);
    vr_pts.CurrentTime = (f-1)/vr_pts.FrameRate;
    F = readFrame(vr_pts);
    C = imcrop(F, roi);
    if size(C,3)==3, G = rgb2gray(C); else, G=C; end
    N = preprocess_for_step1(G, S1);

    % robust click loop (inside ROI image only)
    [pts,~] = prompt_clicks_inside_roi(N, ...
      sprintf('Pick 5 points (frame %d): center, top, right, bottom, left',f), ...
      {'center','top','right','bottom','left'});

    c = pts(1,:); t = pts(2,:); r = pts(3,:); b = pts(4,:); l = pts(5,:);
    rx = max(4, mean([abs(r(1)-c(1)) abs(c(1)-l(1))]));
    ry = max(4, mean([abs(c(2)-t(2)) abs(b(2)-c(2))]));

    [meanIn, meanOut, ringC, perGrad] = ring_contrast_ellipse(N, c, rx, ry, max(3,S1.RING_WIDTH_PX));

    seeds(i).frame=f; seeds(i).center=c; seeds(i).edges=pts(2:5,:);
    seeds(i).rx=rx; seeds(i).ry=ry; seeds(i).meanIn=meanIn; seeds(i).meanOut=meanOut;
    seeds(i).contrast=ringC; seeds(i).perGrad=perGrad;
    imwrite(N, fullfile(DIRS.seedimg, sprintf('frame_%05d.png', f)));
    save(fullfile(DIRS.seedpts, sprintf('points_%05d.mat', f)),'pts','c','rx','ry','meanIn','meanOut','ringC','perGrad');
  end

  % Build REF from the 2 seeds
  cents    = vertcat(seeds.center);
  rxv      = vertcat(seeds.rx); ryv = vertcat(seeds.ry);
  meanInV  = vertcat(seeds.meanIn); meanOutV = vertcat(seeds.meanOut);
  contrV   = vertcat(seeds.contrast); perGV = vertcat(seeds.perGrad);

  REF.center = mean(cents,1);
  REF.centerStd = std(cents,0,1);
  REF.rxMed = median(rxv); REF.ryMed = median(ryv);
  REF.areaMed = pi*REF.rxMed*REF.ryMed;
  REF.circMed = ellipse_circularity(REF.rxMed, REF.ryMed);
  REF.solMed  = 0.95;
  REF.eccMed  = sqrt(1-(min(REF.rxMed,REF.ryMed)/max(REF.rxMed,REF.ryMed))^2);
  REF.meanInMed = median(meanInV); REF.meanOutMed = median(meanOutV);
  REF.ringContrastMed = median(contrV);
  REF.perimGradMed    = median(perGV);

  if isempty(S1.SEARCH_RADIUS_PX)
    S1.SEARCH_RADIUS_PX = round(1.6 * max(REF.rxMed, REF.ryMed));
  end
  S1.MIN_PUPIL_AREA = max(S1.MIN_PUPIL_AREA, round(S1.MIN_AREA_FRAC_OF_SEED * REF.areaMed));
  S1.MAX_PUPIL_AREA = min(S1.MAX_PUPIL_AREA, round(REF.areaMed * S1.MAX_AREA_REL_TO_SEED));

  lockedCenter = REF.center;
  S1.DELTA_I_MIN = max(S1.DELTA_I_MIN, round(0.50*abs(REF.ringContrastMed)));
  REF.centerTolPx = max(6, ceil(2.5*max([REF.centerStd(:);1])));
  REF.areaTolAbs = max(round(0.25*REF.areaMed), round(0.8*(pi*mad(rxv,1)*mad(ryv,1)+1)));
  REF.contrastTol= max(2, round(0.5*max(1,REF.ringContrastMed)));
  REF.edgeTol    = max(1, 0.60*max(1,REF.perimGradMed));
  REF.shapeTol   = struct('circDrop',0.07,'solDrop',0.04,'eccRise',0.12);

  REF.isDarkPupil = REF.meanInMed < REF.meanOutMed;
  REF.thrMid  = 0.5*(REF.meanInMed + REF.meanOutMed);
  REF.thrBand= max(3, 0.40*abs(REF.meanOutMed - REF.meanInMed));
  REF.radLock= max(8, round(0.5*max(REF.rxMed, REF.ryMed)));

  S1.CENTER_TOL_PX       = max(round(0.18*min(roi(3),roi(4))), REF.centerTolPx);
  S1.JUMP_MAX_PX         = max(S1.JUMP_MAX_PX, 40);
  S1.AREA_CHANGE_MAX_FRAC= max(S1.AREA_CHANGE_MAX_FRAC, 0.90);
  S1.DELTA_I_MIN         = max(S1.DELTA_I_MIN, round(0.35*abs(REF.ringContrastMed)));
  if REF.isDarkPupil, S1.USE_BRIGHT_POL=false; else, S1.USE_BRIGHT_POL=false; end

  save(fullfile(RUN_ROOT,'seed_ref_profile.mat'),'REF');
  write_text(fullfile(RUN_ROOT,'seed_ref_profile.txt'), ...
    sprintf(['POINT SEEDS REF:\n center=[%.1f %.1f]±[%.1f %.1f] px\n rx=%.1f, ry=%.1f (area≈%.0f)\n' ...
    ' meanIn=%.1f, meanOut=%.1f, ringContrast=%.1f, perimGrad=%.2f\n' ...
    ' TOL: centerTolPx=%g, areaTolAbs=%g, contrastTol=%g, edgeTol=%.2f\n'], ...
    REF.center, REF.centerStd, REF.rxMed, REF.ryMed, REF.areaMed, ...
    REF.meanInMed, REF.meanOutMed, REF.ringContrastMed, REF.perimGradMed, ...
    REF.centerTolPx, REF.areaTolAbs, REF.contrastTol, REF.edgeTol));

  % start just after the 2nd seed frame
  PTS.frames = sort(PTS.frames);
  lastCenter = seeds(2).center;
  lastArea   = pi*seeds(2).rx*seeds(2).ry;
  vr.CurrentTime = (PTS.frames(2)) / vr.FrameRate;
  startFrameNum  = PTS.frames(2);
end

% ===== third-seed force timing =====
if isempty(S1.THIRD.FORCE_AFTER_SEEN)
  THIRD_FORCE_AFTER_SEEN = max(200, round(0.15 * S1.STOP_AFTER_PAIRS));
else
  THIRD_FORCE_AFTER_SEEN = max(1, S1.THIRD.FORCE_AFTER_SEEN);
end

% =================== STEP-1: PSEUDO-LABELS ===================
vr = VideoReader(VIDEO);
resumeFrom = max(0, (startFrameNum-1)/vr.FrameRate);
vr.CurrentTime = resumeFrom;
frameNum = max(0, startFrameNum-1);
lastSavedGray = [];

seOpen = []; if S1.OPEN_RADIUS>0,  seOpen = strel('disk',S1.OPEN_RADIUS);  end
seClose= []; if S1.CLOSE_RADIUS>0, seClose= strel('disk',S1.CLOSE_RADIUS); end

counts = struct('noComp',0,'shapeFail',0,'contrastFail',0,'temporalFail',0,'borderFail',0,'saved',0);
framesSaved=0; skipLogFile = fullfile(DIRS.audit,'step1_skip_log.txt');
framesSeen = 0;

t0 = tic;
while hasFrame(vr)
  if framesSaved >= S1.STOP_AFTER_PAIRS, break; end
  frameNum = frameNum + 1;
  F = readFrame(vr);
  if mod(frameNum, S1.SAVE_EVERY_N_FRAMES) ~= 0, continue; end
  framesSeen = framesSeen + 1;

  crop = imcrop(F, roi); if isempty(crop), continue; end
  if size(crop,3)==3, G = rgb2gray(crop); else, G=crop; end
  N = preprocess_for_step1(G, S1);

  % global polarity backups (early for scoring)
  pLow = max(5, min(70, S1.PERCENTILE_DARK+3));
  bwDark = N <= uint8(prctile(double(N(:)), pLow));
  pHi = min(99, max(50, 100 - S1.BRIGHT_PCT + 5));
  bwBright= N >= uint8(prctile(double(N(:)), pHi));

  if ~isempty(seOpen),  bwDark = imopen(bwDark, seOpen);   bwBright = imopen(bwBright, seOpen);  end
  if ~isempty(seClose), bwDark = imclose(bwDark, seClose); bwBright = imclose(bwBright, seClose); end
  bwDark = imfill(bwDark,'holes'); bwBright = imfill(bwBright,'holes');

  % local seed window + polarity anchor
  seedMask = false(size(N)); seedWin = false(size(N));
  if ~isempty(REF)
    cx = round(REF.center(1)); cy = round(REF.center(2));
    rad = round(REF.radLock);
    yy = max(1,cy-rad):min(size(N,1), cy+rad);
    xx = max(1,cx-rad):min(size(N,2), cx+rad);
    seedWin(yy,xx) = true;

    if REF.isDarkPupil
      Tlo = max(0, REF.thrMid - REF.thrBand);
      seedMask = (N <= uint8(Tlo)) & seedWin;
    else
      Thi = min(255, REF.thrMid + REF.thrBand);
      seedMask = (N >= uint8(Thi)) & seedWin;
    end
    seedMask = imfill(imopen(seedMask, strel('disk',1)),'holes');
  end

  if ~isempty(REF)
    if     REF.isDarkPupil, bwDark   = bwDark   | seedMask;
    else                    bwBright = bwBright | seedMask; end
  end

  % soft radial mask around lockedCenter
  [W,H] = deal(size(N,2), size(N,1));
  [X,Y] = meshgrid(1:W, 1:H);
  Rmask = ((X - lockedCenter(1)).^2 + (Y - lockedCenter(2)).^2) <= (S1.SEARCH_RADIUS_PX + 3)^2;


  % candidates
  [MW, STW, SCW] = watershed_pupil_candidate(N, REF, S1, lockedCenter);
  [M1,ST1,SC1] = pick_best_by_contrast_scored(bwDark,  true,  N, S1, lockedCenter);
  [M2,ST2,SC2] = pick_best_by_contrast_scored(bwBright,false, N, S1, lockedCenter);

  scoresVec = [ tern(isempty(M1),-inf,SC1), tern(isempty(M2),-inf,SC2), tern(isempty(MW),-inf,SCW) ];
  bestScore = max(scoresVec);
  diffMetric = ifelse(isfinite(bestScore), max(0,-bestScore), 1e9);
  DIFF.frames(end+1)     = frameNum;
  DIFF.bestScore(end+1)  = bestScore;
  DIFF.difficulty(end+1) = diffMetric;

  if DEBUG.ENABLE && DEBUG.SAVE_ALL_CANDS
    if ~isempty(M1), debug_dump(N, M1, ST1, 'cand_dark',     SC1, frameNum, DEBUG.OUT); end
    if ~isempty(M2), debug_dump(N, M2, ST2, 'cand_bright',   SC2, frameNum, DEBUG.OUT); end
    if ~isempty(MW), debug_dump(N, MW, STW, 'cand_watershed',SCW, frameNum, DEBUG.OUT); end
  end

  cand = {M1, ST1, SC1; M2, ST2, SC2; MW, STW, SCW};
  [M, ST, ~] = pick_best_candidate(cand, S1);

  savedThisFrame = false;

  if isempty(M)
    counts.noComp = counts.noComp + 1; log_skip(skipLogFile, frameNum, 'noComp');
    if DEBUG.ENABLE, debug_dump(N, false(size(N)), [], 'noComp', -inf, frameNum, DEBUG.OUT); end
  else
    if ~isempty(seClose) && S1.CLOSE_RADIUS>0, M = imclose(M,seClose); end
    M = imfill(M,'holes');
    [M, ST] = refine_shrink_to_edge(N, M, ST, S1);

    % temporal/diversity gates (unchanged)
    doSkip = false;
    if ~isempty(lastCenter)
      d = hypot(ST.Centroid(1)-lastCenter(1), ST.Centroid(2)-lastCenter(2));
      okJ = d <= S1.JUMP_MAX_PX;
      okA = isnan(lastArea) || abs(ST.Area-lastArea) <= S1.AREA_CHANGE_MAX_FRAC*max(1,lastArea);
      if ~(okJ && okA)
        counts.temporalFail=counts.temporalFail+1; log_skip(skipLogFile, frameNum, 'temporalFail');
        if DEBUG.ENABLE, debug_dump(N, M, ST, sprintf('temporal d=%.1f okJ=%d okA=%d',d,okJ,okA), 0, frameNum, DEBUG.OUT); end
        doSkip = true;
      end
    else
      dctr = hypot(ST.Centroid(1)-REF.center(1), ST.Centroid(2)-REF.center(2));
      hardTol = max(8, min(REF.centerTolPx, round(0.25*min(roi(3),roi(4)))));
      if dctr > hardTol
        counts.shapeFail=counts.shapeFail+1; log_skip(skipLogFile, frameNum, 'firstSave_centerFail');
        if DEBUG.ENABLE, debug_dump(N, M, ST, sprintf('first_center d=%.1f tol=%.1f',dctr,hardTol), 0, frameNum, DEBUG.OUT); end
        doSkip = true;
      end
    end

    if ~doSkip && framesSaved >= 5
      currGray = rgbOrGray(crop);
      if isempty(lastSavedGray), madv = Inf;
      else
        if ~isequal(size(currGray), size(lastSavedGray)), currGray = imresize(currGray, size(lastSavedGray)); end
        madv = median(abs(double(currGray(:)) - double(lastSavedGray(:))));
      end
      if ~isnan(lastArea)
        if abs(ST.Area-lastArea)/max(1,lastArea) < S1.DIVERSITY_MIN_DA && madv < S1.MIN_MEDIAN_ABS_DIFF
          lastCenter = ST.Centroid; lastArea = ST.Area;
          log_skip(skipLogFile, frameNum, 'diversitySkip');
          if DEBUG.ENABLE, debug_dump(N, M, ST, sprintf('diversity madv=%.2f',madv), 0, frameNum, DEBUG.OUT); end
          doSkip = true;
        end
      end
    end

    if ~doSkip
      id = frameNum;
      imwrite(N, fullfile(DIRS.img, sprintf('frame_%05d.png', id)));
      imwrite(uint8(M)*255, fullfile(DIRS.lab, sprintf('frame_%05d_mask.png', id)));
      framesSaved = framesSaved + 1; counts.saved = counts.saved + 1;
      savedThisFrame = true;

      lastSavedGray = rgbOrGray(crop);
      lastCenter = ST.Centroid; lastArea = ST.Area;

      if mod(framesSaved, S1.QC_EVERY_SAVED)==0
        ov = repmat(N,[1 1 3]); ed = bwperim(M);
        ov(:,:,2) = max(ov(:,:,2), uint8(255*M));
        ov(:,:,1) = max(ov(:,:,1), uint8(255*ed));
        c = round(ST.Centroid);
        rr = max(1,c(2)-2):min(size(ov,1),c(2)+2);
        cc = max(1,c(1)-2):min(size(ov,2),c(1)+2);
        ov(rr,cc,1) = 255; ov(rr,cc,2) = 0; ov(rr,cc,3) = 0;
        imwrite(ov, fullfile(DIRS.qc, sprintf('qc_%05d.png', id)));
      end
    end
  end

  % relax if nothing saved by 200 frames
  if framesSaved==0 && frameNum > 200 && frameNum <= 230
    S1.DELTA_I_MIN = 2; S1.MIN_CIRCULARITY = 0.45; S1.MAX_ECCENTRICITY=0.97;
    S1.JUMP_MAX_PX = 60; S1.CENTER_TOL_PX = round(0.35*min(roi(3),roi(4))); S1.BORDER_MARGIN_PX = 1;
  end

  % ---------- MID-STEP1 3rd SEED: auto OR forced (ALWAYS_ASK) ----------
  if S1.THIRD.ENABLE && ~THIRD.clicked
    failsTotal = counts.noComp + counts.shapeFail + counts.contrastFail + counts.temporalFail + counts.borderFail;
    failRate   = failsTotal / max(1, framesSeen);

    autoBad = framesSeen >= S1.THIRD.MIN_PROCESSED && ...
              failsTotal >= S1.THIRD.MIN_FAILS && ...
              failRate   >= S1.THIRD.MIN_FAIL_RATE;

    forceNow = S1.THIRD.ALWAYS_ASK && (framesSeen >= THIRD_FORCE_AFTER_SEEN);

    if autoBad || forceNow
      [~,ixHard] = max(DIFF.difficulty);
      hardFrame = DIFF.frames(max(1,ixHard));
      THIRD.frame = hardFrame;

      try
        vr_mid = VideoReader(VIDEO);
        Nm = local_read_norm_roi(vr_mid, hardFrame, roi, S1);
        if isempty(Nm), error('Could not read ROI for frame %d', hardFrame); end

        [pts,~] = prompt_clicks_inside_roi(Nm, ...
          sprintf('3rd seed @ hardest frame %d: center, top, right, bottom, left',hardFrame), ...
          {'center','top','right','bottom','left'});

        c = pts(1,:); t = pts(2,:); r = pts(3,:); b = pts(4,:); l = pts(5,:);
        rx = max(4, mean([abs(r(1)-c(1)) abs(c(1)-l(1))]));
        ry = max(4, mean([abs(c(2)-t(2)) abs(b(2)-c(2))]));

        [meanIn, meanOut, ringC, perGrad] = ring_contrast_ellipse(Nm, c, rx, ry, max(3,S1.RING_WIDTH_PX));

        k = numel(seeds)+1;
        seeds(k).frame=hardFrame; seeds(k).center=c; seeds(k).edges=pts(2:5,:);
        seeds(k).rx=rx; seeds(k).ry=ry; seeds(k).meanIn=meanIn; seeds(k).meanOut=meanOut;
        seeds(k).contrast=ringC; seeds(k).perGrad=perGrad;

        imwrite(Nm, fullfile(DIRS.seedimg, sprintf('frame_%05d.png', hardFrame)));
        save(fullfile(DIRS.seedpts, sprintf('points_%05d.mat', hardFrame)), ...
             'pts','c','rx','ry','meanIn','meanOut','ringC','perGrad');

        % update REF from all seeds (now 3)
        cents   = vertcat(seeds.center);
        rxv     = vertcat(seeds.rx);      ryv = vertcat(seeds.ry);
        meanInV = vertcat(seeds.meanIn);  meanOutV = vertcat(seeds.meanOut);
        contrV  = vertcat(seeds.contrast); perGV = vertcat(seeds.perGrad);

        REF.center = mean(cents,1);
        REF.centerStd = std(cents,0,1);
        REF.rxMed = median(rxv); REF.ryMed = median(ryv);
        REF.areaMed = pi*REF.rxMed*REF.ryMed;
        REF.circMed = ellipse_circularity(REF.rxMed, REF.ryMed);
        REF.solMed  = 0.95;
        REF.eccMed  = sqrt(1-(min(REF.rxMed,REF.ryMed)/max(REF.rxMed,REF.ryMed))^2);
        REF.meanInMed = median(meanInV); REF.meanOutMed = median(meanOutV);
        REF.ringContrastMed = median(contrV);
        REF.perimGradMed    = median(perGV);

        if isempty(S1.SEARCH_RADIUS_PX), S1.SEARCH_RADIUS_PX = round(1.6 * max(REF.rxMed, REF.ryMed)); end
        S1.MIN_PUPIL_AREA = max(S1.MIN_PUPIL_AREA, round(S1.MIN_AREA_FRAC_OF_SEED * REF.areaMed));
        S1.MAX_PUPIL_AREA = min(S1.MAX_PUPIL_AREA, round(REF.areaMed * S1.MAX_AREA_REL_TO_SEED));
        S1.DELTA_I_MIN    = max(S1.DELTA_I_MIN, round(0.50*abs(REF.ringContrastMed)));

        REF.centerTolPx = max(6, ceil(2.5*max([REF.centerStd(:);1])));
        REF.areaTolAbs  = max(round(0.25*REF.areaMed), round(0.8*(pi*mad(rxv,1)*mad(ryv,1)+1)));
        REF.contrastTol = max(2, round(0.5*max(1,REF.ringContrastMed)));
        REF.edgeTol     = max(1, 0.60*max(1,REF.perimGradMed));
        REF.shapeTol    = struct('circDrop',0.07,'solDrop',0.04,'eccRise',0.12);

        REF.isDarkPupil = REF.meanInMed < REF.meanOutMed;
        REF.thrMid  = 0.5*(REF.meanInMed + REF.meanOutMed);
        REF.thrBand = max(3, 0.40*abs(REF.meanOutMed - REF.meanInMed));
        REF.radLock = max(8, round(0.5*max(REF.rxMed, REF.ryMed)));

        lockedCenter = REF.center;
        S1.CENTER_TOL_PX        = max(round(0.18*min(roi(3),roi(4))), REF.centerTolPx);
        S1.DELTA_I_MIN          = max(S1.DELTA_I_MIN, round(0.35*abs(REF.ringContrastMed)));
        S1.USE_BRIGHT_POL       = false;

        save(fullfile(RUN_ROOT,'seed_ref_profile.mat'),'REF');
        write_text(fullfile(RUN_ROOT,'seed_ref_profile.txt'), 'UPDATED with 3rd seed');

        THIRD.clicked = true;
        fprintf('3rd seed collected at hardest frame %d; REF updated.\n', hardFrame);
      catch ME
        warning('3rd-seed prompt failed: %s', ME.message);
        THIRD.clicked = true; % avoid re-prompt loop
      end
    end
  end

  % ===== post-3rd-seed stats & hardest-since-3rd for the GUIDE =====
  if THIRD.clicked
    STATS.processedSinceThird = STATS.processedSinceThird + 1;
    if ~savedThisFrame, STATS.failsSinceThird = STATS.failsSinceThird + 1; end
    DIFF2.frames(end+1)     = frameNum;
    DIFF2.bestScore(end+1)  = bestScore;
    DIFF2.difficulty(end+1) = diffMetric;
  end

  % ---------- EXTRA GUIDE (inside/outside) ----------
  if S1.GUIDE.ENABLE && THIRD.clicked && ~GUIDE.clicked
    failRate2 = STATS.failsSinceThird / max(1, STATS.processedSinceThird);
    if STATS.processedSinceThird >= S1.GUIDE.MIN_AFTER_THIRD_FRAMES && ...
       STATS.failsSinceThird    >= S1.GUIDE.MIN_AFTER_THIRD_FAILS  && ...
       failRate2                >= S1.GUIDE.MIN_AFTER_THIRD_RATE
      try
        [~,ix2] = max(DIFF2.difficulty);
        hard2 = DIFF2.frames(ix2);
        GUIDE.frame = hard2;

        vr_g = VideoReader(VIDEO);
        Ng = local_read_norm_roi(vr_g, hard2, roi, S1);
        if isempty(Ng), error('Could not read ROI for frame %d', hard2); end

        % 1) inside; 2) just outside
        pIn  = prompt_single_click_inside_roi(Ng, 'Click ONE pixel INSIDE the pupil (inside this ROI view)');
        pOut = prompt_single_click_inside_roi(Ng, 'Click ONE pixel just OUTSIDE the pupil (iris/sclera), inside ROI');

        radS = max(2, S1.GUIDE.DISK_RADIUS_PX);
        mIn  = sample_mean_disk(Ng, pIn,  radS);
        mOut = sample_mean_disk(Ng, pOut, radS);

        REF.isDarkPupil = (mIn < mOut);
        REF.meanInMed   = mIn; REF.meanOutMed = mOut;
        REF.thrMid      = 0.5*(mIn + mOut);
        REF.thrBand     = max(3, S1.GUIDE.BAND_FRAC * abs(mOut - mIn));
        S1.DELTA_I_MIN  = max(S1.DELTA_I_MIN, round(0.50*abs(mOut - mIn)));
        S1.USE_BRIGHT_POL = false;

        save(fullfile(RUN_ROOT,'seed_ref_profile.mat'),'REF');
        write_text(fullfile(RUN_ROOT,'guide_update.txt'), ...
          sprintf('GUIDE @ frame %d: mIn=%.2f, mOut=%.2f\n', hard2, mIn, mOut));
        save(fullfile(DIRS.seedpts,'guide_points.mat'),'hard2','pIn','pOut','mIn','mOut','radS');

        GUIDE.clicked = true;
        fprintf('GUIDE collected at frame %d; USE_BRIGHT_POpolarity/thresholds updated.\n', hard2);
      catch ME
        warning('GUIDE prompt failed: %s', ME.message);
        GUIDE.clicked = true;
      end
    end
  end

  if DEBUG.ENABLE && mod(frameNum, DEBUG.PRINT_EVERY)==0
    fprintf('[dbg] frame=%d saved=%d\n', frameNum, framesSaved);
  end
end
fprintf('STEP-1 main pass: saved %d pairs in %.1f s\n', framesSaved, toc(t0));

% ----- STEP-1 AUDIT -----
try
  ref = []; rf = fullfile(RUN_ROOT,'seed_ref_profile.mat');
  if exist(rf,'file'), S = load(rf); if isfield(S,'REF'), ref = S.REF; end, end
  [T1, summary1] = audit_step1(DIRS.img, DIRS.lab, DIRS.audit, ref);
  if ~exist(DIRS.audit,'dir'), mkdir(DIRS.audit); end
  writetable(T1, fullfile(DIRS.audit,'step1_metrics_full.csv'));
  write_text(fullfile(DIRS.audit,'step1_summary.txt'), summary1);
  fprintf('%s\n', summary1);
catch ME
  warning('audit_step1 failed: %s', ME.message);
end

% ----- FINAL QC GALLERY -----
try, save_qc_gallery(DIRS.img, DIRS.lab, DIRS.qc, S1.QC_GALLERY_MAX);
catch ME, warning('QC gallery failed: %s', ME.message); end

% ----- SAVE STEP-1 CONFIG + COUNTS -----
try
  save(fullfile(RUN_ROOT,'step1_config.mat'), 'S1','roi','RUN_ROOT','DIRS','counts');
  fid = fopen(fullfile(RUN_ROOT,'step1_counts.txt'),'w');
  fprintf(fid,'noComp=%d, shapeFail=%d, contrastFail=%d, temporalFail=%d, borderFail=%d, saved=%d\n', ...
    counts.noComp,counts.shapeFail,counts.contrastFail,counts.temporalFail,counts.borderFail,counts.saved);
  fclose(fid);
catch ME
  warning('Could not save step1 config: %s', ME.message);
end

% ===================== STEP-2: TRAINING =====================
% (unchanged from your version)
imgList = dir(fullfile(DIRS.img, 'frame_*.png'));
labList = dir(fullfile(DIRS.lab, 'frame_*_mask.png'));
getID = @(s) sscanf(s, 'frame_%d');
imgID = arrayfun(@(d) getID(d.name), imgList);
labID = arrayfun(@(d) getID(d.name), labList);
[ids, ia, ib] = intersect(imgID, labID);
[ids, order] = sort(ids); ia = ia(order); ib = ib(order);

if isempty(ids)
  warning('Step-1 yielded 0 pairs. Fallback: synthesize from point seeds.');
  vr_tmp = VideoReader(VIDEO);

  seedFramesAvail = unique(PTS.frames);
  if isfinite(THIRD.frame) && ~isnan(THIRD.frame)
    if exist(fullfile(DIRS.seedpts, sprintf('points_%05d.mat', THIRD.frame)),'file')
      seedFramesAvail = unique([seedFramesAvail, THIRD.frame]);
    end
  end

  for i=1:min(3, numel(seedFramesAvail))
    f = seedFramesAvail(i);
    vr_tmp.CurrentTime = (f-1)/vr_tmp.FrameRate;
    F = readFrame(vr_tmp);
    C = imcrop(F, roi);
    if size(C,3)==3, G = rgb2gray(C); else, G=C; end
    N = preprocess_for_step1(G, S1);
    Sps = load(fullfile(DIRS.seedpts, sprintf('points_%05d.mat', f)));
    M = ellipse_mask(size(N), Sps.c, Sps.rx, Sps.ry);
    imwrite(N, fullfile(DIRS.img, sprintf('frame_%05d.png', f)));
    imwrite(uint8(M)*255, fullfile(DIRS.lab, sprintf('frame_%05d_mask.png', f)));
  end
  imgList = dir(fullfile(DIRS.img, 'frame_*.png'));
  labList = dir(fullfile(DIRS.lab, 'frame_*_mask.png'));
  imgID = arrayfun(@(d) getID(d.name), imgList);
  labID = arrayfun(@(d) getID(d.name), labList);
  [ids, ia, ib] = intersect(imgID, labID);
  [ids, order] = sort(ids); ia = ia(order); ib = ib(order);
end
assert(~isempty(ids), 'No (image,mask) pairs found.');

imgFiles = fullfile({imgList(ia).folder}, {imgList(ia).name});
labFiles = fullfile({labList(ib).folder}, {labList(ib).name});
Npairs = numel(imgFiles);

% ensure masks are {0,255}
for k = 1:Npairs
  M = imread(labFiles{k}); if size(M,3)>1, M = rgb2gray(M); end
  M = uint8(M>0)*255; imwrite(M, labFiles{k});
end

imds = imageDatastore(imgFiles, 'ReadFcn', @(f) im2uint8(imread(f)));
classNames = ["background","pupil"]; labelIDs = [0 255];
pxds = pixelLabelDatastore(labFiles, classNames, labelIDs);

ord = 1:Npairs;
valCount = min( max(1, floor(0.25*Npairs)), max(1, Npairs-2) );
if Npairs <= 6, valCount = max(1, min(2, Npairs-2)); end
va = unique(round(linspace(1, Npairs, valCount)));
guard = (Npairs <= 12) * 2 + (Npairs > 12) * 1;
mask = true(1, Npairs);
for v = va
  lo = max(1, v-guard); hi = min(Npairs, v+guard);
  mask(lo:hi) = false; mask(v) = true;
end
va = find(~mask | ismember(1:Npairs, va));
tr = setdiff(ord, va, 'stable');
while numel(tr) < max(2, round(0.5*Npairs)) && numel(va) > 1
  va = va(1:end-1); tr = setdiff(ord, va, 'stable');
end
tr = unique(tr,'stable'); va = setdiff(unique(va,'stable'), tr,'stable');
fprintf('TRAIN/VAL split -> train=%d, val=%d (N=%d)\n', numel(tr), numel(va), Npairs);

imdsTrain = subset(imds, tr); pxdsTrain = subset(pxds, tr);
imdsVal   = subset(imds, va); pxdsVal   = subset(pxds, va);

Hroi = roi(4); Wroi = roi(3); blk = 2^S2.encoderDepth;
padH = ceil(Hroi / blk)*blk; padW = ceil(Wroi / blk)*blk;
if any([padH padW] ~= [Hroi Wroi])
  fprintf('U-Net input padded to [%d %d] (from ROI [%d %d])\n', padH, padW, Hroi, Wroi);
end

augCfg.padH=padH; augCfg.padW=padW; augCfg.classNames=classNames;
dsTrain = pixelLabelImageDatastore(imdsTrain, pxdsTrain, 'ColorPreprocessing','none');
dsVal   = pixelLabelImageDatastore(imdsVal,   pxdsVal,   'ColorPreprocessing','none');
dsTrain = transform(dsTrain, @(d) pad_and_augment_roi(d, augCfg));
dsVal   = transform(dsVal,   @(d) pad_and_augment_roi(d, setfield(augCfg,'doAug',false)));

imageSize = [padH padW 1];
lgraph = unetLayers(imageSize, numel(classNames), ...
  'EncoderDepth', S2.encoderDepth, 'NumFirstEncoderFilters', S2.numFirstFilters);

isPx = arrayfun(@(L) isa(L,'nnet.cnn.layer.PixelClassificationLayer') || ...
  contains(class(L),'PixelClassification','IgnoreCase',true), lgraph.Layers);
pxIdx = find(isPx,1,'last'); pxName = lgraph.Layers(pxIdx).Name;

useDice=false;
tbl = countEachLabel(pxdsTrain);
freq = tbl.PixelCount / sum(tbl.PixelCount);
w = 1 ./ sqrt(max(freq,1e-6)); w = w / mean(w); w = min(w, S2.weightClipMax);
try
  diceL = dicePixelClassificationLayer('Name','diceLoss','Classes',tbl.Name,'ClassWeights',w);
  lgraph = replaceLayer(lgraph, pxName, diceL); useDice=true;
catch
  wce = pixelClassificationLayer('Name','wce','Classes',tbl.Name,'ClassWeights',w);
  lgraph = replaceLayer(lgraph, pxName, wce); useDice=false;
end

ckDir = fullfile(DIRS.model, ['ckpt_' STAMP]); if ~exist(ckDir,'dir'), mkdir(ckDir); end
opts = trainingOptions('adam', ...
  'InitialLearnRate', S2.baseLR, 'MaxEpochs', S2.maxEpochs, 'MiniBatchSize', S2.miniBatch, ...
  'Shuffle','every-epoch', 'ValidationData', dsVal, 'ValidationFrequency', 20, ...
  'ValidationPatience', 12, 'LearnRateSchedule','piecewise', ...
  'LearnRateDropPeriod', round(S2.maxEpochs/2), 'LearnRateDropFactor', 0.2, ...
  'GradientThresholdMethod','l2norm','GradientThreshold',1.0, ...
  'L2Regularization',5e-4,'ExecutionEnvironment','auto', ...
  'CheckpointPath', ckDir, 'Verbose', true, 'Plots','none');

fprintf('Step-2: depth=%d filters=%d input=[%d %d] epochs=%d (loss=%s)\n', ...
  S2.encoderDepth, S2.numFirstFilters, padH, padW, S2.maxEpochs, tern(useDice,'Dice','WCE'));
[net, info] = trainNetwork(dsTrain, lgraph, opts);

bestNet = net; bestDice = -inf;
if exist(ckDir,'dir')
  files = dir(fullfile(ckDir,'*.mat'));
  for k = 1:numel(files)
    S = load(fullfile(ckDir, files(k).name));
    if isfield(S,'net')
      try
        [Ttmp, ~] = quick_validate_dice(S.net, imdsVal, pxdsVal, [padH padW], min(20,numel(va)));
        medD = median(Ttmp.dice,'omitnan');
        if medD > bestDice, bestDice = medD; bestNet = S.net; end
      catch
      end
    end
  end
end
net = bestNet;
modelFileTS = fullfile(DIRS.model, sprintf('trainedPupilNet_%s.mat', STAMP));
save(modelFileTS,'net');
copyfile(modelFileTS, fullfile(RUN_ROOT,'trainedPupilNet.mat'));

[valTbl, valSummary] = quick_validate_dice(net, imdsVal, pxdsVal, [padH padW], min(20,numel(va)));
writetable(valTbl, fullfile(RUN_ROOT,'step2_val_metrics.csv'));
write_text(fullfile(RUN_ROOT,'step2_val_summary.txt'), valSummary);
try, Tsum = struct2table(info); writetable(Tsum, fullfile(DIRS.model, ['step2_summary_' STAMP '.csv'])); end

fprintf('\n=== ALL DONE ===\nOutputs under: %s\n', RUN_ROOT);

% ====================== LOCAL FUNCTIONS ======================

function N = preprocess_for_step1(G, S1)
if ~isa(G,'uint8'), G = im2uint8(mat2gray(G)); end
Gd = G;
if S1.SUPPRESS_GLINT
  Gl = Gd >= prctile(double(Gd(:)), 98);
  Gl = imdilate(Gl, strel('disk',2));
  if any(Gl(:))
    medVal = median(double(Gd(~Gl)),'omitnan');
    Gd(Gl) = uint8(medVal);
  end
end
if S1.USE_CLAHE
  N = adapthisteq(Gd);
elseif S1.USE_STRETCHLIM
  N = imadjust(Gd, stretchlim(Gd, S1.ADAPTIVE_CLIP), []);
else
  N = Gd;
end
end

function Iu = rgbOrGray(C)
if size(C,3)==3, C = rgb2gray(C); end
if ~isa(C,'uint8'), Iu = im2uint8(C); else, Iu = C; end
end

function [pts, ax] = prompt_clicks_inside_roi(N, winTitle, labels)
figure('Color','w','Name',winTitle); ax = axes; imshow(N,'Parent',ax,'InitialMagnification','fit'); hold(ax,'on');
W = size(N,2); H = size(N,1);
pts = nan(numel(labels),2);
for k=1:numel(labels)
  title(ax, sprintf('Click %s point (inside this ROI image only)', labels{k}));
  pts(k,:) = get_point_inside(ax, W, H);
  plot(ax, pts(k,1), pts(k,2), 'r+','MarkerSize',10,'LineWidth',1.5);
  text(pts(k,1)+4, pts(k,2), labels{k}, 'Color','y','FontSize',9,'FontWeight','bold','Parent',ax);
end
hold(ax,'off'); drawnow; close(gcf);
end

function p = prompt_single_click_inside_roi(N, msg)
figure('Color','w','Name',msg); ax = axes; imshow(N,'Parent',ax,'InitialMagnification','fit'); hold(ax,'on');
W = size(N,2); H = size(N,1);
title(ax, sprintf('%s', msg));
p = get_point_inside(ax, W, H);
plot(ax, p(1), p(2), 'g+','MarkerSize',10,'LineWidth',1.5);
hold(ax,'off'); drawnow; close(gcf);
end

function p = get_point_inside(ax, W, H)
while true
  if exist('drawpoint','file')
    h = drawpoint(ax); wait(h); p = h.Position; try delete(h); end
  else
    [x,y] = ginput(1); p=[x y];
  end
  if isfinite(p(1)) && isfinite(p(2)) && p(1)>=1 && p(1)<=W && p(2)>=1 && p(2)<=H
    break;
  else
    title(ax, 'Please click INSIDE the ROI image'); drawnow;
  end
end
end

function [meanIn, meanOut, ringC, perGrad] = ring_contrast_ellipse(N, c, rx, ry, ringW)
[X,Y] = meshgrid(1:size(N,2), 1:size(N,1));
E = ((X-c(1)).^2)/(rx^2) + ((Y-c(2)).^2)/(ry^2) <= 1;
D_in = bwdist(~E); D_out = bwdist(E);
rin = E & D_in >= 1 & D_in <= ringW;
rout = ~E & D_out >= 1 & D_out <= ringW;
Ii = double(N(rin)); Io = double(N(rout));
meanIn = mean(Ii,'omitnan'); meanOut = mean(Io,'omitnan');
ringC = meanOut - meanIn; % positive => dark pupil
per = bwperim(E);
[Gmag,~] = imgradient(N);
perGrad = mean(Gmag(per),'omitnan');
end

function mu = sample_mean_disk(I, p, rad)
[h,w] = size(I);
[xg,yg] = meshgrid(1:w,1:h);
mask = (xg - p(1)).^2 + (yg - p(2)).^2 <= rad^2;
vals = double(I(mask));
mu = mean(vals,'omitnan');
end

function M = ellipse_mask(sz, c, rx, ry)
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
M = ((X-c(1)).^2)/(rx^2) + ((Y-c(2)).^2)/(ry^2) <= 1;
end

function circ = ellipse_circularity(rx, ry)
a = max(rx,ry); b = min(rx,ry);
A = pi*a*b;
h = ((a-b)^2)/((a+b)^2);
P = pi*(a+b)*(1 + (3*h)/(10+sqrt(4-3*h)));
circ = 4*pi*A/(P*P);
end

function [bestMask,bestProps,bestScore] = pick_best_by_contrast_scored(BW, isDarkPol, N, S1, lockedCenter)
bestMask=[]; bestProps=[]; bestScore=-inf;
CC=bwconncomp(BW); if CC.NumObjects==0, return; end
st=regionprops(CC,'Area','Perimeter','Solidity','Eccentricity','Centroid','PixelIdxList');
for k=1:numel(st)
  A = st(k).Area; P = max(st(k).Perimeter,1);
  circ = 4*pi*A/(P*P); sld = st(k).Solidity; ecc = st(k).Eccentricity;
  if A<S1.MIN_PUPIL_AREA || A>S1.MAX_PUPIL_AREA, continue; end
  if circ<S1.MIN_CIRCULARITY || sld<S1.MIN_SOLIDITY || ecc>S1.MAX_ECCENTRICITY, continue; end
  Mtmp=false(size(BW)); Mtmp(st(k).PixelIdxList)=true;
  if touches_border_mask(Mtmp,S1.BORDER_MARGIN_PX), continue; end

  cxy = st(k).Centroid;
  dctr = hypot(cxy(1)-lockedCenter(1), cxy(2)-lockedCenter(2));
  if dctr > S1.SEARCH_RADIUS_PX, continue; end

  D_in = bwdist(~Mtmp); D_out = bwdist(Mtmp);
  rin = Mtmp & D_in >= 1 & D_in <= S1.RING_WIDTH_PX;
  rout = ~Mtmp & D_out>= 1 & D_out<= S1.RING_WIDTH_PX;
  if ~any(rin(:)) || ~any(rout(:)), continue; end
  Iin = double(N(rin)); Iout = double(N(rout));
  if isDarkPol
    dI = mean(Iout) - mean(Iin);
    if dI < S1.DELTA_I_MIN, continue; end
    thrBright = prctile(double(N(:)), 95);
    glintFrac = mean(Iin > thrBright);
    if glintFrac > S1.GLINT_PCT_MAX, continue; end
  else
    dI = mean(Iin) - mean(Iout);
    if dI < S1.DELTA_I_MIN, continue; end
  end

  tol = max(1, scalar_tol(S1.CENTER_TOL_PX));
  distPenalty = (dctr / tol);
  score = 3.5*dI + 45*circ + 12*sld - 15*distPenalty - 8*ecc;
  if score > bestScore, bestScore = score; bestMask = Mtmp; bestProps = st(k); end
end
end

function [MW, STW, score] = watershed_pupil_candidate(N, REF, S1, lockedCenter)
MW = []; STW = []; score = -inf;
if isempty(REF) || isempty(REF.center) || ~all(isfinite(REF.center)), return; end
Ng = imgaussfilt(N, 1.2);
[Gmag,~] = imgradient(Ng);

H = size(N,1); W = size(N,2);
c = REF.center; c = max([1 1], min([W H], c));
radIn  = max(8, round(0.7*max(REF.rxMed, REF.ryMed)));
radOut = min(max(H,W), round(1.6*max(REF.rxMed, REF.ryMed)));

[X,Y] = meshgrid(1:W, 1:H);
M_in   = ((X-c(1)).^2 + (Y-c(2)).^2) <= radIn^2;
M_ring = ((X-c(1)).^2 + (Y-c(2)).^2) <= radOut^2 & ~M_in;

markers = false(H,W);
markers(M_in) = true;
markers = markers | bwperim(M_ring);

G2 = imimposemin(Gmag, markers);
L2 = watershed(G2);

cxi = max(1,min(W,round(c(1)))); cyi = max(1,min(H,round(c(2))));
labAtCtr = L2(cyi,cxi);
if labAtCtr==0
  nb = L2(max(1,cyi-2):min(H,cyi+2), max(1,cxi-2):min(W,cxi+2));
  labs = nb(nb>0);
  if isempty(labs), return; end
  labAtCtr = mode(labs(:));
end
MW = L2==labAtCtr;
MW = imfill(imopen(MW, strel('disk',1)), 'holes');

CC=bwconncomp(MW); if CC.NumObjects==0, MW=[]; return; end
st=regionprops(CC,'Area','Perimeter','Solidity','Eccentricity','Centroid','PixelIdxList');
[~,kk]=max([st.Area]); ST=st(kk);
A = ST.Area; P = max(ST.Perimeter,1);
circ = 4*pi*A/(P*P); sld = ST.Solidity; ecc = ST.Eccentricity;

if A<S1.MIN_PUPIL_AREA || A>S1.MAX_PUPIL_AREA, MW=[]; return; end
if circ<S1.MIN_CIRCULARITY || sld<S1.MIN_SOLIDITY || ecc>S1.MAX_ECCENTRICITY, MW=[]; return; end
if touches_border_mask(MW,S1.BORDER_MARGIN_PX), MW=[]; return; end

D_in = bwdist(~MW); D_out = bwdist(MW);
rin = MW & D_in >= 1 & D_in <= max(3,S1.RING_WIDTH_PX);
rout = ~MW & D_out >= 1 & D_out <= max(3,S1.RING_WIDTH_PX);
if ~any(rin(:)) || ~any(rout(:)), MW=[]; return; end
Iin = double(N(rin)); Iout = double(N(rout));

isDark = true; if ~isempty(REF), isDark = REF.isDarkPupil; end
if isDark
  dI = mean(Iout) - mean(Iin);
  if dI < S1.DELTA_I_MIN, MW=[]; return; end
  thrBright = prctile(double(N(:)), 95);
  glintFrac = mean(Iin > thrBright);
  if glintFrac > S1.GLINT_PCT_MAX, MW=[]; return; end
else
  dI = mean(Iin) - mean(Iout);
  if dI < S1.DELTA_I_MIN, MW=[]; return; end
end

cxy = ST.Centroid;
dctr = hypot(cxy(1)-lockedCenter(1), cxy(2)-lockedCenter(2));
if dctr > S1.SEARCH_RADIUS_PX, MW = []; STW = []; score = -inf; return; end

tol = max(1, scalar_tol(S1.CENTER_TOL_PX));
distPenalty = (dctr / tol);
score = 3.5*dI + 45*circ + 12*sld - 15*distPenalty - 8*ecc;
end

function [M, ST, pickedFrom] = pick_best_candidate(cand, S1)
scores = cellfun(@(x) tern(isempty(x), -inf, x), cand(:,3));
[~,ix] = max(scores);
M = cand{ix,1}; ST = cand{ix,2}; pickedFrom = ix;
if ~S1.USE_BRIGHT_POL
  s1 = scores(1); s2 = scores(2); s3 = scores(3);
  if ix==2 && (s2 < max(s1,s3)*1.05)
    if ~isempty(cand{1,1}) && s1>=s3
      M = cand{1,1}; ST = cand{1,2}; pickedFrom = 1;
    elseif ~isempty(cand{3,1})
      M = cand{3,1}; ST = cand{3,2}; pickedFrom = 3;
    end
  end
end
end

function touch = touches_border_mask(M, mpx)
[h,w]=size(M); B=false(h,w);
if mpx>0, B(1:mpx,:)=true; B(end-mpx+1:end,:)=true; B(:,1:mpx)=true; B(:,end-mpx+1:end)=true; end
per=bwperim(M); touch=any(per & B,'all');
end

function [Mbest, STbest] = refine_shrink_to_edge(N, M, ST, S1)
maxS = 0; if isfield(S1,'SHRINK_MAX_PX'), maxS = max(0, round(S1.SHRINK_MAX_PX)); end
if maxS==0, Mbest = M; STbest = ST; return; end
[Gmag,~] = imgradient(N);
bestScore = -inf; Mbest = M; STbest = ST;
for s = 0:maxS
  if s==0, Mi = M; else, Mi = imerode(M, strel('disk', s)); if ~any(Mi(:)), break; end, end
  CC = bwconncomp(Mi); if CC.NumObjects==0, continue; end
  st = regionprops(CC,'Area','Perimeter','Centroid','PixelIdxList'); [~,kk] = max([st.Area]); Si = st(kk);
  Din  = bwdist(~Mi); Dout = bwdist(Mi); rW = max(3, S1.RING_WIDTH_PX);
  rin  = Mi  & Din  >= 1 & Din  <= rW; rout = ~Mi & Dout >= 1 & Dout <= rW;
  if ~any(rin(:)) || ~any(rout(:)), continue; end
  Ii = double(N(rin)); Io = double(N(rout));
  dI = mean(Io,'omitnan') - mean(Ii,'omitnan'); per = bwperim(Mi); perG = mean(Gmag(per),'omitnan');
  score = 3.0*dI + 0.05*perG - 0.6*s;
  if score > bestScore, bestScore = score; Mbest = Mi; STbest = Si; end
end
end

function [T, summary] = audit_step1(imgDir, labDir, outDir, ref)
if nargin < 4, ref = []; end
dI = dir(fullfile(imgDir,'frame_*.png'));
dL = dir(fullfile(labDir,'frame_*_mask.png'));
getID = @(s) sscanf(s,'frame_%d');
imgID = arrayfun(@(d) getID(d.name), dI);
labID = arrayfun(@(d) getID(d.name), dL);
[ids, ia, ib] = intersect(imgID, labID); if isempty(ids), T=table(); summary='STEP-1 AUDIT: no pairs.'; return; end
[ids, ord] = sort(ids); ia=ia(ord); ib=ib(ord);

n=numel(ids);
circ=nan(n,1); solid=nan(n,1); ecc=nan(n,1); area=nan(n,1);
cx=nan(n,1); cy=nan(n,1); cdist=nan(n,1);
meanIn=nan(n,1); meanOut=nan(n,1); ringC=nan(n,1); perimG=nan(n,1);
pass_center=false(n,1); pass_dark=false(n,1); pass_size=false(n,1);
pass_edge=false(n,1); pass_shape=false(n,1);

if ~exist(outDir,'dir'), mkdir(outDir); end
qcCSV = fullfile(outDir,'step1_qc_vs_seeds.csv');

for i=1:n
  I = imread(fullfile(dI(ia(i)).folder, dI(ia(i)).name));
  M = imread(fullfile(dL(ib(i)).folder, dL(ib(i)).name)); if size(M,3)>1, M=rgb2gray(M); end, M=M>0;
  CC=bwconncomp(M); if CC.NumObjects==0, continue; end
  st=regionprops(CC,'Area','Perimeter','Solidity','Eccentricity','Centroid');
  [~,k]=max([st.Area]); S=st(k); P=max(S.Perimeter,1);
  area(i) = S.Area; circ(i) = 4*pi*S.Area/(P*P); solid(i) = S.Solidity; ecc(i) = S.Eccentricity;
  cx(i) = S.Centroid(1); cy(i) = S.Centroid(2);

  D_in = bwdist(~M); D_out = bwdist(M);
  rin = M & D_in >= 1 & D_in <= 6; rout = ~M & D_out >= 1 & D_out <= 6;
  Ii = double(I(rin)); Io = double(I(rout));
  meanIn(i) = mean(Ii,'omitnan'); meanOut(i) = mean(Io,'omitnan'); ringC(i) = meanOut(i) - meanIn(i);
  per = bwperim(M); [Gmag,~] = imgradient(I); perimG(i) = mean(Gmag(per),'omitnan');

  if ~isempty(ref) && all(isfield(ref,{'center','centerTolPx','areaMed','areaTolAbs','ringContrastMed','contrastTol','perimGradMed','edgeTol','shapeTol'}))
    cdist(i) = hypot(cx(i)-ref.center(1), cy(i)-ref.center(2));
    pass_center(i)= cdist(i) <= ref.centerTolPx;
    pass_size(i)  = abs(area(i) - ref.areaMed) <= ref.areaTolAbs;
    pass_dark(i)  = ringC(i) >= max(1, ref.ringContrastMed - ref.contrastTol);
    pass_edge(i)  = perimG(i) >= ref.edgeTol;
    pass_shape(i) = (circ(i) >= ref.circMed - ref.shapeTol.circDrop) && ...
                    (solid(i) >= ref.solMed - ref.shapeTol.solDrop) && ...
                    (ecc(i) <= ref.eccMed + ref.shapeTol.eccRise);
  else
    cdist(i)=NaN; pass_center(i)=true; pass_size(i)=true; pass_dark(i)=true; pass_edge(i)=true; pass_shape(i)=true;
  end
end

qcScore = (double(pass_center)+double(pass_dark)+double(pass_size)+double(pass_edge)+double(pass_shape))/5;

T = table(ids(:), area, circ, solid, ecc, cx, cy, cdist, meanIn, meanOut, ringC, perimG, ...
  pass_center, pass_dark, pass_size, pass_edge, pass_shape, qcScore, ...
  'VariableNames',{'frame_id','area','circularity','solidity','eccentricity', ...
  'cx','cy','center_dist','mean_in','mean_out','ring_contrast','perim_grad', ...
  'pass_center','pass_dark','pass_size','pass_edge','pass_shape','qc_score'});

writetable(T, qcCSV);

medCirc = median(circ,'omitnan'); medSol = median(solid,'omitnan'); medEcc = median(ecc,'omitnan');
passRate = mean(qcScore>=0.8,'omitnan');
summary = sprintf(['STEP-1 AUDIT (N=%d): medCirc=%.3f, medSol=%.3f, medEcc=%.3f\n' ...
  'QC vs seeds: pass>=80%% on %.1f%% of frames; med centerDist=%.1f px, ' ...
  'med ringContrast=%.1f, med edge=%.2f'], ...
  height(T), medCirc, medSol, medEcc, 100*passRate, ...
  median(T.center_dist,'omitnan'), median(T.ring_contrast,'omitnan'), median(T.perim_grad,'omitnan'));
write_text(fullfile(outDir,'step1_qc_vs_seeds.txt'), summary);
end

function write_text(fname,s)
fid = fopen(fname,'w'); if fid==-1, return; end
fprintf(fid,'%s\n',s); fclose(fid);
end

function save_qc_gallery(imgDir, labDir, outDir, K)
dI = dir(fullfile(imgDir,'frame_*.png'));
dL = dir(fullfile(labDir,'frame_*_mask.png'));
getID = @(s) sscanf(s,'frame_%d');
imgID = arrayfun(@(d) getID(d.name), dI);
labID = arrayfun(@(d) getID(d.name), dL);
[ids, ia, ib] = intersect(imgID, labID); if isempty(ids), return; end
[ids, ord] = sort(ids); ia=ia(ord); ib=ib(ord);
take = unique(round(linspace(1, numel(ids), min(K, numel(ids)))));
if ~exist(outDir,'dir'), mkdir(outDir); end
for t = 1:numel(take)
  i = take(t);
  I = imread(fullfile(dI(ia(i)).folder, dI(ia(i)).name));
  M = imread(fullfile(dL(ib(i)).folder, dL(ib(i)).name)); if size(M,3)>1, M=rgb2gray(M); end; M=M>0;
  if size(I,3)==1, O=repmat(I,[1 1 3]); else, O=I; end
  ed = bwperim(M);
  O(:,:,2) = max(O(:,:,2), uint8(255*M));
  O(:,:,1) = max(O(:,:,1), uint8(255*ed));
  imwrite(O, fullfile(outDir, sprintf('qc_gallery_%05d.png', ids(i))));
end
end

function dataOut = pad_and_augment_roi(dataIn, cfg)
if istable(dataIn)
  I = dataIn.inputImage{1};
  Lraw = dataIn.pixelLabelImage{1};
  wrapAsTable = true;
else
  I = dataIn{1};
  Lraw = dataIn{2};
  wrapAsTable = false;
end
if size(I,3)>1, I = rgb2gray(I); end
if iscategorical(Lraw)
  Lbin = (Lraw == cfg.classNames(2));
else
  if size(Lraw,3)>1, Lraw = rgb2gray(Lraw); end
  Lbin = Lraw > 127;
end
I  = center_fit_numeric(I,  cfg.padH, cfg.padW, 0);
Lb = center_fit_numeric(uint8(Lbin), cfg.padH, cfg.padW, 0) > 0;
Lcat = categorical(Lb, [false true], cellstr(cfg.classNames));
if wrapAsTable, dataOut = table({I},{Lcat},'VariableNames',{'inputImage','pixelLabelImage'});
else, dataOut = {I, Lcat}; end
end

function d = dice_binary(A,B)
A = logical(A); B = logical(B); inter = nnz(A & B);
d = (2*inter) / max(1, nnz(A) + nnz(B));
end

function [T, summary] = quick_validate_dice(net, imdsVal, pxdsVal, targetSize, maxEval)
if isempty(imdsVal.Files), T=table(); summary='VAL: no validation files.'; return; end
N = numel(imdsVal.Files); nEval = min(N, maxEval); idx = randperm(N, nEval);
diceVals = nan(nEval,1);
[valDir,~] = fileparts(imdsVal.Files{1}); dumpDir = fullfile(valDir, '..', 'val_pred_overlays');
try, if ~exist(dumpDir,'dir'), mkdir(dumpDir); end, catch, end
for i=1:nEval
  I = im2uint8(imread(imdsVal.Files{idx(i)})); if size(I,3)>1, I = rgb2gray(I); end
  if ~isequal([size(I,1) size(I,2)], targetSize(1:2)), I = center_pad_numeric(I, targetSize(1), targetSize(2), 0); end
  L = readimage(pxdsVal, idx(i));
  if iscategorical(L), Lbin = (L=='pupil'); else, if size(L,3)>1, L = rgb2gray(L); end, Lbin = L>127; end
  if ~isequal([size(Lbin,1) size(Lbin,2)], targetSize(1:2)), Lbin = center_pad_numeric(uint8(Lbin), targetSize(1), targetSize(2), 0)>0; end
  try, C = semanticseg(I, net); catch, C = semanticseg(repmat(I,[1 1 3]), net); end
  Mpred = (C=='pupil'); Mgt = Lbin;
  diceVals(i) = dice_binary(Mpred, Mgt);
  try
    O = repmat(I,[1 1 3]);
    O(:,:,2) = max(O(:,:,2), uint8(255*Mgt));
    O(:,:,1) = max(O(:,:,1), uint8(255*bwperim(Mgt)));
    O(:,:,3) = max(O(:,:,3), uint8(255*bwperim(Mpred)));
    imwrite(O, fullfile(dumpDir, sprintf('val_%03d.png', i)));
  catch, end
end
medDice = median(diceVals,'omitnan');
summary = sprintf('STEP-2 validation: median Dice on %d pairs = %.3f', nEval, medDice);
T = table((1:nEval)', diceVals, 'VariableNames',{'case','dice'});
end

function debug_dump(I, M, ST, why, score, id, outdir)
try
  if isempty(M), return; end
  if size(I,3)==1, O=repmat(I,[1 1 3]); else, O=I; end
  Mb = logical(M); ed = bwperim(Mb);
  O(:,:,2) = max(O(:,:,2), uint8(255*Mb));
  O(:,:,1) = max(O(:,:,1), uint8(255*ed));
  if ~isempty(ST) && isfield(ST,'Centroid')
    c = round(ST.Centroid);
    rr = max(1,c(2)-2):min(size(O,1),c(2)+2);
    cc = max(1,c(1)-2):min(size(O,2),c(1)+2);
    O(rr,cc,1) = 255; O(rr,cc,2) = 0; O(rr,cc,3) = 0;
  end
  fname = sprintf('dbg_%05d_%s.png', id, regexprep(why,'[^a-zA-Z0-9_\.]','_'));
  imwrite(O, fullfile(outdir, fname));
  fid = fopen(fullfile(outdir, 'reasons.txt'),'a');
  if fid>0
    if isnumeric(score) && isfinite(score)
      fprintf(fid,'frame=%d %s score=%.3f\n', id, why, score);
    else
      fprintf(fid,'frame=%d %s\n', id, why);
    end
    fclose(fid);
  end
catch
end
end

function log_skip(fpath, frameNum, reason)
try
  if ~exist(fileparts(fpath),'dir'), mkdir(fileparts(fpath)); end
  fid = fopen(fpath,'a'); if fid>0, fprintf(fid,'frame=%d reason=%s\n', frameNum, reason); fclose(fid); end
catch
end
end

function out = tern(cond, a, b)
if cond, out=a; else, out=b; end
end

function y = ifelse(cond,a,b), if cond, y=a; else, y=b; end, end

function J = center_fit_numeric(I, H, W, padval)
if nargin<4, padval = 0; end
h = size(I,1); w = size(I,2);
topCrop = max(0, floor((h - H)/2)); bottomCrop = max(0, h - H - topCrop);
leftCrop= max(0, floor((w - W)/2)); rightCrop = max(0, w - W - leftCrop);
if topCrop>0 || bottomCrop>0 || leftCrop>0 || rightCrop>0
  rs = (1+topCrop):(h-bottomCrop); cs = (1+leftCrop):(w-rightCrop);
  I = I(rs, cs, :); h = size(I,1); w = size(I,2);
end
top = max(0, floor((H-h)/2)); bottom = max(0, H-h-top);
left= max(0, floor((W-w)/2)); right = max(0, W-w-left);
if ndims(I)==2
  J = padarray(I, [top left], padval, 'pre');
  J = padarray(J,[bottom right], padval, 'post');
else
  J = padarray(I, [top left 0], padval, 'pre');
  J = padarray(J,[bottom right 0], padval, 'post');
end
end

function J = center_pad_numeric(varargin), J = center_fit_numeric(varargin{:}); end

function s = scalar_tol(x)
if isempty(x), s = 1; return; end
if numel(x) == 1, s = double(x); return; end
s = double(max(x(:)));
end

function N = local_read_norm_roi(vrScan, frameIdx, roi, S1)
try
  vrScan.CurrentTime = max(0, (frameIdx-1)/vrScan.FrameRate);
  F = readFrame(vrScan);
  C = imcrop(F, roi);
  if isempty(C), N=[]; return; end
  if size(C,3)==3, G = rgb2gray(C); else, G=C; end
  N = preprocess_for_step1(G, S1);
  if ~isa(N,'uint8'), N = im2uint8(mat2gray(N)); end
catch
  N = [];
end
end
