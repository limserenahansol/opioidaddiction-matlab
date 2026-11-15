function [m, info] = fullpupil_mask(C)
% fullpupil_mask  The whole dark area inside the eyelids on one eye crop (uint8 or double, the ROI crop that
% pupil_mask_frame makes). This is the pupil when it is fully dilated (rule from the user, 2026-11-15: "when the whole
% eye is dark, the pupil is the whole dark area inside the lids").
%   1) smooth (Gaussian sigma 1.5); glints = brightest 1 % of the crop, dilated 3 px
%   2) Otsu threshold between the dark eye and the skin (glint pixels left out); dark = below it
%   3) open with a disk of 4 px (cuts thin dark lines: lid margin hairs, the spout's shadow); keep the dark component
%      with the largest area in the central half of the crop; fill holes (glints inside the eye)
%   4) glints on the edge of the eye: their pixels are added where a closing (disk 6) of the mask covers them; a glint
%      outside the eye (below the lid) is not added
% info: thr (Otsu level), diam (equivalent diameter, px), skin (median outside the mask), S (smoothed crop), G (glints)
I = double(C); [h, w] = size(I);
S = imgaussfilt(I, 1.5);
G = imdilate(S >= prctile(S(:), 99), strel('disk', 3));
t = multithresh(S(~G), 1);
A = S < t;
A = imopen(A, strel('disk', 4));
CC = bwconncomp(A); best = 0; bestArea = 0;
cen = false(h, w); cen(round(h / 4):round(3 * h / 4), round(w / 4):round(3 * w / 4)) = true;
for c = 1:CC.NumObjects
    a = nnz(cen(CC.PixelIdxList{c}));
    if a > bestArea, bestArea = a; best = c; end
end
m = false(h, w);
if best > 0
    m(CC.PixelIdxList{best}) = true;
    m = imfill(m, 'holes');
    mc = imfill(imclose(m, strel('disk', 6)), 'holes');
    m = m | (mc & G);
    m = imfill(m, 'holes');
end
out = ~imdilate(m, strel('disk', 3)) & ~G;
info = struct('thr', t, 'diam', 2 * sqrt(nnz(m) / pi), 'skin', median(S(out)), 'S', S, 'G', G);
end
