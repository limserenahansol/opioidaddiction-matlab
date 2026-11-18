function J = eye_prep(I)
% eye_prep  Input of the whole-eye model (2026-11-17, 9670f black): the grey eye crop with adaptive histogram
% equalisation, as the package does before its U-Net, but WITHOUT the glint erasing: a cloudy cornea is bright and
% must stay. The same function is used for training (eye_finetune) and tracking (eye_track_session).
if size(I, 3) > 1, I = rgb2gray(I); end
J = adapthisteq(im2uint8(I));
end
