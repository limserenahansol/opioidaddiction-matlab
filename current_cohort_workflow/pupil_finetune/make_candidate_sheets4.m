function make_candidate_sheets4(diagDir, outDir, modelFile, pkgPupilDir)
% Labelling sheets: for every session 2 sheets x 12 sampled frames. Each frame twice: left = enlarged eye crop
% (contrast stretched), right = the same with ALL pupil candidates (U-Net blobs + dark-region chains) drawn as
% numbered contours. Candidates are numbered left to right (no ranking). Features of every candidate -> candidates.csv
addpath(pkgPupilDir);
if ~exist(outDir, 'dir'), mkdir(outDir); end
M = load(modelFile);
meta = jsondecode(fileread(fullfile(diagDir, 'meta.json')));
cols = [1 0 0; 0 1 0; 0.25 0.55 1; 1 1 0; 1 0 1; 0 1 1; 1 0.55 0; 1 1 1];
rows = {}; maxC = 8;
for k = 1:numel(meta)
    key = meta(k).key;
    S = load(fullfile(diagDir, [key '.mat']));
    X = S.samples; nF = size(X, 3);
    r0 = double(S.roi_in_crop(:)'); rr = [r0(1) + 1, r0(2) + 1, r0(3), r0(4)];
    F = arrayfun(@(j) X(:, :, j), 1:nF, 'UniformOutput', false)';
    [masks, ~, crops] = pupil_mask_frame(F, rr, M, 'IrisCorrection', true);
    allc = cell(nF, 1);
    for sh = 1:ceil(nF / 4)
        fr = (sh - 1) * 4 + (1:4); fr = fr(fr <= nF);
        fig = figure('Visible', 'off', 'Position', [0 0 2000 860]);
        tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        for j = fr
            C = crops{j};
            cand = pupil_hybrid_candidates(C, masks{j});
            if numel(cand) > maxC   % keep the strongest by edge contrast, then renumber left to right
                [~, o] = sort([cand.con], 'descend'); cand = cand(sort(o(1:maxC)));
            end
            [~, o] = sort([cand.cx]); cand = cand(o); allc{j} = rmfield(cand, {'sec'});
            lo = prctile(double(C(:)), 1); hi = prctile(double(C(:)), 75);   % dark-range stretch: the eye fills the grey range
            V = imresize(uint8(255 * min(max((double(C) - lo) / max(1, hi - lo), 0), 1)), 3, 'bicubic');
            nexttile(tl); imshow(V); title(sprintf('frame %d', j - 1), 'FontSize', 13);
            nexttile(tl); imshow(V); hold on
            placed = zeros(0, 2);
            for c = 1:numel(cand)
                B = bwboundaries(imresize(cand(c).m, 3, 'nearest'));
                for q = 1:numel(B), plot(B{q}(:,2), B{q}(:,1), 'Color', cols(c, :), 'LineWidth', 1.4); end
                [yy, xx] = find(imresize(cand(c).m, 3, 'nearest'));
                [~, ii] = min(yy - 0.3 * xx); tx = xx(ii) + 6; ty = yy(ii) - 8;
                while any(hypot(placed(:,1) - tx, placed(:,2) - ty) < 16), ty = ty + 16; end
                placed(end+1, :) = [tx ty]; %#ok<AGROW>
                text(tx, ty, sprintf('%d', c), 'Color', cols(c, :), 'FontWeight', 'bold', 'FontSize', 15, ...
                    'BackgroundColor', 'k', 'Margin', 1);
                rows(end+1, :) = {key, j - 1, c, cand(c).src, cand(c).area, cand(c).diam, cand(c).cx, cand(c).cy, ...
                    cand(c).circ, cand(c).sol, cand(c).ecc, cand(c).con, cand(c).sec(1), cand(c).sec(2), cand(c).sec(3), ...
                    cand(c).sec(4), cand(c).dark, cand(c).stab, cand(c).dglint, cand(c).dctr, cand(c).t, cand(c).unet_ov}; %#ok<AGROW>
            end
            if isempty(cand)
                title(sprintf('frame %d: no candidate', j - 1), 'FontSize', 10, 'Color', [0.8 0 0]);
            else
                title(sprintf('frame %d: candidates 1-%d', j - 1, numel(cand)), 'FontSize', 10);
            end
        end
        title(tl, sprintf('%s   sheet %s (frames %d-%d)   left: eye crop (dark range stretched)   right: numbered candidate contours', key, ...
            char('A' + sh - 1), fr(1) - 1, fr(end) - 1), 'Interpreter', 'none', 'FontSize', 15);
        exportgraphics(fig, fullfile(outDir, sprintf('%s__cand%s.png', key, char('A' + sh - 1))), 'Resolution', 72);
        close(fig);
    end
    save(fullfile(outDir, [key '_cands.mat']), 'allc', '-v7');
    fprintf('%s done\n', key);
end
T = cell2table(rows, 'VariableNames', {'session', 'frame', 'cand', 'src', 'area', 'diam', 'cx', 'cy', 'circ', 'sol', 'ecc', ...
    'con', 'sec1', 'sec2', 'sec3', 'sec4', 'dark', 'stab', 'dglint', 'dctr', 't', 'unet_ov'});
writetable(T, fullfile(outDir, 'candidates.csv'));
end
