function eye_track_all(cropsJson, modelFile, outRoot, varargin)
% eye_track_all  eye_track_session for every session in analysis/wholeeye/crops.json (the mice whose pupil cannot be
% seen: 9670f black and 9669m black, 2026-11-17), one after another. A session whose wholeeye_frames.csv exists is
% skipped. Output: <outRoot>/day<N>/<cage>/<mouse>/
C = jsondecode(fileread(cropsJson));
names = fieldnames(C);
for i = 1:numel(names)
    c = C.(names{i});
    out = fullfile(outRoot, sprintf('day%d', c.day), c.cage, c.mouse);
    if exist(fullfile(out, 'wholeeye_frames.csv'), 'file')
        fprintf('[eye] %s done before - skipped\n', names{i}); continue
    end
    if ~exist(out, 'dir'), mkdir(out); end
    try
        eye_track_session(c.video, c.crop(:)', modelFile, out, varargin{:});
    catch err
        fprintf(2, '[eye] %s FAILED: %s\n', names{i}, err.message);
    end
end
end
