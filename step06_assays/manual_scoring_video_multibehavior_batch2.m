function manual_scoring_video_multibehavior_batch2(folderPath)
    clc; close all;

    if nargin < 1
        folderPath = 'H:\addiction_dec29\day18\hot';
    end
    baseFolder = folderPath;
    outputFolder = fullfile(baseFolder, 'output');

    if ~exist(baseFolder, 'dir')
        error('Base folder does not exist: %s', baseFolder);
    end

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Check write access
    testFile = fullfile(outputFolder, 'test_write.txt');
    fid = fopen(testFile, 'w');
    if fid == -1
        error('Cannot write to %s. Check permissions or drive status.', outputFolder);
    else
        fclose(fid); delete(testFile);
    end

    % Find video files
    videoFiles = [dir(fullfile(baseFolder, '*.mp4')); ...
                  dir(fullfile(baseFolder, '*.avi')); ...
                  dir(fullfile(baseFolder, '*.mov'))];
    if isempty(videoFiles)
        error('No video files found in %s', baseFolder);
    end
    videoNames = {videoFiles.name};

    processAnother = true;
    while processAnother
        [indx, tf] = listdlg('PromptString','Select a video to process:', ...
                             'SelectionMode','single', 'ListString', videoNames);
        if ~tf
            disp('No video selected. Exiting...');
            break;
        end
        selectedVideo = videoFiles(indx);
        videoPath = fullfile(selectedVideo.folder, selectedVideo.name);
        fprintf('\n========== Processing Video: %s ==========\n', videoPath);

        scoreSingleVideo(videoPath, outputFolder);
        fprintf('Done with video: %s\n', selectedVideo.name);

        choice = questdlg('Do you want to process another video?', ...
            'Continue Processing', 'Yes','No','Yes');
        if strcmpi(choice, 'No')
            processAnother = false;
        end
    end
    disp('Video processing completed!');
end

%% Subfunction
function scoreSingleVideo(videoFile, outputFolder)
    global currentBehavior isPaused isStopped videoStarted;
    currentBehavior = 0;
    isPaused = false;
    isStopped = false;
    videoStarted = false;

    vidObj = VideoReader(videoFile);
    frameRate = vidObj.FrameRate;
    numFrames = floor(vidObj.Duration * frameRate);
    score = zeros(numFrames, 1);

    hFig = figure('Name','Manual Scoring','NumberTitle','off');
    set(hFig, 'WindowKeyPressFcn', @keyDownFcn, 'WindowKeyReleaseFcn', @keyUpFcn);
    % Display an initial placeholder
    imshow(zeros(480,640));
    text(320,240,'Press ENTER to start video playback', ...
        'Color','w','FontSize',14,'HorizontalAlignment','center');
    drawnow;

    disp('Waiting for ENTER key to start video playback...');
    while ~videoStarted
        pause(0.1);
        if ~ishandle(hFig)
            disp('Figure closed before starting. Skipping...');
            return;
        end
    end

    frameIdx = 1;
    playbackSpeed = 2;  % Set desired speed here
    frameDelay = 1 / (frameRate * playbackSpeed);

    % Read the first frame and create the image object once
    frame = readFrame(vidObj);
    hImg = imshow(frame);
    title(sprintf('Frame %d of %d (%.1fx speed)', frameIdx, numFrames, playbackSpeed));
    drawnow;

    % Main processing loop
    while hasFrame(vidObj) && ishandle(hFig) && ~isStopped && frameIdx <= numFrames
        if ~isPaused
            tStart = tic;
            frame = readFrame(vidObj);
            score(frameIdx) = currentBehavior;
            
            % Update the image data rather than calling imshow repeatedly
            set(hImg, 'CData', frame);
            title(sprintf('Frame %d of %d (%.1fx speed)', frameIdx, numFrames, playbackSpeed));
            drawnow;
            
            frameIdx = frameIdx + 1;
            elapsed = toc(tStart);
            pause(max(0, frameDelay - elapsed));  % Ensures intended pacing
        else
            pause(0.05);
        end
    end

    %% Post-processing
    timeAxis = (0:numFrames-1) / frameRate;
    % Define behavior names (with 'Jump' added as behavior 6)
    behaviorNames = {'Non-moving','Licking','Rearing','Flinching','Hindlimb Licking','Jump'};
    behaviorColors = lines(6);
    
    % Create jump variable (logical vector: true when score equals jump code 5)
    jump = (score == 5);

    [~, vidName, ~] = fileparts(videoFile);

    % 1. Time Series Plot
    figTimeSeries = figure('Name','Behavior Time Series','NumberTitle','off');
    hold on;
    for b = 0:5
        idx = (score == b);
        if any(idx)
            scatter(timeAxis(idx), score(idx), 36, behaviorColors(b+1,:), 'filled', ...
                    'DisplayName', behaviorNames{b+1});
        end
    end
    xlabel('Time (s)'); ylabel('Behavior Code');
    title('Behavior Time Series'); legend('show'); grid on;
    saveas(figTimeSeries, fullfile(outputFolder, ['BehaviorTimeSeries_', vidName, '.png']));

    % 2. Raster Plot
    figRaster = figure('Name','Raster Plots','NumberTitle','off');
    for b = 0:5
        subplot(6,1,b+1);
        binaryScore = (score == b);
        diffScore = diff([0; binaryScore(:); 0]);
        startIdx = find(diffScore == 1);
        stopIdx  = find(diffScore == -1);

        % Clip indices to within bounds
        startIdx(startIdx > numFrames) = numFrames;
        stopIdx(stopIdx > numFrames) = numFrames;

        % Ensure equal length
        minLen = min(length(startIdx), length(stopIdx));
        if length(startIdx) ~= length(stopIdx)
            warning('⚠️ Unmatched start/stop pairs for %s. Truncating.', behaviorNames{b+1});
        end
        startIdx = startIdx(1:minLen);
        stopIdx  = stopIdx(1:minLen);

        startTimes = timeAxis(startIdx);
        stopTimes  = timeAxis(stopIdx);

        for j = 1:minLen
            line([startTimes(j), stopTimes(j)], [j j], 'LineWidth', 2, ...
                 'Color', behaviorColors(b+1,:));
            hold on;
        end
        xlabel('Time (s)'); ylabel('Episode');
        title(behaviorNames{b+1}); grid on;
    end

    saveas(figRaster, fullfile(outputFolder, ['RasterPlot_', vidName, '.png']));

    % 3. Save .mat file (including jump variable)
    matName = ['ManualScoringResults_', vidName, '.mat'];
    matPath = fullfile(outputFolder, matName);
    save(matPath, 'score', 'timeAxis', 'behaviorNames', 'frameRate', 'numFrames', 'jump');
    fprintf('✅ Scoring results saved to: %s\n', matPath);
end

%% Key Callbacks
function keyDownFcn(~, event)
    global currentBehavior isPaused isStopped videoStarted;
    switch lower(event.Key)
        case 'return'
            if ~videoStarted
                videoStarted = true;
                disp('Starting video playback...');
            end
        case 'a', currentBehavior = 1;
        case 's', currentBehavior = 2;
        case 'd', currentBehavior = 3;
        case 'f', currentBehavior = 4;
        case 'v', currentBehavior = 5;  % Added key for "jump"
        case 'space', isPaused = ~isPaused;
        case 'q', isStopped = true; disp('Playback stopped by user.');
    end
end

function keyUpFcn(~, event)
    global currentBehavior;
    switch lower(event.Key)
        case {'a','s','d','f','v'}  % Reset behavior when key is released
            currentBehavior = 0;
    end
end
