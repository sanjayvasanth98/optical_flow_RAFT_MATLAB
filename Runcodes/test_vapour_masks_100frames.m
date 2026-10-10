%% Standalone dark-vapour segmentation test (no RAFT calculation)
% Requires MATLAB and Image Processing Toolbox.
% Edit the inputs below, then run this script. Existing RAFT code is untouched.
% Frame ranges are inclusive, using original 1-based video frame numbers.
% Threshold is in grayscale intensity units on a 0..255 scale.

%% Inputs
videoPaths = {
    "E:\Sept 2026 Flowfield data\P10S30\P10S30_3_472_40lpm.avi"
    % 'D:\path\to\another_video.avi'
};
frameRanges = [2000 2100];             % One [firstFrame lastFrame] row per video.

% Supply a logical/binary H-by-W mask shared by all videos, or a cell array
% with one mask per video. [] explicitly selects the whole image.
% Example with different ROIs: roiMask = {maskForVideo1, maskForVideo2};
roiFile = 'E:\Sept 2026 Flowfield data\P10S30\P10S30_3_472_40lpm\P10S30_3_472_40lpm_ROI.mat';
assert(isfile(roiFile), 'ROI file not found: %s', roiFile);
roiData = load(roiFile, 'maskROI');
assert(isfield(roiData, 'maskROI'), ...
    'ROI file must contain the variable maskROI: %s', roiFile);
roiMask = roiData.maskROI;

threshold = 20;                   % Starting value; adjust for your videos.
minArea = 20;                     % Keep blobs with at least this many pixels.
nBgFrames = 50;                   % Evenly sampled across the WHOLE video.
outputRoot = 'D:\Working codes\RAFT code\optical_flow_RAFT_MATLAB\test100frames';

%% Validate configuration and run each video
videoPaths = cellstr(string(videoPaths(:)));
nVideos = numel(videoPaths);
assert(nVideos > 0, 'Provide at least one .avi path in videoPaths.');
validateattributes(frameRanges, {'numeric'}, ...
    {'2d', 'size', [nVideos 2], 'integer', 'positive', 'finite'}, ...
    mfilename, 'frameRanges');
assert(all(frameRanges(:, 2) >= frameRanges(:, 1)), ...
    'Each frame range must have lastFrame >= firstFrame.');
validateattributes(threshold, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'threshold');
validateattributes(minArea, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'minArea');
validateattributes(nBgFrames, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'nBgFrames');
if iscell(roiMask)
    assert(numel(roiMask) == nVideos, 'Provide one ROI mask per video.');
end

videoStems = cell(nVideos, 1);
for videoIndex = 1:nVideos
    assert(isfile(videoPaths{videoIndex}), ...
        'Video file not found: %s', videoPaths{videoIndex});
    [~, videoStems{videoIndex}, extension] = fileparts(videoPaths{videoIndex});
    assert(strcmpi(extension, '.avi'), ...
        'Expected an .avi file: %s', videoPaths{videoIndex});
end
if ~isfolder(outputRoot), mkdir(outputRoot); end

for videoIndex = 1:nVideos
    if iscell(roiMask)
        videoROI = roiMask{videoIndex};
    else
        videoROI = roiMask;
    end
    % Prefix duplicate stems so distinct input videos cannot share outputs.
    folderName = videoStems{videoIndex};
    if nnz(strcmpi(videoStems, folderName)) > 1
        folderName = sprintf('%03d_%s', videoIndex, folderName);
    end
    processVapourVideo(videoPaths{videoIndex}, frameRanges(videoIndex, :), ...
        videoROI, threshold, minArea, nBgFrames, ...
        fullfile(outputRoot, folderName));
end

%% Local helpers
function processVapourVideo(videoPath, frameRange, roiMask, ...
        threshold, minArea, nBgFrames, outputFolder)
    videoTimer = tic;
    video = VideoReader(videoPath);
    totalFrames = countVideoFrames(video);
    assert(frameRange(2) <= totalFrames, ...
        'Requested frame %d, but %s has only %d frames.', ...
        frameRange(2), videoPath, totalFrames);

    firstRawFrame = read(video, 1);
    imageSize = [size(firstRawFrame, 1), size(firstRawFrame, 2)];
    if isempty(roiMask)
        roiMask = true(imageSize);
    else
        validateattributes(roiMask, {'logical', 'numeric'}, ...
            {'2d', 'size', imageSize, 'real', 'finite'}, ...
            mfilename, 'roiMask');
        assert(all(roiMask(:) == 0 | roiMask(:) == 1), ...
            'ROI must be a logical or binary (0/1) mask: %s', videoPath);
        roiMask = logical(roiMask);
    end
    roiAreaPixels = nnz(roiMask);
    assert(roiAreaPixels > 0, 'ROI must contain at least one pixel: %s', videoPath);

    effectiveBgFrames = min(nBgFrames, totalFrames);
    if effectiveBgFrames < nBgFrames
        warning('vapourTest:ShortVideo', ...
            '%s has %d frames; using all frames instead of %d background samples.', ...
            videoPath, totalFrames, nBgFrames);
    end
    [background, backgroundFrameNumbers] = computeVapourBackground(video, nBgFrames);

    if ~isfolder(outputFolder), mkdir(outputFolder); end
    qcFolder = fullfile(outputFolder, 'qc');
    if ~isfolder(qcFolder), mkdir(qcFolder); end
    imwrite(uint8(round(background)), fullfile(outputFolder, 'background.png'));

    frameNumbers = (frameRange(1):frameRange(2)).';
    framesProcessed = numel(frameNumbers);
    maskPixelIndices = cell(framesProcessed, 1);
    vapourAreaPixels = zeros(framesProcessed, 1);
    vapourAreaFraction = zeros(framesProcessed, 1);
    qcFigure = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1000 700], 'NumberTitle', 'off');
    figureCleanup = onCleanup(@() close(qcFigure));
    qcAxes = axes('Parent', qcFigure);
    fprintf('\nVideo: %s\n', videoPath);
    fprintf('Background: %d samples over frames 1-%d; test range: %d-%d.\n', ...
        effectiveBgFrames, totalFrames, frameRange(1), frameRange(2));

    for frameIndex = 1:framesProcessed
        frameNumber = frameNumbers(frameIndex);
        rawFrame = read(video, frameNumber);
        coreMask = detectVapourCore(rawFrame, background, roiMask, threshold, minArea);
        % Store the core only: no dilation or expansion is applied.
        maskPixelIndices{frameIndex} = find(coreMask);
        vapourAreaPixels(frameIndex) = numel(maskPixelIndices{frameIndex});
        vapourAreaFraction(frameIndex) = vapourAreaPixels(frameIndex) / roiAreaPixels;

        cla(qcAxes);
        imshow(rawFrame, 'Parent', qcAxes);
        hold(qcAxes, 'on');
        boundaries = bwboundaries(coreMask, 8, 'holes');
        for boundaryIndex = 1:numel(boundaries)
            boundary = boundaries{boundaryIndex};
            plot(qcAxes, boundary(:, 2), boundary(:, 1), ...
                'r-', 'LineWidth', 1);
        end
        hold(qcAxes, 'off');
        title(qcAxes, sprintf('Frame %d | Vapour area / ROI = %.4f', ...
            frameNumber, vapourAreaFraction(frameIndex)), 'Interpreter', 'none');
        exportgraphics(qcAxes, fullfile(qcFolder, ...
            sprintf('frame_%06d.png', frameNumber)), 'Resolution', 150);
    end

    meanVapourAreaFraction = mean(vapourAreaFraction);
    parameters = struct('threshold', threshold, 'minArea', minArea, ...
        'nBgFrames', nBgFrames, 'effectiveBgFrames', effectiveBgFrames, ...
        'frameRange', frameRange, 'totalVideoFrames', totalFrames, ...
        'backgroundFrameNumbers', backgroundFrameNumbers, ...
        'intensityScale', [0 255], 'blobConnectivity', 8, ...
        'holeFillConnectivity', 4, 'fillHoles', true, 'dilationApplied', false, ...
        'areaFractionDenominator', 'ROI pixels', ...
        'frameRate', video.FrameRate, 'videoDurationSeconds', video.Duration, ...
        'outputFolder', outputFolder, 'qcResolutionDPI', 150);
    % Row k in frameNumbers corresponds to cell k in maskPixelIndices.
    % Reconstruct: mask = false(imageSize); mask(maskPixelIndices{k}) = true;
    resultFile = fullfile(outputFolder, 'vapour_core_masks.mat');
    save(resultFile, 'videoPath', 'frameNumbers', 'maskPixelIndices', ...
        'imageSize', 'roiMask', 'background', 'backgroundFrameNumbers', ...
        'parameters', 'framesProcessed', 'roiAreaPixels', 'vapourAreaPixels', ...
        'vapourAreaFraction', 'meanVapourAreaFraction', '-v7.3');
    runtimeSeconds = toc(videoTimer);
    save(resultFile, 'runtimeSeconds', '-append');
    fprintf('Frames processed: %d | Mean vapour-area fraction: %.6f (%.3f%% of ROI) | Runtime: %.2f s\n', ...
        framesProcessed, meanVapourAreaFraction, ...
        100 * meanVapourAreaFraction, runtimeSeconds);
    fprintf('Saved: %s\n', resultFile);
end

function totalFrames = countVideoFrames(video)
    % Use the video frame count, not a rounded Duration * FrameRate estimate.
    totalFrames = [];
    try
        totalFrames = double(video.NumFrames);
    catch
        % Some readers cannot supply a frame count; count decoded frames.
    end
    if isempty(totalFrames) || ~isscalar(totalFrames) || ...
            ~isfinite(totalFrames) || totalFrames < 1 || totalFrames ~= fix(totalFrames)
        video.CurrentTime = 0;
        totalFrames = 0;
        while hasFrame(video)
            readFrame(video);
            totalFrames = totalFrames + 1;
        end
        video.CurrentTime = 0;
    end
    assert(totalFrames > 0, 'Video contains no readable frames.');
end
