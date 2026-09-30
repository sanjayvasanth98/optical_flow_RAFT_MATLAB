%% ------------------------------------------------------------------------
% Optical Flow RAFT -- local desktop version
% Based on raftmatlabsideview_ARC.m, but intended for interactive PC tests.
% This file does not modify the cluster ARC runner.
%% ------------------------------------------------------------------------
clear; clc; close all;

%% ------------------------- USER SETTINGS -------------------------------
% Set this to the total number of video frames to read (including the first
% reference frame).  Use Inf to process the entire video.
numFramesToProcess = 20;

showFigures = true;
savePlots = true;
raftIters = 8;
raftTolerance = 1e-6;
accelMode = "auto";

% Optional: set the full path here to skip the file-selection dialog.
% Leave this as "" to choose a video interactively when the script runs.
videoPath = "E:\Sept 2026 Flowfield data\P10S20\P10S20_3_475_40lpm.avi";

% Cropped videos may not need an ROI or a throat reference.
useFullFrameROI = true;   % true: use every pixel and skip ROI selection
useThroat = false;        % false: skip throat selection and omit its plot line

% Contrast preprocessing applied to every frame sent to RAFT.
% Run once with "none" and once with "clahe" to compare the results.
contrastMode = "clahe";   % "none" or "clahe"

% Compare a few evenly spaced frames within numFramesToProcess.
showCLAHEPreview = false;
previewOnly = false;      % false: continue with RAFT after the preview
numPreviewFrames = 4;

% Calibration parameters
mm_per_pixel = 0.00828164;          % [mm/pixel]
fps = 102247;                       % [frames/second]
m_per_pixel = mm_per_pixel / 1000;  % [m/pixel]

validContrastModes = ["none", "clahe"];
contrastMode = lower(string(contrastMode));
if ~ismember(contrastMode, validContrastModes)
    error('contrastMode must be "none" or "clahe".');
end
contrastLabel = char(contrastMode);

%% ------------------------ SELECT VIDEO ---------------------------------
if strlength(videoPath) == 0
    [videoName, videoFolder] = uigetfile( ...
        {'*.avi;*.mp4;*.mov;*.mkv', 'Video files (*.avi, *.mp4, *.mov, *.mkv)'}, ...
        'Select the video to process');
    if isequal(videoName, 0)
        error('No video selected.');
    end
    videoPath = fullfile(videoFolder, videoName);
else
    videoPath = char(videoPath);
    if ~isfile(videoPath)
        error('Video file not found: %s', videoPath);
    end
end

[filepath, filename] = fileparts(videoPath);
v = VideoReader(videoPath);

availableFrames = max(2, floor(v.Duration * v.FrameRate));
if isinf(numFramesToProcess)
    numFrames = availableFrames;
else
    validateattributes(numFramesToProcess, {'numeric'}, ...
        {'scalar', 'integer', '>=', 2}, mfilename, 'numFramesToProcess');
    numFrames = min(numFramesToProcess, availableFrames);
end
numFlowFrames = numFrames - 1;

fprintf('Video: %s\n', videoPath);
%% ------------------------- CLAHE PREVIEW -------------------------------
% Use the same preprocessing helper as RAFT, regardless of contrastMode.
if showCLAHEPreview || previewOnly
    validateattributes(numPreviewFrames, {'numeric'}, ...
        {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'numPreviewFrames');
    previewIndices = unique(round(linspace(1, numFrames, ...
        min(numPreviewFrames, numFrames))));
    if savePlots
        previewDir = fullfile(filepath, 'RAFT_results', 'clahe', 'plots');
        if ~isfolder(previewDir), mkdir(previewDir); end
    end
    for frameIndex = previewIndices
        rawFrame = read(v, frameIndex);
        original = prepareFrame(rawFrame, "none");
        enhanced = prepareFrame(rawFrame, "clahe");
        previewFig = figure('Name', sprintf('CLAHE preview - frame %d', frameIndex), ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Visible', tern(showFigures, 'on', 'off'));
        tiledlayout(previewFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        % Fixed [0, 1] display range makes the contrast comparison faithful.
        nexttile; imshow(im2double(original), [0 1]);
        title(sprintf('Original grayscale - frame %d', frameIndex));
        nexttile; imshow(im2double(enhanced), [0 1]);
        title(sprintf('CLAHE - frame %d', frameIndex));
        drawnow;
        if savePlots
            previewFile = fullfile(previewDir, sprintf('%s_CLAHE_preview_frame%04d.png', ...
                filename, frameIndex));
            exportgraphics(previewFig, previewFile, 'Resolution', 200);
            fprintf('Saved CLAHE preview: %s\n', previewFile);
        end
        if ~showFigures, close(previewFig); end
    end
    v.CurrentTime = 0;
end
if previewOnly
    fprintf('Preview complete. Set previewOnly = false to run RAFT.\n');
    return;
end

if canUseGPU
    executionEnv = "gpu";
else
    executionEnv = "cpu";
end

fprintf('Processing %d of %d video frames (%d flow fields).\n', ...
    numFrames, availableFrames, numFlowFrames);
fprintf('Contrast preprocessing: %s\n', contrastLabel);
fprintf('Execution environment: %s\n', executionEnv);

%% -------------------------- OUTPUTS ------------------------------------
% Keep the two contrast runs separate for a direct comparison.
resultsDir = fullfile(filepath, 'RAFT_results', contrastLabel);
plotsDir = fullfile(resultsDir, 'plots');
matDir = fullfile(resultsDir, 'mat_files');
if ~isfolder(plotsDir), mkdir(plotsDir); end
if ~isfolder(matDir), mkdir(matDir); end
fprintf('Results folder: %s\n', resultsDir);

%% -------------------------- ROI ----------------------------------------
v.CurrentTime = 0;
framePrev = im2gray(readFrame(v));
[H, W] = size(framePrev);
roiFile = fullfile(filepath, [filename '_ROI.mat']);

if useFullFrameROI
    maskROI = true(H, W);
elseif isfile(roiFile)
    S = load(roiFile, 'maskROI');
    maskROI = logical(S.maskROI);
    if any(size(maskROI) ~= [H W])
        error('ROI mask size does not match the selected video.');
    end
    useSavedROI = questdlg('Use the saved ROI for this video?', ...
        'ROI', 'Use saved ROI', 'Draw a new ROI', 'Use saved ROI');
else
    useSavedROI = 'Draw a new ROI';
end

if ~useFullFrameROI && strcmp(useSavedROI, 'Draw a new ROI')
    figure('Name', 'ROI Selection', 'NumberTitle', 'off');
    imshow(framePrev, []);
    title('Draw the ROI polygon; double-click to finish');
    maskROI = roipoly;
    if isempty(maskROI), error('No ROI selected.'); end
    save(roiFile, 'maskROI');
    close(gcf);
end

throatFile = fullfile(filepath, [filename '_throat.mat']);
if ~useThroat
    x_throat_pixel = NaN;
    y_throat_pixel = NaN;
elseif isfile(throatFile)
    Sth = load(throatFile, 'x_throat', 'y_throat');
    x_throat_pixel = double(Sth.x_throat);
    if isfield(Sth, 'y_throat'), y_throat_pixel = double(Sth.y_throat); else, y_throat_pixel = NaN; end
else
    figure('Name', 'Throat Selection', 'NumberTitle', 'off');
    imshow(framePrev, []); title('Click the throat location');
    [x_throat_pixel, y_throat_pixel] = ginput(1);
    save(throatFile, 'x_throat_pixel', 'y_throat_pixel');
    % Keep the variable names compatible with ROI_and_throatloader.m.
    x_throat = x_throat_pixel; y_throat = y_throat_pixel;
    save(throatFile, 'x_throat', 'y_throat');
    close(gcf);
end

%% ---------------------- RAFT PROCESSING --------------------------------
opticalFlowObj = opticalFlowRAFT;
sumU = zeros(H, W, 'double');
sumV = zeros(H, W, 'double');
sumMag = zeros(H, W, 'double');

% This local test version retains the instantaneous fields in memory.
% Use a modest numFramesToProcess value; the ARC script is the version for
% large, disk-streamed cluster processing.
u_all = zeros(H, W, numFlowFrames, 'single');
v_all = zeros(H, W, numFlowFrames, 'single');

v.CurrentTime = 0;
framePrev = prepareFrame(readFrame(v), contrastMode);
fprintf('Running RAFT...\n');
tic;
for k = 1:numFlowFrames
    frameCurr = prepareFrame(readFrame(v), contrastMode);
    flow = estimateFlow(opticalFlowObj, frameCurr, ...
        ExecutionEnvironment=executionEnv, Acceleration=accelMode, ...
        MaxIterations=raftIters, Tolerance=raftTolerance);

    u_phys = flow.Vx .* maskROI * m_per_pixel * fps;
    v_phys = flow.Vy .* maskROI * m_per_pixel * fps;
    u_all(:, :, k) = single(u_phys);
    v_all(:, :, k) = single(v_phys);
    sumU = sumU + u_phys;
    sumV = sumV + v_phys;
    sumMag = sumMag + hypot(u_phys, v_phys);

    elapsed = toc;
    eta = elapsed / k * (numFlowFrames - k);
    fprintf('Flow frame %d/%d (%.1f%%), ETA %.1f s\n', ...
        k, numFlowFrames, 100 * k / numFlowFrames, eta);
    framePrev = frameCurr;
end

%% ------------------------ SAVE RESULTS ---------------------------------
u_mean = single(sumU / numFlowFrames);
v_mean = single(sumV / numFlowFrames);
velMean = single(sumMag / numFlowFrames);
u_mean(~maskROI) = NaN;
v_mean(~maskROI) = NaN;
velMean(~maskROI) = NaN;

u_mean_phys = flipud(u_mean);
v_mean_phys = -flipud(v_mean);
velMean_phys = flipud(velMean);
x_mm = (0:W-1) * mm_per_pixel;
y_mm = (0:H-1) * mm_per_pixel;
x_throat_mm = x_throat_pixel * mm_per_pixel;

outFile = fullfile(matDir, sprintf('%s_velocity_local_%s_%dframes.mat', ...
    filename, contrastLabel, numFrames));
save(outFile, 'u_all', 'v_all', 'u_mean', 'v_mean', 'velMean', ...
    'u_mean_phys', 'v_mean_phys', 'velMean_phys', 'maskROI', ...
    'x_mm', 'y_mm', 'x_throat_pixel', 'y_throat_pixel', 'x_throat_mm', ...
    'mm_per_pixel', 'm_per_pixel', 'fps', 'numFrames', 'numFlowFrames', '-v7.3');
fprintf('Saved local results: %s\n', outFile);

if savePlots
    fig = figure('Visible', tern(showFigures, 'on', 'off'), 'Color', 'w');
    imagesc(x_mm, y_mm, velMean_phys); axis image; set(gca, 'YDir', 'normal');
    colormap(turbo); colorbar;
    if useThroat
        hold on; plot([x_throat_mm x_throat_mm], [y_mm(1) y_mm(end)], 'w:', 'LineWidth', 1.5);
    end
    xlabel('X (mm)'); ylabel('Y (mm)'); title('Time-Averaged Velocity Magnitude');
    plotFile = fullfile(plotsDir, sprintf('%s_TimeAvgVelMag_local_%s_%dframes.png', ...
        filename, contrastLabel, numFrames));
    exportgraphics(fig, plotFile, 'Resolution', 300);
    fprintf('Saved local plot: %s\n', plotFile);

    figU = figure('Visible', tern(showFigures, 'on', 'off'), 'Color', 'w');
    imagesc(x_mm, y_mm, u_mean_phys); axis image; set(gca, 'YDir', 'normal');
    colormap(turbo); colorbar;
    xlabel('X (mm)'); ylabel('Y (mm)'); title('Time-Averaged Horizontal Velocity (Vx)');
    plotUFile = fullfile(plotsDir, sprintf('%s_TimeAvgUVel_local_%s_%dframes.png', ...
        filename, contrastLabel, numFrames));
    exportgraphics(figU, plotUFile, 'Resolution', 300);
    fprintf('Saved mean U plot: %s\n', plotUFile);

    figV = figure('Visible', tern(showFigures, 'on', 'off'), 'Color', 'w');
    imagesc(x_mm, y_mm, v_mean_phys); axis image; set(gca, 'YDir', 'normal');
    colormap(turbo); colorbar;
    xlabel('X (mm)'); ylabel('Y (mm)'); title('Time-Averaged Vertical Velocity (Vy)');
    plotVFile = fullfile(plotsDir, sprintf('%s_TimeAvgVVel_local_%s_%dframes.png', ...
        filename, contrastLabel, numFrames));
    exportgraphics(figV, plotVFile, 'Resolution', 300);
    fprintf('Saved mean V plot: %s\n', plotVFile);
end

function out = tern(condition, whenTrue, whenFalse)
if condition, out = whenTrue; else, out = whenFalse; end
end

function frame = prepareFrame(frame, contrastMode)
frame = im2gray(frame);
if contrastMode == "clahe"
    if canUseGPU
        frame = gpuClahe(frame);
    else
        frame = adapthisteq(frame);
    end
end
end

function frame = gpuClahe(frame)
frameGPU = gpuArray(im2single(frame));
[height, width] = size(frameGPU);
numTileRows = min(8, height);
numTileCols = min(8, width);
outputGPU = zeros(height, width, 'like', frameGPU);

for tileRow = 1:numTileRows
    rowStart = floor((tileRow - 1) * height / numTileRows) + 1;
    rowEnd = floor(tileRow * height / numTileRows);
    for tileCol = 1:numTileCols
        colStart = floor((tileCol - 1) * width / numTileCols) + 1;
        colEnd = floor(tileCol * width / numTileCols);
        tileGPU = frameGPU(rowStart:rowEnd, colStart:colEnd);

        [histogramGPU, ~] = histcounts(tileGPU, 256, 'BinLimits', [0 1]);
        clipLimit = max(1, 0.01 * numel(tileGPU) / 256);
        excess = sum(max(histogramGPU - clipLimit, 0));
        histogramGPU = min(histogramGPU, clipLimit) + excess / 256;
        cdfGPU = cumsum(histogramGPU) / numel(tileGPU);

        binGPU = min(floor(tileGPU * 256) + 1, 256);
        outputGPU(rowStart:rowEnd, colStart:colEnd) = cdfGPU(binGPU);
    end
end

frame = gather(im2uint8(outputGPU));
end
