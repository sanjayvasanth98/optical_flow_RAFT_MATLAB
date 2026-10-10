%% ------------------------------------------------------------------------
% Optical Flow RAFT -- local desktop version
% Based on raftmatlabsideview_ARC.m, but intended for interactive PC tests.
% This file does not modify the cluster ARC runner.
%% ------------------------------------------------------------------------
clear; clc; close all;

%% ------------------------- USER SETTINGS -------------------------------
% Set this to the total number of video frames to read (including the first
% reference frame). Use Inf to process from startFrame to the video end.
startFrame = 2000;               % Original image number; current masks cover 2000-2100.
numFramesToProcess = 10;

showFigures = true;
savePlots = true;
raftIters = 8;
raftTolerance = 1e-6;
accelMode = "auto";

% Optional: set the full path here to skip the file-selection dialog.
% Leave this as "" to choose a video interactively when the script runs.
videoPath = "E:\Sept 2026 Flowfield data\P10S30\P10S30_3_472_40lpm.avi";

% Core masks from test_vapour_masks_100frames.m for this same video.
vapourMaskFile = "D:\Working codes\RAFT code\optical_flow_RAFT_MATLAB\test100frames\P10S30_3_472_40lpm\vapour_core_masks.mat";
vapourQCFolder = "D:\Working codes\RAFT code\optical_flow_RAFT_MATLAB\test100frames\stage2";
dilatePx = 5;                    % Disk radius in pixels; 0 keeps the union.
vapourQCColorLimits = [-10 10];   % Fixed u display range (m/s) for all QC frames.

% Cropped videos may not need an ROI or a throat reference.
useFullFrameROI = true;   % true: use every pixel and skip ROI selection
useThroat = false;        % false: skip throat selection and omit its plot line

% Contrast preprocessing applied to every frame sent to RAFT.
% Run once with "none" and once with "clahe" to compare the results.
contrastMode = "none";   % "none" or "clahe"

% Compare a few evenly spaced frames within numFramesToProcess.
showCLAHEPreview = false;
previewOnly = false;      % false: continue with RAFT after the preview
numPreviewFrames = 4;

% Interactive vapor inspection after RAFT finishes (independent of showFigures).
vaporTestingMode = "off";  % "on": select an inspection rectangle on each sample
numVaporSamples = 4;       % number of evenly spaced source frames to inspect
vaporFrameNumbers = [];   % explicit flow indices, e.g. [1 3 7]; overrides sample count
% These inspection numbers are flow indices, from 1 to numFrames-1.
% Index k shows image startFrame+k-1 with forward flow to startFrame+k.
% Draw/resize the rectangle, then press Enter to save it and advance.
% Inspect vaporInspection.vapor_1_frame_1.u (m/s) in the workspace afterward.

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

vaporTestingMode = lower(string(vaporTestingMode));
if ~isscalar(vaporTestingMode) || ~ismember(vaporTestingMode, ["on", "off"])
    error('vaporTestingMode must be "on" or "off".');
end

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

availableFrames = double(v.NumFrames);
validateattributes(startFrame, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'startFrame');
if startFrame >= availableFrames
    error('startFrame must leave at least two images to calculate flow (video has %d frames).', ...
        availableFrames);
end
remainingFrames = availableFrames - startFrame + 1;
if isequal(numFramesToProcess, Inf)
    numFrames = remainingFrames;
else
    validateattributes(numFramesToProcess, {'numeric'}, ...
        {'scalar', 'integer', 'finite', '>=', 2}, mfilename, 'numFramesToProcess');
    numFrames = min(numFramesToProcess, remainingFrames);
end
numFlowFrames = numFrames - 1;
imageFrameNumbers = (startFrame:startFrame + numFrames - 1).';
% Row k records the original source and destination images of u_all(:,:,k).
flowFrameNumbers = [imageFrameNumbers(1:end-1), imageFrameNumbers(2:end)];

% Validate inspection settings before running the expensive RAFT calculation.
vaporSampleFrameNumbers = [];
if vaporTestingMode == "on"
    if isempty(vaporFrameNumbers)
        validateattributes(numVaporSamples, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'numVaporSamples');
        sampleCount = min(numVaporSamples, numFlowFrames);
        if numVaporSamples > numFlowFrames
            warning('Only %d source frames have flow; inspecting all of them.', numFlowFrames);
        end
        if sampleCount == 1
            vaporSampleFrameNumbers = 1;
        else
            vaporSampleFrameNumbers = unique(round(linspace(1, numFlowFrames, sampleCount)));
        end
    else
        validateattributes(vaporFrameNumbers, {'numeric'}, ...
            {'vector', 'integer', 'positive', 'finite', '<=', numFlowFrames}, ...
            mfilename, 'vaporFrameNumbers');
        vaporSampleFrameNumbers = reshape(unique(vaporFrameNumbers, 'stable'), 1, []);
    end
end

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
        frameNumber = imageFrameNumbers(frameIndex);
        rawFrame = read(v, frameNumber);
        original = prepareFrame(rawFrame, "none");
        enhanced = prepareFrame(rawFrame, "clahe");
        previewFig = figure('Name', sprintf('CLAHE preview - frame %d', frameNumber), ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Visible', tern(showFigures, 'on', 'off'));
        tiledlayout(previewFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        % Fixed [0, 1] display range makes the contrast comparison faithful.
        nexttile; imshow(im2double(original), [0 1]);
        title(sprintf('Original grayscale - frame %d', frameNumber));
        nexttile; imshow(im2double(enhanced), [0 1]);
        title(sprintf('CLAHE - frame %d', frameNumber));
        drawnow;
        if savePlots
            previewFile = fullfile(previewDir, sprintf('%s_CLAHE_preview_frame%04d.png', ...
                filename, frameNumber));
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

% Catch mask-file/range mismatches before calculating or saving any flow.
loadVapourMasksForFrames(vapourMaskFile, [v.Height v.Width], imageFrameNumbers);

if canUseGPU
    executionEnv = "gpu";
else
    executionEnv = "cpu";
end

fprintf('Processing video frames %d-%d of %d (%d flow fields).\n', ...
    imageFrameNumbers(1), imageFrameNumbers(end), availableFrames, numFlowFrames);
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
framePrev = im2gray(read(v, startFrame));
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

framePrev = prepareFrame(read(v, startFrame), contrastMode);
fprintf('Running RAFT...\n');
tic;
% RAFT caches the previous image internally. Its first estimate is zero;
% Seed it with startFrame so flow k is startFrame+k-1 -> startFrame+k.
estimateFlow(opticalFlowObj, framePrev, ...
    ExecutionEnvironment=executionEnv, Acceleration=accelMode, ...
    MaxIterations=raftIters, Tolerance=raftTolerance);
for k = 1:numFlowFrames
    frameCurr = prepareFrame(read(v, startFrame + k), contrastMode);
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
completedFlowFrames = numFlowFrames;

outFile = fullfile(matDir, sprintf('%s_velocity_local_%s_%dframes.mat', ...
    filename, contrastLabel, numFrames));
save(outFile, 'u_all', 'v_all', 'u_mean', 'v_mean', 'velMean', ...
    'u_mean_phys', 'v_mean_phys', 'velMean_phys', 'maskROI', ...
    'x_mm', 'y_mm', 'x_throat_pixel', 'y_throat_pixel', 'x_throat_mm', ...
    'mm_per_pixel', 'm_per_pixel', 'fps', 'numFrames', 'numFlowFrames', ...
    'completedFlowFrames', 'startFrame', 'imageFrameNumbers', 'flowFrameNumbers', '-v7.3');
fprintf('Saved local results: %s\n', outFile);

%% ------------------------ VAPOUR MASK QC -------------------------------
% Match core masks by ORIGINAL image number, never by their cell position.
[vapourCore, coreLocations] = loadVapourMasksForFrames( ...
    vapourMaskFile, [H W], imageFrameNumbers);
validateattributes(dilatePx, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'dilatePx');
validateattributes(vapourQCColorLimits, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite'}, mfilename, 'vapourQCColorLimits');
assert(vapourQCColorLimits(1) < vapourQCColorLimits(2), ...
    'vapourQCColorLimits must be [lower upper] with lower < upper.');
if dilatePx > 0
    vapourDisk = strel('disk', dilatePx, 0);
end
flowMaskPixelIndices = cell(numFlowFrames, 1);
flowCoreAreaPixels = zeros(numFlowFrames, 1);
flowMaskedAreaPixels = zeros(numFlowFrames, 1);
for k = 1:numFlowFrames
    % Store core(n) OR core(n+1) before dilation, using linear indices.
    flowCoreMask = false(H, W);
    flowCoreMask(vapourCore.maskPixelIndices{coreLocations(k)}) = true;
    flowCoreMask(vapourCore.maskPixelIndices{coreLocations(k + 1)}) = true;
    flowMaskPixelIndices{k} = find(flowCoreMask);
    flowCoreAreaPixels(k) = nnz(flowCoreMask);
    flowDilatedMask = flowCoreMask;
    if dilatePx > 0
        flowDilatedMask = imdilate(flowCoreMask, vapourDisk);
    end
    flowMaskedAreaPixels(k) = nnz(flowDilatedMask);
end
% Rank by the area of the dilated flow-frame mask; ties retain flow order.
[~, vapourAreaOrder] = sort(flowMaskedAreaPixels, 'descend');
qcFlowIndices = vapourAreaOrder(1:min(10, numFlowFrames));
if ~isfolder(vapourQCFolder), mkdir(vapourQCFolder); end
imageSize = [H W];
flowMaskFile = fullfile(vapourQCFolder, ...
    sprintf('%s_vapour_flow_masks_local_%s_%dframes.mat', filename, contrastLabel, numFrames));
% Reconstruct the undilated union: mask = false(imageSize);
% mask(flowMaskPixelIndices{k}) = true. Row k of flowFrameNumbers gives its pair.
save(flowMaskFile, 'flowMaskPixelIndices', 'flowFrameNumbers', 'imageSize', ...
    'flowCoreAreaPixels', 'flowMaskedAreaPixels', 'qcFlowIndices', 'dilatePx', ...
    'vapourQCColorLimits', 'vapourMaskFile', 'videoPath', 'startFrame', '-v7.3');
for qcIndex = 1:numel(qcFlowIndices)
    k = qcFlowIndices(qcIndex);
    flowDilatedMask = false(H, W);
    flowDilatedMask(flowMaskPixelIndices{k}) = true;
    if dilatePx > 0
        flowDilatedMask = imdilate(flowDilatedMask, vapourDisk);
    end
    originalGray = prepareFrame(read(v, flowFrameNumbers(k, 1)), "none");
    qcFile = fullfile(vapourQCFolder, ...
        sprintf('%s_%s_vapour_QC_%02d_flow%06d_frames%06d_%06d.png', ...
        filename, contrastLabel, qcIndex, k, flowFrameNumbers(k, 1), flowFrameNumbers(k, 2)));
    saveVapourQCFrame(originalGray, u_all(:, :, k), flowDilatedMask, ...
        flowFrameNumbers(k, :), k, flowMaskedAreaPixels(k), vapourQCColorLimits, qcFile);
end
fprintf('Saved vapour flow masks: %s\n', flowMaskFile);
fprintf('Saved %d vapour QC images: %s\n', numel(qcFlowIndices), vapourQCFolder);

%% ------------------- INTERACTIVE VAPOR INSPECTION -----------------------
% Entries stay in the base workspace and are saved after every confirmed ROI.
% Each entry contains the cropped u matrix, original grayscale pixels, pixel
% coordinates, and frame-pair metadata. Image rows increase downward; no flip
% or interpolation is applied. Pixels outside maskROI are NaN in the u crop.
vaporInspection = struct();
if vaporTestingMode == "on"
    vaporInspectionFile = fullfile(matDir, ...
        sprintf('%s_vapor_inspection_local_%s_%dframes.mat', ...
        filename, contrastLabel, numFrames));

    % Use one color range for all samples, so color changes reflect velocity.
    vaporColorLimits = [Inf -Inf];
    for frameNumber = vaporSampleFrameNumbers
        sampleU = u_all(:, :, frameNumber);
        sampleValues = sampleU(maskROI & isfinite(sampleU));
        if ~isempty(sampleValues)
            vaporColorLimits(1) = min(vaporColorLimits(1), double(min(sampleValues)));
            vaporColorLimits(2) = max(vaporColorLimits(2), double(max(sampleValues)));
        end
    end
    if any(~isfinite(vaporColorLimits))
        vaporColorLimits = [-1 1];
    elseif vaporColorLimits(1) == vaporColorLimits(2)
        padding = max(1e-6, abs(vaporColorLimits(1)) * 0.01);
        vaporColorLimits = vaporColorLimits + [-padding padding];
    end

    fprintf('\nVapor inspection: %d samples. Draw a rectangle on the top image,\n', ...
        numel(vaporSampleFrameNumbers));
    fprintf('adjust it, and press Enter to save and advance. Close the window to stop.\n');
    for sampleIndex = 1:numel(vaporSampleFrameNumbers)
        flowIndex = vaporSampleFrameNumbers(sampleIndex);
        frameNumber = flowFrameNumbers(flowIndex, 1);
        originalGray = prepareFrame(read(v, frameNumber), "none");
        [inspection, confirmed] = inspectVaporFrame(originalGray, ...
            u_all(:, :, flowIndex), maskROI, frameNumber, sampleIndex, ...
            numel(vaporSampleFrameNumbers), vaporColorLimits, mm_per_pixel);
        if ~confirmed
            fprintf('Vapor inspection stopped; previously saved areas are retained.\n');
            break;
        end
        inspection.flow_index = flowIndex;
        fieldName = sprintf('vapor_%d_frame_%d', sampleIndex, frameNumber);
        vaporInspection.(fieldName) = inspection;
        save(vaporInspectionFile, 'vaporInspection', 'vaporSampleFrameNumbers', ...
            'vaporTestingMode', 'videoPath', 'contrastMode', 'mm_per_pixel', ...
            'm_per_pixel', 'fps', '-v7.3');
        fprintf('Saved vaporInspection.%s.u: %d rows x %d columns (m/s).\n', ...
            fieldName, size(inspection.u, 1), size(inspection.u, 2));
    end
    % Also keep the inspection alongside the full velocity results.
    save(outFile, 'vaporInspection', 'vaporSampleFrameNumbers', 'vaporTestingMode', '-append');
    if ~isempty(fieldnames(vaporInspection))
        fprintf('Saved vapor inspection results: %s\n', vaporInspectionFile);
    end
end

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
    frame = adapthisteq(frame);
end
end

function [vapourCore, coreLocations] = loadVapourMasksForFrames(maskFile, imageSize, imageFrameNumbers)
assert(isfile(maskFile), ['Vapour core mask file not found: %s\n' ...
    'Set vapourMaskFile to masks generated for this same video.'], maskFile);
vapourCore = load(maskFile, 'frameNumbers', 'maskPixelIndices', 'imageSize');
assert(all(isfield(vapourCore, {'frameNumbers', 'maskPixelIndices', 'imageSize'})), ...
    'Vapour mask file must contain frameNumbers, maskPixelIndices and imageSize.');
assert(isequal(double(vapourCore.imageSize(:).'), imageSize), ...
    'Vapour core mask image size does not match the flow images.');
coreFrameNumbers = vapourCore.frameNumbers(:);
assert(iscell(vapourCore.maskPixelIndices) && ~isempty(coreFrameNumbers) && ...
    numel(coreFrameNumbers) == numel(vapourCore.maskPixelIndices), ...
    'Vapour frameNumbers and maskPixelIndices must have one entry per image.');
assert(numel(unique(coreFrameNumbers)) == numel(coreFrameNumbers), ...
    'Vapour frameNumbers must be unique for unambiguous mask matching.');
[haveCoreMasks, coreLocations] = ismember(imageFrameNumbers, coreFrameNumbers);
if any(~haveCoreMasks)
    error('vapourQC:MissingFrames', ...
        ['Missing vapour core masks for original image frame(s): %s\n' ...
        'Requested image range: %d-%d. Available mask range: %d-%d.\n' ...
        'Set startFrame/numFramesToProcess to a covered range, or regenerate masks for the requested images.'], ...
        mat2str(imageFrameNumbers(~haveCoreMasks).'), ...
        imageFrameNumbers(1), imageFrameNumbers(end), min(coreFrameNumbers), max(coreFrameNumbers));
end
end

function saveVapourQCFrame(originalGray, uFrame, dilatedMask, framePair, ...
    flowIndex, maskedArea, colorLimits, outputFile)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 800]);
figureCleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
grayAxes = nexttile(layout);
imshow(im2double(originalGray), [0 1], 'Parent', grayAxes);
colormap(grayAxes, gray(256));
axis(grayAxes, 'image'); axis(grayAxes, 'ij');
title(grayAxes, sprintf('Grayscale image %d | dilated mask area: %d pixels', ...
    framePair(1), maskedArea));
flowAxes = nexttile(layout);
imagesc(flowAxes, uFrame, colorLimits);
axis(flowAxes, 'image'); axis(flowAxes, 'ij');
colormap(flowAxes, turbo(256));
flowColorbar = colorbar(flowAxes);
flowColorbar.Label.String = 'u (m/s)';
title(flowAxes, sprintf('Flow %d | u: image %d -> %d', ...
    flowIndex, framePair(1), framePair(2)));
boundaries = bwboundaries(dilatedMask, 8, 'holes');
for ax = [grayAxes flowAxes]
    hold(ax, 'on');
    for boundaryIndex = 1:numel(boundaries)
        boundary = boundaries{boundaryIndex};
        plot(ax, boundary(:, 2), boundary(:, 1), 'r-', 'LineWidth', 1);
    end
    hold(ax, 'off');
    xlabel(ax, 'X (pixel column)'); ylabel(ax, 'Y (pixel row)');
end
linkaxes([grayAxes flowAxes], 'xy');
exportgraphics(fig, outputFile, 'Resolution', 150);
end

function [inspection, confirmed] = inspectVaporFrame(originalGray, uFrame, ...
    maskROI, frameNumber, sampleIndex, sampleCount, colorLimits, mmPerPixel)
inspection = struct();
confirmed = false;
[height, width] = size(originalGray);
uFrame(~maskROI) = NaN;

fig = figure('Name', sprintf('Vapor inspection %d/%d - frame %d', ...
    sampleIndex, sampleCount, frameNumber), 'NumberTitle', 'off', ...
    'Color', 'w', 'Visible', 'on', 'Position', [100 100 1000 800]);
% Clean up this window on completion or error without closing other plots.
figureCleanup = onCleanup(@() deleteInspectionFigure(fig));
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(layout, sprintf('Sample %d/%d: draw/resize an area above, then press Enter', ...
    sampleIndex, sampleCount));

grayAxes = nexttile(layout);
imshow(im2double(originalGray), [0 1], 'Parent', grayAxes);
colormap(grayAxes, gray(256));
axis(grayAxes, 'on');
title(grayAxes, sprintf('Original grayscale - source frame %d', frameNumber));
xlabel(grayAxes, 'X (pixel column)'); ylabel(grayAxes, 'Y (pixel row)');

flowAxes = nexttile(layout);
flowImage = imagesc(flowAxes, uFrame, colorLimits);
flowImage.AlphaData = isfinite(uFrame);
axis(flowAxes, 'image');
set(flowAxes, 'YDir', 'reverse', 'Color', [0.85 0.85 0.85]);
colormap(flowAxes, turbo(256));
flowColorbar = colorbar(flowAxes);
flowColorbar.Label.String = 'u (m/s)';
title(flowAxes, sprintf('u component - forward flow %d -> %d', frameNumber, frameNumber + 1));
xlabel(flowAxes, 'X (pixel column)'); ylabel(flowAxes, 'Y (pixel row)');
linkaxes([grayAxes flowAxes], 'xy');
drawnow;

% Releasing the mouse finishes drawing; the rectangle remains editable until
% Enter. Its outline is mirrored below while it is moved or resized.
try
    inspectionROI = drawrectangle(grayAxes, 'Color', 'r', 'FaceAlpha', 0, ...
        'Rotatable', false, 'Deletable', false, 'DrawingArea', [0.5 0.5 width height]);
catch exception
    if ~isgraphics(fig), return; end
    rethrow(exception);
end
if ~isgraphics(fig) || isempty(inspectionROI) || ~isvalid(inspectionROI) || ...
        isempty(inspectionROI.Position)
    return;
end
flowOutline = rectangle(flowAxes, 'Position', inspectionROI.Position, ...
    'EdgeColor', 'r', 'LineWidth', 1.5);
addlistener(inspectionROI, 'MovingROI', ...
    @(~, event) set(flowOutline, 'Position', event.CurrentPosition));
addlistener(inspectionROI, 'ROIMoved', ...
    @(~, event) set(flowOutline, 'Position', event.CurrentPosition));
fig.WindowKeyPressFcn = @(source, event) confirmVaporROI(source, event, inspectionROI);

% Only Enter accepts the ROI; other keys and double-clicks do not advance.
uiwait(fig);
if ~isgraphics(fig) || ~isvalid(inspectionROI)
    return;
end
position = inspectionROI.Position;
% Include pixel centers inside the drawn rectangle, clipped to image bounds.
columns = max(1, ceil(position(1))):min(width, floor(position(1) + position(3)));
rows = max(1, ceil(position(2))):min(height, floor(position(2) + position(4)));
inspection.frame_number = frameNumber;
inspection.next_frame_number = frameNumber + 1;
inspection.flow_index = frameNumber;
inspection.roi_position_pixels = position; % [x y width height]
inspection.row_indices = rows;
inspection.column_indices = columns;
inspection.x_mm = (columns - 1) * mmPerPixel;
inspection.y_mm = (rows - 1) * mmPerPixel; % image coordinates, downward positive
inspection.u = uFrame(rows, columns);     % signed horizontal velocity, m/s
inspection.grayscale = originalGray(rows, columns);
inspection.maskROI = maskROI(rows, columns);
inspection.velocity_units = 'm/s';
confirmed = true;
end

function confirmVaporROI(fig, event, inspectionROI)
if ~ismember(event.Key, {'return', 'enter'}) || ~isvalid(inspectionROI)
    return;
end
position = inspectionROI.Position;
% Do not accept an empty or subpixel rectangle containing no pixel centers.
if numel(position) ~= 4 || any(~isfinite(position)) || ...
        any(position(3:4) <= 0) || ...
        ceil(position(1)) > floor(position(1) + position(3)) || ...
        ceil(position(2)) > floor(position(2) + position(4))
    title(inspectionROI.Parent, 'Select an area containing at least one pixel, then press Enter');
    return;
end
uiresume(fig);
end

function deleteInspectionFigure(fig)
if isgraphics(fig), delete(fig); end
end
