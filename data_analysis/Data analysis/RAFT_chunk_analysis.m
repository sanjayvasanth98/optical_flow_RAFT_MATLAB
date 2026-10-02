%% RAFT_chunk_analysis.m
% Analyze selected frames from RAFT velocity MAT files.
% Read frames in small blocks to limit memory use.
% Flow frame k describes the velocity between video frames k and k+1.

clearvars; clc;

%% 1. Choose velocity files and case labels
% Add one MAT file per case. Use the labels in profile legends and map titles.
matPaths = { ...
    'E:\Sept 2026 Flowfield data\Processed data\P10S20\mat files\P10S20_3_475_40lpm_velocity.mat';
    % 'E:\path\to\case2_velocity.mat';
};

% Keep labels in the same order as the MAT files.
caseLabels = { ...
    'Rough5';  % Change to the roughness label you want shown.
    % 'Case 2 roughness';
};

%% 2. Choose the phase and frames to analyze
% Set the phase label for this run: pre-inception, inception, or desinence.
analysisPhase = "Inception";

% Give each MAT file one [firstFrame lastFrame] row, in the same order.
% Both endpoints are included; [2000 5000] processes frames 2000:5000.
frameRanges = [ ...
    2000, 5000;  % case 1
    % 7000, 9000;  % case 2 (uncomment when case 2 path is added)
];

%% 3. Set reading and output options
% Limit frames per disk read, choose whether to make a frame summary,
% and set where plots are saved.
maxFramesPerRead = 100;  % Upper bound per disk read; lower this for less RAM.
makeFrameSummary = false;  % Set true to also scan the full ROI for a time trace.
outputDir = fullfile(fileparts(mfilename('fullpath')), 'results');
% Older MAT files may lack a completion counter. Set true only when you
% know that every allocated frame in those files was written successfully.
assumeCompleteWhenNoCounter = false;

%% 4. Check the inputs and prepare file information
% Confirm that the files, labels, and frame ranges are usable.
% Read file metadata needed for the analysis without loading all frames.
assert(~isempty(matPaths), 'Add at least one MAT-file path.');
assert(numel(caseLabels) == numel(matPaths), ...
    'Add one case label for each MAT-file path.');
assert(all(strlength(string(caseLabels)) > 0), ...
    'Case labels must not be empty.');
analysisPhase = lower(string(analysisPhase));
assert(isscalar(analysisPhase) && ismember(analysisPhase, ...
    ["pre-inception", "inception", "desinence"]), ...
    'analysisPhase must be pre-inception, inception, or desinence.');
assert(isscalar(maxFramesPerRead) && isnumeric(maxFramesPerRead) && ...
    isfinite(maxFramesPerRead) && maxFramesPerRead >= 1 && ...
    maxFramesPerRead == fix(maxFramesPerRead), ...
    'maxFramesPerRead must be a positive integer.');
assert(isnumeric(frameRanges) && size(frameRanges,2) == 2 && ...
    size(frameRanges,1) == numel(matPaths), ...
    'frameRanges needs one [firstFrame lastFrame] row per MAT-file path.');
fileInfo = cell(numel(matPaths), 1);
for fileIndex = 1:numel(matPaths)
    fileInfo{fileIndex} = inspectVelocityFile(matPaths{fileIndex}, ...
        assumeCompleteWhenNoCounter);
    firstFrame = frameRanges(fileIndex,1);
    lastFrame = frameRanges(fileIndex,2);
    info = fileInfo{fileIndex};
    assert(isfinite(firstFrame) && isfinite(lastFrame) && ...
        firstFrame == fix(firstFrame) && lastFrame == fix(lastFrame) && ...
        firstFrame >= 1 && lastFrame >= firstFrame && ...
        lastFrame <= info.completedFrames, ...
        'Case %d range must lie within completed flow frames 1:%d.', ...
        fileIndex, info.completedFrames);
end

%% 5. Optionally summarize each frame across the ROI
% If makeFrameSummary is true, calculate mean U, V, and speed for every
% selected frame using valid pixels inside the ROI. Print the first 10 rows.
% U points right and V points down in the source image.
frameSummary = table();
if makeFrameSummary
nSelected = sum(frameRanges(:,2) - frameRanges(:,1) + 1);
caseIndex = zeros(nSelected,1);
phase = repmat(analysisPhase,nSelected,1);
flowFrame = zeros(nSelected,1);
time_s = zeros(nSelected,1);
meanU_mps = nan(nSelected,1);
meanV_mps = nan(nSelected,1);
meanSpeed_mps = nan(nSelected,1);
validROIPixels = zeros(nSelected,1);
outRow = 0;

for i = 1:numel(matPaths)
    info = fileInfo{i};
    fprintf('Case %d (%s, %s): %s, flow frames %d:%d\n', ...
        i, char(string(caseLabels{i})), char(analysisPhase), char(string(matPaths{i})), ...
        frameRanges(i,1), frameRanges(i,2));

    for k0 = frameRanges(i,1):maxFramesPerRead:frameRanges(i,2)
        k1 = min(k0 + maxFramesPerRead - 1, frameRanges(i,2));
        [U, V] = readVelocityBlock(info, k0, k1);

        % Add future inception/desinence calculations here while this case's
        % selected U/V block is available. Columns correspond to k0:k1.

        for j = 1:(k1-k0+1)
            outRow = outRow + 1;
            u = double(U(:,j));
            v = double(V(:,j));
            good = isfinite(u) & isfinite(v);

            caseIndex(outRow) = i;
            flowFrame(outRow) = k0+j-1;
            time_s(outRow) = (flowFrame(outRow)-0.5) / info.fps;
            validROIPixels(outRow) = nnz(good);
            if any(good)
                meanU_mps(outRow) = mean(u(good));
                meanV_mps(outRow) = mean(v(good));
                meanSpeed_mps(outRow) = mean(hypot(u(good),v(good)));
            end
        end
        clear U V
    end
end

frameSummary = table(caseIndex, phase, flowFrame, time_s, ...
    meanU_mps, meanV_mps, meanSpeed_mps, validROIPixels);
disp(frameSummary(1:min(10,height(frameSummary)),:));
fprintf('Summarized %d selected flow frames. Full velocity blocks were released after each read.\n', ...
    height(frameSummary));
end

%% 6. Plot mean speed profiles at selected stations
% Choose distances downstream from the throat. At each station, average
% speed over the selected frames at every height. Plot one figure per
% station, with all cases overlaid. The plot's y axis points upward.

profileOffset_mm = [1.5 3.0 4.5];

assert(isnumeric(profileOffset_mm) && isvector(profileOffset_mm) && ...
    ~isempty(profileOffset_mm) && all(isfinite(profileOffset_mm)) && ...
    all(profileOffset_mm >= 0), ...
    'profileOffset_mm must be a nonempty vector of nonnegative distances in mm.');

profileInfo = cell(numel(matPaths),numel(profileOffset_mm));
for i = 1:numel(matPaths)
    for stationIndex = 1:numel(profileOffset_mm)
        profileInfo{i,stationIndex} = prepareProfileColumn( ...
            fileInfo{i},profileOffset_mm(stationIndex),matPaths{i});
    end
end

% Rows of meanSpeedProfiles are cases; columns are stations.
meanSpeedProfiles = cell(numel(matPaths),numel(profileOffset_mm));
maxProfileSpeed = 0;
minProfileY = inf;
maxProfileY = -inf;
for stationIndex = 1:numel(profileOffset_mm)
    for i = 1:numel(matPaths)
        info = profileInfo{i,stationIndex};
        profileSum = zeros(numel(info.profileImageRows),1);
        profileCount = zeros(numel(info.profileImageRows),1);
        fprintf('Profile for %s: flow frames %d:%d, x = %.4f mm (%.4f mm after throat).\n', ...
            char(string(caseLabels{i})),frameRanges(i,1),frameRanges(i,2), ...
            info.profileX_mm,info.actualOffset_mm);
        for k0 = frameRanges(i,1):maxFramesPerRead:frameRanges(i,2)
            k1 = min(k0+maxFramesPerRead-1,frameRanges(i,2));
            [uColumn,vColumn] = readProfileBlock(info,k0,k1);
            uColumn = double(uColumn);
            vColumn = double(vColumn);
            validColumn = isfinite(uColumn) & isfinite(vColumn);
            speedColumn = hypot(uColumn,vColumn);
            speedColumn(~validColumn) = 0;
            profileSum = profileSum + sum(speedColumn,2);
            profileCount = profileCount + sum(validColumn,2);
            clear uColumn vColumn validColumn speedColumn
        end
        speedByImageRow = nan(size(info.maskROI,1),1);
        hasSamples = profileCount > 0;
        assert(any(hasSamples), ...
            'Case %d has no finite U/V samples at station %.3f mm.', ...
            i,profileOffset_mm(stationIndex));
        speedByImageRow(info.profileImageRows(hasSamples)) = ...
            profileSum(hasSamples) ./ profileCount(hasSamples);
        speedPhysical = flipud(speedByImageRow);
        meanSpeedProfiles{i,stationIndex} = struct( ...
            'label',string(caseLabels{i}), ...
            'frameRange',frameRanges(i,:), ...
            'x_mm',info.profileX_mm, ...
            'offsetFromThroat_mm',info.actualOffset_mm, ...
            'y_mm',info.y_mm(:), ...
            'meanSpeed_mps',speedPhysical);
        validProfileRows = isfinite(speedPhysical);
        maxProfileSpeed = max(maxProfileSpeed,max(speedPhysical(validProfileRows)));
        minProfileY = min(minProfileY,min(info.y_mm(validProfileRows)));
        maxProfileY = max(maxProfileY,max(info.y_mm(validProfileRows)));
    end
end

% A consistent speed scale makes profiles at different stations comparable.
if maxProfileSpeed > 0
    speedAxisMax = 1.05 * maxProfileSpeed;
else
    speedAxisMax = 1;
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Reference order: open black circles, green crosses, solid blue, dashed red.
% Keep each case's style the same at every station.
for stationIndex = 1:numel(profileOffset_mm)
    fig = figure('Color','w','Position',[100 100 760 650]);
    ax = axes(fig);
    hold(ax,'on');
    for i = 1:numel(matPaths)
        profile = meanSpeedProfiles{i,stationIndex};
        plotProfileSeries(ax,profile,mod(i-1,4)+1,char(string(caseLabels{i})));
    end
    styleProfileAxes(ax,sprintf('Mean speed | throat + %.2f mm | %s', ...
        profileOffset_mm(stationIndex),char(analysisPhase)), ...
        speedAxisMax,[minProfileY maxProfileY]);
    lgd = legend(ax,'show','Location','best','Interpreter','none','Box','off');
    set(lgd,'FontName','Times New Roman','FontSize',11);

    plotName = sprintf('mean_speed_profile_xplus_%.2fmm_%s', ...
        profileOffset_mm(stationIndex),char(analysisPhase));
    if numel(profileOffset_mm) > 1
        plotName = sprintf('%s_station%02d',plotName,stationIndex);
    end
    profilePlotBase = fullfile(outputDir,plotName);
    profilePlotPath = [profilePlotBase '.png'];
    profileFigPath = [profilePlotBase '.fig'];
    exportgraphics(fig,profilePlotPath,'Resolution',600);
    savefig(fig,profileFigPath);
    fprintf('Saved comparison plots: %s and %s\n',profilePlotPath,profileFigPath);
end

%% 7. Plot mean speed maps with station lines
% Average speed over the selected frames at every valid ROI pixel.
% Save one map per case with vertical lines showing the Section 6 stations.
% Read a few full-ROI frames at a time to limit memory use.
mapFramesPerRead = min(maxFramesPerRead,5);
meanSpeedMaps = cell(numel(matPaths),1);
maxMapSpeed = 0;
for i = 1:numel(matPaths)
    info = fileInfo{i};
    nROIPixels = nnz(info.maskROI);
    speedSum = zeros(nROIPixels,1);
    validCount = zeros(nROIPixels,1);
    fprintf('Mean speed map for %s: flow frames %d:%d.\n', ...
        char(string(caseLabels{i})),frameRanges(i,1),frameRanges(i,2));
    for k0 = frameRanges(i,1):mapFramesPerRead:frameRanges(i,2)
        k1 = min(k0+mapFramesPerRead-1,frameRanges(i,2));
        [U,V] = readVelocityBlock(info,k0,k1);
        validBlock = isfinite(U) & isfinite(V);
        speedBlock = hypot(U,V);
        speedBlock(~validBlock) = 0;
        speedSum = speedSum + sum(double(speedBlock),2);
        validCount = validCount + sum(validBlock,2);
        clear U V validBlock speedBlock
    end
    hasSamples = validCount > 0;
    assert(any(hasSamples), ...
        'Case %d has no finite U/V samples in the selected frame range.',i);
    meanSpeedROI = nan(nROIPixels,1,'single');
    meanSpeedROI(hasSamples) = single(speedSum(hasSamples) ./ validCount(hasSamples));
    speedImage = nan(size(info.maskROI),'single');
    speedImage(info.maskROI) = meanSpeedROI;
    meanSpeedMaps{i} = flipud(speedImage);
    maxMapSpeed = max(maxMapSpeed,double(max(meanSpeedROI(hasSamples))));
    clear speedSum validCount meanSpeedROI speedImage
end

% Use the same color scale across cases for direct comparison.
if maxMapSpeed > 0
    mapColorMax = double(maxMapSpeed);
else
    mapColorMax = 1;
end
stationColors = [ ...
    0.95 0.70 0.00;  % gold
    0.95 0.22 0.45;  % pink
    0.00 0.65 0.78;  % cyan
    0.68 0.28 0.88;  % purple
    0.97 0.46 0.08]; % orange
for i = 1:numel(matPaths)
    info = fileInfo{i};
    mapPhysical = meanSpeedMaps{i};
    fig = figure('Color','w','Position',[100 100 1200 800]);
    ax = axes(fig);
    imageHandle = imagesc(ax,info.x_mm,info.y_mm,mapPhysical);
    set(imageHandle,'AlphaData',single(isfinite(mapPhysical)), ...
        'HandleVisibility','off');
    set(ax,'YDir','normal','FontName','Times New Roman','FontSize',12, ...
        'LineWidth',1,'TickDir','in','Box','on','Layer','top', ...
        'XColor',[0.18 0.18 0.18],'YColor',[0.18 0.18 0.18]);
    axis(ax,'image');
    colormap(ax,turbo(256));
    caxis(ax,[0 mapColorMax]);
    hold(ax,'on');
    for stationIndex = 1:numel(profileOffset_mm)
        stationX = profileInfo{i,stationIndex}.profileX_mm;
        stationColor = stationColors(mod(stationIndex-1,size(stationColors,1))+1,:);
        plot(ax,[stationX stationX],[info.y_mm(1) info.y_mm(end)],'-', ...
            'Color',[0.08 0.08 0.08],'LineWidth',3.4, ...
            'HandleVisibility','off');
        plot(ax,[stationX stationX],[info.y_mm(1) info.y_mm(end)],'--', ...
            'Color',stationColor,'LineWidth',2, ...
            'DisplayName',sprintf('Station %d: +%.2f mm', ...
                stationIndex,profileOffset_mm(stationIndex)));
    end
    xlabel(ax,'x (mm)','FontName','Times New Roman','FontSize',14);
    ylabel(ax,'y (mm)','FontName','Times New Roman','FontSize',14);
    title(ax,sprintf('%s | Mean speed | frames %d:%d | %s', ...
        char(string(caseLabels{i})),frameRanges(i,1),frameRanges(i,2), ...
        char(analysisPhase)), ...
        'FontName','Times New Roman','FontSize',15, ...
        'FontWeight','normal','Interpreter','none');
    cb = colorbar(ax,'eastoutside');
    set(cb,'FontName','Times New Roman','FontSize',11,'TickDirection','in');
    cb.Label.String = 'Mean speed (m/s)';
    cb.Label.FontName = 'Times New Roman';
    cb.Label.FontSize = 13;
    lgd = legend(ax,'show','Location','southoutside', ...
        'Orientation','horizontal','Interpreter','none','Box','off');
    set(lgd,'FontName','Times New Roman','FontSize',11);

    safeLabel = regexprep(char(string(caseLabels{i})),'[^A-Za-z0-9_-]','_');
    mapPlotBase = fullfile(outputDir,sprintf( ...
        'case%02d_%s_mean_speed_map_%s_frames%d-%d', ...
        i,safeLabel,char(analysisPhase),frameRanges(i,1),frameRanges(i,2)));
    mapPngPath = [mapPlotBase '.png'];
    mapFigPath = [mapPlotBase '.fig'];
    exportgraphics(fig,mapPngPath,'Resolution',600);
    savefig(fig,mapFigPath);
    fprintf('Saved mean speed map: %s and %s\n',mapPngPath,mapFigPath);
end

%% 8. Notes for future calculations
% This section only shows how to rebuild an image from ROI data; it runs
% no analysis. U/V rows follow find(info.maskROI), with one frame per column.
% To rebuild one image-coordinate frame:
%   uImage = nan(size(info.maskROI), 'single');
%   uImage(info.maskROI) = U(:,1);
% To get physical coordinates (origin at lower left):
%   uPhysical = flipud(uImage);  vPhysical = -flipud(vImage);

%% Helper functions used by the sections above
function info = inspectVelocityFile(matPath, assumeCompleteWhenNoCounter)
    assert(isfile(matPath), 'MAT file not found: %s', char(string(matPath)));
    % Partial matfile reads require v7.3/HDF5 storage. Reject older MAT
    % formats before any velocity access, since they may load a whole array.
    try
        h5info(matPath, '/u_all');
    catch
        error('Expected a v7.3 MAT file with /u_all: %s', char(string(matPath)));
    end
    M = matfile(matPath);  % Read-only handle; does not load the velocity arrays.
    vars = who(M);
    required = {'u_all','v_all','maskROI','fps','mm_per_pixel'};
    assert(all(ismember(required,vars)), ...
        'MAT file is missing velocity or calibration metadata: %s', ...
        char(string(matPath)));

    maskROI = logical(M.maskROI);
    assert(ismatrix(maskROI) && any(maskROI(:)), ...
        'maskROI must be a nonempty 2-D mask: %s', char(string(matPath)));
    fps = double(M.fps);
    assert(isscalar(fps) && isfinite(fps) && fps > 0, ...
        'fps must be positive: %s', char(string(matPath)));
    mmPerPixel = double(M.mm_per_pixel);
    assert(isscalar(mmPerPixel) && isfinite(mmPerPixel) && mmPerPixel > 0, ...
        'mm_per_pixel must be positive: %s', char(string(matPath)));

    [heightPx,widthPx] = size(maskROI);
    if ismember('x_mm',vars)
        x_mm = double(M.x_mm);
        x_mm = x_mm(:);
    else
        x_mm = (0:widthPx-1)' * mmPerPixel;
    end
    if ismember('y_mm',vars)
        y_mm = double(M.y_mm);
        y_mm = y_mm(:);
    else
        y_mm = (0:heightPx-1)' * mmPerPixel;
    end
    assert(numel(x_mm) == widthPx && all(isfinite(x_mm)) && ...
        all(diff(x_mm) > 0) && numel(y_mm) == heightPx && ...
        all(isfinite(y_mm)) && all(diff(y_mm) > 0), ...
        'x_mm or y_mm does not match the velocity grid: %s', char(string(matPath)));
    if ismember('x_throat_mm',vars)
        xThroat_mm = double(M.x_throat_mm);
    elseif ismember('x_throat_pixel',vars)
        xThroat_mm = double(M.x_throat_pixel) * mmPerPixel;
    else
        error('x_throat_mm is missing; the profile location cannot be set: %s', ...
            char(string(matPath)));
    end
    assert(isscalar(xThroat_mm) && isfinite(xThroat_mm), ...
        'x_throat_mm must be finite: %s', char(string(matPath)));

    uSize = size(M,'u_all');
    vSize = size(M,'v_all');
    assert(isequal(uSize,vSize), 'u_all and v_all have different sizes: %s', ...
        char(string(matPath)));
    packed = ismember('instantaneousStorage',vars);
    if packed
        assert(strcmp(M.instantaneousStorage,'roi_pixels_by_frame_v1'), ...
            'Unknown instantaneousStorage: %s', char(string(matPath)));
        assert(numel(uSize) == 2 && uSize(1) == nnz(maskROI), ...
            'Packed velocity dimensions do not match maskROI: %s', char(string(matPath)));
        allocatedFrames = uSize(2);
    else
        assert(uSize(1) == size(maskROI,1) && uSize(2) == size(maskROI,2), ...
            'Full-frame velocity dimensions do not match maskROI: %s', ...
            char(string(matPath)));
        if numel(uSize) < 3
            allocatedFrames = 1;
        else
            allocatedFrames = uSize(3);
        end
    end

    if ismember('completedFlowFrames',vars)
        completedFrames = double(M.completedFlowFrames);
        assert(isscalar(completedFrames) && isfinite(completedFrames) && ...
            completedFrames == fix(completedFrames) && completedFrames >= 0 && ...
            completedFrames <= allocatedFrames, ...
            'Invalid completedFlowFrames in %s', char(string(matPath)));
    else
        assert(assumeCompleteWhenNoCounter, ...
            ['No completedFlowFrames in %s. If this file is complete, set ' ...
             'assumeCompleteWhenNoCounter = true.'], char(string(matPath)));
        completedFrames = allocatedFrames;
    end

    info = struct('M',M, 'maskROI',maskROI, 'packed',packed, ...
        'fps',fps, 'completedFrames',completedFrames, ...
        'x_mm',x_mm, 'y_mm',y_mm, 'xThroat_mm',xThroat_mm);
end

function info = prepareProfileColumn(info,offset_mm,matPath)
    xTarget_mm = info.xThroat_mm + offset_mm;
    assert(xTarget_mm >= info.x_mm(1) && xTarget_mm <= info.x_mm(end), ...
        'Throat + %.3f mm is outside the x grid: %s', ...
        offset_mm,char(string(matPath)));
    [~,xIndex] = min(abs(info.x_mm-xTarget_mm));
    imageRows = find(info.maskROI(:,xIndex));
    assert(~isempty(imageRows), ...
        'No ROI pixels at throat + %.3f mm in %s', ...
        offset_mm,char(string(matPath)));
    roiPixels = find(info.maskROI);
    columnPixels = (xIndex-1)*size(info.maskROI,1) + imageRows;
    [present,roiRows] = ismember(columnPixels,roiPixels);
    assert(all(present), 'Could not map the profile column to ROI storage.');
    assert(all(diff(roiRows) == 1), ...
        'Profile pixels are not contiguous in packed ROI storage.');
    info.profileImageRows = imageRows;
    info.profileROIRows = roiRows;
    info.profileXIndex = xIndex;
    info.profileX_mm = info.x_mm(xIndex);
    info.actualOffset_mm = info.profileX_mm - info.xThroat_mm;
end

function [U,V] = readProfileBlock(info,k0,k1)
    assert(k0 >= 1 && k1 >= k0 && k1 <= info.completedFrames, ...
        'Requested profile frames exceed the committed range.');
    M = info.M;
    nFrames = k1-k0+1;
    if info.packed
        roiRange = info.profileROIRows(1):info.profileROIRows(end);
        U = M.u_all(roiRange,k0:k1);
        V = M.v_all(roiRange,k0:k1);
    else
        % Legacy layout: read just one image column per selected frame.
        if size(M,'u_all',3) == 1
            uColumn = M.u_all(:,info.profileXIndex);
            vColumn = M.v_all(:,info.profileXIndex);
        else
            uColumn = M.u_all(:,info.profileXIndex,k0:k1);
            vColumn = M.v_all(:,info.profileXIndex,k0:k1);
        end
        uColumn = reshape(uColumn,[],nFrames);
        vColumn = reshape(vColumn,[],nFrames);
        U = uColumn(info.profileImageRows,:);
        V = vColumn(info.profileImageRows,:);
    end
end

function [U,V] = readVelocityBlock(info,k0,k1)
    assert(k0 >= 1 && k1 >= k0 && k1 <= info.completedFrames, ...
        'Requested flow frames exceed the committed range.');
    M = info.M;
    if info.packed
        % Stored rows already follow find(maskROI).
        U = M.u_all(:,k0:k1);
        V = M.v_all(:,k0:k1);
    else
        % Legacy files store complete images. Read only these time slices,
        % then retain pixels inside the ROI in the same order as packed files.
        if size(M,'u_all',3) == 1
            uFull = M.u_all(:,:);
            vFull = M.v_all(:,:);
        else
            uFull = M.u_all(:,:,k0:k1);
            vFull = M.v_all(:,:,k0:k1);
        end
        U = reshape(uFull,[],k1-k0+1);
        V = reshape(vFull,[],k1-k0+1);
        U = U(info.maskROI(:),:);
        V = V(info.maskROI(:),:);
    end
end

function plotProfileSeries(ax,profile,styleIndex,seriesLabel)
    validRows = find(isfinite(profile.meanSpeed_mps) & isfinite(profile.y_mm));
    assert(~isempty(validRows), 'Cannot plot a profile without finite values.');
    switch styleIndex
        case 1  % open black circles
            lineColor = [0 0 0];
            lineStyle = '-';
            markerSymbol = 'o';
            markerCount = 25;
            lineWidth = 1.3;
            markerSize = 6;
        case 2  % dense green crosses
            lineColor = [0 0.55 0];
            lineStyle = '-';
            markerSymbol = 'x';
            markerCount = 70;
            lineWidth = 1.4;
            markerSize = 5;
        case 3  % solid blue curve
            lineColor = [0 0 1];
            lineStyle = '-';
            markerSymbol = 'none';
            markerCount = 0;
            lineWidth = 2.5;
            markerSize = 6;
        case 4  % dashed red curve
            lineColor = [1 0 0];
            lineStyle = '--';
            markerSymbol = 'none';
            markerCount = 0;
            lineWidth = 2;
            markerSize = 6;
    end
    markerRows = [];
    if markerCount > 0
        nMarkers = min(markerCount,numel(validRows));
        markerRows = unique(validRows(round(linspace(1,numel(validRows),nMarkers))));
    end
    plot(ax,profile.meanSpeed_mps,profile.y_mm, ...
        'Color',lineColor,'LineStyle',lineStyle,'LineWidth',lineWidth, ...
        'Marker',markerSymbol,'MarkerIndices',markerRows, ...
        'MarkerSize',markerSize,'MarkerFaceColor','w', ...
        'MarkerEdgeColor',lineColor, ...
        'DisplayName',seriesLabel);
end

function styleProfileAxes(ax,plotTitle,speedAxisMax,yLimits)
    fontName = 'Times New Roman';
    set(ax,'FontName',fontName,'FontSize',12,'LineWidth',1, ...
        'TickDir','in','Box','on','Color','w', ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XGrid','on','YGrid','on','GridColor',[0.67 0.67 0.67], ...
        'GridAlpha',0.6,'Layer','top');
    xlim(ax,[0 speedAxisMax]);
    if yLimits(2) > yLimits(1)
        ylim(ax,yLimits);
    end
    xlabel(ax,'Mean speed (m/s)','FontName',fontName,'FontSize',14);
    ylabel(ax,'y (mm)','FontName',fontName,'FontSize',14);
    title(ax,plotTitle,'FontName',fontName,'FontSize',15, ...
        'FontWeight','normal','Interpreter','none');
end
