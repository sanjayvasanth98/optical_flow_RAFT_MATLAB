%% RAFT_chunk_analysis.m
% Inspect one selected flow-frame range per MAT file from
% raftmatlabsideview_ARC.m output, for one phase per run.
% No full u_all/v_all arrays are loaded. Each read is restricted to one
% case's frame range and at most maxFramesPerRead consecutive frames.
% Flow frame k is the velocity between video frames k and k+1 (1-based).

clearvars; clc;

%% 1. MAT-file paths (edit these)
matPaths = { ...
    'E:\Sept 2026 Flowfield data\Processed data\P10S20\mat files\P10S20_3_475_40lpm_velocity.mat';
    % 'E:\path\to\case2_velocity.mat';
};

% Labels appear in the comparison plot legend, in the same order as matPaths.
caseLabels = { ...
    'Rough5';  % Change to the roughness label you want shown.
    % 'Case 2 roughness';
};

%% 2. Phase and frame range for each MAT file (edit these)
% Run one phase at a time. Set this label to "pre-inception", "inception",
% or "desinence" when you change the ranges for the next run.
analysisPhase = "Inception";

% One row per path in matPaths, in the same order. A comma separates the
% start/end frames within a case; a semicolon separates cases.
% Endpoints are inclusive. Example: case 1 uses flow frames 2000:5000.
frameRanges = [ ...
    2000, 5000;  % case 1
    % 7000, 9000;  % case 2 (uncomment when case 2 path is added)
];

%% 3. Read settings
maxFramesPerRead = 100;  % Upper bound per disk read; lower this for less RAM.
makeFrameSummary = false;  % Set true to also scan the full ROI for a time trace.
outputDir = fullfile(fileparts(mfilename('fullpath')), 'results');
% Older MAT files may lack a completion counter. Set true only when you
% know that every allocated frame in those files was written successfully.
assumeCompleteWhenNoCounter = false;

%% 4. Validate paths, frame ranges, and small metadata
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

%% 5. Analysis 1: time trace of ROI-averaged velocity in each case range
% Output contains one small row per selected flow frame. U and V retain
% image coordinates: +U is right and +V is down. Speed is hypot(U,V).
% Only finite U/V pairs inside maskROI contribute to each average.
% This optional section reads the full ROI for the selected frames.
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

%% 6. Mean speed profiles at selected stations downstream of the throat
% Mean speed means mean(hypot(U,V)), matching the writer's velMean field,
% calculated from each case's selected frame range. The plotted y axis is
% physical orientation: y increases upward from the bottom of the image.
% This section reads only the selected station columns and frame ranges.
% Specify one or more distances downstream from each case's throat (mm).
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
palette = [ ...
    0.16 0.34 0.58;  % blue
    0.72 0.36 0.20;  % burnt orange
    0.22 0.49 0.39;  % green
    0.47 0.39 0.61;  % purple
    0.67 0.53 0.24;  % ochre
    0.19 0.52 0.61]; % teal
markers = {'o','s','^','d','v','>','<','p','h'};
caseLineStyles = {'-','--','-.',':'};
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Keep an individual comparison figure for each station.
for stationIndex = 1:numel(profileOffset_mm)
    fig = figure('Color','w','Position',[100 100 1050 720]);
    ax = axes(fig);
    hold(ax,'on');
    for i = 1:numel(matPaths)
        profile = meanSpeedProfiles{i,stationIndex};
        caseColor = palette(mod(i-1,size(palette,1))+1,:);
        caseMarker = markers{mod(i-1,numel(markers))+1};
        plotProfileSeries(ax,profile,caseColor,caseMarker,'-', ...
            char(string(caseLabels{i})));
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

% The combined figure uses color and marker for station, line style for case.
combinedFig = figure('Color','w','Position',[100 100 1200 760]);
combinedAx = axes(combinedFig);
hold(combinedAx,'on');
for stationIndex = 1:numel(profileOffset_mm)
    stationColor = palette(mod(stationIndex-1,size(palette,1))+1,:);
    stationMarker = markers{mod(stationIndex-1,numel(markers))+1};
    for i = 1:numel(matPaths)
        profile = meanSpeedProfiles{i,stationIndex};
        caseLineStyle = caseLineStyles{mod(i-1,numel(caseLineStyles))+1};
        seriesLabel = sprintf('%s | +%.2f mm', ...
            char(string(caseLabels{i})),profileOffset_mm(stationIndex));
        plotProfileSeries(combinedAx,profile,stationColor,stationMarker, ...
            caseLineStyle,seriesLabel);
    end
end
styleProfileAxes(combinedAx,sprintf('Mean speed profiles | %s', ...
    char(analysisPhase)),speedAxisMax,[minProfileY maxProfileY]);
lgd = legend(combinedAx,'show','Location','eastoutside', ...
    'Interpreter','none','Box','off');
set(lgd,'FontName','Times New Roman','FontSize',11);
combinedPlotBase = fullfile(outputDir,sprintf( ...
    'mean_speed_profiles_all_stations_%s',char(analysisPhase)));
combinedPngPath = [combinedPlotBase '.png'];
combinedFigPath = [combinedPlotBase '.fig'];
exportgraphics(combinedFig,combinedPngPath,'Resolution',600);
savefig(combinedFig,combinedFigPath);
fprintf('Saved combined plots: %s and %s\n',combinedPngPath,combinedFigPath);

%% 7. Later analyses
% For future whole-ROI calculations, set makeFrameSummary=true and add them
% inside Section 5's read loop, where U/V are available before being
% cleared. U/V have one ROI pixel per row, ordered as find(info.maskROI),
% and one selected frame per column. Section 6 reads only the station columns.
% To reconstruct one image-coordinate frame:
%   uImage = nan(size(info.maskROI), 'single');
%   uImage(info.maskROI) = U(:,1);
% To get physical coordinates (origin at lower left):
%   uPhysical = flipud(uImage);  vPhysical = -flipud(vImage);

%% Local helpers
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

function plotProfileSeries(ax,profile,lineColor,markerSymbol,lineStyle,seriesLabel)
    validRows = find(isfinite(profile.meanSpeed_mps) & isfinite(profile.y_mm));
    assert(~isempty(validRows), 'Cannot plot a profile without finite values.');
    % Show a small, even sample of markers so dense pixel rows stay legible.
    nMarkers = min(12,numel(validRows));
    markerRows = unique(validRows(round(linspace(1,numel(validRows),nMarkers))));
    plot(ax,profile.meanSpeed_mps,profile.y_mm, ...
        'Color',lineColor,'LineStyle',lineStyle,'LineWidth',1.8, ...
        'Marker',markerSymbol,'MarkerIndices',markerRows, ...
        'MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',lineColor, ...
        'DisplayName',seriesLabel);
end

function styleProfileAxes(ax,plotTitle,speedAxisMax,yLimits)
    fontName = 'Times New Roman';
    set(ax,'FontName',fontName,'FontSize',12,'LineWidth',1, ...
        'TickDir','out','Box','off','Color','w', ...
        'XColor',[0.18 0.18 0.18],'YColor',[0.18 0.18 0.18], ...
        'XGrid','on','YGrid','off','GridColor',[0.84 0.86 0.88], ...
        'GridAlpha',0.35,'Layer','top');
    xlim(ax,[0 speedAxisMax]);
    if yLimits(2) > yLimits(1)
        ylim(ax,yLimits);
    end
    xlabel(ax,'Mean speed (m/s)','FontName',fontName,'FontSize',14);
    ylabel(ax,'y (mm)','FontName',fontName,'FontSize',14);
    title(ax,plotTitle,'FontName',fontName,'FontSize',15, ...
        'FontWeight','normal','Interpreter','none');
end
