%% RAFT_chunk_analysis.m
% Analyze selected frames from RAFT velocity MAT files.
% Read frames in small blocks to limit memory use.
% Flow frame k describes the velocity between video frames k and k+1.

clearvars; clc;
% Analysis and export sections log local start/end timestamps and elapsed seconds.

%% Global parameters
throatHeight_mm = 10;  % H, used to scale spatial coordinates.
throatWidth_mm = 5;     % Channel width into the image plane.
flowRate_Lpm = 40;      % Shared flow rate, or one value per MAT file in order.


%% Analysis
% Choose one or more phases for this run; other configured ranges are ignored.
% Possible options: "pre-inception", "inception", "desinence", "fullcase".
% For all phases: ["pre-inception", "inception", "desinence", "fullcase"].
analysisPhases = ["inception"];

% Set true to run an analysis or false to skip it; no commenting is needed.
analysis1_verticalProfiles = true;
analysis2_instantaneousFrameSaving = false;
analysis3_profileStationVisualization = false;
analysis4_frameSummary = false;

% A true switch runs only in the phases listed for that analysis.
% Example: run all phases but exclude fullcase from vertical profiles.
analysis1_phases = ["pre-inception", "inception", "desinence"];
analysis2_phases = ["pre-inception", "inception", "desinence", "fullcase"];
analysis3_phases = ["pre-inception", "inception", "desinence", "fullcase"];
analysis4_phases = ["pre-inception", "inception", "desinence", "fullcase"];

% Optional case filters, using the original matPaths order. [] means all cases.
% Example: analysis1_caseIndices = [1 3 6] runs profiles for those cases only.
analysis1_caseIndices = [];
analysis2_caseIndices = [];
analysis3_caseIndices = [];
analysis4_caseIndices = [];

% Shared station locations for analyses 1 and 3.
profileOffset_mm = [1.5 3.0 4.5];
% Number of frames per case to save when analysis 2 is enabled.
instantaneousFramesToSave = 1;


%% Choose velocity files and case labels
% Add one MAT file per case. Use the labels in profile legends and map titles.

matPaths = { ...
    "E:\Sept 2026 Flowfield data\Processed data\P10S20\mat files\P10S20_3_475_40lpm_velocity.mat";
    "E:\Sept 2026 Flowfield data\Processed data\P10S30\mat files\P10S30_3_472_40lpm_velocity.mat";
    "E:\Sept 2026 Flowfield data\Processed data\P10S50\mat files\P10S50_3_456_40lpm_velocity.mat";
    "E:\Sept 2026 Flowfield data\Processed data\P10S70\mat files\P10S70_3_453_40lpm_velocity.mat";
    "E:\Sept 2026 Flowfield data\Processed data\P10S100\mat files\P10S100_3_440_40lpm_velocity.mat";
    "E:\Sept 2026 Flowfield data\Processed data\Smooth\mat files\Smooth_3_444_40lpm_velocity.mat";
    % 'E:\path\to\case6_velocity.mat';
};

% Keep labels in the same order as the MAT files.
caseLabels = { ...
    'Rough5';  % Change to the roughness label you want shown.
    'Rough4';  % Change to the roughness label you want shown.
    'Rough3';  % Change to the roughness label you want shown.
    'Rough2';  % Change to the roughness label you want shown.
    'Rough1';  % Change to the roughness label you want shown.
    'Smooth';  % Change to the roughness label you want shown.
    % 'Case 3 roughness';
    % 'Case 4 roughness';
    % 'Case 5 roughness';
    % 'Case 6 roughness';
};


%% Choose frame ranges for each phase
% Each phase has one [firstFrame lastFrame] row per original MAT file.
% Both endpoints are included; [2000 5000] processes flow frames 2000:5000.
% Pre-inception/desinence ranges below are dummy examples. Replace them with
% your measured ranges before selecting those phases.
% Unused phase/case ranges are not checked and may be left as [NaN NaN].
frameRangesByPhase.pre_inception = [ ...
    1, 1999;  % case 1
    1, 1999;  % case 2
    1, 1999;  % case 3
    1, 1999;  % case 4
    1, 1999;  % case 5
    1, 1999;  % case 6
];

frameRangesByPhase.inception = [ ...
    2000, 5000;  % case 1
    2000, 5000;  % case 2
    2000, 5000;  % case 3
    2000, 5000;  % case 4
    2000, 5000;  % case 5
    2000, 5000;  % case 6
];

frameRangesByPhase.desinence = [ ...
    5001, 6000;  % case 1
    5001, 6000;  % case 2
    5001, 6000;  % case 3
    5001, 6000;  % case 4
    5001, 6000;  % case 5
    5001, 6000;  % case 6
];

% fullcase is always [1 completedFlowFrames] for each case, set automatically.
% Do not specify a manual fullcase range.
frameRangesByPhase.fullcase = [];  % Populated from the file metadata below.


%% Set reading and output options
% Limit frames per disk read and set where plots and their MAT data are saved.

maxFramesPerRead = 100;  % Upper bound per disk read; lower this for less RAM.
% Use the local date when this run starts (for example, results/2026-10-01).
runDate = datestr(now,'yyyy-mm-dd');
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outputDir = fullfile(fileparts(mfilename('fullpath')), 'results', runDate, ...
    ['run_' runStamp]);
% Keep each run together so standalone plotting never mixes section data.
folderNumber = 2;
outputBase = outputDir;
while exist(outputDir,'dir') || exist(outputDir,'file')
    outputDir = sprintf('%s_%02d',outputBase,folderNumber);
    folderNumber = folderNumber + 1;
end
% Analyze older MAT files without a completion counter as complete.
% This assumes every allocated frame was written successfully; set false
% when analyzing files from interrupted runs whose completion is unknown.
assumeCompleteWhenNoCounter = true;


%% Check the inputs and prepare file information
% Confirm that the files, labels, and frame ranges are usable.
% Read file metadata needed for the analysis without loading all frames.

analysisSwitches = {analysis1_verticalProfiles, ...
    analysis2_instantaneousFrameSaving, ...
    analysis3_profileStationVisualization, analysis4_frameSummary};
assert(all(cellfun(@(value) islogical(value) && isscalar(value), ...
    analysisSwitches)), 'Analysis switches must be scalar true or false values.');
% Retain these names in run metadata for existing plotting workflows.
makeFrameSummary = analysis4_frameSummary;
makeMeanSpeedMaps = analysis3_profileStationVisualization;

assert(~isempty(matPaths), 'Add at least one MAT-file path.');
assert(numel(matPaths) <= 6, 'This plot defines distinct styles for up to six cases.');
assert(numel(caseLabels) == numel(matPaths), ...
    'Add one case label for each MAT-file path.');
assert(all(strlength(string(caseLabels)) > 0), ...
    'Case labels must not be empty.');
analysisPhases = normalizePhases(analysisPhases,'analysisPhases',false);
analysisPhaseFilters = {analysis1_phases,analysis2_phases, ...
    analysis3_phases,analysis4_phases};
analysisCaseIndices = {analysis1_caseIndices,analysis2_caseIndices, ...
    analysis3_caseIndices,analysis4_caseIndices};
analysisNames = {'vertical profiles','instantaneous frames', ...
    'profile station visualization','frame summary'};
phaseAnalysisEnabled = false(numel(analysisPhases),numel(analysisSwitches));
for analysisIndex = 1:numel(analysisSwitches)
    analysisPhaseFilters{analysisIndex} = normalizePhases( ...
        analysisPhaseFilters{analysisIndex}, ...
        sprintf('analysis%d_phases',analysisIndex),true);
    indices = analysisCaseIndices{analysisIndex};
    if isempty(indices)
        indices = 1:numel(matPaths);
    end
    assert(isnumeric(indices) && isreal(indices) && isvector(indices) && ...
        all(isfinite(indices)) && all(indices == fix(indices)) && ...
        all(indices >= 1 & indices <= numel(matPaths)) && ...
        numel(unique(indices)) == numel(indices), ...
        'analysis%d_caseIndices must contain distinct case numbers from 1:%d.', ...
        analysisIndex,numel(matPaths));
    analysisCaseIndices{analysisIndex} = indices(:)';
    phaseAnalysisEnabled(:,analysisIndex) = analysisSwitches{analysisIndex} & ...
        ismember(analysisPhases(:),analysisPhaseFilters{analysisIndex});
end
assert(any(phaseAnalysisEnabled(:)), ...
    'No analyses are enabled for analysisPhases. Check the switches and analysisN_phases.');
assert(isscalar(maxFramesPerRead) && isnumeric(maxFramesPerRead) && ...
    isfinite(maxFramesPerRead) && maxFramesPerRead >= 1 && ...
    maxFramesPerRead == fix(maxFramesPerRead), ...
    'maxFramesPerRead must be a positive integer.');
% Inspect only cases used by an enabled analysis in a requested phase.
activeCases = unique([analysisCaseIndices{any(phaseAnalysisEnabled,1)}]);
fileInfo = cell(numel(matPaths), 1);
for fileIndex = activeCases
    fileInfo{fileIndex} = inspectVelocityFile(matPaths{fileIndex}, ...
        assumeCompleteWhenNoCounter);
end
frameRangesByPhase.fullcase = nan(numel(matPaths),2);
for fileIndex = activeCases
    frameRangesByPhase.fullcase(fileIndex,:) = [1 fileInfo{fileIndex}.completedFrames];
end
for phaseIndex = 1:numel(analysisPhases)
    enabledAnalyses = phaseAnalysisEnabled(phaseIndex,:);
    if ~any(enabledAnalyses), continue; end
    phaseKey = char(strrep(analysisPhases(phaseIndex),'-','_'));
    assert(isfield(frameRangesByPhase,phaseKey), ...
        'Specify frameRangesByPhase.%s for the selected phase.',phaseKey);
    ranges = frameRangesByPhase.(phaseKey);
    assert(isnumeric(ranges) && isreal(ranges) && ismatrix(ranges) && ...
        size(ranges,2) == 2 && size(ranges,1) == numel(matPaths), ...
        'frameRangesByPhase.%s needs one [firstFrame lastFrame] row per MAT file.',phaseKey);
    phaseCases = unique([analysisCaseIndices{enabledAnalyses}]);
    for fileIndex = phaseCases
        firstFrame = ranges(fileIndex,1);
        lastFrame = ranges(fileIndex,2);
        info = fileInfo{fileIndex};
        assert(isfinite(firstFrame) && isfinite(lastFrame) && ...
            firstFrame == fix(firstFrame) && lastFrame == fix(lastFrame) && ...
            firstFrame >= 1 && lastFrame >= firstFrame && ...
            lastFrame <= info.completedFrames, ...
            'Phase %s, case %d range must lie within completed flow frames 1:%d.', ...
            char(analysisPhases(phaseIndex)),fileIndex,info.completedFrames);
    end
end

% Save provenance before computing outputs; plotting uses this run folder.
[folderCreated,folderMessage] = mkdir(outputDir);
assert(folderCreated,'Cannot create run folder: %s',folderMessage);
runOutputDir = outputDir;
save(fullfile(runOutputDir,'run_metadata.mat'), ...
    'runStamp','matPaths','caseLabels','frameRangesByPhase','analysisPhases', ...
    'analysisPhaseFilters','analysisCaseIndices','phaseAnalysisEnabled','analysisNames', ...
    'throatHeight_mm','throatWidth_mm','flowRate_Lpm', ...
    'makeFrameSummary','makeMeanSpeedMaps', ...
    'analysis1_verticalProfiles','analysis2_instantaneousFrameSaving', ...
    'analysis3_profileStationVisualization','analysis4_frameSummary', ...
    'profileOffset_mm','instantaneousFramesToSave');
fprintf('Analysis output: %s\n',outputDir);


%% Prepare shared profile stations
profileAnalyses = any(phaseAnalysisEnabled,1) & [true false true false];
profileInfo = cell(numel(matPaths),numel(profileOffset_mm));
if any(profileAnalyses)
assert(isscalar(throatHeight_mm) && isfinite(throatHeight_mm) && ...
    throatHeight_mm > 0, 'Throat height must be positive and finite in mm.');
assert(isnumeric(profileOffset_mm) && isvector(profileOffset_mm) && ...
    ~isempty(profileOffset_mm) && all(isfinite(profileOffset_mm)) && ...
    all(profileOffset_mm >= 0), ...
    'profileOffset_mm must be a nonempty vector of nonnegative distances in mm.');

profileCases = unique([analysisCaseIndices{profileAnalyses}]);
for i = profileCases
    for stationIndex = 1:numel(profileOffset_mm)
        profileInfo{i,stationIndex} = prepareProfileColumn( ...
            fileInfo{i},profileOffset_mm(stationIndex),matPaths{i});
    end
end

end

% Keep the original case order available while each analysis uses its subset.
allInputs = struct('matPaths',{matPaths},'caseLabels',{caseLabels}, ...
    'fileInfo',{fileInfo},'profileInfo',{profileInfo});

%% Run each requested phase, with independent analysis and case filters
for phaseIndex = 1:numel(analysisPhases)
analysisPhase = analysisPhases(phaseIndex);
enabledAnalyses = phaseAnalysisEnabled(phaseIndex,:);
if ~any(enabledAnalyses)
    fprintf('Skipping phase %s: no enabled analyses.\n',char(analysisPhase));
    continue
end
phaseTimer = logSectionStart(sprintf('Phase %s',char(analysisPhase)));
phaseFrameRanges = frameRangesByPhase.(char(strrep(analysisPhase,'-','_')));
phaseOutputDir = fullfile(runOutputDir,char(analysisPhase));
[folderCreated,folderMessage] = mkdir(phaseOutputDir);
assert(folderCreated,'Cannot create phase folder: %s',folderMessage);
% Phase metadata keeps all cases; each analysis MAT records its selected cases.
matPaths = allInputs.matPaths;
caseLabels = allInputs.caseLabels;
frameRanges = phaseFrameRanges;
save(fullfile(phaseOutputDir,'run_metadata.mat'), ...
    'runStamp','analysisPhase','matPaths','caseLabels','frameRanges', ...
    'enabledAnalyses','analysisNames','analysisCaseIndices');

%% Analysis 1. Vertical profiles of mean velocity and velocity fluctuations
if enabledAnalyses(1)
[matPaths,caseLabels,fileInfo,frameRanges,caseIndices,profileInfo] = ...
    selectAnalysisInputs(allInputs,analysisCaseIndices{1},phaseFrameRanges);
[outputDir,plotMatDir] = prepareAnalysisOutput(phaseOutputDir,analysisNames{1});
sectionTimer = logSectionStart('Analysis 1. Vertical profiles of mean velocity and velocity fluctuations');

% Choose distances downstream from the throat. At each station, average
% instantaneous U and V at each ROI height. V is positive downward in the
% MAT file, so reverse its sign for upward-positive physical coordinates.
% Each figure overlays the cases at one station; y = 0 is the lower ROI wall.

% Marker controls for the black and green profile styles.
markerOptions.blackCount = 48;       % Markers spread over the full profile.
markerOptions.blackSize = 6;
markerOptions.blackSymbol = 'o';
markerOptions.greenCount = 80;
markerOptions.greenSize = 5;
markerOptions.greenSymbol = 'p';

assert(isscalar(throatHeight_mm) && isfinite(throatHeight_mm) && throatHeight_mm > 0 && ...
    isscalar(throatWidth_mm) && isfinite(throatWidth_mm) && throatWidth_mm > 0, ...
    'Throat height and width must be positive finite values in mm.');
assert(isnumeric(flowRate_Lpm) && ...
    (isscalar(flowRate_Lpm) || numel(flowRate_Lpm) == numel(allInputs.matPaths)) && ...
    all(isfinite(flowRate_Lpm)) && all(flowRate_Lpm > 0), ...
    'flowRate_Lpm must be one shared positive value or one per MAT file.');
markerCounts = [markerOptions.blackCount markerOptions.greenCount];
assert(all(isfinite(markerCounts)) && all(markerCounts >= 0) && ...
    all(markerCounts == fix(markerCounts)), ...
    'Marker counts must be nonnegative integers.');
throatArea_m2 = throatHeight_mm * throatWidth_mm * 1e-6;
if isscalar(flowRate_Lpm)
    caseFlowRate_Lpm = repmat(flowRate_Lpm,numel(matPaths),1);
else
    caseFlowRate_Lpm = flowRate_Lpm(caseIndices);
    caseFlowRate_Lpm = caseFlowRate_Lpm(:);
end
bulkThroatSpeed_mps = (caseFlowRate_Lpm * 1e-3 / 60) / throatArea_m2;

% Rows are cases; columns are stations. Each statistic uses time samples
% where both U and V are finite, so covariance and component RMS agree.
verticalProfiles = cell(numel(matPaths),numel(profileOffset_mm));
maxProfileY = 0;
for stationIndex = 1:numel(profileOffset_mm)
    for i = 1:numel(matPaths)
        info = profileInfo{i,stationIndex};
        nRows = numel(info.profileImageRows);
        sampleCount = zeros(nRows,1);
        sumU = zeros(nRows,1);
        sumV = zeros(nRows,1);
        sumU2 = zeros(nRows,1);
        sumV2 = zeros(nRows,1);
        sumUV = zeros(nRows,1);
        fprintf('Profile for %s: flow frames %d:%d, x = %.4f mm (%.4f mm after throat).\n', ...
            char(string(caseLabels{i})),frameRanges(i,1),frameRanges(i,2), ...
            info.profileX_mm,info.actualOffset_mm);
        for k0 = frameRanges(i,1):maxFramesPerRead:frameRanges(i,2)
            k1 = min(k0+maxFramesPerRead-1,frameRanges(i,2));
            [uColumn,vColumn] = readUVProfileBlock(info,k0,k1);
            uColumn = double(uColumn);
            vColumn = -double(vColumn);
            validColumn = isfinite(uColumn) & isfinite(vColumn);
            uColumn(~validColumn) = 0;
            vColumn(~validColumn) = 0;
            sampleCount = sampleCount + sum(validColumn,2);
            sumU = sumU + sum(uColumn,2);
            sumV = sumV + sum(vColumn,2);
            sumU2 = sumU2 + sum(uColumn.^2,2);
            sumV2 = sumV2 + sum(vColumn.^2,2);
            sumUV = sumUV + sum(uColumn.*vColumn,2);
            clear uColumn vColumn validColumn
        end
        assert(any(sampleCount >= 2), ...
            'Case %d needs at least two valid U/V samples at station %.3f mm.', ...
            i,profileOffset_mm(stationIndex));
        hasSamples = sampleCount >= 2;
        meanU = nan(nRows,1);
        meanV = nan(nRows,1);
        varU = nan(nRows,1);
        varV = nan(nRows,1);
        covUV = nan(nRows,1);
        meanU(hasSamples) = sumU(hasSamples) ./ sampleCount(hasSamples);
        meanV(hasSamples) = sumV(hasSamples) ./ sampleCount(hasSamples);
        % Time-mean fluctuation moments (population convention, 1/N).
        varU(hasSamples) = max(0, sumU2(hasSamples) ./ ...
            sampleCount(hasSamples) - meanU(hasSamples).^2);
        varV(hasSamples) = max(0, sumV2(hasSamples) ./ ...
            sampleCount(hasSamples) - meanV(hasSamples).^2);
        covUV(hasSamples) = sumUV(hasSamples) ./ ...
            sampleCount(hasSamples) - meanU(hasSamples).*meanV(hasSamples);
        meanU = profileToPhysicalRows(meanU,info);
        sigmaU = profileToPhysicalRows(sqrt(varU),info);
        sigmaV = profileToPhysicalRows(sqrt(varV),info);
        reynoldsShear = profileToPhysicalRows(-covUV,info);
        % The reference uses sqrt(-<u'v'>)/U_b. This is undefined where
        % -<u'v'> is negative; retain the signed stress separately below.
        sqrtReynoldsShear = nan(size(reynoldsShear));
        nonnegativeShear = isfinite(reynoldsShear) & reynoldsShear >= 0;
        sqrtReynoldsShear(nonnegativeShear) = ...
            sqrt(reynoldsShear(nonnegativeShear));
        tkeInPlane = 0.5 * (sigmaU.^2 + sigmaV.^2);
        % Set the lower ROI boundary (the wall) to y = 0 at this station.
        profileY_mm = info.y_mm(:) - info.wallY_mm;
        verticalProfiles{i,stationIndex} = struct( ...
            'label',string(caseLabels{i}), ...
            'frameRange',frameRanges(i,:), ...
            'x_mm',info.profileX_mm, ...
            'offsetFromThroat_mm',info.actualOffset_mm, ...
            'y_mm',profileY_mm, ...
            'y_over_H',profileY_mm / throatHeight_mm, ...
            'sampleCount',profileToPhysicalRows(sampleCount,info), ...
            'meanU_mps',meanU, ...
            'bulkThroatSpeed_mps',bulkThroatSpeed_mps(i), ...
            'meanU_over_Ub',meanU / bulkThroatSpeed_mps(i), ...
            'sigmaU_mps',sigmaU, ...
            'sigmaV_mps',sigmaV, ...
            'sigmaU_over_Ub',sigmaU / bulkThroatSpeed_mps(i), ...
            'sigmaV_over_Ub',sigmaV / bulkThroatSpeed_mps(i), ...
            'reynoldsShear_m2ps2',reynoldsShear, ...
            'reynoldsShear_over_Ub2',reynoldsShear / bulkThroatSpeed_mps(i)^2, ...
            'sqrtReynoldsShear_over_Ub',sqrtReynoldsShear / bulkThroatSpeed_mps(i), ...
            'tkeInPlane_m2ps2',tkeInPlane, ...
            'tkeInPlane_over_Ub2',tkeInPlane / bulkThroatSpeed_mps(i)^2);
        validProfileRows = isfinite(meanU);
        maxProfileY = max(maxProfileY,max(profileY_mm(validProfileRows) / throatHeight_mm));
    end
end

verticalProfilesDir = outputDir;
save(fullfile(plotMatDir,'vertical_profile_data.mat'), ...
    'verticalProfiles','profileOffset_mm','matPaths','caseLabels', ...
    'frameRanges','analysisPhase','caseIndices','throatHeight_mm', ...
    'caseFlowRate_Lpm','bulkThroatSpeed_mps','markerOptions','-v7.3');

% Keep all five profile metrics for each station in its own plot folder.
stationDirs = cell(1,numel(profileOffset_mm));
for stationIndex = 1:numel(profileOffset_mm)
    stationDirs{stationIndex} = fullfile(verticalProfilesDir, ...
        sprintf('station%02d_xplus_%.2fmm',stationIndex,profileOffset_mm(stationIndex)));
    [folderCreated,folderMessage] = mkdir(stationDirs{stationIndex});
    assert(folderCreated,'Cannot create station folder: %s',folderMessage);
end

% Each metric uses one x scale across all stations and cases. Reference
% order: black circles, green stars, blue, red, orange, purple.
metricFields = {'meanU_over_Ub','sqrtReynoldsShear_over_Ub', ...
    'sigmaU_over_Ub','sigmaV_over_Ub','tkeInPlane_over_Ub2'};
metricNames = {'Mean streamwise U','Reynolds shear stress', ...
    'Streamwise velocity standard deviation', ...
    'Wall-normal velocity standard deviation', ...
    'In-plane turbulent kinetic energy'};
metricLabels = {'$\overline{u}/U_b$', ...
    '$\sqrt{-\overline{u''v''}}/U_b$', ...
    '$\sigma_u/U_b$','$\sigma_v/U_b$', ...
    '$k_{2D}/U_b^2$'};
metricFileTags = {'mean_streamwise_U','sqrt_reynolds_shear', ...
    'sigma_streamwise_U','sigma_wall_normal_V','tke_in_plane'};
for metricIndex = 1:numel(metricFields)
    metricMin = inf;
    metricMax = -inf;
    for profileIndex = 1:numel(verticalProfiles)
        values = verticalProfiles{profileIndex}.(metricFields{metricIndex});
        values = values(isfinite(values));
        if ~isempty(values)
            metricMin = min(metricMin,min(values));
            metricMax = max(metricMax,max(values));
        end
    end
    if ~isfinite(metricMin)
        warning('No finite values for %s; no figures saved.',metricNames{metricIndex});
        continue
    end
    xAxisLimits = [min(0,metricMin) max(0,metricMax)];
    for stationIndex = 1:numel(profileOffset_mm)
        fig = figure('Color','w','Position',[100 100 760 650]);
        ax = axes(fig);
        hold(ax,'on');
        plottedCases = 0;
        for i = 1:numel(matPaths)
            profile = verticalProfiles{i,stationIndex};
            values = profile.(metricFields{metricIndex});
            if ~any(isfinite(values) & isfinite(profile.y_over_H))
                warning('%s has no finite %s values at station %d.', ...
                    char(string(caseLabels{i})),metricNames{metricIndex},stationIndex);
                continue
            end
            plotProfileSeries(ax,profile,metricFields{metricIndex},caseIndices(i), ...
                char(string(caseLabels{i})),markerOptions);
            plottedCases = plottedCases + 1;
        end
        if plottedCases == 0
            close(fig);
            continue
        end
        styleProfileAxes(ax,sprintf('%s | throat + %.3f H | %s', ...
            metricNames{metricIndex}, ...
            profileOffset_mm(stationIndex) / throatHeight_mm, ...
            char(analysisPhase)),xAxisLimits,[0 maxProfileY], ...
            metricLabels{metricIndex});
        lgd = legend(ax,'show','Location','best', ...
            'Interpreter','none','Box','off');
        set(lgd,'FontName','Times New Roman','FontSize',11);
        plotName = sprintf('%s_profile_xplus_%.2fmm_%s_station%02d', ...
            metricFileTags{metricIndex},profileOffset_mm(stationIndex), ...
            char(analysisPhase),stationIndex);
        profilePlotBase = fullfile(stationDirs{stationIndex},plotName);
        profilePlotPath = [profilePlotBase '.png'];
        profileFigPath = [profilePlotBase '.fig'];
        exportgraphics(fig,profilePlotPath,'Resolution',600);
        savefig(fig,profileFigPath);
        close(fig);
        fprintf('Saved comparison plots: %s and %s\n', ...
            profilePlotPath,profileFigPath);
    end
end

logSectionEnd('Analysis 1. Vertical profiles of mean velocity and velocity fluctuations',sectionTimer);

end

%% Analysis 2. Instantaneous frame saving and MAT subset
if enabledAnalyses(2)
[matPaths,caseLabels,fileInfo,frameRanges,caseIndices] = ...
    selectAnalysisInputs(allInputs,analysisCaseIndices{2},phaseFrameRanges);
[outputDir,plotMatDir] = prepareAnalysisOutput(phaseOutputDir,analysisNames{2});
sectionTimer = logSectionStart('Analysis 2. Instantaneous frame saving and MAT subset');
% Save instantaneousFramesToSave frames from the start of each case's
% selected frame range, with a subfolder for each case's PNG frames.
% MAT subsets go in this phase/analysis's plot mat files folder.
% The MAT subset uses physical coordinates: x right, y up, V positive up.
streamlineDensity = 1.6;  % Denser than the mean plot to show small instantaneous turns.
negativeColorMin_mps = -1.2;  % U at or below this value uses the lowest color.

assert(isnumeric(instantaneousFramesToSave) && ...
    isscalar(instantaneousFramesToSave) && ...
    isfinite(instantaneousFramesToSave) && ...
    instantaneousFramesToSave >= 1 && ...
    instantaneousFramesToSave == fix(instantaneousFramesToSave), ...
    'instantaneousFramesToSave must be a positive integer when analysis 2 is enabled.');
assert(isnumeric(streamlineDensity) && isscalar(streamlineDensity) && ...
    isfinite(streamlineDensity) && streamlineDensity > 0, ...
    'streamlineDensity must be positive.');
assert(isnumeric(negativeColorMin_mps) && ...
    isscalar(negativeColorMin_mps) && ...
    isfinite(negativeColorMin_mps) && negativeColorMin_mps < 0, ...
    'negativeColorMin_mps must be negative.');

    selectedCounts = frameRanges(:,2) - frameRanges(:,1) + 1;
    assert(all(selectedCounts >= instantaneousFramesToSave), ...
        ['Each frame range must contain at least %d frames for the ' ...
         'instantaneous export.'], instantaneousFramesToSave);

    % Scan the frames being exported for one shared, linear U color scale.
    % This keeps positive velocities evenly spaced from zero to their max.
    minimumObservedU = inf;
    maximumObservedU = -inf;
    negativePixelCount = 0;
    finitePixelCount = 0;
    for i = 1:numel(fileInfo)
        info = fileInfo{i};
        lastExportFrame = frameRanges(i,1) + instantaneousFramesToSave - 1;
        for sourceFrame = frameRanges(i,1):lastExportFrame
            [uROI,~] = readVelocityBlock(info,sourceFrame,sourceFrame);
            finiteROI = uROI(isfinite(uROI));
            if ~isempty(finiteROI)
                minimumObservedU = min(minimumObservedU,min(double(finiteROI)));
                maximumObservedU = max(maximumObservedU,max(double(finiteROI)));
                negativePixelCount = negativePixelCount + nnz(finiteROI < 0);
                finitePixelCount = finitePixelCount + numel(finiteROI);
            end
        end
    end
    assert(finitePixelCount > 0, ...
        'No finite U was found for the instantaneous color scale.');
    assert(maximumObservedU > 0, ...
        'The exported frames need positive U for the 0-to-max color scale.');
    uColorLimits = [negativeColorMin_mps maximumObservedU];
    positiveMax = maximumObservedU;
    uColormap = interp1( ...
        [negativeColorMin_mps 0 0.25*positiveMax 0.5*positiveMax ...
         0.75*positiveMax positiveMax], ...
        [0.05 0.10 0.48; 0.18 0.64 0.87; 0.25 0.79 0.40; ...
         0.98 0.91 0.20; 0.97 0.46 0.08; 0.68 0.03 0.09], ...
        linspace(negativeColorMin_mps,positiveMax,256));
    uTickValues = [negativeColorMin_mps linspace(0,positiveMax,7)];
    fprintf('Instantaneous U color limits: %.3f to %.3f m/s\n',uColorLimits);
    fprintf('Observed minimum U: %.3f m/s; %d negative pixels among %d finite pixels.\n', ...
        minimumObservedU,negativePixelCount,finitePixelCount);

    instantaneousDir = outputDir;
    fprintf('Instantaneous output: %s\n',instantaneousDir);

    for i = 1:numel(matPaths)
        info = fileInfo{i};
        safeLabel = regexprep(char(string(caseLabels{i})), ...
            '[^A-Za-z0-9_-]','_');
        caseBase = sprintf('case%02d_%s',caseIndices(i),safeLabel);
        frameDir = fullfile(instantaneousDir,[caseBase '_frames']);
        mkdir(frameDir);
        matPath = fullfile(plotMatDir,[caseBase '.mat']);

        % MAT arrays are written one frame at a time to keep RAM bounded.
        x_mm = info.x_mm; %#ok<NASGU>
        y_mm = info.y_mm; %#ok<NASGU>
        fps_case = info.fps; %#ok<NASGU>
        sourceFrames = frameRanges(i,1) + (0:instantaneousFramesToSave-1); %#ok<NASGU>
        caseLabel = string(caseLabels{i}); %#ok<NASGU>
        sourceMatPath = string(matPaths{i}); %#ok<NASGU>
        sourceCaseIndex = caseIndices(i); %#ok<NASGU>
        save(matPath,'x_mm','y_mm','fps_case','sourceFrames', ...
            'caseLabel','sourceMatPath','sourceCaseIndex','analysisPhase','uColorLimits', ...
            'streamlineDensity','negativeColorMin_mps','-v7.3');
        subset = matfile(matPath,'Writable',true);
        [heightPx,widthPx] = size(info.maskROI);
        if instantaneousFramesToSave == 1
            subset.u_firstN(heightPx,widthPx) = single(0);
            subset.v_firstN(heightPx,widthPx) = single(0);
            subset.velMag_firstN(heightPx,widthPx) = single(0);
        else
            subset.u_firstN(heightPx,widthPx,instantaneousFramesToSave) = single(0);
            subset.v_firstN(heightPx,widthPx,instantaneousFramesToSave) = single(0);
            subset.velMag_firstN(heightPx,widthPx,instantaneousFramesToSave) = single(0);
        end

        fig = figure('Color','w','Units','pixels', ...
            'Position',[100 100 900 700], ...
            'PaperPositionMode','auto','Visible','off');
        ax = axes(fig);
        hold(ax,'on');
        maskPhysical = flipud(info.maskROI);
        [streamX,streamY] = meshgrid(info.x_mm,info.y_mm);
        streamlineHandles = gobjects(0);

        for j = 1:instantaneousFramesToSave
            sourceFrame = sourceFrames(j);
            [uROI,vROI] = readVelocityBlock(info,sourceFrame,sourceFrame);
            uImage = nan(heightPx,widthPx,'single');
            vImage = nan(heightPx,widthPx,'single');
            uImage(info.maskROI) = single(uROI);
            vImage(info.maskROI) = single(vROI);
            uPhysical = flipud(uImage);
            vPhysical = -flipud(vImage);
            speedPhysical = hypot(uPhysical,vPhysical);

            if instantaneousFramesToSave == 1
                subset.u_firstN(:,:) = uPhysical;
                subset.v_firstN(:,:) = vPhysical;
                subset.velMag_firstN(:,:) = speedPhysical;
            else
                subset.u_firstN(:,:,j) = uPhysical;
                subset.v_firstN(:,:,j) = vPhysical;
                subset.velMag_firstN(:,:,j) = speedPhysical;
            end

            if j == 1
                imageHandle = imagesc(ax,info.x_mm,info.y_mm,uPhysical);
                set(imageHandle,'AlphaData',isfinite(uPhysical));
                set(ax,'YDir','normal');
                axis(ax,'image');
                colormap(ax,uColormap);
                caxis(ax,uColorLimits);
                colorbarHandle = colorbar(ax);
                colorbarHandle.Label.String = 'Streamwise U (m/s)';
                colorbarHandle.Ticks = uTickValues;
                colorbarHandle.TickLabels = arrayfun( ...
                    @(v) sprintf('%.3g',v),uTickValues,'UniformOutput',false);
                box(ax,'on');
                set(ax,'FontName','Times New Roman','FontSize',12, ...
                    'LineWidth',1,'XMinorTick','on','YMinorTick','on', ...
                    'TickDir','out');
                xlabel(ax,'x (mm)');
                ylabel(ax,'y (mm)');
            else
                set(imageHandle,'CData',uPhysical, ...
                    'AlphaData',isfinite(uPhysical));
            end
            if ~isempty(streamlineHandles)
                delete(streamlineHandles(isgraphics(streamlineHandles)));
            end
            streamlineHandles = drawArrowedStreamlines(ax,streamX,streamY, ...
                double(uPhysical),double(vPhysical),maskPhysical, ...
                streamlineDensity);
            title(ax,sprintf('%s: Instantaneous streamwise U (flow frame %d)', ...
                char(string(caseLabels{i})),sourceFrame),'Interpreter','none');
            pngPath = fullfile(frameDir,sprintf('frame_%06d.png',sourceFrame));
            exportgraphics(fig,pngPath,'Resolution',300);
            if mod(j,50) == 0 || j == instantaneousFramesToSave
                fprintf('  %s: saved %d/%d frames\n',caseBase,j,instantaneousFramesToSave);
            end
        end
        close(fig);
        fprintf('Saved %s and PNG frames in %s\n',matPath,frameDir);
    end

logSectionEnd('Analysis 2. Instantaneous frame saving and MAT subset',sectionTimer);
end

%% Analysis 3. Profile station visualization
if enabledAnalyses(3)
[matPaths,caseLabels,fileInfo,frameRanges,caseIndices,profileInfo] = ...
    selectAnalysisInputs(allInputs,analysisCaseIndices{3},phaseFrameRanges);
[outputDir,plotMatDir] = prepareAnalysisOutput(phaseOutputDir,analysisNames{3});
sectionTimer = logSectionStart('Analysis 3. Profile station visualization');
% Average speed over the selected frames at every valid ROI pixel.
% Save one map per case with vertical lines showing the analysis 1 stations.
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

% Save calibrated coordinates and station positions without matfile handles.
mapCoordinates = cell(numel(fileInfo),1);
for i = 1:numel(fileInfo)
    info = fileInfo{i};
    actualOffsets_mm = cellfun(@(p) p.actualOffset_mm,profileInfo(i,:));
    mapCoordinates{i} = struct('x_mm',info.x_mm,'y_mm',info.y_mm, ...
        'xThroat_mm',info.xThroat_mm,'actualOffsets_mm',actualOffsets_mm);
end
save(fullfile(plotMatDir,'mean_speed_map_data.mat'), ...
    'meanSpeedMaps','mapCoordinates','caseLabels','frameRanges', ...
    'analysisPhase','caseIndices','throatHeight_mm','profileOffset_mm','matPaths','-v7.3');

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
    imageHandle = imagesc(ax,(info.x_mm-info.xThroat_mm)/throatHeight_mm, ...
        info.y_mm/throatHeight_mm,mapPhysical);
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
        stationX = profileInfo{i,stationIndex}.actualOffset_mm/throatHeight_mm;
        stationColor = stationColors(mod(stationIndex-1,size(stationColors,1))+1,:);
        plot(ax,[stationX stationX],info.y_mm([1 end])/throatHeight_mm,'-', ...
            'Color',[0.08 0.08 0.08],'LineWidth',3.4, ...
            'HandleVisibility','off');
        plot(ax,[stationX stationX],info.y_mm([1 end])/throatHeight_mm,'--', ...
            'Color',stationColor,'LineWidth',2, ...
            'DisplayName',sprintf('Station %d: +%.3f H', ...
                stationIndex,profileOffset_mm(stationIndex)/throatHeight_mm));
    end
    xlabel(ax,'(x - x_{throat})/H','FontName','Times New Roman','FontSize',14);
    ylabel(ax,'y/H','FontName','Times New Roman','FontSize',14);
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
        caseIndices(i),safeLabel,char(analysisPhase),frameRanges(i,1),frameRanges(i,2)));
    mapPngPath = [mapPlotBase '.png'];
    mapFigPath = [mapPlotBase '.fig'];
    exportgraphics(fig,mapPngPath,'Resolution',600);
    savefig(fig,mapFigPath);
    fprintf('Saved mean speed map: %s and %s\n',mapPngPath,mapFigPath);
    close(fig);
end

logSectionEnd('Analysis 3. Profile station visualization',sectionTimer);
end

%% Analysis 4. Frame summary across the ROI
if enabledAnalyses(4)
[matPaths,caseLabels,fileInfo,frameRanges,caseIndices] = ...
    selectAnalysisInputs(allInputs,analysisCaseIndices{4},phaseFrameRanges);
[outputDir,plotMatDir] = prepareAnalysisOutput(phaseOutputDir,analysisNames{4});
sectionTimer = logSectionStart('Analysis 4. Frame summary across the ROI');
% Calculate mean U, V, and speed for every
% selected frame using valid pixels inside the ROI. Print the first 10 rows.
% U points right and V points down in the source image.
frameSummary = table();
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
        caseIndices(i), char(string(caseLabels{i})), char(analysisPhase), char(string(matPaths{i})), ...
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

            caseIndex(outRow) = caseIndices(i);
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
save(fullfile(plotMatDir,'frame_summary_data.mat'), ...
    'frameSummary','caseLabels','caseIndices','frameRanges','analysisPhase','matPaths','-v7.3');

logSectionEnd('Analysis 4. Frame summary across the ROI',sectionTimer);
end

logSectionEnd(sprintf('Phase %s',char(analysisPhase)),phaseTimer);
end
fprintf('Completed analysis run: %s\n',runOutputDir);

%% Notes for future calculations

% This section only shows how to rebuild an image from ROI data; it runs
% no analysis. U/V rows follow find(info.maskROI), with one frame per column.
% To rebuild one image-coordinate frame:
%   uImage = nan(size(info.maskROI), 'single');
%   uImage(info.maskROI) = U(:,1);
% To get physical coordinates (origin at lower left):
%   uPhysical = flipud(uImage);  vPhysical = -flipud(vImage);


%% Helper functions used by the sections above

function phases = normalizePhases(values,parameterName,allowEmpty)
    phases = lower(strtrim(string(values)));
    if isempty(phases) && allowEmpty
        phases = strings(1,0);
        return
    end
    validPhases = ["pre-inception", "inception", "desinence", "fullcase"];
    assert(isvector(phases) && ~isempty(phases) && ...
        all(ismember(phases,validPhases)) && ...
        numel(unique(phases)) == numel(phases), ...
        '%s must list distinct phases: pre-inception, inception, desinence, fullcase.', ...
        parameterName);
    phases = phases(:)';
end

function [paths,labels,info,ranges,caseIndices,profiles] = ...
        selectAnalysisInputs(allInputs,caseIndices,phaseRanges)
    paths = allInputs.matPaths(caseIndices);
    labels = allInputs.caseLabels(caseIndices);
    info = allInputs.fileInfo(caseIndices);
    ranges = phaseRanges(caseIndices,:);
    profiles = allInputs.profileInfo(caseIndices,:);
end

function [folder,plotMatDir] = prepareAnalysisOutput(phaseFolder,analysisName)
    folder = fullfile(phaseFolder,analysisName);
    plotMatDir = fullfile(folder,'plot mat files');
    [ok,message] = mkdir(plotMatDir);
    assert(ok,'Cannot create analysis output folder: %s',message);
end

function sectionTimer = logSectionStart(sectionName)
    % Print local wall-clock time; use tic/toc for elapsed duration.
    fprintf('[%s] START | %s\n', ...
        datestr(now,'yyyy-mm-dd HH:MM:SS'),sectionName);
    sectionTimer = tic;
end

function logSectionEnd(sectionName,sectionTimer)
    fprintf('[%s] END   | %s | Elapsed: %.3f s\n', ...
        datestr(now,'yyyy-mm-dd HH:MM:SS'),sectionName,toc(sectionTimer));
end

function handles = drawArrowedStreamlines(ax,X,Y,U,V,roiMask,density)
    % Use the same automatic streamline and arrow placement as the ARC plot.
    % Extra density lets short paths show local rotations and fluctuations.
    U(~roiMask) = NaN;
    V(~roiMask) = NaN;
    handles = streamslice(ax,X,Y,U,V,density);
    set(handles,'Color','k','LineWidth',0.8);
    for k = 1:numel(handles)
        xData = handles(k).XData;
        yData = handles(k).YData;
        if isempty(xData) || isempty(yData), continue; end
        insideROI = interp2(X,Y,single(roiMask),xData,yData, ...
            'nearest',0) >= 0.5;
        xData(~insideROI) = NaN;
        yData(~insideROI) = NaN;
        handles(k).XData = xData;
        handles(k).YData = yData;
    end
end

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
    % The local runner saves NaN when throat selection is disabled. Frame
    % summaries and absolute-coordinate exports do not need a throat.
    xThroat_mm = NaN;
    if ismember('x_throat_mm',vars)
        xThroat_mm = double(M.x_throat_mm);
    end
    if isscalar(xThroat_mm) && isnan(xThroat_mm) && ...
            ismember('x_throat_pixel',vars)
        xThroat_mm = double(M.x_throat_pixel) * mmPerPixel;
    end
    assert(isscalar(xThroat_mm) && isreal(xThroat_mm) && ...
        (isfinite(xThroat_mm) || isnan(xThroat_mm)), ...
        'x_throat_mm must be a finite scalar or NaN: %s', char(string(matPath)));

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
    elseif ~packed && ismember('numFrames',vars) && ...
            ismember('numFlowFrames',vars) && ...
            contains(char(string(matPath)), '_velocity_local_')
        % The local runner writes its full-frame arrays in a single save
        % after processing finishes. Older local files lack a counter.
        nVideoFrames = double(M.numFrames);
        nFlowFrames = double(M.numFlowFrames);
        assert(isscalar(nVideoFrames) && isfinite(nVideoFrames) && ...
            nVideoFrames == fix(nVideoFrames) && nVideoFrames >= 2 && ...
            isscalar(nFlowFrames) && isfinite(nFlowFrames) && ...
            nFlowFrames == allocatedFrames && nFlowFrames == nVideoFrames-1, ...
            'Invalid local frame counts in %s', char(string(matPath)));
        completedFrames = nFlowFrames;
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
    assert(isfinite(info.xThroat_mm), ...
        ['Throat-relative profiles require a throat location in %s. ' ...
         'Rerun RAFT with useThroat = true or supply x_throat_mm.'], ...
        char(string(matPath)));
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
    % Image rows increase downward; y_mm increases upward from the frame
    % bottom. The last ROI row in this column is the lower wall.
    wallPhysicalRow = size(info.maskROI,1) - imageRows(end) + 1;
    info.wallY_mm = info.y_mm(wallPhysicalRow);
end

function [U,V] = readUVProfileBlock(info,k0,k1)
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

function valuesPhysical = profileToPhysicalRows(valuesROI,info)
    valuesImage = nan(size(info.maskROI,1),1);
    valuesImage(info.profileImageRows) = valuesROI;
    valuesPhysical = flipud(valuesImage);
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

function plotProfileSeries(ax,profile,metricField,styleIndex,seriesLabel,markerOptions)
    % Swap the styles assigned to Rough5 (case 1) and Smooth (case 6).
    if styleIndex == 1
        styleIndex = 6;
    elseif styleIndex == 6
        styleIndex = 1;
    end
    values = profile.(metricField);
    validRows = find(isfinite(values) & isfinite(profile.y_over_H));
    assert(~isempty(validRows), 'Cannot plot a profile without finite values.');
    markerRows = [];  % Black, green, and orange series use markers.
    markerCount = 0;
    markerFaceColor = 'none';
    switch styleIndex
        case 1  % open black circles
            lineColor = [0 0 0];
            lineStyle = 'none';
            markerSymbol = markerOptions.blackSymbol;
            markerCount = markerOptions.blackCount;
            markerFaceColor = 'none';
            lineWidth = 1.3;
            markerSize = markerOptions.blackSize;
        case 2  % filled five-point green stars
            lineColor = [0 0.55 0];
            lineStyle = 'none';
            markerSymbol = markerOptions.greenSymbol;
            markerCount = markerOptions.greenCount;
            markerFaceColor = lineColor;
            lineWidth = 1.4;
            markerSize = markerOptions.greenSize;
        case 3  % solid blue curve
            lineColor = [0 0 1];
            lineStyle = '-';
            markerSymbol = 'none';
            lineWidth = 2.5;
            markerSize = 6;
        case 4  % dashed red curve
            lineColor = [1 0 0];
            lineStyle = '--';
            markerSymbol = 'none';
            lineWidth = 2;
            markerSize = 6;
        case 5  % open orange triangles (Rough1)
            lineColor = [0.90 0.45 0];
            lineStyle = 'none';
            markerSymbol = '^';
            markerCount = markerOptions.blackCount;
            markerFaceColor = 'none';
            lineWidth = 2;
            markerSize = markerOptions.greenSize;
        case 6  % dotted purple curve
            lineColor = [0.55 0.20 0.70];
            lineStyle = ':';
            markerSymbol = 'none';
            lineWidth = 2.5;
            markerSize = 6;
    end
    if markerCount > 0
        % Spread markers evenly over the valid profile rows.
        nMarkers = min(markerCount,numel(validRows));
        markerRows = unique(validRows(round(linspace(1,numel(validRows),nMarkers))));
    end
    plot(ax,values,profile.y_over_H, ...
        'Color',lineColor,'LineStyle',lineStyle,'LineWidth',lineWidth, ...
        'Marker',markerSymbol,'MarkerIndices',markerRows, ...
        'MarkerSize',markerSize,'MarkerFaceColor',markerFaceColor, ...
        'MarkerEdgeColor',lineColor, ...
        'DisplayName',seriesLabel);
end

function styleProfileAxes(ax,plotTitle,xAxisLimits,yLimits,xAxisLabel)
    fontName = 'Times New Roman';
    set(ax,'FontName',fontName,'FontSize',12,'LineWidth',1, ...
        'TickDir','in','Box','on','Color','w', ...
        'XColor',[0 0 0],'YColor',[0 0 0], ...
        'XGrid','on','YGrid','on','GridColor',[0.67 0.67 0.67], ...
        'GridAlpha',0.6,'Layer','top');
    % Scale padding to each metric; keep zero visible for positive profiles.
    if diff(xAxisLimits) > 0
        xPadding = 0.03*diff(xAxisLimits);
    else
        xPadding = max(0.05*abs(xAxisLimits(1)),0.01);
    end
    xLower = xAxisLimits(1) - xPadding;
    if xAxisLimits(1) == 0, xLower = 0; end
    xlim(ax,[xLower xAxisLimits(2)+xPadding]);
    yPadding = max(0.03*(yLimits(2)-yLimits(1)),0.01);
    ylim(ax,[0 yLimits(2)+yPadding]);
    xlabel(ax,xAxisLabel,'Interpreter','latex', ...
        'FontName',fontName,'FontSize',14,'FontAngle','italic');
    ylabel(ax,'y/H','FontName',fontName,'FontSize',14,'FontAngle','italic');
    title(ax,plotTitle,'FontName',fontName,'FontSize',15, ...
        'FontWeight','normal','Interpreter','none');
end
