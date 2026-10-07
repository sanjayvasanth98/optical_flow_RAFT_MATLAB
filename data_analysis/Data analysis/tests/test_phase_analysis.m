function test_phase_analysis()
% Integration checks with tiny v7.3 fixtures; never reads experiment files.
% Run: addpath('tests'); test_phase_analysis
analysisRoot = fileparts(fileparts(mfilename('fullpath')));
plotRoot = fullfile(analysisRoot,'plotting codes');
addpath(plotRoot);
taskDir = tempname;
mkdir(taskDir);
oldVisibility = get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
cleanup = onCleanup(@() finishChecks(taskDir,oldVisibility)); %#ok<NASGU>

% Case 2 deliberately does not exist: every analysis excludes it.
paths = {fullfile(taskDir,'full.mat'),fullfile(taskDir,'missing.mat'), ...
    fullfile(taskDir,'packed.mat')};
maskROI = true(4,6);
fps = 1000;
mm_per_pixel = 1;
x_throat_mm = 0;
[row,column] = ndgrid(1:4,1:6);
u_all = zeros(4,6,8,'single');
v_all = u_all;
for k = 1:8
    u_all(:,:,k) = k + row/10 + column/100;
    v_all(:,:,k) = 0.2*k;
end
completedFlowFrames = 6;
u_all(:,:,7:8) = 10000; % Allocated but uncommitted frames must be excluded.
save(paths{1},'u_all','v_all','maskROI','fps','mm_per_pixel', ...
    'x_throat_mm','completedFlowFrames','-v7.3');
u_all = zeros(nnz(maskROI),9,'single');
v_all = u_all;
for k = 1:9
    image = single(10 + k + row/10 + column/100);
    u_all(:,k) = image(:);
    v_all(:,k) = 0.3*k;
end
completedFlowFrames = 7;
u_all(:,8:9) = 10000;
instantaneousStorage = 'roi_pixels_by_frame_v1';
save(paths{3},'u_all','v_all','maskROI','fps','mm_per_pixel', ...
    'x_throat_mm','completedFlowFrames','instantaneousStorage','-v7.3');

sourceText = fileread(fullfile(analysisRoot,'RAFT_chunk_analysis.m'));
body = sourceText(strfind(sourceText,'%% Set reading and output options'):end);
body = strrep(body,'maxFramesPerRead = 100;','maxFramesPerRead = 2;');
body = strrep(body,'''Resolution'',600','''Resolution'',72');
body = strrep(body,'''Resolution'',300','''Resolution'',72');
config = strjoin({ ...
    'clearvars;', ...
    'throatHeight_mm = 10; throatWidth_mm = 5; flowRate_Lpm = [10 20 30];', ...
    'analysis1_verticalProfiles = true; analysis2_instantaneousFrameSaving = true;', ...
    'analysis3_profileStationVisualization = true; analysis4_frameSummary = true;', ...
    'analysis1_phases = ["pre-inception", "inception", "desinence"];', ...
    'analysis2_phases = "fullcase";', ...
    'analysis3_phases = ["desinence", "fullcase"];', ...
    'analysis4_phases = ["pre-inception", "inception", "desinence", "fullcase"];', ...
    'analysis1_caseIndices = [3 1]; analysis2_caseIndices = 1;', ...
    'analysis3_caseIndices = 3; analysis4_caseIndices = [1 3];', ...
    'profileOffset_mm = 1; instantaneousFramesToSave = 1;', ...
    sprintf('matPaths = {''%s''; ''%s''; ''%s''};',paths{:}), ...
    'caseLabels = {''Full''; ''Unused''; ''Packed''};', ...
    'analysisPhases = ["pre-inception", "inception", "desinence", "fullcase"];', ...
    'frameRangesByPhase.pre_inception = [1 2; NaN NaN; 2 3];', ...
    'frameRangesByPhase.inception = [3 4; NaN NaN; 4 5];', ...
    'frameRangesByPhase.desinence = [5 6; NaN NaN; 6 7];'},newline);
scriptPath = fullfile(taskDir,'fixture_analysis.m');
runDir = executeFixture(scriptPath,config,body);

phases = ["pre-inception", "inception", "desinence", "fullcase"];
expectedRanges = {[1 2;2 3],[3 4;4 5],[5 6;6 7],[1 6;1 7]};
for phase = 1:numel(phases)
    phaseDir = fullfile(runDir,char(phases(phase)));
    summaryPath = fullfile(phaseDir,'frame summary','plot mat files','frame_summary_data.mat');
    data = load(summaryPath);
    assert(isequal(data.caseIndices,[1 3]));
    assert(isequal(data.frameRanges,expectedRanges{phase}));
    assert(all(data.frameSummary.phase == phases(phase)));
    assert(height(data.frameSummary) == sum(diff(expectedRanges{phase},1,2)+1));
    for index = [1 3]
        rows = data.frameSummary.caseIndex == index;
        expectedU = double(data.frameSummary.flowFrame(rows)) + 0.285 + 10*(index==3);
        assert(max(abs(data.frameSummary.meanU_mps(rows)-expectedU)) < 1e-5);
    end
    profilePath = fullfile(phaseDir,'vertical profiles','plot mat files','vertical_profile_data.mat');
    if phases(phase) == "fullcase"
        assert(~isfolder(fullfile(phaseDir,'vertical profiles')));
    else
        profile = load(profilePath);
        assert(isequal(profile.caseIndices,[3 1]));
        assert(isequal(profile.caseFlowRate_Lpm,[30;10]));
        assert(isequal(profile.frameRanges,expectedRanges{phase}([2 1],:)));
        expectedMean = mean(expectedRanges{phase}(2,:)) + 10 + flipud(row(:,2)/10) + 0.02;
        assert(max(abs(profile.verticalProfiles{1}.meanU_mps-expectedMean)) < 1e-5);
    end
end
instantPath = fullfile(runDir,'fullcase','instantaneous frames','plot mat files','case01_Full.mat');
instant = load(instantPath);
assert(instant.analysisPhase == "fullcase" && instant.sourceCaseIndex == 1);
assert(isequal(instant.sourceFrames,1));
assert(max(abs(double(instant.u_firstN(:))-reshape(flipud(single(1+row/10+column/100)),[],1))) < 1e-5);
assert(~isfolder(fullfile(runDir,'inception','instantaneous frames')));
maps = load(fullfile(runDir,'fullcase','profile station visualization', ...
    'plot mat files','mean_speed_map_data.mat'));
assert(isequal(maps.caseIndices,3) && isequal(maps.frameRanges,[1 7]));
expectedSpeed = zeros(4,6);
for k = 1:7
    expectedSpeed = expectedSpeed + hypot(10+k+row/10+column/100,0.3*k)/7;
end
assert(max(abs(double(maps.meanSpeedMaps{1}(:))-reshape(flipud(expectedSpeed),[],1))) < 1e-5);

% Replot routing must select one phase, preserve case identity, and support old data.
options = raft_plot_options();
options.runDir = runDir;
mustFail(@() raftplot.source(options,'summary'),'Multiple phases');
options.analysisPhase = 'inception';
[files,selected] = raftplot.source(options,'summary');
assert(contains(selected,fullfile('inception','frame summary')));
assert(numel(files) == 1);
options.visible = 'off'; options.closeAfterSave = true;
options.savePNG = false; options.saveFIG = false;
raftplot.summary(options);
raftplot.profiles(options,1);
options.analysisPhase = 'fullcase';
mustFail(@() raftplot.source(options,'profiles'),'No profiles MAT data');
raftplot.maps(options);
options.analysisPhase = '';
[files,~] = raftplot.source(options,'instantaneous');
assert(numel(files) == 1 && strcmp(files{1},instantPath));
raftplot.instantaneous(options);
options.runDir = fullfile(runDir,'inception');
[files,~] = raftplot.source(options,'profiles');
assert(contains(files{1},fullfile('inception','vertical profiles')));
options.runDir = fileparts(fileparts(files{1}));
raftplot.source(options,'profiles');
options.analysisPhase = 'fullcase';
mustFail(@() raftplot.source(options,'profiles'),'conflicts');
legacyDir = fullfile(taskDir,'legacy','plot mat files'); mkdir(legacyDir);
copyfile(summaryPath,fullfile(legacyDir,'frame_summary_data.mat'));
options.runDir = fileparts(legacyDir);
raftplot.source(options,'summary');
options.analysisPhase = 'inception';
mustFail(@() raftplot.source(options,'summary'),'Saved phase does not match');
options.dataFile = instantPath;
raftplot.source(options,'instantaneous');

% A single phase ignores all other ranges and never prepares excluded profiles.
singleConfig = strrep(config, ...
    'analysisPhases = ["pre-inception", "inception", "desinence", "fullcase"];', ...
    'analysisPhases = "inception";');
singleConfig = strrep(singleConfig, ...
    'frameRangesByPhase.desinence = [5 6; NaN NaN; 6 7];', ...
    'frameRangesByPhase.desinence = NaN;');
singleDir = executeFixture(scriptPath,singleConfig,body);
assert(isfolder(fullfile(singleDir,'inception')));
assert(~isfolder(fullfile(singleDir,'fullcase')) && ~isfolder(fullfile(singleDir,'desinence')));
fullConfig = strrep(singleConfig,'analysisPhases = "inception";', ...
    'analysisPhases = ["fullcase", "desinence"];');
fullConfig = strrep(fullConfig,'analysis1_verticalProfiles = true;', ...
    'analysis1_verticalProfiles = false;');
fullConfig = strrep(fullConfig,'analysis2_instantaneousFrameSaving = true;', ...
    'analysis2_instantaneousFrameSaving = false;');
fullConfig = strrep(fullConfig,'analysis3_profileStationVisualization = true;', ...
    'analysis3_profileStationVisualization = false;');
fullConfig = strrep(fullConfig, ...
    'analysis4_phases = ["pre-inception", "inception", "desinence", "fullcase"];', ...
    'analysis4_phases = "fullcase";');
fullDir = executeFixture(scriptPath,fullConfig,body);
assert(isfolder(fullfile(fullDir,'fullcase')) && ~isfolder(fullfile(fullDir,'desinence')));

badConfig = strrep(singleConfig,'analysisPhases = "inception";', ...
    'analysisPhases = ["inception", "inception"];');
mustFail(@() executeFixture(scriptPath,badConfig,body),'distinct phases');
badConfig = strrep(singleConfig,'[3 4; NaN NaN; 4 5]','[3 7; NaN NaN; 4 5]');
mustFail(@() executeFixture(scriptPath,badConfig,body),'completed flow frames');
badConfig = strrep(singleConfig,'analysis1_caseIndices = [3 1];', ...
    'analysis1_caseIndices = [3 3];');
mustFail(@() executeFixture(scriptPath,badConfig,body),'distinct case numbers');
fprintf('PASS: phase routing, case filters, completed-frame limits, saved data, replots, and validation.\n');
end

function runDir = executeFixture(scriptPath,config,body)
fid = fopen(scriptPath,'w');
assert(fid >= 0);
fwrite(fid,sprintf('%s\n%s',config,body),'char');
fclose(fid);
evalin('base',sprintf('run(''%s'');',strrep(scriptPath,'''','''''')));
runDir = evalin('base','runOutputDir');
end

function mustFail(action,messagePart)
try
    action();
catch exception
    assert(contains(exception.message,messagePart), ...
        'Unexpected error: %s',exception.message);
    return
end
error('Expected an error containing: %s',messagePart);
end

function finishChecks(taskDir,oldVisibility)
set(groot,'defaultFigureVisible',oldVisibility);
evalin('base','clearvars'); % Release fixture matfile handles before removal.
if isfolder(taskDir), rmdir(taskDir,'s'); end
end
