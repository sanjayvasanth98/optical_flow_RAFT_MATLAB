%% Section 9: instantaneous streamwise U with arrowed streamlines.
% Run independently; no variables from RAFT_chunk_analysis are needed.
addpath(fileparts(mfilename('fullpath')));
options = raft_plot_options();
options.showGrid = false;
options.figurePosition = [100 100 900 700];
options.resolution = 300;

%% Refinement controls (uncomment or edit as needed)
% options.runDir = fullfile(options.resultsRoot,'2026-10-07');
% options.analysisPhase = 'inception';
% options.dataFile = 'D:\path\to\saved_section_data.mat';
% options.xLimits = [0 1];
% options.yLimits = [0 0.5];
% options.caseIndices = [1 3 6];
% options.stationIndices = [1 2];
% options.useSavedMarkerOptions = false;
% options.markerOptions.blackCount = 30;
% options.savePDF = true;
% options.closeAfterSave = true;

%% Recreate plots from saved analysis data
raftplot.instantaneous(options);
