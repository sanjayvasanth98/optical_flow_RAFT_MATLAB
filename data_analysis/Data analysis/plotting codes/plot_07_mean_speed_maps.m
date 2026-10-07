%% Section 7: mean speed maps with profile station lines.
% Run independently; no variables from RAFT_chunk_analysis are needed.
addpath(fileparts(mfilename('fullpath')));
options = raft_plot_options();
options.showGrid = false;
options.legendLocation = 'southoutside';
options.figurePosition = [100 100 1200 800];

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
raftplot.maps(options);
