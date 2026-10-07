function options = raft_plot_options()
% Shared refinement controls. Individual plotting scripts can override these.
root = fileparts(mfilename('fullpath'));
options.resultsRoot = fullfile(fileparts(root),'results');
options.runDir = ''; % Empty: latest run. Or set an older run/date folder.
options.analysisPhase = ''; % One phase; empty auto-selects only if unambiguous.
options.dataFile = ''; % Optional exact MAT path (instantaneous: one case).
options.outputDir = ''; % Empty: unique replots/<section>_<timestamp> folder.
options.visible = 'on';
options.closeAfterSave = false; % Keep figures open for interactive refinement.
options.savePNG = true;
options.saveFIG = true;
options.savePDF = false;
options.resolution = 600;
options.fontName = 'Times New Roman';
options.fontSize = 12;
options.labelFontSize = 14;
options.titleFontSize = 15;
options.figurePosition = [100 100 760 650];
options.legendLocation = 'best';
options.showGrid = true;
options.xLimits = []; % Empty: automatic. Profiles use one scale per metric.
options.yLimits = [];
options.colorLimits = []; % Empty: shared scale across all selected cases/frames.
options.caseIndices = []; % Empty: all saved cases.
options.stationIndices = []; % Empty: all saved profile stations.
options.markerOptions = struct('blackCount',48,'blackSize',6, ...
    'blackSymbol','o','greenCount',80,'greenSize',5,'greenSymbol','p');
options.useSavedMarkerOptions = true;
options.streamlineDensity = 1.6;
options.showStreamlines = true;
options.negativeColorMin_mps = -1.2;
options.frameIndices = []; % Instantaneous subset positions, not source frame IDs.
options.useSavedColorLimits = true;
end
