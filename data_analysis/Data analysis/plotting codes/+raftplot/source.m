function [files,runDir] = source(options,kind)
% Resolve one coherent run; never fall back to older section data silently.
if ~isempty(options.dataFile)
    assert(isfile(options.dataFile),'MAT file not found: %s',options.dataFile);
    files = {char(options.dataFile)};
    runDir = fileparts(files{1});
    [~,folderName] = fileparts(runDir);
    if strcmp(folderName,'plot mat files'), runDir = fileparts(runDir); end
    fprintf('Plot source: %s\n',files{1});
    return
end
runDir = char(options.runDir);
if isempty(runDir)
    dates = dir(options.resultsRoot);
    dates = dates([dates.isdir]);
    dates = dates(~cellfun('isempty',regexp({dates.name}, ...
        '^\d{4}-\d{2}-\d{2}$','once')));
    [~,order] = sort(string({dates.name}),'descend');
    for index = order
        dateDir = fullfile(dates(index).folder,dates(index).name);
        runs = dir(fullfile(dateDir,'run_*'));
        runs = runs([runs.isdir]);
        if ~isempty(runs)
            [~,newest] = sort(string({runs.name}),'descend');
            runDir = fullfile(runs(newest(1)).folder,runs(newest(1)).name);
            break
        end
        % Legacy date folders may include empty directories from failed runs.
        data = dir(fullfile(dateDir,'**','*.mat'));
        if ~isempty(data)
            runDir = dateDir;
            break
        end
    end
end
assert(~isempty(runDir) && isfolder(runDir), ...
    'No saved analysis run found. Run RAFT_chunk_analysis or set options.runDir.');
[~,folderName] = fileparts(runDir);
if strcmp(folderName,'plot mat files'), runDir = fileparts(runDir); end
switch kind
    case 'profiles'
        pattern = 'vertical_profile_data.mat';
        analysisFolder = 'vertical profiles';
    case 'summary'
        pattern = 'frame_summary_data.mat';
        analysisFolder = 'frame summary';
    case 'maps'
        pattern = 'mean_speed_map_data.mat';
        analysisFolder = 'profile station visualization';
    case 'instantaneous'
        pattern = 'case*.mat';
        analysisFolder = 'instantaneous frames';
    otherwise
        error('Unknown plot section: %s',kind);
end
% An explicit date folder also selects one run, rather than mixing its runs.
runs = dir(fullfile(runDir,'run_*'));
runs = runs([runs.isdir]);
if ~isempty(runs)
    [~,order] = sort(string({runs.name}),'descend');
    runDir = fullfile(runs(order(1)).folder,runs(order(1)).name);
end

validPhases = ["pre-inception", "inception", "desinence", "fullcase"];
requestedPhase = "";
if isfield(options,'analysisPhase') && ~isempty(options.analysisPhase)
    requestedPhase = lower(strtrim(string(options.analysisPhase)));
    assert(isscalar(requestedPhase) && ...
        (strlength(requestedPhase) == 0 || ismember(requestedPhase,validPhases)), ...
        'options.analysisPhase must be pre-inception, inception, desinence, or fullcase.');
end
[~,folderName] = fileparts(runDir);
[~,parentName] = fileparts(fileparts(runDir));
phaseLayout = false;
% Accept an exact phase folder or one of its analysis folders as runDir.
if ismember(string(folderName),validPhases)
    phaseLayout = true;
    selectedPhase = string(folderName);
    assert(strlength(requestedPhase) == 0 || requestedPhase == selectedPhase, ...
        'options.analysisPhase conflicts with the selected phase folder: %s.',runDir);
    runDir = fullfile(runDir,analysisFolder);
elseif ismember(string(parentName),validPhases)
    phaseLayout = true;
    selectedPhase = string(parentName);
    assert(strcmp(folderName,analysisFolder), ...
        'Selected analysis folder does not contain %s data: %s.',kind,runDir);
    assert(strlength(requestedPhase) == 0 || requestedPhase == selectedPhase, ...
        'options.analysisPhase conflicts with the selected analysis folder: %s.',runDir);
else
    phaseFolders = validPhases(arrayfun(@(p) ...
        isfolder(fullfile(runDir,char(p))),validPhases));
    if ~isempty(phaseFolders)
        phaseLayout = true;
        if strlength(requestedPhase) > 0
            selectedPhase = requestedPhase;
            assert(ismember(selectedPhase,phaseFolders), ...
                'Phase %s was not saved in selected run: %s.',char(selectedPhase),runDir);
        else
            available = false(size(phaseFolders));
            for index = 1:numel(phaseFolders)
                candidates = dir(fullfile(runDir,char(phaseFolders(index)), ...
                    analysisFolder,'plot mat files',pattern));
                available(index) = ~isempty(candidates);
            end
            candidates = phaseFolders(available);
            assert(~isempty(candidates), ...
                'No %s data in selected run: %s. Enable that analysis for a requested phase.',kind,runDir);
            assert(numel(candidates) == 1, ...
                'Multiple phases contain %s data (%s). Set options.analysisPhase.', ...
                kind,char(strjoin(candidates,', ')));
            selectedPhase = candidates(1);
        end
        runDir = fullfile(runDir,char(selectedPhase),analysisFolder);
    elseif strlength(requestedPhase) > 0
        % Legacy runs carry the phase in their saved MAT files.
        selectedPhase = requestedPhase;
    end
end
plotMatDir = fullfile(runDir,'plot mat files');
if phaseLayout || isfolder(plotMatDir)
    matches = dir(fullfile(plotMatDir,pattern));
else
    % Older analyses stored MAT data beside plots or in export subfolders.
    matches = dir(fullfile(runDir,'**',pattern));
end
assert(~isempty(matches), ...
    ['No %s MAT data in latest/selected run: %s. Enable the corresponding ' ...
     'section in RAFT_chunk_analysis, or explicitly choose an older ' ...
     'options.runDir/options.dataFile.'],kind,runDir);
% Multiple legacy exports on one date: choose the latest section export only.
[~,newest] = max([matches.datenum]);
if strcmp(kind,'instantaneous')
    matches = matches(strcmp({matches.folder},matches(newest).folder));
    [~,order] = sort({matches.name});
    matches = matches(order);
else
    matches = matches(newest);
end
files = arrayfun(@(f) fullfile(f.folder,f.name),matches,'UniformOutput',false);
if strlength(requestedPhase) > 0
    for index = 1:numel(files)
        metadata = load(files{index},'analysisPhase');
        assert(isfield(metadata,'analysisPhase') && ...
            strcmpi(string(metadata.analysisPhase),requestedPhase), ...
            'Saved phase does not match options.analysisPhase in %s.',files{index});
    end
end
fprintf('Selected analysis run: %s\n',runDir);
for index = 1:numel(files), fprintf('Plot source: %s\n',files{index}); end
end
