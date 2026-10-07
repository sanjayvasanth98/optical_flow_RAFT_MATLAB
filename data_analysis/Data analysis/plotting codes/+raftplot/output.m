function [folder,plotMatDir] = output(options,runDir,section)
folder = char(options.outputDir);
if isempty(folder)
    base = fullfile(runDir,'replots', ...
        [section '_' datestr(now,'yyyymmdd_HHMMSS')]);
    folder = base;
    index = 2;
    while exist(folder,'dir') || exist(folder,'file')
        folder = sprintf('%s_%02d',base,index);
        index = index+1;
    end
end
if ~isfolder(folder)
    [ok,message] = mkdir(folder);
    assert(ok,'Cannot create plot folder: %s',message);
end
plotMatDir = fullfile(folder,'plot mat files');
[ok,message] = mkdir(plotMatDir);
assert(ok,'Cannot create plot MAT folder: %s',message);
fprintf('Refined plots: %s\n',folder);
end
