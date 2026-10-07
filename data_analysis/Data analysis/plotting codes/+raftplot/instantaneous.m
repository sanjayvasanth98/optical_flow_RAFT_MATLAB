function instantaneous(options)
% Read one physical-coordinate saved frame at a time, including during scaling.
[files,runDir] = raftplot.source(options,'instantaneous');
cases = raftplot.indices(options.caseIndices,numel(files),'caseIndices');
metadata = cell(numel(files),1); frames = cell(numel(files),1);
maximum = -inf;
for index = cases
    M = matfile(files{index}); vars = who(M);
    assert(all(ismember({'u_firstN','v_firstN','x_mm','y_mm','sourceFrames','fps_case'},vars)), ...
        'Not an instantaneous velocity subset: %s',files{index});
    info = struct('x_mm',M.x_mm,'y_mm',M.y_mm, ...
        'sourceFrames',M.sourceFrames,'fps_case',M.fps_case);
    info.singleFrame = numel(size(M,'u_firstN')) < 3;
    [~,label] = fileparts(files{index});
    if ismember('caseLabel',vars), label = char(string(M.caseLabel)); end
    info.label = label;
    info.savedLimits = [];
    if ismember('uColorLimits',vars), info.savedLimits = M.uColorLimits; end
    metadata{index} = info;
    frames{index} = raftplot.indices(options.frameIndices,numel(info.sourceFrames),'frameIndices');
    if isempty(options.colorLimits)
        for frame = frames{index}
            if info.singleFrame
                u = M.u_firstN(:,:);
            else
                u = M.u_firstN(:,:,frame);
            end
            values = u(isfinite(u));
            if ~isempty(values), maximum = max(maximum,double(max(values))); end
        end
    end
end
limits = options.colorLimits;
if isempty(limits)
    assert(isfinite(maximum),'No finite instantaneous U values in selected frames.');
    limits = [options.negativeColorMin_mps maximum];
    if limits(2) <= limits(1), limits(2) = limits(1)+1; end
    savedLimits = cellfun(@(info) info.savedLimits,metadata(cases),'UniformOutput',false);
    if options.useSavedColorLimits && ~isempty(savedLimits{1}) && ...
            all(cellfun(@(value) isequal(value,savedLimits{1}),savedLimits))
        limits = savedLimits{1};
    end
end
assert(numel(limits)==2 && all(isfinite(limits)) && limits(2)>limits(1), ...
    'colorLimits must be [minimum maximum] with finite increasing values.');
if limits(1) < 0 && limits(2) > 0
    anchors = [limits(1) 0 0.25*limits(2) 0.5*limits(2) 0.75*limits(2) limits(2)];
    colors = [0.05 0.10 0.48;0.18 0.64 0.87;0.25 0.79 0.40; ...
        0.98 0.91 0.20;0.97 0.46 0.08;0.68 0.03 0.09];
    cmap = interp1(anchors,colors,linspace(limits(1),limits(2),256));
    ticks = [limits(1) linspace(0,limits(2),7)];
else
    cmap = turbo(256); ticks = linspace(limits(1),limits(2),8);
end
[folder,plotMatDir] = raftplot.output(options,runDir,'instantaneous');
sourceFiles = files;
save(fullfile(plotMatDir,'plot_settings.mat'),'options','sourceFiles','limits');
for index = cases
    info = metadata{index}; M = matfile(files{index});
    [~,caseBase] = fileparts(files{index});
    caseFolder = fullfile(folder,caseBase); mkdir(caseFolder);
    [X,Y] = meshgrid(info.x_mm,info.y_mm);
    for frame = frames{index}
        % Already x-right/y-up/V-up in the subset: no extra flip or sign change.
        if info.singleFrame
            u = M.u_firstN(:,:); v = M.v_firstN(:,:);
        else
            u = M.u_firstN(:,:,frame); v = M.v_firstN(:,:,frame);
        end
        fig = figure('Color','w','Position',options.figurePosition,'Visible',options.visible);
        ax = axes(fig); h = imagesc(ax,info.x_mm,info.y_mm,u);
        set(h,'AlphaData',isfinite(u)); set(ax,'YDir','normal');
        axis(ax,'image'); hold(ax,'on'); colormap(ax,cmap); clim(ax,limits);
        if options.showStreamlines
            raftplot.streamlines(ax,X,Y,double(u),double(v), ...
                isfinite(u) & isfinite(v),options.streamlineDensity);
        end
        raftplot.style(ax,options);
        xlabel(ax,'x (mm)','FontSize',options.labelFontSize);
        ylabel(ax,'y (mm)','FontSize',options.labelFontSize);
        sourceFrame = info.sourceFrames(frame);
        title(ax,sprintf('%s: Instantaneous streamwise U (flow frame %d)', ...
            info.label,sourceFrame),'Interpreter','none', ...
            'FontSize',options.titleFontSize,'FontWeight','normal');
        cb = colorbar(ax); cb.Label.String = 'Streamwise U (m/s)';
        cb.Ticks = ticks;
        cb.TickLabels = arrayfun(@(value) sprintf('%.3g',value),ticks,'UniformOutput',false);
        set(cb,'FontName',options.fontName,'FontSize',options.fontSize-1);
        raftplot.export(fig,caseFolder,sprintf('frame_%06d',sourceFrame),options);
    end
end
end
