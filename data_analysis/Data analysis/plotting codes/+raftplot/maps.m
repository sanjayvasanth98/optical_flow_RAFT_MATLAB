function maps(options)
[files,runDir] = raftplot.source(options,'maps');
data = load(files{1});
cases = raftplot.indices(options.caseIndices,numel(data.meanSpeedMaps),'caseIndices');
sourceCaseIndices = 1:numel(data.meanSpeedMaps);
if isfield(data,'caseIndices'), sourceCaseIndices = data.caseIndices; end
stations = raftplot.indices(options.stationIndices,numel(data.profileOffset_mm),'stationIndices');
limits = options.colorLimits;
if isempty(limits)
    maximum = 0;
    for index = 1:numel(data.meanSpeedMaps)
        values = data.meanSpeedMaps{index}; values = values(isfinite(values));
        if ~isempty(values), maximum = max(maximum,double(max(values))); end
    end
    if maximum == 0, maximum = 1; end
    limits = [0 maximum];
end
[folder,plotMatDir] = raftplot.output(options,runDir,'mean_speed_maps');
sourceFiles = files;
save(fullfile(plotMatDir,'plot_settings.mat'),'options','sourceFiles','limits');
stationColors = [0.95 0.70 0;0.95 0.22 0.45;0 0.65 0.78;0.68 0.28 0.88;0.97 0.46 0.08];
for index = cases
    coordinates = data.mapCoordinates{index};
    speed = data.meanSpeedMaps{index};
    fig = figure('Color','w','Position',options.figurePosition,'Visible',options.visible);
    ax = axes(fig);
    h = imagesc(ax,(coordinates.x_mm-coordinates.xThroat_mm)/data.throatHeight_mm, ...
        coordinates.y_mm/data.throatHeight_mm,speed);
    set(h,'AlphaData',isfinite(speed),'HandleVisibility','off');
    set(ax,'YDir','normal'); axis(ax,'image'); hold(ax,'on');
    colormap(ax,turbo(256)); clim(ax,limits);
    for station = stations
        x = coordinates.actualOffsets_mm(station)/data.throatHeight_mm;
        y = coordinates.y_mm([1 end])/data.throatHeight_mm;
        color = stationColors(mod(station-1,size(stationColors,1))+1,:);
        plot(ax,[x x],y,'-','Color',[0.08 0.08 0.08],'LineWidth',3.4,'HandleVisibility','off');
        plot(ax,[x x],y,'--','Color',color,'LineWidth',2, ...
            'DisplayName',sprintf('Station %d: +%.3f H',station, ...
            data.profileOffset_mm(station)/data.throatHeight_mm));
    end
    raftplot.style(ax,options);
    xlabel(ax,'(x - x_{throat})/H','FontSize',options.labelFontSize);
    ylabel(ax,'y/H','FontSize',options.labelFontSize);
    title(ax,sprintf('%s | Mean speed | frames %d:%d | %s', ...
        char(string(data.caseLabels{index})),data.frameRanges(index,1), ...
        data.frameRanges(index,2),char(data.analysisPhase)), ...
        'Interpreter','none','FontSize',options.titleFontSize,'FontWeight','normal');
    cb = colorbar(ax); cb.Label.String = 'Mean speed (m/s)';
    set(cb,'FontName',options.fontName,'FontSize',options.fontSize-1);
    legend(ax,'show','Location',options.legendLocation,'Orientation','horizontal', ...
        'Interpreter','none','Box','off','FontName',options.fontName,'FontSize',options.fontSize-1);
    safeLabel = regexprep(char(string(data.caseLabels{index})),'[^A-Za-z0-9_-]','_');
    name = sprintf('case%02d_%s_mean_speed_map_%s_frames%d-%d',sourceCaseIndices(index),safeLabel, ...
        char(data.analysisPhase),data.frameRanges(index,1),data.frameRanges(index,2));
    raftplot.export(fig,folder,name,options);
end
end
