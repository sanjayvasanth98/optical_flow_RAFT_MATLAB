function style(ax,options)
set(ax,'FontName',options.fontName,'FontSize',options.fontSize, ...
    'LineWidth',1,'TickDir','in','Box','on','Color','w','Layer','top');
if options.showGrid, grid(ax,'on'); else, grid(ax,'off'); end
if ~isempty(options.xLimits), xlim(ax,options.xLimits); end
if ~isempty(options.yLimits), ylim(ax,options.yLimits); end
end
