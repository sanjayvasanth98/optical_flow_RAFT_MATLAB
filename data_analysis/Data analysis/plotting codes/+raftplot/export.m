function export(fig,folder,name,options)
base = fullfile(folder,name);
if options.savePNG, exportgraphics(fig,[base '.png'],'Resolution',options.resolution); end
if options.saveFIG, savefig(fig,[base '.fig']); end
if options.savePDF, exportgraphics(fig,[base '.pdf'],'ContentType','vector'); end
if options.closeAfterSave, close(fig); end
end
