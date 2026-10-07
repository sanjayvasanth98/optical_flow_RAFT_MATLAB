function summary(options)
[files,runDir] = raftplot.source(options,'summary');
data = load(files{1});
t = data.frameSummary;
cases = raftplot.indices(options.caseIndices,numel(data.caseLabels),'caseIndices');
sourceCaseIndices = 1:numel(data.caseLabels);
if isfield(data,'caseIndices'), sourceCaseIndices = data.caseIndices; end
[folder,plotMatDir] = raftplot.output(options,runDir,'frame_summary');
sourceFiles = files;
save(fullfile(plotMatDir,'plot_settings.mat'),'options','sourceFiles');
fields = {'meanU_mps','meanV_mps','meanSpeed_mps'};
labels = {'Mean U (m/s)','Mean V, image-down positive (m/s)','Mean speed (m/s)'};
colors = [0 0 0;0 0.55 0;0 0 1;1 0 0;0.90 0.45 0;0.55 0.20 0.70];
styles = {'-','--','-','--','-.',':'};
fig = figure('Color','w','Visible',options.visible,'Position',options.figurePosition);
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
for metric = 1:numel(fields)
    ax = nexttile(layout); hold(ax,'on');
    for index = cases
        sourceCaseIndex = sourceCaseIndices(index);
        rows = t.caseIndex == sourceCaseIndex;
        [time,order] = sort(t.time_s(rows));
        values = t.(fields{metric})(rows); values = values(order);
        plot(ax,time,values,'Color',colors(mod(sourceCaseIndex-1,6)+1,:), ...
            'LineStyle',styles{mod(sourceCaseIndex-1,6)+1},'LineWidth',1.5, ...
            'DisplayName',char(string(data.caseLabels{index})));
    end
    raftplot.style(ax,options);
    ylabel(ax,labels{metric},'FontSize',options.labelFontSize);
    xlabel(ax,'Time (s)','FontSize',options.labelFontSize);
    if metric == 1
        title(ax,sprintf('ROI frame means | %s',char(data.analysisPhase)), ...
            'Interpreter','none','FontSize',options.titleFontSize,'FontWeight','normal');
        legend(ax,'show','Location',options.legendLocation,'Interpreter','none','Box','off');
    end
end
raftplot.export(fig,folder,'roi_frame_summary',options);
end
