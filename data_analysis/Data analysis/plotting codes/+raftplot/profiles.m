function profiles(options,metricIndices)
[files,runDir] = raftplot.source(options,'profiles');
data = load(files{1});
fields = {'meanU_over_Ub','sqrtReynoldsShear_over_Ub', ...
    'sigmaU_over_Ub','sigmaV_over_Ub','tkeInPlane_over_Ub2'};
names = {'Mean streamwise U','Reynolds shear stress', ...
    'Streamwise velocity standard deviation', ...
    'Wall-normal velocity standard deviation','In-plane turbulent kinetic energy'};
labels = {'$\overline{u}/U_b$','$\sqrt{-\overline{u''v''}}/U_b$', ...
    '$\sigma_u/U_b$','$\sigma_v/U_b$','$k_{2D}/U_b^2$'};
tags = {'mean_streamwise_U','sqrt_reynolds_shear', ...
    'sigma_streamwise_U','sigma_wall_normal_V','tke_in_plane'};
cases = raftplot.indices(options.caseIndices,size(data.verticalProfiles,1),'caseIndices');
sourceCaseIndices = 1:size(data.verticalProfiles,1);
if isfield(data,'caseIndices'), sourceCaseIndices = data.caseIndices; end
stations = raftplot.indices(options.stationIndices,size(data.verticalProfiles,2),'stationIndices');
assert(max(cases) <= 6,'Profile styles support up to six saved cases.');
if options.useSavedMarkerOptions && isfield(data,'markerOptions')
    options.markerOptions = data.markerOptions;
end
[folder,plotMatDir] = raftplot.output(options,runDir,'vertical_profiles');
sourceFiles = files;
save(fullfile(plotMatDir,'plot_settings.mat'),'options','sourceFiles','metricIndices');
stationDirs = cell(1,size(data.verticalProfiles,2));
for station = stations
    stationDirs{station} = fullfile(folder, ...
        sprintf('station%02d_xplus_%.2fmm',station,data.profileOffset_mm(station)));
    [ok,message] = mkdir(stationDirs{station});
    assert(ok,'Cannot create station folder: %s',message);
end
% Retain original y extent and metric scales across ALL saved cases/stations.
maxY = 0;
for index = 1:numel(data.verticalProfiles)
    p = data.verticalProfiles{index};
    validY = p.y_over_H(isfinite(p.meanU_mps) & isfinite(p.y_over_H));
    if ~isempty(validY), maxY = max(maxY,max(validY)); end
end
for metric = metricIndices
    minimum = inf; maximum = -inf;
    for index = 1:numel(data.verticalProfiles)
        values = data.verticalProfiles{index}.(fields{metric});
        values = values(isfinite(values));
        if ~isempty(values)
            minimum = min(minimum,min(values)); maximum = max(maximum,max(values));
        end
    end
    if ~isfinite(minimum)
        warning('No finite values for %s; skipping.',names{metric}); continue
    end
    limits = [min(0,minimum) max(0,maximum)];
    padding = max(0.03*diff(limits),0.01*(diff(limits)==0));
    if limits(1) ~= 0, limits(1) = limits(1)-padding; end
    limits(2) = limits(2)+padding;
    for station = stations
        fig = figure('Color','w','Position',options.figurePosition,'Visible',options.visible);
        ax = axes(fig); hold(ax,'on'); count = 0;
        for index = cases
            p = data.verticalProfiles{index,station};
            if ~any(isfinite(p.(fields{metric})) & isfinite(p.y_over_H)), continue; end
            raftplot.profileSeries(ax,p,fields{metric},sourceCaseIndices(index), ...
                char(string(data.caseLabels{index})),options.markerOptions);
            count = count+1;
        end
        if count == 0, close(fig); continue; end
        xlim(ax,limits); ylim(ax,[0 maxY+max(0.03*maxY,0.01)]);
        raftplot.style(ax,options);
        xlabel(ax,labels{metric},'Interpreter','latex','FontSize',options.labelFontSize);
        ylabel(ax,'y/H','FontSize',options.labelFontSize);
        title(ax,sprintf('%s | throat + %.3f H | %s',names{metric}, ...
            data.profileOffset_mm(station)/data.throatHeight_mm,char(data.analysisPhase)), ...
            'Interpreter','none','FontSize',options.titleFontSize,'FontWeight','normal');
        legend(ax,'show','Location',options.legendLocation,'Interpreter','none', ...
            'Box','off','FontName',options.fontName,'FontSize',options.fontSize-1);
        name = sprintf('%s_profile_xplus_%.2fmm_%s_station%02d',tags{metric}, ...
            data.profileOffset_mm(station),char(data.analysisPhase),station);
        raftplot.export(fig,stationDirs{station},name,options);
    end
end
end
