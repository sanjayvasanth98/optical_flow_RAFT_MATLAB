function profileSeries(ax,profile,field,index,label,markers)
% Swap the styles assigned to Rough5 (case 1) and Smooth (case 6).
if index == 1
    index = 6;
elseif index == 6
    index = 1;
end
values = profile.(field);
valid = find(isfinite(values) & isfinite(profile.y_over_H));
colors = [0 0 0;0 0.55 0;0 0 1;1 0 0;0.90 0.45 0;0.55 0.20 0.70];
styles = {'none','none','-','--','none',':'};
widths = [1.3 1.4 2.5 2 2 2.5];
symbol = 'none'; face = 'none'; markerRows = []; markerSize = 6;
if index <= 2 || index == 5
    if index == 1
        count = markers.blackCount; symbol = markers.blackSymbol;
        markerSize = markers.blackSize; face = 'none';
    elseif index == 2
        count = markers.greenCount; symbol = markers.greenSymbol;
        markerSize = markers.greenSize; face = colors(index,:);
    else
        count = markers.blackCount; symbol = '^'; face = 'none';
        markerSize = markers.greenSize;
    end
    assert(isscalar(count) && isfinite(count) && count >= 0 && count == fix(count), ...
        'Marker counts must be nonnegative integers.');
    if count > 0
        markerRows = unique(valid(round(linspace(1,numel(valid),min(count,numel(valid))))));
    end
end
plot(ax,values,profile.y_over_H,'Color',colors(index,:), ...
    'LineStyle',styles{index},'LineWidth',widths(index), ...
    'Marker',symbol,'MarkerIndices',markerRows,'MarkerSize',markerSize, ...
    'MarkerFaceColor',face,'MarkerEdgeColor',colors(index,:),'DisplayName',label);
end
