function handles = streamlines(ax,X,Y,U,V,mask,density)
U(~mask) = NaN; V(~mask) = NaN;
handles = streamslice(ax,X,Y,U,V,density);
set(handles,'Color','k','LineWidth',0.8,'HandleVisibility','off');
for index = 1:numel(handles)
    x = handles(index).XData; y = handles(index).YData;
    if isempty(x) || isempty(y), continue; end
    inside = interp2(X,Y,single(mask),x,y,'nearest',0) >= 0.5;
    x(~inside) = NaN; y(~inside) = NaN;
    handles(index).XData = x; handles(index).YData = y;
end
end
