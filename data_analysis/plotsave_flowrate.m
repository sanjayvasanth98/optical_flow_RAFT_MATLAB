%% Roughness vs Flow Rate Plot (Publication Quality)

% -------------------------------------------------------
% Input Data
% -------------------------------------------------------
roughness_um = [5 12 20 30 53 80];
flowrate_lpm = [47.5 50 51 52.3 54.5 55.5];

labels   = 'Flow Rate';
savePath = 'G:\My Drive\Research\cavitation\My research papers\Pre-Cavitation dynamics of the nuclei\inkscape';
plotTitle = 'Roughness_vs_FlowRate';
xlab = 'Surface Roughness (micrometer)';
ylab = 'Flow Rate (L/min)';

% -------------------------------------------------------
% Create Figure (normal visible figure)
% -------------------------------------------------------
F = figure('Color','white','Position',[200 200 900 600]);
hold on;

% Scientific deep blue
col = [0 0.2 0.8];

% Plot
plot(roughness_um, flowrate_lpm, ...
    '-o', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', col, ...
    'Color', col, ...
    'DisplayName', labels);

% Axis Padding
xlim([0, max(roughness_um)*1.1]);
ylim([min(flowrate_lpm)*0.9, max(flowrate_lpm)*1.1]);

% Aesthetics
set(gca, 'FontName','Times New Roman', ...
         'FontSize',18, ...
         'LineWidth',1.4, ...
         'Box','on', ...
         'TickDir','out', ...
         'XMinorTick','on', ...
         'YMinorTick','on');

xlabel(xlab, 'FontSize',20,'FontName','Times New Roman', 'Interpreter','none');
ylabel(ylab, 'FontSize',20,'FontName','Times New Roman', 'Interpreter','none');
title(strrep(plotTitle,'_',' '), ...
      'FontSize',22,'FontName','Times New Roman', 'Interpreter','none');

legend('Location','northwest','FontSize',16,'Interpreter','none');

hold off;

% -------------------------------------------------------
% Save Figure (simple + robust)
% -------------------------------------------------------
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

pngFile = fullfile(savePath, 'flowrate.png');
epsFile = fullfile(savePath, 'flowrate.eps');

saveas(F, pngFile);              % saves PNG
saveas(F, epsFile, 'epsc');      % saves color EPS

fprintf('Saved:\n%s\n%s\n', pngFile, epsFile);
