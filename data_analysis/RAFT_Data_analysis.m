%% main_postprocess_raft.m
% Cluster-safe postprocessing for RAFT instantaneous velocity MAT-files
% - Tiled (subplot-style) profiles: throat + downstream stations
% - Profiles for: |V|_mean, U_mean (Vx), V_mean (Vy)
% - Multi-case overlay (legend at bottom of whole figure)
% - Reynolds stresses: uu, vv, uv (turbo colormap)
% - Low RAM: uses matfile + chunked accumulation (no full 3D load)
%
% OUTPUT (all saved under user-chosen folder):
%   <userOutDir>/
%     Velocity profiles/
%       VelProfile_Vmag_dxXmm_NstationsY.png
%       VelProfile_U_dxXmm_NstationsY.png
%       VelProfile_V_dxXmm_NstationsY.png
%     Reynolds Stresses/
%       <caseName>_ReStress_uu.png, _vv.png, _uv.png
%       <caseName>_ReynoldsStresses.mat
%
% NOTE:
% - No temp folders created.
% - Uses *_velocity.mat produced by your RAFT script:
%   u_all, v_all, mm_per_pixel, fps, m_per_pixel, maskROI,
%   (often) x_mm, y_mm, maskROI_phys, u_mean_phys, v_mean_phys, velMean_phys, x_throat_mm
%
% Author: Sanjay Vasanth | last edit: 2/17/2026
%----------------------------------------------
clear; clc;

% ---------------- USER INPUTS ----------------
% 1) List your cases here
matPaths = { ...
    "/home/kbsanjayvasanth/Inception_raft_test/smooth_w_particles/mat files/Smooth_48lpm_inc_velocity.mat"    % <-- EDIT THIS
    % "/path/to/case2_velocity.mat"
};

% 2) Labels (same order as matPaths) for legend
caseLabels = { ...
    "Smooth 48 lpm"  % <-- EDIT THIS
    % "Case2"
};

% 3) Output folder (ALL results go under this directory)
userOutDir = "/home/kbsanjayvasanth/Inception_raft_test/data_analysis/results";  % <-- EDIT THIS

% 4) Settings
dx_mm       = 0.25;     % profile spacing downstream from throat   % <-- EDIT THIS if needed
nStations   = 10;       % number of downstream profile locations to plot (excluding throat)   % <-- EDIT THIS
chunkFrames = 100;      % chunk size for stress computation (IO/RAM tradeoff)
dpiExport   = 600;      % export DPI

makeProfiles = true;
makeStresses = true;

% ---------------- RUN ----------------
tAll = tic;
fprintf('\n=== RAFT POSTPROCESS START ===\n');
fprintf('Cases: %d\n', numel(matPaths));
assert(numel(matPaths) == numel(caseLabels), "caseLabels must match matPaths length.");

% Make user output subfolders
profilesDir = fullfile(userOutDir, "Velocity profiles");
stressesDir = fullfile(userOutDir, "Reynolds Stresses");
if ~exist(profilesDir,'dir'), mkdir(profilesDir); end
if ~exist(stressesDir,'dir'), mkdir(stressesDir); end

fprintf('Output root: %s\n', userOutDir);
fprintf('  Profiles:  %s\n', profilesDir);
fprintf('  Stresses:  %s\n', stressesDir);

% Read info handles
infos = cell(numel(matPaths),1);
for i = 1:numel(matPaths)
    infos{i} = read_case_info(matPaths{i});
end

if makeProfiles
    plot_profiles_tiled_from_instantaneous(infos, caseLabels, dx_mm, nStations, profilesDir, dpiExport);
end

if makeStresses
    plot_reynolds_stresses_from_instantaneous(infos, chunkFrames, stressesDir, dpiExport);
end

fprintf('=== DONE. Total elapsed: %.1f s ===\n\n', toc(tAll));

%% =========================================================================
% 1) TILED PROFILES (OVERLAY): throat + N downstream stations
%    One figure per field, overlays all cases in each tile
% =========================================================================
function plot_profiles_tiled_from_instantaneous(infos, caseLabels, dx_mm, nStations, profilesDir, dpiExport)
if nargin < 3 || isempty(dx_mm), dx_mm = 0.25; end
if nargin < 4 || isempty(nStations), nStations = 10; end

nC = numel(infos);
assert(nC >= 1);

% Check throats exist
for c = 1:nC
    if isnan(infos{c}.x_throat_mm)
        error("x_throat_mm missing in %s. Your RAFT script must save x_throat_mm into *_velocity.mat.", infos{c}.matPath);
    end
end

% Use stations based on reference case (case 1)
ref = infos{1};
xStations = ref.x_throat_mm + (0:nStations)*dx_mm;  % includes throat
xStations = xStations(xStations <= ref.x_mm(end));

% Map stations per case
ixStationsAll = cell(nC,1);
for c = 1:nC
    ixStationsAll{c} = nearest_index_vec(infos{c}.x_mm, xStations);
end

% y-range intersection across cases
yMin = max(cellfun(@(ci) ci.y_min_roi_mm, infos));
yMax = min(cellfun(@(ci) ci.y_max_roi_mm, infos));

% Load/compute means (2D only)
uMean = cell(nC,1); vMean = cell(nC,1); qMean = cell(nC,1);
for c = 1:nC
    info = infos{c};
    fprintf('[Means] %s\n', info.caseName);
    uMean{c} = get_mean_field(info.M, "u_mean_phys", info, 'u');
    vMean{c} = get_mean_field(info.M, "v_mean_phys", info, 'v');
    qMean{c} = get_mean_field(info.M, "velMean_phys", info, 'q');
end

% Make overlay figures (legend bottom of whole plot)
make_profile_tiled_overlay(infos, caseLabels, qMean, ixStationsAll, xStations, yMin, yMax, "Vmag", dx_mm, profilesDir, dpiExport);
make_profile_tiled_overlay(infos, caseLabels, uMean, ixStationsAll, xStations, yMin, yMax, "U",    dx_mm, profilesDir, dpiExport);
make_profile_tiled_overlay(infos, caseLabels, vMean, ixStationsAll, xStations, yMin, yMax, "V",    dx_mm, profilesDir, dpiExport);
end

function make_profile_tiled_overlay(infos, caseLabels, Fmeans, ixStationsAll, xStations, yMin, yMax, fieldTag, dx_mm, profilesDir, dpiExport)
nC = numel(infos);
nP = numel(xStations);         % tiles = throat + downstream

% tile geometry: rectangles (wider than tall)
nCols = min(8, nP);
nRows = ceil(nP / nCols);

% Figure
fig = figure('Visible','off','Color','w','Position',[100 100 1700 900]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

% Big title + xlabel text
switch fieldTag
    case "Vmag"
        bigTitle = sprintf('Mean |V| profiles (dx=%.2f mm)', dx_mm);
        xLab = '$\overline{|V|}\;(\mathrm{m/s})$';
    case "U"
        bigTitle = sprintf('Mean U=V_x profiles (dx=%.2f mm)', dx_mm);
        xLab = '$\overline{u}\;(\mathrm{m/s})$';
    case "V"
        bigTitle = sprintf('Mean V=V_y profiles (dx=%.2f mm)', dx_mm);
        xLab = '$\overline{v}\;(\mathrm{m/s})$';
    otherwise
        bigTitle = sprintf('Mean profiles (dx=%.2f mm)', dx_mm);
        xLab = '$\overline{(\cdot)}$';
end
title(tl, bigTitle, 'Interpreter','none', 'FontName','Times New Roman', 'FontSize', 20);

% colors: MATLAB chooses distinct colors
cols = lines(max(nC,7));

% Create a hidden axes ONLY for the legend (MATLAB cannot legend(tiledlayout))
legendAx = axes(fig, 'Visible','off', 'Position',[0 0 1 1], 'HitTest','off');
hold(legendAx,'on');
legendLineHandles = gobjects(nC,1);
for c = 1:nC
    legendLineHandles(c) = plot(legendAx, nan, nan, '-', ...
        'Color', cols(c,:), 'LineWidth', 1.25);
end

% Plot tiles
for p = 1:nP
    ax = nexttile(tl);
    hold(ax,'on'); box(ax,'on');

    % title per tile
    if p == 1
        tstr = "Throat";
    else
        dx_here = xStations(p) - xStations(1);
        tstr = sprintf('x_t + %.2f mm', dx_here);
    end
    title(ax, tstr, 'Interpreter','none', 'FontName','Times New Roman', 'FontSize', 12);

    % Each case overlay
    for c = 1:nC
        info = infos{c};
        Fmean = Fmeans{c};

        xi = ixStationsAll{c}(p);

        % y mask intersection
        yMask = info.y_mm >= yMin & info.y_mm <= yMax;
        yPlot = info.y_mm(yMask);

        prof = Fmean(:,xi);
        prof = prof(yMask);

        plot(ax, prof, yPlot, '-', 'Color', cols(c,:), 'LineWidth', 1.25);
    end

    % Style: Times New Roman, dotted minor grid
    set(ax,'FontName','Times New Roman','FontSize',11,'LineWidth',1.1,'Box','on');
    grid(ax,'on'); grid(ax,'minor');
    ax.GridLineStyle      = '-';
    ax.MinorGridLineStyle = ':';
    ax.GridAlpha          = 0.20;
    ax.MinorGridAlpha     = 0.15;
    ax.XMinorTick         = 'on';
    ax.YMinorTick         = 'on';

    ylim(ax, [yMin, yMax]);

    % X limit: extend +1.5 beyond max value (user request)
    xl = xlim(ax);
    xlim(ax, [xl(1), xl(2)+1.5]);

    % labels on EVERY tile (your request)
    xlabel(ax, xLab, 'Interpreter','latex', 'FontSize', 12);
    ylabel(ax, '$y\;(\mathrm{mm})$', 'Interpreter','latex', 'FontSize', 12);
end

% ---- Legend at the bottom of the whole plot ----
lg = legend(legendAx, legendLineHandles, caseLabels, ...
    'Interpreter','none', 'Box','off', 'Orientation','horizontal');
lg.FontName = 'Times New Roman';
lg.FontSize = 12;
lg.Units = 'normalized';
lg.Location = 'none';
lg.Position = [0.10 0.015 0.80 0.04];   % [left bottom width height]

% Give space at bottom so legend doesn't overlap x-labels
tl.OuterPosition = [0 0.06 1 0.94];

% Save (600 dpi)
outPng = fullfile(profilesDir, sprintf('VelProfile_%s_dx%.3fmm_N%d.png', fieldTag, dx_mm, nP-1));
exportgraphics(fig, outPng, 'Resolution', dpiExport);
close(fig);

fprintf('Saved: %s\n', outPng);
end

%% =========================================================================
% 2) Reynolds stresses (uu, vv, uv) per case from instantaneous u_all/v_all
% =========================================================================
function plot_reynolds_stresses_from_instantaneous(infos, chunkFrames, stressesDir, dpiExport)
if nargin < 2 || isempty(chunkFrames), chunkFrames = 100; end

for c = 1:numel(infos)
    tCase = tic;
    info = infos{c};
    fprintf('\n[STRESS %d/%d] %s\n', c, numel(infos), info.caseName);

    uMean = get_mean_field(info.M, "u_mean_phys", info, 'u');
    vMean = get_mean_field(info.M, "v_mean_phys", info, 'v');

    sumUU = zeros(info.H, info.W, 'double');
    sumVV = zeros(info.H, info.W, 'double');
    sumUV = zeros(info.H, info.W, 'double');
    nTot  = 0;

    for k0 = 1:chunkFrames:info.K
        k1 = min(info.K, k0+chunkFrames-1);

        U = info.M.u_all(:,:,k0:k1);
        V = info.M.v_all(:,:,k0:k1);

        U = flipud(U);
        V = -flipud(V);

        for kk = 1:size(U,3)
            u = double(U(:,:,kk));
            v = double(V(:,:,kk));
            sumUU = sumUU + u.^2;
            sumVV = sumVV + v.^2;
            sumUV = sumUV + (u.*v);
        end
        nTot = nTot + size(U,3);

        if k1==info.K || mod(k0-1, chunkFrames*10)==0
            fprintf('  Frames: %d/%d\n', k1, info.K);
        end
    end

    EUU = single(sumUU/max(nTot,1) - double(uMean).^2);
    EVV = single(sumVV/max(nTot,1) - double(vMean).^2);
    EUV = single(sumUV/max(nTot,1) - double(uMean).*double(vMean));

    EUU(~info.maskROI_phys) = NaN;
    EVV(~info.maskROI_phys) = NaN;
    EUV(~info.maskROI_phys) = NaN;

    % Save plots and MAT into user folder (use caseName prefix)
    plot_map_pretty(EUU, info, stressesDir, sprintf('%s_ReStress_uu.png', info.caseName), '$\langle u''u''\rangle\;(\mathrm{m^2/s^2})$', dpiExport);
    plot_map_pretty(EVV, info, stressesDir, sprintf('%s_ReStress_vv.png', info.caseName), '$\langle v''v''\rangle\;(\mathrm{m^2/s^2})$', dpiExport);
    plot_map_pretty(EUV, info, stressesDir, sprintf('%s_ReStress_uv.png', info.caseName), '$\langle u''v''\rangle\;(\mathrm{m^2/s^2})$', dpiExport);

    save(fullfile(stressesDir, sprintf('%s_ReynoldsStresses.mat', info.caseName)), 'EUU','EVV','EUV','-v7.3');
    fprintf('  [OK] Saved stresses + plots (user folder). Elapsed: %.1f s\n', toc(tCase));
end
end

function plot_map_pretty(F, info, outDir, outName, cbarLabelLatex, dpiExport)
fig = figure('Visible','off','Color','w','Position',[100 100 1200 800]);
ax = axes(fig);

yMask = info.y_mm >= info.y_min_roi_mm & info.y_mm <= info.y_max_roi_mm;
F2 = F(yMask,:);
yPlot = info.y_mm(yMask);

imagesc(ax, info.x_mm, yPlot, F2);
set(ax,'YDir','normal'); axis(ax,'image');
colormap(ax, turbo);

cb = colorbar(ax, 'eastoutside');
cb.TickDirection = 'out';
cb.Label.String  = cbarLabelLatex;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;
cb.Label.Rotation = 90;
cb.Label.Units = 'normalized';
cb.Label.Position = [3 0.5 0];

ax.Position = [0.08 0.10 0.78 0.82];

xlabel(ax, 'X (mm)', 'FontSize', 16, 'Interpreter','latex');
ylabel(ax, 'Y (mm)', 'FontSize', 16, 'Interpreter','latex');
title(ax, sprintf('%s — Reynolds stress', info.caseName), 'FontSize', 18, 'Interpreter','none');
set(ax,'FontSize',13,'LineWidth',1.2,'Box','on','FontName','Times New Roman');

outPng = fullfile(outDir, outName);
exportgraphics(fig, outPng, 'Resolution', dpiExport);
close(fig);

fprintf('    Saved: %s\n', outPng);
end

%% =========================================================================
% Helpers
% =========================================================================
function info = read_case_info(matPath)
if ~isfile(matPath), error("MAT file not found: %s", matPath); end

[~, baseName, ~] = fileparts(matPath);
caseName = erase(baseName, "_velocity");

M = matfile(matPath);

wu = whos(M, 'u_all');
if isempty(wu), error('u_all not found in %s', matPath); end
sz = wu.size;
H = sz(1); W = sz(2);
if numel(sz) >= 3, K = sz(3); else, K = 1; end

mm_per_pixel = M.mm_per_pixel;

maskROI = logical(M.maskROI);
if has_var(M, "maskROI_phys")
    maskROI_phys = logical(M.maskROI_phys);
else
    maskROI_phys = flipud(maskROI);
end

if has_var(M, "x_mm"), x_mm = M.x_mm; else, x_mm = (0:W-1)*mm_per_pixel; end
if has_var(M, "y_mm"), y_mm = M.y_mm; else, y_mm = (0:H-1)*mm_per_pixel; end

rowsIn = find(any(maskROI_phys,2));
if isempty(rowsIn)
    y_min_roi_mm = y_mm(1); y_max_roi_mm = y_mm(end);
else
    y_min_roi_mm = y_mm(min(rowsIn)); y_max_roi_mm = y_mm(max(rowsIn));
end

x_throat_mm = NaN;
if has_var(M, "x_throat_mm"), x_throat_mm = M.x_throat_mm; end

% Store the original MAT file folder only as metadata (NOT used for saving)
[baseDir,~,~] = fileparts(matPath);

info = struct( ...
    'matPath',matPath, ...
    'baseDir',baseDir, ...
    'caseName',caseName, ...
    'M',M, ...
    'H',H,'W',W,'K',K, ...
    'mm_per_pixel',mm_per_pixel, ...
    'maskROI_phys',maskROI_phys, ...
    'x_mm',x_mm,'y_mm',y_mm, ...
    'y_min_roi_mm',y_min_roi_mm,'y_max_roi_mm',y_max_roi_mm, ...
    'x_throat_mm',x_throat_mm);
end

function tf = has_var(M, varName)
try
    tf = ~isempty(whos(M, varName));
catch
    tf = false;
end
end

function ix = nearest_index_vec(xVec, xQueries)
ix = zeros(size(xQueries));
for i = 1:numel(xQueries)
    [~,ix(i)] = min(abs(xVec - xQueries(i)));
end
end

function Fmean = get_mean_field(M, varName, info, which)
if has_var(M, varName)
    Fmean = M.(varName);
else
    Fmean = compute_mean_from_instant(M, info, which, 100);
end
Fmean(~info.maskROI_phys) = NaN;
end

function meanField = compute_mean_from_instant(M, info, which, chunkFrames)
sumF = zeros(info.H, info.W, 'double');
nTot = 0;

for k0 = 1:chunkFrames:info.K
    k1 = min(info.K, k0+chunkFrames-1);

    switch which
        case 'u'
            F = M.u_all(:,:,k0:k1);
            F = flipud(F);
        case 'v'
            F = M.v_all(:,:,k0:k1);
            F = -flipud(F);
        case 'q'
            U = M.u_all(:,:,k0:k1); U = flipud(U);
            V = M.v_all(:,:,k0:k1); V = -flipud(V);
            F = zeros(info.H, info.W, size(U,3), 'like', U);
            for kk = 1:size(U,3)
                F(:,:,kk) = hypot(U(:,:,kk), V(:,:,kk));
            end
        otherwise
            error('Unknown mean request: %s', which);
    end

    for kk = 1:size(F,3)
        sumF = sumF + double(F(:,:,kk));
    end
    nTot = nTot + size(F,3);
end

meanField = single(sumF / max(nTot,1));
end

