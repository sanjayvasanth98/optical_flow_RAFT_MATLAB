%% ============================================================
% Plot 6 velocity magnitude cases with throat location overlaid
% - Reads 6 velocity .mat files
% - Reads 6 throat .mat files
% - Flips velMean_phys upside down
% - Uses user-input calibration (mm/pixel)
% - Overlays throat location as red '+'
% - One figure with 6 subplots
% ============================================================

clear; clc; close all;

%% ---------------- USER INPUT ----------------

% Folder paths
velPath    = "E:\VelocitydataRAFT\velocity mag 5x\";
throatPath = 'E:\VelocitydataRAFT\velocity mag 5x\throat\';

% Velocity files
velFiles = {
    'smooth_TimeAvgVelField.mat'
    'p10s20_TimeAvgVelField.mat'
    'p10s30_TimeAvgVelField.mat'
    'p10s50_TimeAvgVelField.mat'
    'p10s70_TimeAvgVelField.mat'
    'p10s100_TimeAvgVelField.mat'
};


% Throat files (same order as velocity files)
throatFiles = {
    'smooth_5x_throat.mat'
    'p10s20_5x_throat.mat'
    'p10s30_5x_throat.mat'
    'p10s50_5x_throat.mat'
    'p10s70_5x_throat.mat'
    'p10s100_5x_throat.mat'
};

%calibration
mm_per_pixel = 0.00375;   % mm per pixel
throatInPixels = true;   % true = throat coordinates in pixels

savePath = 'E:\VelocitydataRAFT\velocity mag 5x\extracted profiles\';
if ~exist(savePath, 'dir')
    mkdir(savePath);
end


% % Ask whether throat coordinates are in pixels or mm
% choice = questdlg('Are throat coordinates stored in PIXELS or MM?', ...
%                   'Throat Coordinate Units', ...
%                   'Pixels', 'MM', 'Pixels');
% 
% if isempty(choice)
%     error('No throat coordinate unit selected.');
% end

% throatInPixels = strcmpi(choice, 'Pixels');



%% ---------------- PRELOAD DATA ----------------
velData  = cell(1,6);
xData    = cell(1,6);
yData    = cell(1,6);
throatXY = nan(6,2);
caseNames = cell(1,6);

globalMin = inf;
globalMax = -inf;

for i = 1:6
    
    %% ---- Load velocity file ----
    velFileFull = fullfile(velPath, velFiles{i});
    S = load(velFileFull);
    
    % Check velocity variable
    if ~isfield(S, 'velMean_phys')
        error('File "%s" does not contain variable "velMean_phys".', velFiles{i});
    end
    
    vel = S.velMean_phys;
    
    if ~ismatrix(vel)
        error('velMean_phys in file "%s" is not a 2D matrix.', velFiles{i});
    end
    
    % Flip upside down
    vel = flipud(vel);
    
    [nRows, nCols] = size(vel);
    
    % Use x_mm and y_mm if available, else create from calibration
    if isfield(S, 'x_mm') && isfield(S, 'y_mm')
        x_mm = S.x_mm;
        y_mm = S.y_mm;
        
        % Make sure vectors are row/column as expected
        x_mm = x_mm(:).';   % force row
        y_mm = y_mm(:);     % force column
        
        % Since vel was flipped, y-axis also needs flipping
        y_mm = flipud(y_mm);
    else
        x_mm = (0:nCols-1) * mm_per_pixel;
        y_mm = (0:nRows-1)' * mm_per_pixel;
    end
    
    velData{i} = vel;
    xData{i}   = x_mm;
    yData{i}   = y_mm;
    
    globalMin = min(globalMin, min(vel(:), [], 'omitnan'));
    globalMax = max(globalMax, max(vel(:), [], 'omitnan'));
    
    [~, nameOnly, ~] = fileparts(velFiles{i});
    caseNames{i} = strrep(nameOnly, '_', '\_');
    
    %% ---- Load throat file ----
    throatFileFull = fullfile(throatPath, throatFiles{i});
    T = load(throatFileFull);
    
    % Try to extract throat coordinates from common variable names
    [xt, yt] = getThroatCoordinates(T);
    
    if isempty(xt) || isempty(yt)
        warning('Could not find throat coordinates in "%s". Skipping throat overlay for this case.', throatFiles{i});
        throatXY(i,:) = [NaN, NaN];
    else
        % Convert throat location to mm if stored in pixels
        if throatInPixels
            % x conversion
            xt_mm = (xt - 1) * mm_per_pixel;
            
            % y must also be flipped because velMean_phys was flipud-ed
            yt_flipped_pix = nRows - yt + 1;
            yt_mm = (yt_flipped_pix - 1) * mm_per_pixel;
        else
            % already in mm, but y still must be mirrored after flipud
            xt_mm = xt;
            y_min = min(y_mm);
            y_max = max(y_mm);
            yt_mm = y_min + y_max - yt;
        end
        
        throatXY(i,:) = [xt_mm, yt_mm];
    end
end

%% ---------------- PLOT ----------------
figure('Color', 'w', 'Position', [100 80 1400 800]);
tiledlayout(2,3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:6
    nexttile;
    
    vel = velData{i};
    x_mm = xData{i};
    y_mm = yData{i};
    
    imagesc(x_mm, y_mm, vel);
    set(gca, 'YDir', 'normal');
    axis image;
    hold on;
    
    if all(~isnan(throatXY(i,:)))
        plot(throatXY(i,1), throatXY(i,2), 'r+', ...
            'MarkerSize', 12, 'LineWidth', 2);
    end
    
    caxis([globalMin globalMax]); % same scale for all cases
    colormap(turbo);
    colorbar;
    
    xlabel('x [mm]');
    ylabel('y [mm]');
    title(caseNames{i}, 'Interpreter', 'tex');
    box on;
end

sgtitle('Velocity Magnitude of All Cases with Throat Location', ...
        'FontWeight', 'bold', 'FontSize', 14);

%% ---------------- HELPER FUNCTION ----------------
function [xt, yt] = getThroatCoordinates(T)
% Tries to find throat coordinates in a loaded throat MAT structure.
% Modify this if your throat files use different variable names.

    xt = [];
    yt = [];
    
    % Case 1: separate variables
    xCandidates = {'throat_x', 'x_throat', 'xt', 'xT', 'throatX', 'x'};
    yCandidates = {'throat_y', 'y_throat', 'yt', 'yT', 'throatY', 'y'};
    
    fn = fieldnames(T);
    
    xFound = '';
    yFound = '';
    
    for k = 1:numel(xCandidates)
        if ismember(xCandidates{k}, fn)
            xFound = xCandidates{k};
            break;
        end
    end
    
    for k = 1:numel(yCandidates)
        if ismember(yCandidates{k}, fn)
            yFound = yCandidates{k};
            break;
        end
    end
    
    if ~isempty(xFound) && ~isempty(yFound)
        xt = T.(xFound);
        yt = T.(yFound);
        xt = xt(1);
        yt = yt(1);
        return;
    end
    
    % Case 2: one variable containing [x y]
    pairCandidates = {'throat', 'throat_xy', 'throatLoc', 'throat_location', 'loc'};
    
    for k = 1:numel(pairCandidates)
        if ismember(pairCandidates{k}, fn)
            val = T.(pairCandidates{k});
            if isnumeric(val)
                if numel(val) >= 2
                    xt = val(1);
                    yt = val(2);
                    return;
                elseif ismatrix(val) && size(val,2) >= 2
                    xt = val(1,1);
                    yt = val(1,2);
                    return;
                end
            end
        end
    end
end

%% ============================================================
% EXTRACT VERTICAL VELOCITY PROFILES FROM THROAT TO END
% - Profiles at throat x, throat+0.25 mm, throat+0.50 mm, ...
% - For each profile, y starts at throat y and goes to end
% - Saves extracted data to .mat
% - Shows extraction lines on plots
% ============================================================

dx_extract_mm = 0.25;   % spacing between profile locations in mm

profileResults = struct();

figure('Color','w','Position',[100 80 1400 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for i = 1:6
    
    vel = velData{i};
    x_mm = xData{i};
    y_mm = yData{i};
    caseName = caseNames{i};
    
    xt_mm = throatXY(i,1);
    yt_mm = throatXY(i,2);
    
    if isnan(xt_mm) || isnan(yt_mm)
        warning('Skipping case %s because throat location is NaN.', caseName);
        continue;
    end
    
    % ---- Find nearest throat indices ----
    [~, ix_throat] = min(abs(x_mm - xt_mm));
    [~, iy_throat] = min(abs(y_mm - yt_mm));
    
    % ---- Build extraction x-locations from throat to end ----
    x_extract_mm = xt_mm : dx_extract_mm : x_mm(end);
    
    % Snap requested x-locations to nearest grid points
    ix_extract = zeros(size(x_extract_mm));
    x_extract_actual_mm = zeros(size(x_extract_mm));
    
    for k = 1:numel(x_extract_mm)
        [~, ix_extract(k)] = min(abs(x_mm - x_extract_mm(k)));
        x_extract_actual_mm(k) = x_mm(ix_extract(k));
    end

    % Remove duplicates in case grid spacing is coarser than 0.25 mm
    [ix_extract_unique, ia] = unique(ix_extract, 'stable');
    x_extract_actual_mm = x_extract_actual_mm(ia);
    x_extract_req_mm    = x_extract_mm(ia);

    % ---- Extract y-range from throat y to top of image ----
    y_profile_mm = y_mm(1:iy_throat);

    % Velocity profiles matrix:
    % rows = y points
    % cols = x extraction stations
    vel_profiles = vel(1:iy_throat, ix_extract_unique);

    % ---- Store results ----
    profileResults(i).caseName = caseName;
    profileResults(i).velFile = velFiles{i};
    profileResults(i).throatFile = throatFiles{i};
    
    profileResults(i).throat_x_mm = xt_mm;
    profileResults(i).throat_y_mm = yt_mm;
    profileResults(i).throat_ix = ix_throat;
    profileResults(i).throat_iy = iy_throat;
    
    profileResults(i).dx_extract_mm = dx_extract_mm;
    
    profileResults(i).x_extract_requested_mm = x_extract_req_mm;
    profileResults(i).x_extract_actual_mm = x_extract_actual_mm;
    profileResults(i).x_extract_indices = ix_extract_unique;
    
    profileResults(i).y_profile_mm = y_profile_mm;
    profileResults(i).y_indices = 1:iy_throat;
    
    profileResults(i).vel_profiles = vel_profiles;  
    % size = [num_y_points x num_x_stations]
    % each column is one vertical profile
    
    profileResults(i).full_x_mm = x_mm;
    profileResults(i).full_y_mm = y_mm;
    
    % ---- Plot field with extraction lines ----
    nexttile;
    imagesc(x_mm, y_mm, vel);
    set(gca,'YDir','normal');
    axis image;
    hold on;
    
    % Throat marker
    plot(xt_mm, yt_mm, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    
    % Horizontal marker from throat y to right edge
    plot([xt_mm x_mm(end)], [yt_mm yt_mm], 'w--', 'LineWidth', 1.2);
    
    % Vertical extraction lines starting at throat y
    for k = 1:numel(x_extract_actual_mm)
        xline_k = x_extract_actual_mm(k);
        plot([xline_k xline_k], [y_mm(1) yt_mm], 'w-', 'LineWidth', 1);
    end
    
    xlabel('x [mm]');
    ylabel('y [mm]');
    title(sprintf('%s', caseName), 'Interpreter','tex');
    colormap(turbo);
    caxis([globalMin globalMax]);
    colorbar;
    box on;
end

sgtitle('Velocity Maps with Vertical Profile Extraction Lines', ...
    'FontWeight','bold','FontSize',14);

%% ---------------- SAVE EXTRACTED PROFILES ----------------
saveFile = fullfile(savePath, 'vertical_velocity_profiles_from_throat.mat');
save(saveFile, 'profileResults', 'dx_extract_mm', 'mm_per_pixel', '-v7.3');

fprintf('Saved extracted profiles to:\n%s\n', saveFile);

%%
%% ============================================================
% COMPARE TWO PROFILE DATASETS (e.g. base magnification vs 5x)
% - Loads two saved profileResults .mat files
% - Matches cases
% - Matches x stations by nearest physical location
% - Interpolates one profile onto the other's y-grid
% - Plots overlays
% - Computes RMSE / mean absolute difference
% ============================================================

clearvars -except profileResults
clc;

%% ---------------- USER INPUT ----------------
fileA = "E:\December 2025- laminar separation bubble\Vel mag mat files\extracted profiles\vertical_velocity_profiles_from_throat.mat";
fileB = "E:\VelocitydataRAFT\velocity mag 5x\extracted profiles\vertical_velocity_profiles_from_throat.mat";

saveComparePath = 'E:\December 2025- laminar separation bubble\Comparison results\';
if ~exist(saveComparePath, 'dir')
    mkdir(saveComparePath);
end

labelA = '2x mag';
labelB = '5x mag';

% tolerance for deciding whether two x-locations are the "same"
xMatchTol_mm = 0.05;

%% ---------------- LOAD ----------------
A = load(fileA);
B = load(fileB);

if ~isfield(A, 'profileResults')
    error('File A does not contain profileResults.');
end
if ~isfield(B, 'profileResults')
    error('File B does not contain profileResults.');
end

profilesA = A.profileResults;
profilesB = B.profileResults;

nCases = min(numel(profilesA), numel(profilesB));
if nCases ~= 6
    warning('Expected 6 cases, but found %d in common comparison loop.', nCases);
end

comparisonResults = struct();

%% ============================================================
% CASE-BY-CASE COMPARISON
% ============================================================

for i = 1:nCases
    
    caseNameA = profilesA(i).caseName;
    caseNameB = profilesB(i).caseName;
    
    fprintf('\nComparing case %d:\n', i);
    fprintf('  A: %s\n', caseNameA);
    fprintf('  B: %s\n', caseNameB);
    
    % -------- Extract data --------
    xA = profilesA(i).x_extract_actual_mm(:);
    yA = profilesA(i).y_profile_mm(:);
    VA = profilesA(i).vel_profiles;   % rows = y, cols = x stations
    
    xB = profilesB(i).x_extract_actual_mm(:);
    yB = profilesB(i).y_profile_mm(:);
    VB = profilesB(i).vel_profiles;
    
    % Make sure y is increasing for interpolation
    if numel(yA) > 1 && any(diff(yA) < 0)
        yA = flipud(yA);
        VA = flipud(VA);
    end
    
    if numel(yB) > 1 && any(diff(yB) < 0)
        yB = flipud(yB);
        VB = flipud(VB);
    end
    
    % -------- Match x stations --------
    matchedXA = [];
    matchedXB = [];
    idxA_list = [];
    idxB_list = [];
    
    for k = 1:numel(xA)
        [dxMin, idxB] = min(abs(xB - xA(k)));
        if dxMin <= xMatchTol_mm
            matchedXA(end+1,1) = xA(k); %#ok<SAGROW>
            matchedXB(end+1,1) = xB(idxB); %#ok<SAGROW>
            idxA_list(end+1,1) = k; %#ok<SAGROW>
            idxB_list(end+1,1) = idxB; %#ok<SAGROW>
        end
    end
    
    if isempty(idxA_list)
        warning('No matched x stations found for case %s.', caseNameA);
        continue;
    end
    
    % -------- Compare profiles at matched x stations --------
    nMatch = numel(idxA_list);
    
    rmseVals   = nan(nMatch,1);
    madVals    = nan(nMatch,1);
    maxErrVals = nan(nMatch,1);
    
    % store interpolated profiles
    profileCompare = struct([]);
    
    fig = figure('Color','w','Position',[100 80 1400 800]);
    tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    
    nTiles = min(nMatch, 6);  % first 6 matched stations shown
    
    for k = 1:nMatch
        
        ia = idxA_list(k);
        ib = idxB_list(k);
        
        profA = VA(:, ia);
        profB = VB(:, ib);
        
        % common y-range
        yMin = max(min(yA), min(yB));
        yMax = min(max(yA), max(yB));
        
        if yMax <= yMin
            warning('No overlapping y-range for case %s at station %d.', caseNameA, k);
            continue;
        end
        
        % choose finer grid as comparison grid
        if numel(yA) >= numel(yB)
            yCommon = yA(yA >= yMin & yA <= yMax);
        else
            yCommon = yB(yB >= yMin & yB <= yMax);
        end
        
        if numel(yCommon) < 2
            warning('Too few overlapping y points for case %s at station %d.', caseNameA, k);
            continue;
        end
        
        profA_i = interp1(yA, profA, yCommon, 'linear');
        profB_i = interp1(yB, profB, yCommon, 'linear');
        
        diffProf = profA_i - profB_i;
        
        rmseVals(k)   = sqrt(mean(diffProf.^2, 'omitnan'));
        madVals(k)    = mean(abs(diffProf), 'omitnan');
        maxErrVals(k) = max(abs(diffProf), [], 'omitnan');
        
        profileCompare(k).xA_mm = xA(ia);
        profileCompare(k).xB_mm = xB(ib);
        profileCompare(k).yCommon_mm = yCommon;
        profileCompare(k).profileA = profA_i;
        profileCompare(k).profileB = profB_i;
        profileCompare(k).difference = diffProf;
        profileCompare(k).rmse = rmseVals(k);
        profileCompare(k).mad = madVals(k);
        profileCompare(k).maxAbsErr = maxErrVals(k);
        
        % ----- Plot first 6 matched stations -----
        if k <= nTiles
            nexttile;
            hold on;
            plot(profA_i, yCommon, 'LineWidth', 1.8, 'DisplayName', labelA);
            plot(profB_i, yCommon, '--', 'LineWidth', 1.8, 'DisplayName', labelB);
            xlabel('|V| [m/s]');
            ylabel('y [mm]');
            title(sprintf('x = %.2f / %.2f mm', xA(ia), xB(ib)));
            legend('Location','best');
            box on;
            grid on;
        end
    end
    
    sgtitle(sprintf('Profile Comparison: %s', strrep(caseNameA, '\_', '_')), ...
        'FontWeight','bold','FontSize',14);
    
    % save figure
    figName = sprintf('%s_profile_comparison.png', sanitize_filename(caseNameA));
    exportgraphics(fig, fullfile(saveComparePath, figName), 'Resolution', 300);
    
    % -------- Summary plot of error vs x --------
    fig2 = figure('Color','w','Position',[120 120 1000 500]);
    hold on;
    plot(matchedXA, rmseVals, '-o', 'LineWidth', 1.5, 'DisplayName', 'RMSE');
    plot(matchedXA, madVals, '-s', 'LineWidth', 1.5, 'DisplayName', 'Mean abs diff');
    plot(matchedXA, maxErrVals, '-^', 'LineWidth', 1.5, 'DisplayName', 'Max abs diff');
    xlabel('x location [mm]');
    ylabel('Error in |V| [m/s]');
    title(sprintf('Error Metrics vs x: %s', strrep(caseNameA, '\_', '_')));
    legend('Location','best');
    grid on;
    box on;
    
    figName2 = sprintf('%s_error_metrics.png', sanitize_filename(caseNameA));
    exportgraphics(fig2, fullfile(saveComparePath, figName2), 'Resolution', 300);
    
    % -------- Save results --------
    comparisonResults(i).caseNameA = caseNameA;
    comparisonResults(i).caseNameB = caseNameB;
    comparisonResults(i).labelA = labelA;
    comparisonResults(i).labelB = labelB;
    comparisonResults(i).matched_xA_mm = matchedXA;
    comparisonResults(i).matched_xB_mm = matchedXB;
    comparisonResults(i).idxA = idxA_list;
    comparisonResults(i).idxB = idxB_list;
    comparisonResults(i).rmse = rmseVals;
    comparisonResults(i).mad = madVals;
    comparisonResults(i).maxAbsErr = maxErrVals;
    comparisonResults(i).profiles = profileCompare;
end

%% ---------------- SAVE COMPARISON STRUCT ----------------
save(fullfile(saveComparePath, 'comparison_two_magnifications.mat'), ...
    'comparisonResults', 'fileA', 'fileB', 'labelA', 'labelB', 'xMatchTol_mm', '-v7.3');

fprintf('\nSaved comparison results to:\n%s\n', ...
    fullfile(saveComparePath, 'comparison_two_magnifications.mat'));

%% ============================================================
% HELPER FUNCTION
% ============================================================
function out = sanitize_filename(str)
    out = strrep(str, '\_', '_');
    out = regexprep(out, '[^\w\- ]', '');
    out = strrep(out, ' ', '_');
end
