%% ------------------------------------------------------------
% Optical Flow RAFT + ROI (roipoly) + Calibrated Velocity in m/s
% Time-averaged velocity + vertical profiles with axes in mm (origin at lower-left)
%
% OUTPUT ORGANIZATION:
%   - All plots / PNG / GIF saved to:   <filepath>/plots/
%   - All .mat files saved to:         <filepath>/.mat files/
% ------------------------------------------------------------
clear all; clc; close all;

%% ------------------------------------------------------------
% MATLAB PATH SETUP (BATCH-SAFE)
%  - In interactive mode: restores full custom path (pathdef.m)
%  - In batch mode: adds ONLY user code folders (no system toolbox paths)
% ------------------------------------------------------------

userPathFile = fullfile(getenv('HOME'), 'matlab', 'pathdef.m');
userCodeRoot = fullfile(getenv('HOME'), 'matlab');   % your personal code location

if usejava('desktop')
    % ----- INTERACTIVE MATLAB SESSION -----
    if isfile(userPathFile)
        run(userPathFile);
        fprintf('Interactive mode: loaded custom MATLAB path:\n  %s\n', userPathFile);
    else
        warning('Interactive mode: pathdef.m not found. Using default MATLAB path.');
    end
else
    % ----- BATCH / NO-DESKTOP SESSION -----
    if isfolder(userCodeRoot)
        addpath(genpath(userCodeRoot));
        fprintf('Batch mode: added personal MATLAB code root:\n  %s\n', userCodeRoot);
    else
        warning('Batch mode: personal MATLAB code folder not found.');
    end
end


try

    %% --- Specify video path ---
    videoPath = '/home/kbsanjayvasanth/Velocityraft2/smooth_52lpm/Smooth_52lpm.avi';  % <-- Adjust if needed
    [filepath, filename, ~] = fileparts(char(videoPath));
    
    %% --- Calibration parameters ---
    mm_per_pixel = 0.0083333;      % [mm/pixel]
    fps          = 102247;         % [frames per second]
    m_per_pixel  = mm_per_pixel / 1000;  % [m/pixel]
    fprintf('Calibration: %.9f m/pixel | Frame rate: %.1f fps\n', m_per_pixel, fps);

    if ~isfile(videoPath)
        error('Video file not found at: %s', videoPath);
    end

    v = VideoReader(videoPath);
    fprintf('Loaded video: %s\n', filename);

    %% --- Output folders (inside existing save directory) ---
    plotsDir = fullfile(filepath, 'plots');
    matDir   = fullfile(filepath, '.mat files');   % exactly as requested

    if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end
    if ~exist(matDir,   'dir'); mkdir(matDir);   end

    fprintf('Outputs will be saved to:\n  Plots: %s\n  MAT:   %s\n', plotsDir, matDir);


    %% ------------------------------------------------------------
    % Auto-load ROI file (prefer MAT folder, otherwise base folder)
    % ------------------------------------------------------------
    fprintf('Searching for ROI file in:\n  %s\n  %s\n', matDir, filepath);

    roiList = dir(fullfile(matDir, '*_ROI.mat'));
    if isempty(roiList)
        roiList = dir(fullfile(filepath, '*_ROI.mat'));
    end

    if ~isempty(roiList)
        if numel(roiList) > 1
            warning('Multiple ROI files found. Using first: %s', roiList(1).name);
        end

        roiFile = fullfile(roiList(1).folder, roiList(1).name);
        S = load(roiFile);

        if ~isfield(S, 'maskROI')
            error('ROI file "%s" does not contain variable "maskROI".', roiList(1).name);
        end

        maskROI = S.maskROI;
        fprintf('Loaded ROI mask: %s\n', roiFile);

    else
        % ROI DOES NOT EXIST → ask user to draw one
        fprintf('No ROI file found — user must draw a new ROI.\n');

        v.CurrentTime = 0;
        frame0 = im2gray(readFrame(v));

        figROI = figure('Name', 'Select ROI');
        imshow(frame0, []);
        title('Click polygon points (double-click to close)');

        try
            maskROI = roipoly;
            if isempty(maskROI)
                warning('No ROI selected. Using full frame.');
                maskROI = true(size(frame0));
            end

            % Save new ROI into MAT folder
            roiFile = fullfile(matDir, [filename '_ROI.mat']);
            save(roiFile, 'maskROI');
            fprintf('Saved new ROI file: %s\n', roiFile);

        catch MEroi
            warning('ROI selection failed (%s). Using full frame.', MEroi.message);
            maskROI = true(size(frame0));
        end

        close(figROI);
    end

    %% ------------------------------------------------------------
    % Optical Flow Setup (RAFT) with GPU auto-selection
    % ------------------------------------------------------------
    fprintf('Initializing RAFT Optical Flow...\n');

    opticalFlowObj = opticalFlowRAFT;   % create the model once

    % Choose automatic device selection
    executionEnv = "auto";  % "auto" uses GPU if available, CPU otherwise
    accelMode    = "auto";  % "auto" = faster, more memory use

    % NOTE: VideoReader.NumFrames can be unreliable for some codecs; keeping your usage.
    numFrames = v.NumFrames;
    fprintf('Running RAFT Optical Flow on %d frames using %s mode...\n', numFrames, executionEnv);

    %% --- Initialize storage ---
    u_all = [];
    v_all = [];
    velMag_all = [];

    %% --- Rewind and process video ---
    v.CurrentTime = 0;
    framePrev = im2gray(readFrame(v)); %#ok<NASGU>

    tic;  % start timing
    for i = 2:numFrames
        frameCurr = im2gray(readFrame(v));

        % --- Estimate optical flow using RAFT ---
        flow = estimateFlow(opticalFlowObj, frameCurr, ...
            ExecutionEnvironment=executionEnv, ...
            Acceleration=accelMode, ...
            MaxIterations=20, Tolerance=1e-6);

        % --- Apply ROI mask ---
        u  = flow.Vx .* maskROI;
        v_ = flow.Vy .* maskROI;

        % --- Convert to physical velocity (m/s) ---
        u_phys = u  * m_per_pixel * fps;
        v_phys = v_ * m_per_pixel * fps;

        % --- Store instantaneous data ---
        u_all(:,:,i-1)      = u_phys;
        v_all(:,:,i-1)      = v_phys;
        velMag_all(:,:,i-1) = sqrt(u_phys.^2 + v_phys.^2);

        % --- Progress ---
        tElapsed = toc;
        estTotal = tElapsed / (i-1) * numFrames;
        pct = 100 * i / numFrames;
        fprintf('Processed %d/%d (%.1f%%) | Elapsed %.1fs | ETA %.1fs\n', ...
            i, numFrames, pct, tElapsed, max(estTotal - tElapsed, 0));
    end

    %% --- Time-averaged fields (image coords, +y down) ---
    u_mean  = mean(u_all,       3, 'omitnan');
    v_mean  = mean(v_all,       3, 'omitnan');
    velMean = mean(velMag_all,  3, 'omitnan');

    % Mask outside ROI
    u_mean(~maskROI)  = NaN;
    v_mean(~maskROI)  = NaN;
    velMean(~maskROI) = NaN;

    fprintf('Time-averaged fields computed.\n');

    %% ------------------------------------------------------------
    % Convert to physical orientation (origin bottom-left, +Y up):
    % ------------------------------------------------------------
    u_mean_phys   = flipud(u_mean);
    v_mean_phys   = -flipud(v_mean);    % image +y down -> physical +y up
    velMean_phys  = flipud(velMean);
    maskROI_phys  = flipud(maskROI);

    % Instantaneous versions in physical orientation
    u_all_phys   = flip(u_all, 1);
    v_all_phys   = -flip(v_all, 1);
    vel_all_phys = flip(velMag_all, 1);

    %% --- Axes in mm (origin at lower-left) ---
    [H, W] = size(velMean_phys);
    x_mm = (0:W-1) * mm_per_pixel;
    y_mm = (0:H-1) * mm_per_pixel;

    %% ------------------------------------------------------------
    % Auto-load throat location file (prefer MAT folder, otherwise base folder)
    % ------------------------------------------------------------
    fprintf('Searching for throat file in:\n  %s\n  %s\n', matDir, filepath);

    throatList = dir(fullfile(matDir, '*_throat.mat'));
    if isempty(throatList)
        throatList = dir(fullfile(filepath, '*_throat.mat'));
    end

    if isempty(throatList)
        error('No throat file found. Expecting "*_throat.mat" in base folder or MAT folder.');
    end

    if numel(throatList) > 1
        warning('Multiple throat files found. Using first: %s', throatList(1).name);
    end

    throatFile = fullfile(throatList(1).folder, throatList(1).name);
    Sth = load(throatFile);

    if ~isfield(Sth, 'x_throat')
        error('Throat file %s does not contain x_throat.', throatList(1).name);
    end

    x_throat_pixel = double(Sth.x_throat);
    y_throat_pixel = NaN;
    if isfield(Sth, 'y_throat')
        y_throat_pixel = double(Sth.y_throat);
    end

    x_throat_mm = x_throat_pixel * mm_per_pixel;
    y_throat_mm = y_throat_pixel * mm_per_pixel;

    fprintf('Loaded throat: x = %.3f mm (%.1f px), y = %.3f mm (%.1f px)\n', ...
        x_throat_mm, x_throat_pixel, y_throat_mm, y_throat_pixel);

    % --- Find nearest x index to throat location ---
    [~, ix_throat] = min(abs(x_mm - x_throat_mm));
    fprintf('Nearest grid x to throat: %.3f mm (index %d)\n', x_mm(ix_throat), ix_throat);

    %% ------------------------------------------------------------
    % Extract velocity profiles from throat downstream every 0.25 mm
    % ------------------------------------------------------------
    x_spacing = 0.25;  % mm
    x_sample_mm = x_mm(ix_throat):x_spacing:x_mm(end);
    ix_samples = arrayfun(@(x) find(abs(x_mm - x) == min(abs(x_mm - x)), 1), x_sample_mm);

    numFrames_stack = size(vel_all_phys, 3);
    numY = numel(y_mm);
    numX = numel(ix_samples);

    u_profiles   = zeros(numY, numX, numFrames_stack);
    v_profiles   = zeros(numY, numX, numFrames_stack);
    vel_profiles = zeros(numY, numX, numFrames_stack);

    fprintf('Extracting %d x-locations every %.2f mm from throat downstream...\n', numX, x_spacing);

    for ii = 1:numX
        xi = ix_samples(ii);
        u_profiles(:, ii, :)   = squeeze(u_all_phys(:, xi, :));
        v_profiles(:, ii, :)   = squeeze(v_all_phys(:, xi, :));
        vel_profiles(:, ii, :) = squeeze(vel_all_phys(:, xi, :));
    end

    % --- Compute time-averaged profiles ---
    u_profiles_mean   = mean(u_profiles,   3, 'omitnan');
    v_profiles_mean   = mean(v_profiles,   3, 'omitnan');
    vel_profiles_mean = mean(vel_profiles, 3, 'omitnan');

    %% --- Save results (.mat -> matDir) ---
    profileFile = fullfile(matDir, [filename '_ThroatProfiles.mat']);
    save(profileFile, ...
        'x_sample_mm', 'y_mm', ...
        'u_profiles', 'v_profiles', 'vel_profiles', ...
        'u_profiles_mean', 'v_profiles_mean', 'vel_profiles_mean', ...
        'x_throat_mm', 'y_throat_mm', ...
        '-v7.3');
    fprintf('Saved throat and downstream profiles to: %s\n', profileFile);

    avgVelFile = fullfile(matDir, [filename '_TimeAvgVelField.mat']);
    save(avgVelFile, 'velMean_phys', 'x_mm', 'y_mm', '-v7.3');
    fprintf('Saved time-averaged velocity field to: %s\n', avgVelFile);

    %% --- Plot: Time-Averaged Velocity Magnitude (m/s) ---
    fig1 = figure('Name', 'Time-Averaged Velocity Magnitude (ROI, mm axes, physical orientation)');
    imagesc(x_mm, y_mm, velMean_phys);
    set(gca, 'YDir', 'normal');
    axis image; colorbar;
    title('Time-Averaged Velocity Magnitude (m/s)');
    xlabel('X (mm)'); ylabel('Y (mm)');

    saveas(fig1, fullfile(plotsDir, [filename '_TimeAvgVelMag.png']));

    %% ------------------------------------------------------------
    % Time-Averaged Magnitude Plot with Throat Line (publication-style)
    % ------------------------------------------------------------
    fprintf('Creating research-grade time-averaged magnitude plot...\n');

    figT = figure('Name','Time-Averaged Velocity (Research Grade)',...
        'Color','w','Position',[200 200 900 700]);

    imagesc(x_mm, y_mm, velMean_phys);
    set(gca, 'YDir','normal');
    axis image;

    colormap(turbo);
    cb = colorbar;
    cb.Label.String = 'Velocity Magnitude (m/s)';
    cb.Label.FontSize = 12;
    cb.LineWidth = 1.2;

    hold on;
    plot([x_throat_mm x_throat_mm], [y_mm(1) y_mm(end)], 'w:', 'LineWidth', 1.25);

    title('Time-Averaged Velocity Magnitude', 'FontSize', 18, 'Interpreter','latex');
    xlabel('X (mm)', 'FontSize', 16, 'Interpreter','latex');
    ylabel('Y (mm)', 'FontSize', 16, 'Interpreter','latex');

    set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'Box','on', 'FontName','Times');

    timeAvgPlotFile = fullfile(plotsDir, [filename '_TimeAvgVelMag_ThroatLine.png']);
    exportgraphics(figT, timeAvgPlotFile, 'Resolution', 350);

    fprintf('Saved time-averaged magnitude plot with throat line: %s\n', timeAvgPlotFile);

    %% ------------------------------------------------------------
    % Vertical Velocity Profiles at User-Specified Stations
    % ------------------------------------------------------------
    fprintf('\n--- Plotting vertical velocity profiles ---\n');

    %numStations = input('Enter number of X-stations to plot vertical profiles (e.g., 5): ');
    numStations = 12;

    if isempty(numStations) || numStations < 1
        numStations = 5;
        fprintf('Using default of %d stations.\n', numStations);
    end

    x_spacing = 0.25; % mm
    x_sample_mm = x_mm(ix_throat):x_spacing:x_mm(end);

    station_idx = round(linspace(1, numel(x_sample_mm), numStations));
    stations_mm  = x_sample_mm(station_idx);

    profiles_to_plot = vel_profiles_mean(:, station_idx);

    figP = figure('Name','Vertical Profiles','Color','w',...
        'Position',[100 100 900 700]);

    tiledlayout(1,1);
    ax = nexttile;
    hold(ax, 'on');

    colors = lines(numStations);
    for s = 1:numStations
        plot(ax, profiles_to_plot(:, s), y_mm, ...
            'LineWidth', 2.2, ...
            'Color', colors(s,:), ...
            'DisplayName', sprintf('x = %.2f mm', stations_mm(s)));
    end

    set(ax, 'YDir','normal', 'FontSize', 14, 'LineWidth', 1.2);
    xlabel(ax, 'Velocity Magnitude (m/s)', 'FontSize', 16, 'Interpreter','latex');
    ylabel(ax, 'Y (mm)', 'FontSize', 16, 'Interpreter','latex');
    title(ax, 'Vertical Velocity Profiles at Selected Streamwise Stations',...
        'FontSize', 18, 'Interpreter','latex');

    grid(ax, 'on');
    legend(ax, 'Location','best', 'Box','off');

    ax.Box = 'on';
    ax.LineWidth = 1.3;
    ax.FontName = 'Times';

    profilePlotFile = fullfile(plotsDir, [filename '_velprofiles.png']);
    exportgraphics(figP, profilePlotFile, 'Resolution', 350);
    fprintf('Saved velocity profile plot: %s\n', profilePlotFile);

    %% ------------------------------------------------------------
    % Section 1: Display Quiver Animation (No Saving)
    % ------------------------------------------------------------
    showAnimation = false;   % <-- set to true to display frames

    if showAnimation
        fprintf('Displaying quiver animation...\n');

        skip = 25;
        [X, Y] = meshgrid(x_mm(1:skip:end), y_mm(1:skip:end));

        startFrame = 2;
        endFrame   = 30;
        frameRange = startFrame:endFrame;
        fprintf('Animating frames %d to %d...\n', startFrame, endFrame);

        figQ = figure('Name', 'Velocity Field Animation (Physical Orientation)');
        axis equal;
        set(gca, 'YDir', 'normal');
        xlabel('X (mm)'); ylabel('Y (mm)');
        title('Instantaneous Velocity Field (m/s)');
        colormap turbo;

        for k = frameRange
            if ~ishghandle(figQ)
                fprintf('Animation interrupted by user (figure closed).\n');
                break;
            end

            u_plot = u_all_phys(1:skip:end, 1:skip:end, k);
            v_plot = v_all_phys(1:skip:end, 1:skip:end, k);

            imagesc(x_mm, y_mm, vel_all_phys(:,:,k));
            hold on;
            set(gca, 'YDir', 'normal');
            axis image;
            caxis([0, max(velMean_phys(:), [], 'omitnan')]);
            colorbar;

            q = quiver(X, Y, u_plot, v_plot, 'k');
            q.AutoScale = 'on';
            q.AutoScaleFactor = 0.8;

            text(x_mm(end)*0.9, y_mm(end)*0.05, sprintf('Frame: %d', k), ...
                'Color','r','FontWeight','bold','FontSize',12, ...
                'HorizontalAlignment','right');

            title(sprintf('Instantaneous Velocity Field (Frame %d)', k));
            drawnow;
            hold off;
        end

        fprintf('Display-only animation complete.\n\n');
    end

    %% ------------------------------------------------------------
    % Section 2: Save Quiver Animation as GIF (No Display Needed)
    % ------------------------------------------------------------
    saveGIF = true;

    if saveGIF
        fprintf('Saving HIGH-QUALITY quiver animation as GIF...\n');

        skip = 15;
        [X, Y] = meshgrid(x_mm(1:skip:end), y_mm(1:skip:end));

        startFrame = 200;
        endFrame   = 300;
        frameRange = startFrame:endFrame;

        gifFile = fullfile(plotsDir, ...
            [filename sprintf('_VelField_%dto%d_highQ.gif', startFrame, endFrame)]);

        delayTime = 0.4;

        figG = figure('Visible','off', ...
            'Position',[100 100 1400 1050], ...
            'Renderer','painters');

        colormap turbo;

        firstFrame = true;

        for k = frameRange
            imagesc(x_mm, y_mm, vel_all_phys(:,:,k));
            axis image;
            set(gca,'YDir','normal');
            caxis([0 max(velMean_phys(:), [], 'omitnan')]);
            colorbar;
            hold on;

            u_plot = u_all_phys(1:skip:end, 1:skip:end, k);
            v_plot = v_all_phys(1:skip:end, 1:skip:end, k);
            q = quiver(X, Y, u_plot, v_plot, 'k', 'LineWidth', 1);
            q.AutoScale = 'on';
            q.AutoScaleFactor = 0.8;

            text(x_mm(end)*0.9, y_mm(end)*0.05, sprintf('Frame: %d', k), ...
                'Color','r','FontWeight','bold','FontSize',16, ...
                'HorizontalAlignment','right');

            drawnow;

            frame = getframe(figG);
            RGB = frame2im(frame);

            scaleFactor = 1.5;
            bigRGB = imresize(RGB, scaleFactor, 'bicubic');

            [imind, cm] = rgb2ind(bigRGB, 256, 'nodither');

            if firstFrame
                imwrite(imind, cm, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
                firstFrame = false;
            else
                imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
            end
            hold off;
        end

        close(figG);
        fprintf('High-quality GIF saved:\n%s\n\n', gifFile);
    end

    %% ------------------------------------------------------------
    % Save velocity components and calibration (.mat -> matDir)
    % ------------------------------------------------------------
    fprintf('Saving velocity components and calibration...\n');
    uvFile = fullfile(matDir, [filename '_velocity.mat']);
    save(uvFile, ...
        'u_all', 'v_all', 'velMag_all', ...                 % image-coord instantaneous (m/s)
        'u_mean', 'v_mean', 'velMean', ...                  % image-coord mean (m/s)
        'u_mean_phys', 'v_mean_phys', 'velMean_phys', ...   % physical-orientation mean (m/s)
        'x_mm', 'y_mm', 'mm_per_pixel', 'm_per_pixel', 'fps', ...
        'maskROI', 'maskROI_phys', '-v7.3');
    fprintf('Saved: %s\n', uvFile);

    %% ------------------------------------------------------------
    % Compute and Save Average Background Image (plot file -> plotsDir)
    % ------------------------------------------------------------
    saveBackground = false;   % <-- set to true to compute/save background
    numBgFrames = min(500, v.NumFrames);

    if saveBackground
        fprintf('Computing average background image using first %d frames...\n', numBgFrames);

        % Rewind video to beginning
        v.CurrentTime = 0;

        frame0 = im2gray(readFrame(v));
        [Hb, Wb] = size(frame0);
        bgSum = zeros(Hb, Wb, 'double');

        v.CurrentTime = 0;
        for k = 1:numBgFrames
            frame = im2gray(readFrame(v));
            bgSum = bgSum + double(frame);

            if mod(k, 20) == 0
                fprintf('  Averaging progress: %d / %d frames\n', k, numBgFrames);
            end
        end

        avgBG = uint8(bgSum / numBgFrames);
        fprintf('Average background computed.\n');

        bgFile = fullfile(plotsDir, [filename '_avgBG.png']);
        imwrite(avgBG, bgFile);
        fprintf('Saved average background image to: %s\n', bgFile);

        figure('Name', 'Average Background Image');
        imshow(avgBG, []);
        title(sprintf('Average Background (%d frames)', numBgFrames));

        v.CurrentTime = 0;
    else
        fprintf('Skipping background averaging section.\n');
    end

    fprintf('--- Analysis complete ---\n');

catch ME
    % --- Cluster-safe error handling ---
    try
        homeDir = getenv('HOME');
        if isempty(homeDir)
            homeDir = '/tmp';
        end
        errFile = [homeDir '/RAFT_errorLog.txt'];

        fid = fopen(errFile, 'w');
        if fid ~= -1
            fprintf(fid, 'Error occurred on %s\n\n', datestr(now));
            fprintf(fid, 'Error message:\n%s\n\n', ME.message);
            fprintf(fid, 'Stack trace:\n');
            for k = 1:length(ME.stack)
                fprintf(fid, '  File: %s (line %d)\n', ME.stack(k).file, ME.stack(k).line);
            end
            fclose(fid);
            fprintf(2, 'An error occurred. Details saved to: %s\n', errFile);
        else
            fprintf(2, 'Failed to open error log file. Check permissions.\n');
        end
    catch
        fprintf(2, 'Error logging failed.\n');
    end
end
