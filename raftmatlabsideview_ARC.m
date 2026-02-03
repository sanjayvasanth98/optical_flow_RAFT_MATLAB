%% ------------------------------------------------------------
% Optical Flow RAFT
% Author: Sanjay Vasanth
% Last modified: 2/2/2026
% Time-averaged velocity + vertical profiles with axes in mm (origin at lower-left)
% Incremental saving to avoid huge end-of-run save time / RAM blowups
% Saves:
%   - MAT files into:   <video_folder>/mat files/
%   - Plots into:       <video_folder>/plots/
% ------------------------------------------------------------
clear all; clc; close all;

%% --------------------- USER TOGGLES --------------------------
showFigures    = false;   % show figures on screen
savePlots      = true;    % save plots
skip_animation = true;    % keep GIF optional off by default

% RAFT controls
raftIters     = 8;
raftTolerance = 1e-6;
accelMode     = "auto";

% Execution environment
if canUseGPU
    executionEnv = "gpu";
else
    executionEnv = "cpu";
end

%% --- Calibration parameters ---
mm_per_pixel = 0.0037852979;         % [mm/pixel]
fps          = 100000;               % [frames per second]
m_per_pixel  = mm_per_pixel / 1000;  % [m/pixel]
fprintf('Calibration: %.9f m/pixel | Frame rate: %.1f fps\n', m_per_pixel, fps);

%% --- Load your personal toolbox path safely ---
userPathFile = fullfile(getenv('HOME'), 'matlab', 'pathdef.m');
if isfile(userPathFile)
    addpath(genpath(fileparts(userPathFile)));
    run(userPathFile);
    fprintf('Loaded custom MATLAB path: %s\n', userPathFile);
else
    warning('Custom pathdef.m not found. Using default MATLAB path.');
end

try
    %% --- Specify video path ---
    videoPath = '/home/kbsanjayvasanth/VelocitydataRAFT/Nocavitation_P10S100/P10S100nocavitation_sideview.avi';
    [filepath, filename, ~] = fileparts(char(videoPath));

    if ~isfile(videoPath)
        error('Video file not found at: %s', videoPath);
    end

    v = VideoReader(videoPath);
    fprintf('Loaded video: %s\n', filename);

    %% ------------------------------------------------------------
    % Create output folders
    %% ------------------------------------------------------------
    plotsDir = fullfile(filepath, 'plots');
    if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end

    matDir = fullfile(filepath, 'mat files');
    if ~exist(matDir, 'dir'); mkdir(matDir); end

    %% ------------------------------------------------------------
    % Load existing ROI mask if available, else select and save a new one
    % NOTE: ROI is saved inside "mat files"
    %% ------------------------------------------------------------
    fprintf('Checking for existing ROI...\n');
    roiFile = fullfile(matDir, [filename '_ROI.mat']);

    if isfile(roiFile)
        S = load(roiFile, 'maskROI');
        maskROI = S.maskROI;
        fprintf('Loaded existing ROI mask from: %s\n', roiFile);
    else
        v.CurrentTime = 0;
        frame0 = im2gray(readFrame(v));

        figROI = figure('Name', 'Select ROI');
        imshow(frame0, []);
        title('Click polygon points (double-click to close and confirm ROI)');
        try
            maskROI = roipoly;
            if isempty(maskROI)
                warning('No ROI selected. Using full frame.');
                maskROI = true(size(frame0));
            end
            save(roiFile, 'maskROI');
            fprintf('Saved ROI mask to: %s\n', roiFile);
        catch MEroi
            warning('ROI selection failed (%s). Using full frame.', MEroi.message);
            maskROI = true(size(frame0));
        end
        try close(figROI); end
    end

    %% ------------------------------------------------------------
    % Optical Flow Setup (RAFT)
    %% ------------------------------------------------------------
    fprintf('Initializing RAFT Optical Flow...\n');
    opticalFlowObj = opticalFlowRAFT;

    % Robust frame count estimate
    numFrames     = max(2, floor(v.Duration * v.FrameRate));
    numFlowFrames = numFrames - 1;

    fprintf('Running RAFT Optical Flow on %d frames using %s mode...\n', numFrames, executionEnv);

    %% ------------------------------------------------------------
    % Rewind and read first frame for sizing
    %% ------------------------------------------------------------
    v.CurrentTime = 0;
    framePrev = im2gray(readFrame(v));
    [H, W] = size(framePrev);

    %% ------------------------------------------------------------
    % STREAMING SAVE SETUP (matfile)
    % Saves the BIG file into "mat files"
    %% ------------------------------------------------------------
    fprintf('Preparing streaming MAT-file for incremental saves...\n');
    uvFile = fullfile(matDir, [filename '_velocity.mat']);

    M = matfile(uvFile, 'Writable', true);

    % Preallocate on disk (NOT in RAM) — SINGLE precision
    M.u_all = zeros(H, W, numFlowFrames, 'single');
    M.v_all = zeros(H, W, numFlowFrames, 'single');

    % Save constants once
    M.mm_per_pixel = mm_per_pixel;
    M.m_per_pixel  = m_per_pixel;
    M.fps          = fps;
    M.maskROI      = maskROI;

    fprintf('Streaming file created: %s\n', uvFile);

    %% ------------------------------------------------------------
    % Running sums for means (memory safe)
    %% ------------------------------------------------------------
    sumU   = zeros(H, W, 'double');
    sumV   = zeros(H, W, 'double');
    sumMag = zeros(H, W, 'double');
    count  = 0;

    %% ------------------------------------------------------------
    % Process frames (write each instantaneous frame to disk)
    %% ------------------------------------------------------------
    tic;
    for i = 2:numFrames
        frameCurr = im2gray(readFrame(v));

        flow = estimateFlow(opticalFlowObj, frameCurr, ...
            ExecutionEnvironment=executionEnv, ...
            Acceleration=accelMode, ...
            MaxIterations=raftIters, ...
            Tolerance=raftTolerance);

        % ROI mask (image coords)
        u = flow.Vx .* maskROI;
        v_ = flow.Vy .* maskROI;

        % Convert to physical velocity (m/s)
        u_phys = u * m_per_pixel * fps;
        v_phys = v_ * m_per_pixel * fps;

        % --- Write instantaneous to disk (single) ---
        k = i - 1;
        M.u_all(:,:,k) = single(u_phys);
        M.v_all(:,:,k) = single(v_phys);

        % --- Update running sums for mean (double for accuracy) ---
        sumU   = sumU   + double(u_phys);
        sumV   = sumV   + double(v_phys);
        sumMag = sumMag + hypot(double(u_phys), double(v_phys));
        count  = count + 1;

        % progress
        tElapsed = toc;
        estTotal = tElapsed / max(k,1) * numFlowFrames;
        pct = 100 * i / numFrames;
        fprintf('Processed %d/%d (%.1f%%) | Elapsed %.1fs | ETA %.1fs\n', ...
            i, numFrames, pct, tElapsed, max(estTotal - tElapsed, 0));

        framePrev = frameCurr;
    end
    toc;

    %% ------------------------------------------------------------
    % Time-averaged fields (computed from running sums)
    %% ------------------------------------------------------------
    u_mean  = single(sumU   / max(count,1));
    v_mean  = single(sumV   / max(count,1));
    velMean = single(sumMag / max(count,1));

    % Mask outside ROI
    u_mean(~maskROI)  = NaN;
    v_mean(~maskROI)  = NaN;
    velMean(~maskROI) = NaN;

    fprintf('Time-averaged fields computed.\n');

    %% ------------------------------------------------------------
    % Physical orientation versions (origin bottom-left, +Y up)
    %% ------------------------------------------------------------
    u_mean_phys  = flipud(u_mean);
    v_mean_phys  = -flipud(v_mean);
    velMean_phys = flipud(velMean);
    maskROI_phys = flipud(maskROI);

    % Axes in mm
    x_mm = (0:W-1) * mm_per_pixel;
    y_mm = (0:H-1) * mm_per_pixel;

    % Save mean fields + axes into SAME big MAT-file
    M.u_mean       = u_mean;
    M.v_mean       = v_mean;
    M.velMean      = velMean;
    M.u_mean_phys  = u_mean_phys;
    M.v_mean_phys  = v_mean_phys;
    M.velMean_phys = velMean_phys;
    M.maskROI_phys = maskROI_phys;
    M.x_mm         = x_mm;
    M.y_mm         = y_mm;

    %% ------------------------------------------------------------
    % Throat profiles (MEAN only, memory-safe)
    %% ------------------------------------------------------------
    fprintf('Extracting throat velocity profiles (mean only, memory-safe)...\n');

    throatFile = fullfile(filepath, [filename '_throat.mat']);
    if ~isfile(throatFile)
        error('Throat file not found: %s\nExpected variable x_throat (px).', throatFile);
    end
    Sth = load(throatFile);

    if isfield(Sth, 'x_throat')
        x_throat_pixel = double(Sth.x_throat);
        x_throat_mm = x_throat_pixel * mm_per_pixel;
    else
        error('Throat file must contain variable "x_throat" (in pixels).');
    end

    if isfield(Sth, 'y_throat')
        y_throat_pixel = double(Sth.y_throat);
        y_throat_mm = y_throat_pixel * mm_per_pixel;
    else
        y_throat_pixel = NaN; y_throat_mm = NaN;
    end

    fprintf('Loaded throat position: x = %.3f mm (%.1f px), y = %.3f mm (%.1f px)\n', ...
        x_throat_mm, x_throat_pixel, y_throat_mm, y_throat_pixel);

    [~, ix_throat] = min(abs(x_mm - x_throat_mm));
    fprintf('Nearest grid x to throat: %.3f mm (index %d)\n', x_mm(ix_throat), ix_throat);

    x_spacing = 0.25;  % mm
    x_sample_mm = x_mm(ix_throat):x_spacing:x_mm(end);
    ix_samples = arrayfun(@(x) find(abs(x_mm - x) == min(abs(x_mm - x)), 1), x_sample_mm);

    numX = numel(ix_samples);
    numY = numel(y_mm);

    sumUp   = zeros(numY, numX, 'double');
    sumVp   = zeros(numY, numX, 'double');
    sumMagp = zeros(numY, numX, 'double');

    for k = 1:numFlowFrames
        u_k = M.u_all(:,:,k);
        v_k = M.v_all(:,:,k);

        u_phys_k = flipud(u_k);
        v_phys_k = -flipud(v_k);

        for ii = 1:numX
            xi = ix_samples(ii);
            uk_col = double(u_phys_k(:, xi));
            vk_col = double(v_phys_k(:, xi));

            sumUp(:, ii)   = sumUp(:, ii)   + uk_col;
            sumVp(:, ii)   = sumVp(:, ii)   + vk_col;
            sumMagp(:, ii) = sumMagp(:, ii) + hypot(uk_col, vk_col);
        end
    end

    u_profiles_mean   = single(sumUp   / numFlowFrames);
    v_profiles_mean   = single(sumVp   / numFlowFrames);
    vel_profiles_mean = single(sumMagp / numFlowFrames);

    % Save mean profiles into "mat files"
    profileFile = fullfile(matDir, [filename '_ThroatProfiles_mean.mat']);
    save(profileFile, ...
        'x_sample_mm', 'y_mm', ...
        'u_profiles_mean', 'v_profiles_mean', 'vel_profiles_mean', ...
        'x_throat_mm', 'y_throat_mm', ...
        '-v7.3');
    fprintf('Saved throat/downstream MEAN profiles to: %s\n', profileFile);

    %% ------------------------------------------------------------
    % Save mean velocity magnitude field separately (small) into "mat files"
    %% ------------------------------------------------------------
    avgVelFile = fullfile(matDir, [filename '_TimeAvgVelField.mat']);
    save(avgVelFile, 'velMean_phys', 'x_mm', 'y_mm', '-v7.3');
    fprintf('Saved time-averaged velocity field to: %s\n', avgVelFile);

    %% ------------------------------------------------------------
    % Plot A: magnitude only + throat line (normal colorbar + no overlap)
    %% ------------------------------------------------------------
    if savePlots
        figA = figure('Visible', tern(showFigures,'on','off'), ...
            'Color','w','Position',[100 100 1200 800]);
        ax = axes(figA);

        imagesc(ax, x_mm, y_mm, velMean_phys);
        set(ax,'YDir','normal');
        axis(ax,'image');
        colormap(ax, turbo);

        cb = colorbar(ax, 'eastoutside');
        cb.TickDirection = 'out';
        cb.Label.String  = 'Velocity Magnitude (m/s)';
        cb.Label.FontSize = 14;
        cb.Label.Rotation = 90;

        cb.Label.Units = 'normalized';
        cb.Label.Position = [3 0.5 0];   % push label right

        ax.Position = [0.08 0.10 0.78 0.82];

        hold(ax,'on');
        plot(ax, [x_throat_mm x_throat_mm], [y_mm(1) y_mm(end)], 'w:', 'LineWidth', 1.7);

        title(ax, 'Time-Averaged Velocity Magnitude', 'FontSize', 20, 'Interpreter','latex');
        xlabel(ax, 'X (mm)', 'FontSize', 18, 'Interpreter','latex');
        ylabel(ax, 'Y (mm)', 'FontSize', 18, 'Interpreter','latex');
        set(ax,'FontSize',15,'LineWidth',1.2,'Box','on','FontName','Times');

        outA = fullfile(plotsDir, [filename '_TimeAvgVelMag_ThroatLine.png']);
        exportgraphics(figA, outA, 'Resolution', 350);
        close(figA);
        fprintf('Saved plot A: %s\n', outA);
    end

    %% ------------------------------------------------------------
    % Plot B: Time-averaged |V| + streamlines ONLY inside ROI
    %% ------------------------------------------------------------
    if savePlots
        fprintf('Saving Plot B (Streamlines strictly clipped to ROI)...\n');

        if ~exist('maskROI','var') || isempty(maskROI)
            error('maskROI not found in workspace. Load ROI before Plot B.');
        end

        maskROI_phys = flipud(logical(maskROI));
        maskOut_phys = ~maskROI_phys;

        [Hf, Wf] = size(velMean_phys);
        if any(size(maskROI_phys) ~= [Hf Wf])
            error('maskROI_phys size (%dx%d) does not match velMean_phys (%dx%d).', ...
                size(maskROI_phys,1), size(maskROI_phys,2), Hf, Wf);
        end

        figB = figure('Visible', tern(showFigures,'on','off'), ...
            'Color','w','Position',[100 100 1200 800]);
        ax = axes(figB);

        uPlot = u_mean_phys;
        vPlot = v_mean_phys;
        qPlot = velMean_phys;

        uPlot(maskOut_phys) = NaN;
        vPlot(maskOut_phys) = NaN;
        qPlot(maskOut_phys) = NaN;

        imagesc(ax, x_mm, y_mm, qPlot);
        set(ax,'YDir','normal');
        axis(ax,'image');
        colormap(ax, turbo);

        cb = colorbar(ax);
        cb.TickDirection  = 'out';
        cb.Label.String   = 'Velocity Magnitude (m/s)';
        cb.Label.FontSize = 13;
        cb.Label.Rotation = 90;
        cb.Label.Units    = 'normalized';
        cb.Label.Position = [3 0.5 0];

        hold(ax,'on');

        [XX, YY] = meshgrid(x_mm, y_mm);
        density = 0.8;
        hss = streamslice(XX, YY, uPlot, vPlot, density);
        set(hss, 'Color','k', 'LineWidth', 0.7);

        % STRICT clip streamline vertices to ROI
        for kk = 1:numel(hss)
            xl = hss(kk).XData;
            yl = hss(kk).YData;
            if isempty(xl) || isempty(yl), continue; end

            in_roi = interp2(XX, YY, single(maskROI_phys), xl, yl, 'nearest', 0);
            outside_idx = (in_roi < 0.5);
            xl(outside_idx) = NaN;
            yl(outside_idx) = NaN;
            hss(kk).XData = xl;
            hss(kk).YData = yl;
        end

        % Opaque outside-ROI overlay (visual blocker)
        navyRGB = zeros(size(maskOut_phys,1), size(maskOut_phys,2), 3, 'single');
        navyRGB(:,:,1) = 0.02;
        navyRGB(:,:,2) = 0.05;
        navyRGB(:,:,3) = 0.25;

        hMask = image(ax, x_mm, y_mm, navyRGB);
        set(hMask, 'AlphaData', 1.0 * single(maskOut_phys));
        uistack(hMask, 'top');

        title(ax, 'Time-Averaged Velocity Magnitude with Streamlines (ROI only)', ...
            'FontSize', 18, 'Interpreter','latex');
        xlabel(ax, 'X (mm)', 'FontSize', 16, 'Interpreter','latex');
        ylabel(ax, 'Y (mm)', 'FontSize', 16, 'Interpreter','latex');
        set(ax,'FontSize',14,'LineWidth',1.2,'Box','on','FontName','Times');

        outB = fullfile(plotsDir, [filename '_TimeAvgVelMag_Streamlines_ROI.png']);
        exportgraphics(figB, outB, 'Resolution', 350);
        close(figB);
        fprintf('Saved Plot B: %s\n', outB);
    end

    %% ------------------------------------------------------------
    % OPTIONAL: GIF saving (reads slices from disk)
    %% ------------------------------------------------------------
    if ~skip_animation
        saveGIF = true;
        if saveGIF
            fprintf('Saving quiver animation as GIF (reading slices from MAT-file)...\n');

            skip = 25;
            [Xq, Yq] = meshgrid(x_mm(1:skip:end), y_mm(1:skip:end));
            startFrame = 200;
            endFrame   = min(300, numFlowFrames);
            frameRange = startFrame:endFrame;

            gifFile = fullfile(plotsDir, [filename sprintf('_VelField_%dto%d.gif', startFrame, endFrame)]);
            delayTime = 0.5;

            figG = figure('Visible', 'off');
            axis equal;
            set(gca, 'YDir', 'normal');
            xlabel('X (mm)'); ylabel('Y (mm)');
            title('Instantaneous Velocity Field (m/s)');
            colormap turbo;

            climMax = max(velMean_phys(:), [], 'omitnan');

            for k = frameRange
                u_k = M.u_all(:,:,k);
                v_k = M.v_all(:,:,k);

                u_phys_k = flipud(u_k);
                v_phys_k = -flipud(v_k);
                mag_k = hypot(u_phys_k, v_phys_k);

                imagesc(x_mm, y_mm, mag_k);
                hold on;
                set(gca, 'YDir', 'normal');
                axis image;
                caxis([0, climMax]);
                colorbar;

                u_plot = u_phys_k(1:skip:end, 1:skip:end);
                v_plot = v_phys_k(1:skip:end, 1:skip:end);
                q = quiver(Xq, Yq, u_plot, v_plot, 'k');
                q.AutoScale = 'on';
                q.AutoScaleFactor = 0.8;

                text(x_mm(end)*0.9, y_mm(end)*0.05, sprintf('Frame: %d', k), ...
                    'Color','r','FontWeight','bold','FontSize',12, ...
                    'HorizontalAlignment','right');

                drawnow;

                frame = getframe(figG);
                [imind, cm] = rgb2ind(frame2im(frame), 256);
                if k == frameRange(1)
                    imwrite(imind, cm, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
                else
                    imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
                end
                hold off;
            end
            close(figG);
            fprintf('GIF saved: %s\n', gifFile);
        end
    end

    %% ------------------------------------------------------------
    % Final summary print (counts)
    %% ------------------------------------------------------------
    fprintf('\nDONE.\n');
    fprintf('MAT outputs saved in: %s\n', matDir);
    fprintf('  1) %s\n', uvFile);
    fprintf('  2) %s\n', profileFile);
    fprintf('  3) %s\n', avgVelFile);
    fprintf('  (ROI file also saved: %s)\n', roiFile);

    fprintf('Plots saved in: %s\n', plotsDir);
    fprintf('  A) %s\n', fullfile(plotsDir, [filename '_TimeAvgVelMag_ThroatLine.png']));
    fprintf('  B) %s\n', fullfile(plotsDir, [filename '_TimeAvgVelMag_Streamlines_ROI.png']));

    fprintf('\nCounts:\n');
    fprintf('  MAT files (excluding ROI): 3\n');
    fprintf('  MAT files (including ROI): 4\n');
    fprintf('  Plots saved: 2\n');

catch ME
    % --- Cluster-safe error handling ---
    try
        homeDir = getenv('HOME');
        if isempty(homeDir), homeDir = '/tmp'; end
        errFile = fullfile(homeDir, 'RAFT_errorLog.txt');

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

%% --------------------- helper ---------------------
function out = tern(cond, a, b)
% ternary helper: out = a if cond else b
if cond
    out = a;
else
    out = b;
end
end
