%% ------------------------------------------------------------
% Optical Flow RAFT
% Author: Sanjay Vasanth
% Last modified: 2/3/2026
% Time-averaged velocity + vertical profiles with axes in mm (origin at lower-left)
% Incremental saving to avoid huge end-of-run save time / RAM blowups
% Saves:
%   - MAT files into:   <video_folder>/mat files/
%   - Plots into:       <video_folder>/plots/
%% ------------------------------------------------------------
clear all; clc; close all;

%% --------------------- USER TOGGLES --------------------------
showFigures    = false;   % MUST be false on cluster
savePlots      = true;    % saves Plot A and Plot B
skip_animation = true;    % if true, GIF block is skipped

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
mm_per_pixel = 0.04917372;         % [mm/pixel]   %<--edit
fps          = 130000;               % [frames per second]  %<--edit
m_per_pixel  = mm_per_pixel / 1000;  % [m/pixel]  
fprintf('Calibration: %.9f m/pixel | Frame rate: %.1f fps\n', m_per_pixel, fps);

%% --- Load your personal toolbox path safely ---
fprintf('\nSTEP 0/7: Loading custom MATLAB path (if available)...\n');
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
    fprintf('\nSTEP 1/7: Opening video...\n');
    videoPath = '/home/kbsanjayvasanth/Inception_raft_test/smooth_naga/smooth30.avi';   %<--edit
    [filepath, filename, ~] = fileparts(char(videoPath));

    if ~isfile(videoPath)
        error('Video file not found at: %s', videoPath);
    end

    v = VideoReader(videoPath);
    fprintf('Loaded video: %s\n', filename);

    %% ------------------------------------------------------------
    % Create output folders
    %% ------------------------------------------------------------
    fprintf('\nSTEP 2/7: Creating output folders...\n');
    plotsDir = fullfile(filepath, 'plots');
    if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end

    matDir = fullfile(filepath, 'mat files');
    if ~exist(matDir, 'dir'); mkdir(matDir); end

    %% ------------------------------------------------------------
    % ROI (CLUSTER SAFE): MUST already exist
    %% ------------------------------------------------------------
    fprintf('\nSTEP 3/7: Loading ROI (non-interactive)...\n');
    roiFile = fullfile(filepath, [filename '_ROI.mat']);

    if isfile(roiFile)
        S = load(roiFile, 'maskROI');
        if ~isfield(S,'maskROI') || isempty(S.maskROI)
            error('ROI file exists but maskROI is missing/empty: %s', roiFile);
        end
        maskROI = logical(S.maskROI);
        fprintf('Loaded ROI mask from: %s\n', roiFile);
    else
        error([ ...
            'ROI file not found (cluster cannot use roipoly / UI).\n' ...
            'Expected: %s\n\n' ...
            'Fix:\n' ...
            '  1) Run a local MATLAB session with display.\n' ...
            '  2) Load the first frame and run: maskROI = roipoly;\n' ...
            '  3) Save: save(roiFile,''maskROI'') into the mat files/ folder.\n' ...
            '  4) Re-run this cluster job.\n'], roiFile);
    end

    %% ------------------------------------------------------------
    % Optical Flow Setup (RAFT)
    %% ------------------------------------------------------------
    fprintf('\nSTEP 4/7: Initializing RAFT Optical Flow...\n');
    opticalFlowObj = opticalFlowRAFT;

    % Robust frame count estimate
    numFrames     = max(2, floor(v.Duration * v.FrameRate));
    numFlowFrames = numFrames - 1;

    fprintf('Planned frames: numFrames=%d (flow frames=%d) | Execution=%s | accel=%s\n', ...
        numFrames, numFlowFrames, executionEnv, accelMode);

    %% ------------------------------------------------------------
    % Rewind and read first frame for sizing
    %% ------------------------------------------------------------
    v.CurrentTime = 0;
    framePrev = im2gray(readFrame(v));
    [H, W] = size(framePrev);

    if any(size(maskROI) ~= [H W])
        error('maskROI size (%dx%d) does not match video frame size (%dx%d).', ...
            size(maskROI,1), size(maskROI,2), H, W);
    end

    %% ------------------------------------------------------------
    % STREAMING SAVE SETUP (matfile) -> BIG file in "mat files"
    %% ------------------------------------------------------------
    fprintf('\nSTEP 5/7: Preparing streaming MAT-file (incremental writes)...\n');
    uvFile = fullfile(matDir, [filename '_velocity.mat']);
    M = matfile(uvFile, 'Writable', true);

    % Preallocate on disk (NOT in RAM)
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
    % Process frames (write instantaneous frames to disk)
    %% ------------------------------------------------------------
    fprintf('\nSTEP 6/7: Running RAFT per frame (progress will print)...\n');
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

        % --- MEMORY LOGGING (every 200 frames + first frame) ---
        if mod(k, 200) == 0 || k == 1
            pid = feature('getpid');  % MATLAB PID
            statmFile = sprintf('/proc/%d/status', pid);
            if isfile(statmFile)
                txt = fileread(statmFile);
                rssLine = regexp(txt, 'VmRSS:\s+(\d+)\s+kB', 'tokens','once');
                vmsLine = regexp(txt, 'VmSize:\s+(\d+)\s+kB', 'tokens','once');
                if ~isempty(rssLine)
                    VmRSS_GB  = str2double(rssLine{1})/1024/1024;
                    VmSize_GB = str2double(vmsLine{1})/1024/1024;
                    fprintf('[MEM] k=%d | VmRSS=%.2f GB | VmSize=%.2f GB\n', ...
                        k, VmRSS_GB, VmSize_GB);
                end
            end
        end


        % --- Update running sums for mean (double) ---
        sumU   = sumU   + double(u_phys);
        sumV   = sumV   + double(v_phys);
        sumMag = sumMag + hypot(double(u_phys), double(v_phys));
        count  = count + 1;

        % progress
        tElapsed = toc;
        estTotal = tElapsed / max(k,1) * numFlowFrames;
        pct = 100 * k / max(numFlowFrames,1);
        fprintf('Processed flow frame %d/%d (%.1f%%) | Elapsed %.1fs | ETA %.1fs\n', ...
            k, numFlowFrames, pct, tElapsed, max(estTotal - tElapsed, 0));

        framePrev = frameCurr;
    end
    toc;

    %% ------------------------------------------------------------
    % Means
    %% ------------------------------------------------------------
    fprintf('\nSTEP 7/7: Computing mean fields + saving small outputs + plots...\n');
    u_mean  = single(sumU   / max(count,1));
    v_mean  = single(sumV   / max(count,1));
    velMean = single(sumMag / max(count,1));

    u_mean(~maskROI)  = NaN;
    v_mean(~maskROI)  = NaN;
    velMean(~maskROI) = NaN;

    % Physical orientation
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

    %% ---------------- Throat profiles (mean-only, disk streaming) -----------
    fprintf('Computing throat/downstream MEAN profiles (disk streaming)...\n');

    throatFile = fullfile(filepath, [filename '_throat.mat']);
    if ~isfile(throatFile)
        error('Throat file not found: %s (expected x_throat in px).', throatFile);
    end
    Sth = load(throatFile);

    if ~isfield(Sth,'x_throat')
        error('Throat file must contain variable "x_throat" (in pixels).');
    end
    x_throat_pixel = double(Sth.x_throat);
    x_throat_mm = x_throat_pixel * mm_per_pixel;
    
    M.x_throat_pixel = x_throat_pixel;
    M.x_throat_mm    = x_throat_mm;
    

    if isfield(Sth,'y_throat')
        y_throat_pixel = double(Sth.y_throat);
        y_throat_mm = y_throat_pixel * mm_per_pixel;
    else
        y_throat_pixel = NaN; y_throat_mm = NaN;
    end
    
    M.y_throat_pixel = y_throat_pixel;
    M.y_throat_mm    = y_throat_mm;

    [~, ix_throat] = min(abs(x_mm - x_throat_mm));

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

    profileFile = fullfile(matDir, [filename '_ThroatProfiles_mean.mat']);
    save(profileFile, ...
        'x_sample_mm', 'y_mm', ...
        'u_profiles_mean', 'v_profiles_mean', 'vel_profiles_mean', ...
        'x_throat_mm', 'fps','mm_per_pixel', 'm_per_pixel','x_throat_pixel','y_throat_pixel','y_throat_mm', ...
        '-v7.3');

    avgVelFile = fullfile(matDir, [filename '_TimeAvgVelField.mat']);
    save(avgVelFile, 'velMean_phys', 'x_mm', 'y_mm', '-v7.3');

    %% --------------------- PLOTS --------------------------
    if savePlots
        % Plot A
        figA = figure('Visible', tern(showFigures,'on','off'), ...
            'Color','w','Position',[100 100 1200 800]);
        ax = axes(figA);

        imagesc(ax, x_mm, y_mm, velMean_phys);
        set(ax,'YDir','normal'); axis(ax,'image');
        colormap(ax, turbo);

        cb = colorbar(ax, 'eastoutside');
        cb.TickDirection = 'out';
        cb.Label.String  = 'Velocity Magnitude (m/s)';
        cb.Label.FontSize = 14;
        cb.Label.Rotation = 90;
        cb.Label.Units = 'normalized';
        cb.Label.Position = [3 0.5 0];   % label outward
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
        
        % Plot B
        figA = figure('Visible', tern(showFigures,'on','off'), ...
            'Color','w','Position',[100 100 1200 800]);
        ax = axes(figA);

        imagesc(ax, x_mm, y_mm, u_mean_phys);
        set(ax,'YDir','normal'); axis(ax,'image');
        colormap(ax, turbo);

        cb = colorbar(ax, 'eastoutside');
        cb.TickDirection = 'out';
        cb.Label.String  = 'Vx-mean velocity(m/s)';
        cb.Label.FontSize = 14;
        cb.Label.Rotation = 90;
        cb.Label.Units = 'normalized';
        cb.Label.Position = [3 0.5 0];   % label outward
        ax.Position = [0.08 0.10 0.78 0.82];

        hold(ax,'on');
        plot(ax, [x_throat_mm x_throat_mm], [y_mm(1) y_mm(end)], 'w:', 'LineWidth', 1.7);

        title(ax, 'Time-Averaged Vx', 'FontSize', 20, 'Interpreter','latex');
        xlabel(ax, 'X (mm)', 'FontSize', 18, 'Interpreter','latex');
        ylabel(ax, 'Y (mm)', 'FontSize', 18, 'Interpreter','latex');
        set(ax,'FontSize',15,'LineWidth',1.2,'Box','on','FontName','Times');

        outA = fullfile(plotsDir, [filename '_TimeAvgUVelMag_ThroatLine.png']);
        exportgraphics(figA, outA, 'Resolution', 350);
        close(figA);
        fprintf('Saved plot B: %s\n', outA);
        
        % Plot C
        figA = figure('Visible', tern(showFigures,'on','off'), ...
            'Color','w','Position',[100 100 1200 800]);
        ax = axes(figA);

        imagesc(ax, x_mm, y_mm, v_mean_phys);
        set(ax,'YDir','normal'); axis(ax,'image');
        colormap(ax, turbo);

        cb = colorbar(ax, 'eastoutside');
        cb.TickDirection = 'out';
        cb.Label.String  = 'Vy-mean velocity(m/s)';
        cb.Label.FontSize = 14;
        cb.Label.Rotation = 90;
        cb.Label.Units = 'normalized';
        cb.Label.Position = [3 0.5 0];   % label outward
        ax.Position = [0.08 0.10 0.78 0.82];

        hold(ax,'on');
        plot(ax, [x_throat_mm x_throat_mm], [y_mm(1) y_mm(end)], 'w:', 'LineWidth', 1.7);

        title(ax, 'Time-Averaged Vy', 'FontSize', 20, 'Interpreter','latex');
        xlabel(ax, 'X (mm)', 'FontSize', 18, 'Interpreter','latex');
        ylabel(ax, 'Y (mm)', 'FontSize', 18, 'Interpreter','latex');
        set(ax,'FontSize',15,'LineWidth',1.2,'Box','on','FontName','Times');

        outA = fullfile(plotsDir, [filename '_TimeAvgVVelMag_ThroatLine.png']);
        exportgraphics(figA, outA, 'Resolution', 350);
        close(figA);
        fprintf('Saved plot C: %s\n', outA);
        
        
        % Plot D
        fprintf('Saving Plot D (Streamlines strictly clipped to ROI)...\n');

        maskROI_phys = flipud(logical(maskROI));
        maskOut_phys = ~maskROI_phys;

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
        set(ax,'YDir','normal'); axis(ax,'image');
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
        set(hss, 'Color','k', 'LineWidth', 0.8);

        for kk = 1:numel(hss)
            xl = hss(kk).XData; yl = hss(kk).YData;
            if isempty(xl) || isempty(yl), continue; end
            in_roi = interp2(XX, YY, single(maskROI_phys), xl, yl, 'nearest', 0);
            outside_idx = (in_roi < 0.5);
            xl(outside_idx) = NaN; yl(outside_idx) = NaN;
            hss(kk).XData = xl; hss(kk).YData = yl;
        end

        navyRGB = zeros(size(maskOut_phys,1), size(maskOut_phys,2), 3, 'single');
        navyRGB(:,:,1) = 0.02; navyRGB(:,:,2) = 0.05; navyRGB(:,:,3) = 0.25;
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
        fprintf('Saved Plot D: %s\n', outB);
    end

    %% ------------------ OPTIONAL GIF -------------------
    if ~skip_animation
        fprintf('GIF block enabled (skip_animation=false). Reading slices from MAT file...\n');
        % (Your existing GIF code can live here; kept omitted for brevity.)
    else
        fprintf('GIF block skipped (skip_animation=true).\n');
    end

    %% --------------------- FINAL SUMMARY --------------------------
    fprintf('\nDONE.\n');
    fprintf('Instantaneous frames: SAVED incrementally to %s (u_all, v_all)\n', uvFile);
    fprintf('Folders:\n  MAT:   %s\n  Plots: %s\n', matDir, plotsDir);

    fprintf('\nSaved MAT files (3) + ROI (1):\n');
    fprintf('  1) %s\n', uvFile);
    fprintf('  2) %s\n', profileFile);
    fprintf('  3) %s\n', avgVelFile);
    fprintf('  ROI) %s\n', roiFile);

    fprintf('\nSaved plots (2):\n');
    fprintf('  A) %s\n', fullfile(plotsDir, [filename '_TimeAvgVelMag_ThroatLine.png']));
    fprintf('  B) %s\n', fullfile(plotsDir, [filename '_TimeAvgUVelMag_ThroatLine.png']));
    fprintf('  C) %s\n', fullfile(plotsDir, [filename '_TimeAvgVVelMag_ThroatLine.png']));
    fprintf('  D) %s\n', fullfile(plotsDir, [filename '_TimeAvgVelMag_Streamlines_ROI.png']));

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
if cond, out = a; else, out = b; end
end
