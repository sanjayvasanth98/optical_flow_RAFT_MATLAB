%% ------------------------------------------------------------
% Optical Flow RAFT
% Author: Sanjay Vasanth
% Last modified: 2/3/2026
% Time-averaged velocity + vertical profiles with axes in mm (origin at lower-left)
% Incremental saving to avoid huge end-of-run save time / RAM blowups
% New files store u_all/v_all as [nnz(maskROI) x numFlowFrames] singles.
% Rows follow find(maskROI), in image coordinates; maskROI restores the grid.
% Saves:
%   - MAT files into:   <video_folder>/mat files/
%   - Plots into:       <video_folder>/plots/
%% ------------------------------------------------------------
clear all; clc; close all;

%% --------------------- USER TOGGLES --------------------------
showFigures    = false;   % MUST be false on cluster
savePlots      = true;    % saves Plot A and Plot B
skip_animation = true;    % if true, GIF block is skipped

% Same contrast preprocessing as the local runner, on every RAFT input.
% All results use the standard mat files/ and plots/ folders.
% For a fresh run, move any existing velocity MAT file out of mat files/ first.
% Set "none" to resume the original unenhanced video results.
contrastMode = "clahe";   % "none" or "clahe"
contrastMode = lower(string(contrastMode));
if ~isscalar(contrastMode) || ~ismember(contrastMode, ["none", "clahe"])
    error('contrastMode must be "none" or "clahe".');
end

% Resume controls: existing velocity files are NEVER reinitialized.
% For the first recovery of an OLD run, supply either the exact N from the
% last "Processed flow frame N/..." line OR that run's SLURM output path.
% These legacy settings apply only to contrastMode="none".
% Later restarts use the completion counter automatically.
legacyCompletedFlowFrames = [];   % e.g. 12345 (NOT the percentage)
legacyProgressLog = '';          % Empty for a new run; set only for legacy recovery.

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
mm_per_pixel = 0.00828164;          % [mm/pixel]   %<--edit
fps          = 102247;               % [frames per second]  %<--edit
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
    videoPath = '/home/kbsanjayvasanth/Sept2026_flowfield_RAFT/P10S20/P10S20_3_475_40lpm.avi';   %<--edit
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
    resultsDir = filepath;
    matDir = fullfile(filepath, 'mat files');
    fprintf('Contrast preprocessing: %s | Results: %s\n', contrastMode, resultsDir);
    plotsDir = fullfile(resultsDir, 'plots');
    if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end

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
    % STREAMING SAVE SETUP (matfile) -> ROI-only instantaneous values
    %% ------------------------------------------------------------
    fprintf('\nSTEP 5/7: Preparing streaming MAT-file (incremental writes)...\n');
    uvFile = fullfile(matDir, [filename '_velocity.mat']);
    videoInfo = dir(videoPath);
    resumeConfig = struct('videoBytes', videoInfo.bytes, ...
        'videoModified', videoInfo.datenum, 'height', H, 'width', W, ...
        'numFlowFrames', numFlowFrames, 'raftIters', raftIters, ...
        'raftTolerance', raftTolerance, 'accelMode', accelMode, ...
        'contrastMode', contrastMode);
    numROIPixels = nnz(maskROI);
    if numROIPixels == 0
        error('ROI mask must contain at least one pixel.');
    end
    roiPacked = true;
    completed = 0;
    if isfile(uvFile)
        M = matfile(uvFile, 'Writable', true);
        vars = who(M);
        required = {'u_all','v_all','mm_per_pixel','m_per_pixel','fps','maskROI'};
        if ~all(ismember(required, vars))
            error('Existing velocity file is incomplete. Preserve it and inspect: %s', uvFile);
        end
        % Explicit metadata distinguishes packed arrays, including single-frame
        % and single-pixel cases, from legacy full-frame arrays.
        roiPacked = ismember('instantaneousStorage', vars);
        if roiPacked
            if ~strcmp(M.instantaneousStorage, 'roi_pixels_by_frame_v1')
                error('Unknown instantaneous storage format in: %s', uvFile);
            end
            expectedSize = [numROIPixels numFlowFrames 1];
        else
            expectedSize = [H W numFlowFrames];
        end
        if ~isequal([size(M,'u_all',1) size(M,'u_all',2) size(M,'u_all',3)], expectedSize) || ...
           ~isequal([size(M,'v_all',1) size(M,'v_all',2) size(M,'v_all',3)], expectedSize)
            error('Existing velocity dimensions do not match this video: %s', uvFile);
        end
        if ~isequal(M.mm_per_pixel, mm_per_pixel) || ~isequal(M.fps, fps) || ...
           ~isequal(M.m_per_pixel, m_per_pixel) || ~isequal(M.maskROI, maskROI)
            error('Calibration or ROI changed. Restore the original settings before resuming.');
        end
        if ismember('resumeConfig', vars)
            savedConfig = M.resumeConfig;
            % Previous ARC versions used grayscale without CLAHE.
            if ~isfield(savedConfig, 'contrastMode')
                savedConfig.contrastMode = "none";
            end
            if ~isequaln(savedConfig, resumeConfig)
                error('Video, RAFT, or contrast settings changed. Restore the original settings before resuming.');
            end
        elseif contrastMode ~= "none"
            error('Existing file has no contrast metadata; cannot resume it as a CLAHE run.');
        end
        if ismember('completedFlowFrames', vars)
            completed = M.completedFlowFrames;
        else
            if contrastMode ~= "none"
                error('CLAHE file is missing its completion counter. Inspect it before resuming.');
            end
            completed = legacyCompletedFlowFrames;
            if isempty(completed) && ~isempty(legacyProgressLog)
                logText = fileread(legacyProgressLog);
                tokens = regexp(logText, 'Processed flow frame (\d+)/(\d+)', 'tokens');
                if ~isempty(tokens)
                    completed = str2double(tokens{end}{1});
                    if str2double(tokens{end}{2}) ~= numFlowFrames
                        error('Legacy log frame total does not match this video.');
                    end
                end
            end
            if isempty(completed)
                error(['Existing file has no completion counter. Set legacyCompletedFlowFrames ' ...
                    'to N from the last Processed flow frame N/... log line, or set ' ...
                    'legacyProgressLog to that job log. Existing data was left untouched.']);
            end
            warning(['Recovering an old file using the supplied progress. Ensure it belongs ' ...
                'to this video and used the same RAFT settings. Saved frames are preserved.']);
        end
        validateattributes(completed, {'numeric'}, ...
            {'scalar','real','finite','integer','>=',0,'<=',numFlowFrames});
        if ~roiPacked
            warning(['Resuming legacy full-frame storage. ROI-only storage applies to new ' ...
                'velocity files; existing files are preserved in their original format.']);
        end
        M.resumeConfig = resumeConfig;
        M.completedFlowFrames = completed;
        fprintf('Resuming existing file: %s | %d/%d flow frames saved.\n', ...
            uvFile, completed, numFlowFrames);
    else
        if contrastMode == "none" && (~isempty(legacyCompletedFlowFrames) || ~isempty(legacyProgressLog))
            error('Recovery requested but the existing velocity file was not found: %s', uvFile);
        end
        M = matfile(uvFile, 'Writable', true);
        % Grow arrays on disk without allocating the full video in RAM.
        M.u_all(numROIPixels,numFlowFrames) = single(0);
        M.v_all(numROIPixels,numFlowFrames) = single(0);
        M.instantaneousStorage = 'roi_pixels_by_frame_v1';
        M.mm_per_pixel = mm_per_pixel;
        M.m_per_pixel = m_per_pixel;
        M.fps = fps;
        M.maskROI = maskROI;
        M.resumeConfig = resumeConfig;
        M.completedFlowFrames = 0;
        fprintf('Streaming file created: %s\n', uvFile);
    end

    if roiPacked
        fprintf('ROI-only storage: %d/%d pixels per frame (%.1f%% of full-grid values).\n', ...
            numROIPixels, H*W, 100*numROIPixels/(H*W));
    end

    % Reconstruct sums from committed slices, never from unwritten tail zeros.
    % Using stored singles for both new and recovered frames keeps means consistent.
    sumU = zeros(H, W, 'double');
    sumV = zeros(H, W, 'double');
    sumMag = zeros(H, W, 'double');
    count = completed;
    for k = 1:completed
        if roiPacked
            savedU = double(M.u_all(:,k));
            savedV = double(M.v_all(:,k));
        else
            if numFlowFrames == 1
                savedU = double(M.u_all(:,:));
                savedV = double(M.v_all(:,:));
            else
                savedU = double(M.u_all(:,:,k));
                savedV = double(M.v_all(:,:,k));
            end
            savedU = savedU(maskROI);
            savedV = savedV(maskROI);
        end
        sumU(maskROI) = sumU(maskROI) + savedU;
        sumV(maskROI) = sumV(maskROI) + savedV;
        sumMag(maskROI) = sumMag(maskROI) + hypot(savedU, savedV);
        if mod(k,200) == 0 || k == completed
            fprintf('Rebuilding averages: %d/%d saved frames\n', k, completed);
        end
    end

    if completed < numFlowFrames
        % Read sequentially to preserve exact frame numbering, including videos
        % for which timestamp seeking is imprecise. No RAFT work during skipping.
        for j = 1:completed
            framePrev = im2gray(readFrame(v));
            if mod(j,1000) == 0 || j == completed
                fprintf('Restoring video position: %d/%d frames\n', j, completed);
            end
        end
        % RAFT caches the previous image internally. Seed it with video frame
        % completed+1 so the next output is the pair completed+1 -> completed+2.
        framePrev = prepareFrame(framePrev, contrastMode);
        estimateFlow(opticalFlowObj, framePrev, ...
            ExecutionEnvironment=executionEnv, Acceleration=accelMode, ...
            MaxIterations=raftIters, Tolerance=raftTolerance);
    end

    %% ------------------------------------------------------------
    % Process frames (write instantaneous frames to disk)
    %% ------------------------------------------------------------
    fprintf('\nSTEP 6/7: Running RAFT per frame (progress will print)...\n');
    tic;
    for i = completed+2:numFrames
        frameCurr = prepareFrame(readFrame(v), contrastMode);

        flow = estimateFlow(opticalFlowObj, frameCurr, ...
            ExecutionEnvironment=executionEnv, ...
            Acceleration=accelMode, ...
            MaxIterations=raftIters, ...
            Tolerance=raftTolerance);

        % Extract only ROI pixels in MATLAB linear (column-major) order.
        % Values remain in image coordinates, calibrated to m/s.
        u_roi = single(flow.Vx(maskROI) * m_per_pixel * fps);
        v_roi = single(flow.Vy(maskROI) * m_per_pixel * fps);

        % --- Write instantaneous to disk (single) ---
        k = i - 1;
        if roiPacked
            M.u_all(:,k) = u_roi;
            M.v_all(:,k) = v_roi;
        else
            % Preserve the layout when resuming an existing legacy file.
            u_frame = zeros(H, W, 'single');
            v_frame = zeros(H, W, 'single');
            u_frame(maskROI) = u_roi;
            v_frame(maskROI) = v_roi;
            if numFlowFrames == 1
                M.u_all(:,:) = u_frame;
                M.v_all(:,:) = v_frame;
            else
                M.u_all(:,:,k) = u_frame;
                M.v_all(:,:,k) = v_frame;
            end
        end
        % Commit only AFTER both velocity components have been written.
        % A timeout before this marker causes this pair to be recomputed.
        M.completedFlowFrames = k;

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
        sumU(maskROI) = sumU(maskROI) + double(u_roi);
        sumV(maskROI) = sumV(maskROI) + double(v_roi);
        sumMag(maskROI) = sumMag(maskROI) + hypot(double(u_roi), double(v_roi));
        count  = count + 1;

        % progress
        tElapsed = toc;
        remainingSeconds = tElapsed / (k - completed) * (numFlowFrames - k);
        pct = 100 * k / max(numFlowFrames,1);
        fprintf('Processed flow frame %d/%d (%.1f%%) | Elapsed %.1fs | ETA %.1fs\n', ...
            k, numFlowFrames, pct, tElapsed, max(remainingSeconds, 0));

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
    fprintf('Computing throat/downstream MEAN profiles from accumulated fields...\n');

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

    % The profiles are column samples of the same means. Reuse the running
    % sums instead of rereading the entire large velocity file a second time.
    % Preserve zero values outside the ROI as in the original profile output.
    u_profiles_mean = single(flipud(sumU(:,ix_samples)) / max(count,1));
    v_profiles_mean = single(-flipud(sumV(:,ix_samples)) / max(count,1));
    vel_profiles_mean = single(flipud(sumMag(:,ix_samples)) / max(count,1));

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
    rethrow(ME);  % Let MATLAB batch / SLURM report failure.
end

%% --------------------- helper ---------------------
function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

% Keep this preprocessing identical to raftmatlabsideview_local.m.
function frame = prepareFrame(frame, contrastMode)
frame = im2gray(frame);
if contrastMode == "clahe"
    if canUseGPU
        frame = gpuClahe(frame);
    else
        frame = adapthisteq(frame);
    end
end
end

function frame = gpuClahe(frame)
frameGPU = gpuArray(im2single(frame));
[height, width] = size(frameGPU);
numTileRows = min(8, height);
numTileCols = min(8, width);
outputGPU = zeros(height, width, 'like', frameGPU);

for tileRow = 1:numTileRows
    rowStart = floor((tileRow - 1) * height / numTileRows) + 1;
    rowEnd = floor(tileRow * height / numTileRows);
    for tileCol = 1:numTileCols
        colStart = floor((tileCol - 1) * width / numTileCols) + 1;
        colEnd = floor(tileCol * width / numTileCols);
        tileGPU = frameGPU(rowStart:rowEnd, colStart:colEnd);

        [histogramGPU, ~] = histcounts(tileGPU, 256, 'BinLimits', [0 1]);
        clipLimit = max(1, 0.01 * numel(tileGPU) / 256);
        excess = sum(max(histogramGPU - clipLimit, 0));
        histogramGPU = min(histogramGPU, clipLimit) + excess / 256;
        cdfGPU = cumsum(histogramGPU) / numel(tileGPU);

        binGPU = min(floor(tileGPU * 256) + 1, 256);
        outputGPU(rowStart:rowEnd, colStart:colEnd) = cdfGPU(binGPU);
    end
end

frame = gather(im2uint8(outputGPU));
end
