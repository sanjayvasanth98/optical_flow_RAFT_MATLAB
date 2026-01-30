%% ------------------------------------------------------------
% Optical Flow RAFT — ARC Optimized
% Author:Sanjay Vasanth KB
% Date: 1/29/2026
% ------------------------------------------------------------
clearvars; clc;

%% --------------------- USER TOGGLES --------------------------
showFigures    = false;
savePlots      = true;
saveGIF        = false;
saveMatOutputs = true;

% RAFT controls
raftIters    = 8;
raftTol      = 1e-6;
accelMode    = "auto";

% Default GPU
if canUseGPU
    executionEnv = "gpu";
else
    executionEnv = "cpu";
end

% ROI-respecting neighbor averaging (masked, normalized 5-point)
useMaskedNeighborAverage = true;

% Buffered writes
writeChunk = 100;

% GIF settings
gifStartFrame = 200;
gifEndFrame   = 300;
gifDelayTime  = 0.4;
gifSkip       = 15;

%% --------------------- FIGURE VISIBILITY ---------------------
if ~showFigures
    set(0, 'DefaultFigureVisible', 'off');
end

%% ------------------------------------------------------------
% PATH SETUP (QUIET + BATCH SAFE)
% ------------------------------------------------------------
warnState = warning;
warning('off','MATLAB:dispatcher:nameConflict');
warning('off','MATLAB:mpath:nameNonexistent');
warning('off','MATLAB:general:NameIsNonexistentOrNotADirectory');
warning('off','MATLAB:rmpath:DirNotFound');

userCodeRoot = fullfile(getenv('HOME'), 'matlab');
if ~usejava('desktop')
    if isfolder(userCodeRoot)
        addpath(userCodeRoot);
        cand = { ...
            fullfile(userCodeRoot,'raft_newcode'), ...
            fullfile(userCodeRoot,'utils'), ...
            fullfile(userCodeRoot,'functions') ...
        };
        for k = 1:numel(cand)
            if isfolder(cand{k}), addpath(genpath(cand{k})); end
        end
    end
end
warning(warnState);

try
    %% --- Specify video path ---
    videoPath = '/home/kbsanjayvasanth/raft_newcode/Smooth_48lpm_inc.avi';
    [filepath, filename, ~] = fileparts(char(videoPath));

    %% --- Calibration parameters ---
    mm_per_pixel = 0.0083333;
    fps          = 102247;
    m_per_pixel  = mm_per_pixel / 1000;
    scale        = single(m_per_pixel * fps);

    fprintf('Calibration: %.9f m/pixel | Frame rate: %.1f fps\n', m_per_pixel, fps);
    fprintf('ExecutionEnvironment: %s | Acceleration: %s | RAFT iters: %d\n', executionEnv, accelMode, raftIters);

    if ~isfile(videoPath), error('Video file not found at: %s', videoPath); end
    v = VideoReader(videoPath);
    fprintf('Loaded video: %s\n', filename);

    %% --- Output folders ---
    plotsDir = fullfile(filepath, 'plots');
    matDir   = fullfile(filepath, 'mat files');
    if ~exist(plotsDir, 'dir'); mkdir(plotsDir); end
    if ~exist(matDir,   'dir'); mkdir(matDir);   end

    %% ------------------------------------------------------------
    % PERMANENT FIX: lock reference size to FIRST decoded frame
    % ------------------------------------------------------------
    v.CurrentTime = 0;
    frame0 = im2gray(readFrame(v));
    [H0, W0] = size(frame0);
    x_mm = (0:W0-1) * mm_per_pixel;
    y_mm = (0:H0-1) * mm_per_pixel;
    fprintf('Reference frame size locked to: %d x %d\n', H0, W0);

    %% ------------------------------------------------------------
    % Load ROI mask (force to ref size; persist)
    % ------------------------------------------------------------
    roiFile = fullfile(matDir, [filename '_ROI.mat']);

    if isfile(roiFile)
        S = load(roiFile);
        if ~isfield(S,'maskROI'), error('ROI file missing maskROI: %s', roiFile); end
        maskROI = logical(S.maskROI);
    else
        maskROI = true(H0, W0);
        save(roiFile, 'maskROI');
        fprintf('No ROI found -> saved full-frame ROI to: %s\n', roiFile);
    end

    maskROI = forceMaskToSize(maskROI, H0, W0);
    save(roiFile, 'maskROI');
    fprintf('ROI loaded/saved: %s (size %dx%d)\n', roiFile, size(maskROI,1), size(maskROI,2));

    %% ------------------------------------------------------------
    % Load throat file (plot A throat line + profile start)
    % ------------------------------------------------------------
    throatList = dir(fullfile(matDir, '*_throat.mat'));
    if isempty(throatList), throatList = dir(fullfile(filepath, '*_throat.mat')); end
    if isempty(throatList), error('No throat file found ("*_throat.mat").'); end
    if numel(throatList) > 1
        fprintf('Multiple throat files found; using first: %s\n', throatList(1).name);
    end
    Sth = load(fullfile(throatList(1).folder, throatList(1).name));
    if ~isfield(Sth,'x_throat'), error('Throat file missing x_throat.'); end
    x_throat_mm = double(Sth.x_throat) * mm_per_pixel;
    [~, ix_throat] = min(abs(x_mm - x_throat_mm));
    fprintf('Throat x = %.3f mm | nearest index = %d\n', x_throat_mm, ix_throat);

    %% ------------------------------------------------------------
    % Define profile sample stations (vectorized extraction)
    % ------------------------------------------------------------
    x_spacing   = 0.25;
    x_sample_mm = x_mm(ix_throat):x_spacing:x_mm(end);
    ix_samples  = arrayfun(@(x) find(abs(x_mm-x)==min(abs(x_mm-x)), 1), x_sample_mm);
    numX = numel(ix_samples);
    numY = numel(y_mm);
    fprintf('Profiles: %d stations every %.2f mm from throat\n', numX, x_spacing);

    %% ------------------------------------------------------------
    % Determine numFrames robustly
    % ------------------------------------------------------------
    if ~isprop(v,'NumFrames') || isempty(v.NumFrames) || isnan(v.NumFrames) || v.NumFrames <= 0
        numFrames = max(2, floor(v.Duration * v.FrameRate));
        fprintf('NumFrames unavailable -> using Duration*FrameRate = %d\n', numFrames);
    else
        numFrames = v.NumFrames;
    end
    fprintf('Processing %d frames...\n', numFrames);

    %% ------------------------------------------------------------
    % Initialize RAFT + warmup
    % ------------------------------------------------------------
    opticalFlowObj = opticalFlowRAFT;

    if executionEnv=="gpu"
        try, gpuDevice; catch, executionEnv="cpu"; fprintf('GPU not available -> switching to CPU\n'); end
    end

    v.CurrentTime = 0;
    f2 = im2gray(readFrame(v)); f2 = forceFrameToSize(f2, H0, W0);
    estimateFlow(opticalFlowObj, f2, ...
        ExecutionEnvironment=executionEnv, Acceleration=accelMode, ...
        MaxIterations=2, Tolerance=1e-3);

    %% ------------------------------------------------------------
    % Precompute ROI-respecting 5-point averaging weights (FULL FRAME)
    % ------------------------------------------------------------
    if useMaskedNeighborAverage
        H = H0; W = W0;
        M = single(maskROI);
        M_up    = [M(2:end,:); zeros(1,W,'single')];
        M_down  = [zeros(1,W,'single'); M(1:end-1,:)];
        M_left  = [M(:,2:end) zeros(H,1,'single')];
        M_right = [zeros(H,1,'single') M(:,1:end-1)];
        denom = M + M_up + M_down + M_left + M_right;
        denom(denom==0) = 1;
    end

    %% ------------------------------------------------------------
    % Streamed sums
    % ------------------------------------------------------------
    u_sum   = zeros(H0, W0, 'single');
    v_sum   = zeros(H0, W0, 'single');
    vel_sum = zeros(H0, W0, 'single');
    n_avg   = 0;

    vel_profiles_sum = zeros(numY, numX, 'single');
    n_profile = 0;

    %% ------------------------------------------------------------
    % Chunked MAT writes (every writeChunk frames)
    % ------------------------------------------------------------
    if saveMatOutputs
        incFile = fullfile(matDir, [filename '_incremental.mat']);
        m = matfile(incFile, 'Writable', true);

        m.mm_per_pixel = mm_per_pixel;
        m.m_per_pixel  = m_per_pixel;
        m.fps          = fps;
        m.x_mm         = x_mm;
        m.y_mm         = y_mm;
        m.maskROI      = maskROI;
        m.x_sample_mm  = x_sample_mm;
        m.ix_samples   = ix_samples;
        m.executionEnv = executionEnv;
        m.accelMode    = accelMode;
        m.raftIters    = raftIters;
        m.raftTol      = raftTol;

        m.vel_prof_frame = zeros(numY, numX, numFrames-1, 'single');
        m.meanVelROI     = zeros(numFrames-1, 1, 'single');

        fprintf('Incremental MAT (chunked writes): %s\n', incFile);
    end

    bufCount   = 0;
    bufStartK  = 1;
    buf_vel_prof = zeros(numY, numX, min(writeChunk, numFrames-1), 'single');
    buf_meanVel  = zeros(min(writeChunk, numFrames-1), 1, 'single');

    %% ------------------------------------------------------------
    % GIF setup (optional)
    % ------------------------------------------------------------
    if saveGIF
        [Xq, Yq] = meshgrid(x_mm(1:gifSkip:end), y_mm(1:gifSkip:end));
        gifFile = fullfile(plotsDir, [filename sprintf('_VelField_%dto%d_highQ.gif', gifStartFrame, gifEndFrame)]);
        figG = figure('Visible','off', 'Position',[100 100 1400 1050], 'Renderer','painters');
        colormap(figG, turbo);
        firstGifFrame = true;
        fprintf('GIF: %s\n', gifFile);
    end

    %% ------------------------------------------------------------
    % MAIN LOOP
    % ------------------------------------------------------------
    v.CurrentTime = 0;
    readFrame(v); % discard frame 1
    t0 = tic;

    for i = 2:numFrames
        frameCurr = im2gray(readFrame(v));
        frameCurr = forceFrameToSize(frameCurr, H0, W0);

        flow = estimateFlow(opticalFlowObj, frameCurr, ...
            ExecutionEnvironment=executionEnv, Acceleration=accelMode, ...
            MaxIterations=raftIters, Tolerance=raftTol);

        U = single(flow.Vx);
        V = single(flow.Vy);

        if useMaskedNeighborAverage
            U_up    = [U(2:end,:); zeros(1,W0,'single')];
            U_down  = [zeros(1,W0,'single'); U(1:end-1,:)];
            U_left  = [U(:,2:end) zeros(H0,1,'single')];
            U_right = [zeros(H0,1,'single') U(:,1:end-1)];

            V_up    = [V(2:end,:); zeros(1,W0,'single')];
            V_down  = [zeros(1,W0,'single'); V(1:end-1,:)];
            V_left  = [V(:,2:end) zeros(H0,1,'single')];
            V_right = [zeros(H0,1,'single') V(:,1:end-1)];

            U_num = U.*M + U_up.*M_up + U_down.*M_down + U_left.*M_left + U_right.*M_right;
            V_num = V.*M + V_up.*M_up + V_down.*M_down + V_left.*M_left + V_right.*M_right;

            U = U_num ./ denom;
            V = V_num ./ denom;
        end

        % Hard ROI mask
        U = U .* single(maskROI);
        V = V .* single(maskROI);

        u_phys = U * scale;
        v_phys = V * scale;
        velMag = hypot(u_phys, v_phys);

        u_sum   = u_sum + u_phys;
        v_sum   = v_sum + v_phys;
        vel_sum = vel_sum + velMag;
        n_avg   = n_avg + 1;

        velP = flipud(velMag);
        profMat = velP(:, ix_samples);

        vel_profiles_sum = vel_profiles_sum + profMat;
        n_profile = n_profile + 1;

        bufCount = bufCount + 1;
        k = i-1;
        buf_vel_prof(:,:,bufCount) = profMat;

        roiVals = velMag(maskROI);
        if ~isempty(roiVals)
            buf_meanVel(bufCount,1) = mean(roiVals, 'omitnan');
        else
            buf_meanVel(bufCount,1) = NaN;
        end

        if saveMatOutputs && (bufCount == writeChunk || k == (numFrames-1))
            kStart = bufStartK;
            kEnd   = bufStartK + bufCount - 1;

            m.vel_prof_frame(:,:,kStart:kEnd) = buf_vel_prof(:,:,1:bufCount);
            m.meanVelROI(kStart:kEnd,1)       = buf_meanVel(1:bufCount,1);

            bufStartK = kEnd + 1;
            bufCount  = 0;
        end

        if mod(i, 10) == 0 || i == numFrames
            tElapsed = toc(t0);
            pct = 100 * i / numFrames;
            estTotal = tElapsed / max(i-1,1) * numFrames;
            fprintf('Processed %d/%d (%.1f%%) | Elapsed %.1fs | ETA %.1fs\n', ...
                i, numFrames, pct, tElapsed, max(estTotal - tElapsed, 0));
        end

        if saveGIF && (i >= gifStartFrame) && (i <= gifEndFrame)
            velPhys = flipud(velMag);
            imagesc(x_mm, y_mm, velPhys); axis image; set(gca,'YDir','normal'); colorbar; hold on;
            uP = flipud(u_phys); vP = -flipud(v_phys);
            quiver(Xq, Yq, uP(1:gifSkip:end,1:gifSkip:end), vP(1:gifSkip:end,1:gifSkip:end), ...
                'k', 'LineWidth', 1, 'AutoScale','on', 'AutoScaleFactor', 0.8);
            text(x_mm(end)*0.98, y_mm(end)*0.03, sprintf('Frame: %d', i), ...
                'Color','r','FontWeight','bold','FontSize',16, 'HorizontalAlignment','right');
            drawnow;

            fr = getframe(figG);
            [imind, cm] = rgb2ind(frame2im(fr), 256, 'nodither');
            if firstGifFrame
                imwrite(imind, cm, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', gifDelayTime);
                firstGifFrame = false;
            else
                imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', gifDelayTime);
            end
            hold off;
        end
    end

    if saveGIF
        close(figG);
        fprintf('GIF saved: %s\n', gifFile);
    end

    %% ------------------------------------------------------------
    % Time-averaged fields
    % ------------------------------------------------------------
    u_mean  = u_sum   / max(n_avg,1);
    v_mean  = v_sum   / max(n_avg,1);
    velMean = vel_sum / max(n_avg,1);

    u_mean(~maskROI)  = NaN;
    v_mean(~maskROI)  = NaN;
    velMean(~maskROI) = NaN;

    u_mean_phys   = flipud(u_mean);
    v_mean_phys   = -flipud(v_mean);
    velMean_phys  = flipud(velMean);
    maskROI_phys  = flipud(maskROI);

    vel_profiles_mean = vel_profiles_sum / max(n_profile,1);

    if saveMatOutputs
        meanFile = fullfile(matDir, [filename '_means.mat']);
        save(meanFile, ...
            'u_mean_phys','v_mean_phys','velMean_phys', ...
            'x_mm','y_mm','maskROI','maskROI_phys', ...
            'mm_per_pixel','m_per_pixel','fps', ...
            'x_sample_mm','ix_samples','vel_profiles_mean', ...
            '-v7.3');
        fprintf('Saved means: %s\n', meanFile);
    end

%% ------------------------------------------------------------
% Plot (A): magnitude only + throat line (normal colorbar + no overlap)
% ------------------------------------------------------------
if savePlots
    figA = figure('Visible', tern(showFigures,'on','off'), ...
        'Color','w','Position',[100 100 1200 800]);  % normal size (not extreme)

    ax = axes(figA);

    imagesc(ax, x_mm, y_mm, velMean_phys);
    set(ax,'YDir','normal');
    axis(ax,'image');
    colormap(ax, turbo);

    % Normal colorbar (MATLAB default sizing) + fix label overlap only
    cb = colorbar(ax, 'eastoutside');
    cb.TickDirection = 'out';
    cb.Label.String  = 'Velocity Magnitude (m/s)';
    cb.Label.FontSize = 14;
    cb.Label.Rotation = 90;

    % ---- overlap fix: move LABEL outward without resizing the bar ----
    cb.Label.Units = 'normalized';
    cb.Label.Position = [3 0.5 0];   % push right; try 1.25–1.6 if needed

    % Optional: give a little extra right whitespace without changing bar size
    ax.Position = [0.08 0.10 0.78 0.82];  % small tweak, keeps bar normal

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
% Plot (B): Time-averaged |V| + arrowed streamlines ONLY inside ROI
% ------------------------------------------------------------
if savePlots
    fprintf('Saving Plot B (Streamlines strictly clipped to ROI)... \n');

    % ----- HARD GUARANTEE: use ROI loaded from file -----
    if ~exist('maskROI','var') || isempty(maskROI)
        error('maskROI not found in workspace. Load ROI before Plot B.');
    end

    % Physical-orientation mask consistent with velMean_phys
    maskROI_phys = flipud(logical(maskROI));   % ROI==1 inside
    maskOut_phys = ~maskROI_phys;              % ROI==0 outside

    % Safety: ensure sizes match
    [Hf, Wf] = size(velMean_phys);
    if any(size(maskROI_phys) ~= [Hf Wf])
        error('maskROI_phys size (%dx%d) does not match velMean_phys (%dx%d).', ...
            size(maskROI_phys,1), size(maskROI_phys,2), Hf, Wf);
    end

    figB = figure('Visible', tern(showFigures,'on','off'), ...
        'Color','w','Position',[100 100 1200 800]);
    ax = axes(figB);

    % ----- Masked fields (Data Prep) -----
    uPlot = u_mean_phys;
    vPlot = v_mean_phys;
    qPlot = velMean_phys;

    % Set Outside ROI to NaN (Standard step)
    uPlot(maskOut_phys) = NaN;
    vPlot(maskOut_phys) = NaN;
    qPlot(maskOut_phys) = NaN;

    % ----- 1. Background magnitude -----
    imagesc(ax, x_mm, y_mm, qPlot);
    set(ax,'YDir','normal');
    axis(ax,'image');
    colormap(ax, turbo);

    % ----- Colorbar -----
    cb = colorbar(ax);
    cb.TickDirection  = 'out';
    cb.Label.String   = 'Velocity Magnitude (m/s)';
    cb.Label.FontSize = 13;
    cb.Label.Rotation = 90;
    cb.Label.Units = 'normalized';
    cb.Label.Position = [3 0.5 0]; 

    hold(ax,'on');

    % ----- 2. Draw Streamlines -----
    [XX, YY] = meshgrid(x_mm, y_mm);
    density = 0.8; 
    hss = streamslice(XX, YY, uPlot, vPlot, density);
    set(hss, 'Color','k', 'LineWidth', 0.7);

    % ----- 2b. STRICT FIX: Physically delete streamline points outside ROI -----
    % Iterate over every streamline object and check its vertices against the mask
    for k = 1:numel(hss)
        xl = hss(k).XData;
        yl = hss(k).YData;
        
        % If streamslice produced empty or invalid handles (rare), skip
        if isempty(xl) || isempty(yl), continue; end
        
        % Interpolate the ROI mask at the streamline vertices
        % 'nearest' is safest for binary masks. 0 is extrapolation value.
        in_roi = interp2(XX, YY, single(maskROI_phys), xl, yl, 'nearest', 0);
        
        % Find points where mask is 0 (outside ROI)
        outside_idx = (in_roi < 0.5);
        
        % Set those specific vertices to NaN (breaks the line)
        xl(outside_idx) = NaN;
        yl(outside_idx) = NaN;
        
        % Update the line object
        hss(k).XData = xl;
        hss(k).YData = yl;
    end
    
    % ----- 3. Dark-blue OPAQUE Overlay (Backup Visual Blocker) -----
    navyRGB = zeros(size(maskOut_phys,1), size(maskOut_phys,2), 3, 'single');
    navyRGB(:,:,1) = 0.02;  % R
    navyRGB(:,:,2) = 0.05;  % G
    navyRGB(:,:,3) = 0.25;  % B (navy)

    hMask = image(ax, x_mm, y_mm, navyRGB);
    set(hMask, 'AlphaData', 1.0 * single(maskOut_phys)); % Fully opaque outside
    uistack(hMask, 'top'); 

    % ----- Cosmetics -----
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
% Save ROI mask as binary image (1 = ROI, 0 = outside)
% ------------------------------------------------------------
roiImg = uint8(maskROI_phys) * 255;   % white ROI, black background

roiFilePNG = fullfile(plotsDir, [filename '_ROI_mask.png']);
imwrite(uint8(~maskROI_phys)*255, fullfile(plotsDir,[filename '_ROI_outside.png']));

fprintf('Saved ROI binary mask image: %s\n', roiFilePNG);

    fprintf('--- Analysis complete ---\n');

catch ME
    homeDir = getenv('HOME'); if isempty(homeDir), homeDir = '/tmp'; end
    errFile = fullfile(homeDir, 'RAFT_errorLog.txt');
    fid = fopen(errFile,'w');
    if fid ~= -1
        fprintf(fid, 'Error occurred on %s\n\n', datestr(now));
        fprintf(fid, 'Error message:\n%s\n\n', ME.message);
        fprintf(fid, 'Stack trace:\n');
        for k = 1:numel(ME.stack)
            fprintf(fid, '  File: %s (line %d)\n', ME.stack(k).file, ME.stack(k).line);
        end
        fclose(fid);
    end
    fprintf(2, 'Error: %s\n', ME.message);
end

%% --------------------- helper functions ------------------------
function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

function fOut = forceFrameToSize(fIn, H0, W0)
[h,w] = size(fIn);
if h==H0 && w==W0
    fOut = fIn; return;
end
if ~isa(fIn,'uint8'), fIn = im2uint8(fIn); end

r1 = max(1, floor((h - H0)/2) + 1);
c1 = max(1, floor((w - W0)/2) + 1);
r2 = min(h, r1 + H0 - 1);
c2 = min(w, c1 + W0 - 1);
crop = fIn(r1:r2, c1:c2);

fOut = crop;
[hc,wc] = size(fOut);
if hc < H0 || wc < W0
    top     = floor((H0-hc)/2);
    bottom = H0-hc-top;
    left    = floor((W0-wc)/2);
    right  = W0-wc-left;
    fOut = padarray(fOut, [top left], 'replicate', 'pre');
    fOut = padarray(fOut, [bottom right], 'replicate', 'post');
end
fOut = fOut(1:H0, 1:W0);
end

function maskOut = forceMaskToSize(maskIn, H0, W0)
maskIn = logical(maskIn);
[h,w] = size(maskIn);
if h==H0 && w==W0
    maskOut = maskIn; return;
end

r1 = max(1, floor((h - H0)/2) + 1);
c1 = max(1, floor((w - W0)/2) + 1);
r2 = min(h, r1 + H0 - 1);
c2 = min(w, c1 + W0 - 1);
crop = maskIn(r1:r2, c1:c2);

maskOut = crop;
[hc,wc] = size(maskOut);
if hc < H0 || wc < W0
    top     = floor((H0-hc)/2);
    bottom = H0-hc-top;
    left    = floor((W0-wc)/2);
    right  = W0-wc-left;
    maskOut = padarray(maskOut, [top left], false, 'pre');
    maskOut = padarray(maskOut, [bottom right], false, 'post');
end
maskOut = maskOut(1:H0, 1:W0);
end
