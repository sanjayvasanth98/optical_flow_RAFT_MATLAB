%% ========================================================================
%  analyze_roughness_cavitation.m
%  Cluster-safe postprocessor for roughness vs cavitation inception
%
%  Assumes each velocity .mat (v7.3) has at least:
%    u_all, v_all, velMag_all          [Ny x Nx x Nt]  (m/s, physical orientation)
%    u_mean_phys, v_mean_phys, velMean_phys  [Ny x Nx]
%    x_mm, y_mm                        (monotone axes in mm)
%    fps                               (scalar)
%
%  Assumes each throat .mat has any of:
%    throatx_mm, throaty_mm
%    OR throatx_pix, throaty_pix
%    OR throatx_idx, throaty_idx
%
%  OUTPUT (all into outDir):
%    - Per-case |V| mean-field plots with streamlines + throat (x & y) lines
%    - Shear strength metrics in cavitation ROI vs roughness
%    - Low-speed area fraction in cavitation ROI vs roughness
%    - Mean |V| in cavitation ROI vs roughness
%    - Spectra at a probe point inside cavitation ROI for all cases
%
%  Run on cluster, e.g.:
%    matlab -nodisplay -nosplash -nodesktop -r "run('analyze_roughness_cavitation.m'); exit"
%
%  Author: your ruthless mentor (Nov 2025)
%% ========================================================================
clear; clc; close all;


% --- Load your personal toolbox path safely ---
userPathFile = fullfile(getenv('HOME'), 'matlab', 'pathdef.m');
if isfile(userPathFile)
    addpath(genpath(fileparts(userPathFile)));  % add folder containing pathdef.m
    run(userPathFile);                          % execute your saved pathdef.m
    fprintf('Loaded custom MATLAB path: %s\n', userPathFile);
else
    warning('Custom pathdef.m not found. Using default MATLAB path.');
end

%% ============================ CONFIG ====================================
%  EDIT THIS BLOCK

% Velocity .mat files (full paths)
velMatFiles = {
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_pressurized_smooth/Pressurized_smooth_velocity.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_P10S100_pressurized/bgsub_side_viewP10S100_velocity.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S70_pressurized/nocav_pressurized_sideviewp10s70_velocity.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S50_pressurized/Pressurized_p10s50_velocity.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S30_pressurized/Pressurized p10s30_velocity.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S20_pressurized/pressurized p10s20final_velocity.mat"
};

% Matching throat-location .mat files (same order)
throatMatFiles = {
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_pressurized_smooth/Pressurized_smooth_throat.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_P10S100_pressurized/AVG_side_viewP10S100-1_throat.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S70_pressurized/nocav_pressurized_sideviewp10s70_throat.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S50_pressurized/Pressurized_p10s50_throat.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S30_pressurized/Pressurized p10s30_throat.mat"
    "/home/kbsanjayvasanth/VelocitydataRAFT/Nocav_S20_pressurized/pressurized p10s20final_throat.mat"
};

% Case labels for plots (LaTeX/TeX-friendly)
caseLabels = {
    'Smooth'
    'SA 26.7~\mum'
    'SA 48.3~\mum'
    'SA 55~\mum'
    'SA 65~\mum'
    'SA 80~\mum'
};

% Roughness values in microns (for x-axis in trend plots)
roughness_um = [0, 26.7, 48.3, 55, 65, 80];

% Cavitation region of interest (ROI) in mm (based on your observation)
ROI_x_mm = [1.0, 1.5];   % streamwise
ROI_y_mm = [0.2, 0.6];   % wall-normal

% Probe point (mm) inside cavitation-hotspot region for spectra
probePoint_mm = [1.25, 0.5];   % [x, y]

% Low-speed threshold factor (fraction of GLOBAL max |V| across all cases)
lowVelFrac = 0.20;   % 20% of global max

% Output root folder
outBase = "/home/kbsanjayvasanth/VelocitydataRAFT/Processed_results/plots";
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outBase, ['Results_' timestamp]);
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('Output directory: %s\n', outDir);

% Figure export resolution
dpi = 300;

%% ========================== GLOBAL STYLE ================================
set(groot, 'defaultFigureVisible','off', ...
           'defaultAxesFontName','Times New Roman', ...
           'defaultTextFontName','Times New Roman', ...
           'defaultAxesFontSize',12, ...
           'defaultLineLineWidth',1.5, ...
           'defaultAxesLineWidth',1.0);

%% ====================== BASIC CHECKS ====================================
nC = numel(velMatFiles);
assert(nC == numel(throatMatFiles), 'velMatFiles and throatMatFiles must match length.');
assert(nC == numel(caseLabels),     'caseLabels must match number of files.');
if ~isempty(roughness_um)
    assert(numel(roughness_um) == nC, 'roughness_um length mismatch.');
end

%% ====================== PASS 1: LOAD MEANS, AXES, THROATS ==============
X = cell(nC,1); Y = cell(nC,1);
Umean = cell(nC,1); Vmean = cell(nC,1); Qmean = cell(nC,1);
fps   = zeros(nC,1);
throat = struct('x_mm',nan(nC,1),'y_mm',nan(nC,1));
globMaxVel = -inf;

fprintf('Pass 1: loading mean fields, axes, throat locations...\n');
for i = 1:nC
    mf = matfile(velMatFiles{i});
    mustHave(mf, {'u_mean_phys','v_mean_phys','velMean_phys','x_mm','y_mm','fps'}, velMatFiles{i});

    Umean{i} = mf.u_mean_phys;
    Vmean{i} = mf.v_mean_phys;
    Qmean{i} = mf.velMean_phys;
    x_mm = mf.x_mm; y_mm = mf.y_mm;
    fps(i) = mf.fps;

    [X{i}, Y{i}] = meshgrid(x_mm, y_mm);

    globMaxVel = max(globMaxVel, max(Qmean{i}(:),[],'omitnan'));

    % Resolve throat x,y for this case
    throat(i) = resolveThroat(throatMatFiles{i}, x_mm, y_mm);
end

if numel(unique(fps)) > 1
    warning('Different fps across cases; spectra will use each case''s own fs.');
end

%% ====================== PER-CASE |V| + STREAMLINES ======================
fprintf('Generating per-case mean |V| plots with streamlines only...\n');

for i = 1:nC
    xi = X{i}(1,:); 
    yi = Y{i}(:,1);

    % Force 2-D fields (Ny x Nx) in case they are Ny x Nx x 1 in the MAT files
    Ui = squeeze(Umean{i});
    Vi = squeeze(Vmean{i});
    Qi = squeeze(Qmean{i});

    [XX,YY] = meshgrid(xi, yi);

    f = figure('Position',[100 100 900 700]);
    imagesc(xi, yi, Qi);
    set(gca,'YDir','normal'); axis image;
    colormap(turbo); caxis([0, globMaxVel]);
    cb = colorbar; cb.Label.String = '|V| (m/s)';

    hold on;

    % 2D streamlines: use density scalar
    density = 0.8;  % tweak 0.5–3.0 for fewer/more arrows
    hss = streamslice(XX, YY, Ui, Vi, density);
    set(hss, 'Color','k', 'LineWidth', 0.6);

    title(sprintf('%s: Time-averaged |V|', caseLabels{i}), ...
          'Interpreter','tex');
    xlabel('x (mm)');
    ylabel('y (mm)');
    box on;

    fname = fullfile(outDir, sprintf('%02d_%s_VelMean_Streamlines.png', ...
                    i, safen(caseLabels{i})));
    exportgraphics(f, fname, 'Resolution', dpi);
    close(f);
end


%% ====================== PASS 2: SHEAR + LOW-SPEED METRICS ==============
fprintf('Pass 2: computing shear & low-speed metrics in cavitation ROI...\n');

shearMean_ROI = zeros(nC,1);
shearMax_ROI  = zeros(nC,1);
lowFrac_ROI   = zeros(nC,1);
meanV_ROI     = zeros(nC,1);

lowVelThresh = lowVelFrac * globMaxVel;   % global threshold

for i = 1:nC
    xi = X{i}(1,:); yi = Y{i}(:,1);
    Ui = Umean{i}; Vi = Vmean{i}; Qi = Qmean{i};

    dx_mm = mean(diff(xi));
    dy_mm = mean(diff(yi));
    dx = dx_mm / 1000;   % convert to meters for 1/s units
    dy = dy_mm / 1000;

    [Ux, Uy] = gradient(Ui, dx, dy);
    [Vx, ~ ] = gradient(Vi, dx, dy);

    % shear magnitude (simple proxy)
    S_mag = sqrt(Uy.^2 + Vx.^2);   % [1/s]

    % Cavitation ROI mask in this grid
        xMask = (xi >= ROI_x_mm(1)) & (xi <= ROI_x_mm(2));
    yMask = (yi >= ROI_y_mm(1)) & (yi <= ROI_y_mm(2));
    if ~any(xMask) || ~any(yMask)
        error('ROI outside domain for case %d (%s).', i, velMatFiles{i});
    end

    % Build ROI mask same size as Qi
    ROI = false(size(Qi));      % Ny x Nx
    ROI(yMask, xMask) = true;

    % Shear stats in ROI
    S_ROI = S_mag(ROI);
    shearMean_ROI(i) = mean(S_ROI,'omitnan');
    shearMax_ROI(i)  = max(S_ROI,[],'omitnan');

    % Low-speed area fraction in ROI
    Q_ROI = Qi(ROI);

    lowMask = (Q_ROI <= lowVelThresh);
    lowFrac_ROI(i) = mean(lowMask,'omitnan');

    % Mean |V| in ROI
    meanV_ROI(i) = mean(Q_ROI,'omitnan');
end

%% ====================== PLOT: SHEAR vs ROUGHNESS ========================
fprintf('Generating shear comparison plots...\n');

f = figure('Position',[100 100 1100 450]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

if ~isempty(roughness_um)
    xvals = roughness_um;
    xlab  = 'k_s (\mum)';
else
    xvals = 1:nC;
    xlab  = 'Case index';
end

% Mean shear
nexttile;
plot(xvals, shearMean_ROI, '-o', 'LineWidth',1.8, 'MarkerSize',7);
grid on; box on;
xlabel(xlab, 'Interpreter','tex');
ylabel('<|S|>_{ROI} (1/s)', 'Interpreter','tex');
title('Mean shear in cavitation region', 'Interpreter','tex');
if isempty(roughness_um)
    xticks(1:nC); xticklabels(caseLabels);
end

% Max shear
nexttile;
plot(xvals, shearMax_ROI, '-o', 'LineWidth',1.8, 'MarkerSize',7);
grid on; box on;
xlabel(xlab, 'Interpreter','tex');
ylabel('max(|S|)_{ROI} (1/s)', 'Interpreter','tex');
title('Max shear in cavitation region', 'Interpreter','tex');
if isempty(roughness_um)
    xticks(1:nC); xticklabels(caseLabels);
end

fname = fullfile(outDir, 'Shear_ROI_vs_Roughness.png');
exportgraphics(f, fname, 'Resolution', dpi);
close(f);

%% ====================== PLOT: LOW-SPEED FRACTION vs ROUGHNESS ===========
fprintf('Generating low-speed area fraction plots...\n');

f = figure('Position',[100 100 900 450]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Low-speed fraction
nexttile;
plot(xvals, lowFrac_ROI, '-s', 'LineWidth',1.8, 'MarkerSize',7);
grid on; box on;
xlabel(xlab, 'Interpreter','tex');
ylabel(sprintf('Area frac(|V| \\le %.2f m/s) in ROI', lowVelThresh), 'Interpreter','tex');
title('Low-speed area fraction in cavitation region', 'Interpreter','tex');
if isempty(roughness_um)
    xticks(1:nC); xticklabels(caseLabels);
end

% Mean |V| in ROI
nexttile;
plot(xvals, meanV_ROI, '-d', 'LineWidth',1.8, 'MarkerSize',7);
grid on; box on;
xlabel(xlab, 'Interpreter','tex');
ylabel('<|V|>_{ROI} (m/s)', 'Interpreter','tex');
title('Mean speed in cavitation region', 'Interpreter','tex');
if isempty(roughness_um)
    xticks(1:nC); xticklabels(caseLabels);
end

fname = fullfile(outDir, 'LowSpeed_ROI_vs_Roughness.png');
exportgraphics(f, fname, 'Resolution', dpi);
close(f);

%% -----------PASS 3: spectra at ONE point relative to throat (pixel-aware) ---
fprintf('Pass 3: computing spectra at a single point relative to the throat...\n');

PSD_all = cell(nC,1);
f_all   = cell(nC,1);

% Relative probe location (in mm) from the throat
% +x = downstream, +y = away from wall (adjust dy_rel sign if needed)
dx_rel = 1.25;   % mm downstream from throat
dy_rel = 0.50;   % mm above throat

for i = 1:nC
    xi = X{i}(1,:);      % x in mm
    yi = Y{i}(:,1);      % y in mm

    % ---- Load throat info for this case (in PIXELS) ----
    Sth = load(throatMatFiles{i});
    ix_th = []; iy_th = [];

    if isfield(Sth,'throatx_pix')
        ix_th = round(Sth.throatx_pix);
    elseif isfield(Sth,'throatx_idx')
        ix_th = round(Sth.throatx_idx);
    elseif isfield(Sth,'throatx_mm')
        % fallback: map mm -> index
        [~, ix_th] = min(abs(xi - Sth.throatx_mm));
    end

    if isfield(Sth,'throaty_pix')
        iy_th = round(Sth.throaty_pix);
    elseif isfield(Sth,'throaty_idx')
        iy_th = round(Sth.throaty_idx);
    elseif isfield(Sth,'throaty_mm')
        [~, iy_th] = min(abs(yi - Sth.throaty_mm));
    end

    % If still empty, fall back to center of domain (shouldn't happen ideally)
    if isempty(ix_th), ix_th = round(numel(xi)/2); end
    if isempty(iy_th), iy_th = round(numel(yi)/2); end

    % Clamp to valid indices
    ix_th = max(1, min(numel(xi), ix_th));
    iy_th = max(1, min(numel(yi), iy_th));

    % ---- Convert throat pixel to physical location (mm) ----
    throat_x_mm = xi(ix_th);
    throat_y_mm = yi(iy_th);

    % ---- Target probe location in physical coordinates ----
    x_target = throat_x_mm + dx_rel;
    y_target = throat_y_mm + dy_rel;

    % ---- Find nearest grid indices to this physical target ----
    [~, ix] = min(abs(xi - x_target));
    [~, iy] = min(abs(yi - y_target));

    % Final safety clamp
    Nx = numel(xi);
    Ny = numel(yi);
    ix = max(1, min(Nx, ix));
    iy = max(1, min(Ny, iy));

    fprintf(['  Case %d (%s): throat at (pix) = (%d,%d) -> (mm) = (%.4f, %.4f);\n' ...
             '                         probe at (mm) = (%.4f, %.4f) -> indices (ix,iy)=(%d,%d)\n'], ...
            i, caseLabels{i}, ix_th, iy_th, throat_x_mm, throat_y_mm, ...
            xi(ix), yi(iy), ix, iy);

    % ---- Load time-resolved |V| ----
    S = load(velMatFiles{i}, 'velMag_all');
    velMag_all = S.velMag_all;   % Ny x Nx x Nt
    [Ny_v, Nx_v, Nt] = size(velMag_all);

    if ix < 1 || ix > Nx_v || iy < 1 || iy > Ny_v
        warning('Probe index out of bounds for case %d (%s). Skipping.', ...
                i, velMatFiles{i});
        continue;
    end

    % ---- Extract time series at this grid point ----
    ts = squeeze(velMag_all(iy, ix, :));
    ts = ts(:);

    % ---- Basic sanity on the time series ----
    finiteMask = isfinite(ts);
    if nnz(finiteMask) < 32
        warning('Too few finite samples at probe for case %d (%s). Skipping.', ...
                i, velMatFiles{i});
        continue;
    end

    ts = ts(finiteMask);

    % Remove DC component
    ts = ts - mean(ts, 'omitnan');

    sig_ts = std(ts, 'omitnan');
    fprintf('    std(|V|(t)) at probe for %s: %.3e m/s\n', caseLabels{i}, sig_ts);

    fs = fps(i);           % sampling frequency for this case
    N  = numel(ts);

    % Welch parameters
    segLen   = min(N, 4096);
    if segLen < 64
        warning('    Too few samples for good PSD in %s (N=%d). Skipping.', caseLabels{i}, N);
        continue;
    end

    win      = hanning(segLen);
    noverlap = floor(segLen/2);
    nfft     = 2^nextpow2(segLen);

    [Pxx, fvec] = pwelch(ts, win, noverlap, nfft, fs, 'onesided');

    % Enforce strictly positive PSD for log plotting
    Pxx = max(Pxx, realmin);

    PSD_all{i} = Pxx;
    f_all{i}   = fvec;
end

%% -----------Plot comparison (log-y, all cases)------------------------
f = figure('Position',[100 100 900 550], 'Color','w');
ax = axes(f); hold(ax,'on');

specColors = lines(nC);
colororder(ax, specColors);

% Determine global frequency range
fMaxGlobal = 0;
for i = 1:nC
    if isempty(f_all{i}), continue; end
    fMaxGlobal = max(fMaxGlobal, max(f_all{i}));
end
if fMaxGlobal <= 0
    warning('No valid spectra computed – check probe location / data.');
    fMaxGlobal = 1;
end
fMinPlot = 10;                       % remove very low freq/near-DC
fMaxPlot = min(fMaxGlobal, 2e4);     % cap at 20 kHz

for i = 1:nC
    fi = f_all{i};
    Pi = PSD_all{i};
    if isempty(fi) || isempty(Pi), continue; end

    idx = (fi >= fMinPlot) & (fi <= fMaxPlot);
    if ~any(idx), continue; end

    semilogy(ax, fi(idx), Pi(idx), ...
             'LineWidth', 1.5, 'Color', specColors(i,:));
end

set(ax, 'YScale','log', ...
        'FontName','Times New Roman', ...
        'FontSize',12, ...
        'LineWidth',1.2, ...
        'Box','on', ...
        'XMinorTick','on', ...
        'YMinorTick','on');
grid(ax,'on'); grid(ax,'minor');

xlabel(ax, '$f$ (Hz)', 'Interpreter','latex');
ylabel(ax, '$\Phi_{|V|}(f)\;[(\mathrm{m/s})^2/\mathrm{Hz}]$', ...
           'Interpreter','latex');

title(ax, sprintf(['Spectra at probe: $x - x_t = %.2f$ mm, ' ...
                   '$y - y_t = %.2f$ mm'], dx_rel, dy_rel), ...
           'Interpreter','latex');

lg = legend(ax, caseLabels, 'Location','best', ...
            'Interpreter','tex', 'Box','off');
lg.ItemTokenSize = [18 9];

fname = fullfile(outDir, 'Spectra_SinglePoint_vs_Roughness.png');
exportgraphics(f, fname, 'Resolution', dpi);
close(f);



%% ====================== DMD ANALYSIS (velMag_all) =======================
fprintf('=== Performing DMD on velMag_all for all cases ===\n');

% ------------------ user-tunable DMD parameters --------------------------
timeDecim    = 5;      % use every Nth frame in time (reduce size / noise)
maxRank      = 50;     % max DMD rank
energyThresh = 0.999;  % keep enough modes to capture this fraction of energy
nModesSave   = 2;      % <-- ONLY first 2 dominant modes (by amplitude)
nModesEnergy = 20;     % number of modes for energy comparison plot

% ------------------------------------------------------------------------
% Preallocate DMD struct
% ------------------------------------------------------------------------
DMD = struct('lambda',[],'omega',[],'freq',[],'Phi',[], ...
             'b',[],'amp',[],'energy',[], ...
             't',[],'a_dom',[], ...
             'freq_pos',[],'amp_pos',[]);

domTime_t   = cell(nC,1);   % time vector for dominant mode
domTime_amp = cell(nC,1);   % |a_dom(t)| for each case

modeColors = lines(nC);

for i = 1:nC
    fprintf('  DMD for case %d/%d: %s\n', i, nC, caseLabels{i});
    S = load(velMatFiles{i}, 'velMag_all');
    velMag_all = S.velMag_all;        % Ny x Nx x Nt

    [Ny, Nx, Nt] = size(velMag_all);
    if Nt < 3
        warning('Not enough time samples for DMD in case %s. Skipping.', caseLabels{i});
        continue;
    end

    % ----- temporal decimation -----
    t_idx = 1:timeDecim:Nt;
    Nt_sub = numel(t_idx);
    if Nt_sub < 3
        warning('Too few decimated frames for DMD in case %s. Skipping.', caseLabels{i});
        continue;
    end

    % ----- build snapshot matrices X1, X2 (flattened fields) -----
    X1 = zeros(Ny*Nx, Nt_sub-1);
    X2 = zeros(Ny*Nx, Nt_sub-1);
    for k = 1:Nt_sub-1
        frame1 = velMag_all(:,:,t_idx(k));
        frame2 = velMag_all(:,:,t_idx(k+1));
        X1(:,k) = frame1(:);
        X2(:,k) = frame2(:);
    end

    % ----- subtract temporal mean to emphasize dynamics -----
    x_mean = mean(X1, 2);
    X1 = X1 - x_mean;
    X2 = X2 - x_mean;

    % ----- SVD and truncation -----
    [U,Sig,V] = svd(X1, 'econ');
    singVals  = diag(Sig);
    cumEnergy = cumsum(singVals.^2) / sum(singVals.^2);
    r_e = find(cumEnergy >= energyThresh, 1, 'first');
    if isempty(r_e)
        r = min(maxRank, numel(singVals));
    else
        r = min(maxRank, r_e);
    end
    U_r   = U(:,1:r);
    Sig_r = Sig(1:r,1:r);
    V_r   = V(:,1:r);

    % ----- reduced operator and eigen-decomposition -----
    A_tilde = U_r' * X2 * V_r / Sig_r;
    [W,D]   = eig(A_tilde);
    lambda  = diag(D);

    % exact DMD modes
    Phi = X2 * V_r / Sig_r * W;   % Ny*Nx x r

    % ----- continuous-time eigenvalues and frequencies -----
    dt_eff = timeDecim / fps(i);        % effective time-step
    omega  = log(lambda) / dt_eff;      % complex growth rates
    freq   = imag(omega) / (2*pi);      % Hz

    % ----- mode amplitudes -----
    x1 = X1(:,1);               % first snapshot (mean-subtracted)
    b  = Phi \ x1;              % modal amplitudes (least-squares)
    amp = abs(b);               % |amplitude|
    energy = amp.^2 / sum(amp.^2 + eps);  % normalized "energy" per mode

    % store in struct
    DMD(i).lambda = lambda;
    DMD(i).omega  = omega;
    DMD(i).freq   = freq;
    DMD(i).Phi    = Phi;
    DMD(i).b      = b;
    DMD(i).amp    = amp;
    DMD(i).energy = energy;

    % ----- TIME EVOLUTION of dominant mode amplitude -----
    [~, idx_sorted] = sort(amp, 'descend');
    dom_idx = idx_sorted(1);

    t_vec = (0:Nt_sub-1).' * dt_eff;            % time vector (s)
    a_dom = b(dom_idx) * exp(omega(dom_idx) * t_vec);  % time evolution
    domTime_t{i}   = t_vec;
    domTime_amp{i} = abs(a_dom);

    DMD(i).t     = t_vec;
    DMD(i).a_dom = a_dom;

    % ================= 1) EIGENVALUE SPECTRUM (PER CASE) =================
    f1 = figure('Position',[100 100 700 650],'Color','w');
    ax1 = axes(f1); hold(ax1,'on');
    theta = linspace(0,2*pi,400);
    plot(ax1, cos(theta), sin(theta), 'k--', 'LineWidth',1.0);  % unit circle

    % color by amplitude (stronger = darker)
    scatter(ax1, real(lambda), imag(lambda), 40, amp, 'filled');
    colormap(ax1, turbo); cb = colorbar(ax1);
    cb.Label.String = '|b_k| (mode amplitude)';

    axis(ax1,'equal');
    xlim(ax1,[-1.2 1.2]); ylim(ax1,[-1.2 1.2]);
    grid(ax1,'on'); box(ax1,'on');
    set(ax1, 'FontName','Times New Roman', 'FontSize',12, ...
             'XMinorTick','on', 'YMinorTick','on');

    xlabel(ax1, '$\Re(\lambda_k)$', 'Interpreter','latex');
    ylabel(ax1, '$\Im(\lambda_k)$', 'Interpreter','latex');
    title(ax1, sprintf('DMD Eigenvalues: %s', caseLabels{i}), 'Interpreter','tex');

    fname = fullfile(outDir, sprintf('%02d_%s_DMD_EigSpectrum.png', ...
                                     i, safen(caseLabels{i})));
    exportgraphics(f1, fname, 'Resolution', 300);
    close(f1);

    % ========== STORE DMD "SPECTRUM": AMPLITUDE vs FREQUENCY (GLOBAL) ====
    % keep non-negative frequencies
    posIdx  = freq >= 0;
    f_pos   = freq(posIdx);
    amp_pos = amp(posIdx);

    % sort by frequency
    [f_pos, sortIdx] = sort(f_pos);
    amp_pos = amp_pos(sortIdx);

    % store for global comparison plot
    DMD(i).freq_pos = f_pos;
    DMD(i).amp_pos  = amp_pos;

    % ========= 3) SPATIAL STRUCTURE OF TOP 2 MODES (PER CASE) ============
    [~, idx_sorted] = sort(amp, 'descend');
    nShow = min(nModesSave, numel(idx_sorted));

    xi = X{i}(1,:);
    yi = Y{i}(:,1);
    Ui_mean = Umean{i};
    Vi_mean = Vmean{i};

    f3 = figure('Position',[100 100 1200 550],'Color','w');
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    for m = 1:nShow
        km = idx_sorted(m);

        % mode field (real part) reshaped back to Ny x Nx
        modeField = reshape(real(Phi(:,km)), Ny, Nx);

        nexttile;
        imagesc(xi, yi, modeField);
        set(gca,'YDir','normal'); axis image;
        colormap(redblue);  % your red-blue colormap

        cmax = max(abs(modeField(:)));
        if cmax == 0, cmax = 1; end
        caxis([-cmax cmax]);
        colorbar;

        hold on; box on;
        set(gca,'FontName','Times New Roman','FontSize',11, ...
                 'XMinorTick','on','YMinorTick','on');

        % Overlay coarsened mean velocity vectors
        stepX = max(1, round(Nx/25));
        stepY = max(1, round(Ny/25));
        xs = 1:stepX:Nx;
        ys = 1:stepY:Ny;
        [XS,YS] = meshgrid(xs, ys);
        quiver(xi(xs), yi(ys), ...
               Ui_mean(ys,xs), Vi_mean(ys,xs), ...
               1.5, 'k');  % scale factor 1.5; adjust if needed

        title(sprintf('Mode %d (rank %d): f = %.1f Hz, |b|=%.2g', ...
              m, km, freq(km), amp(km)), 'Interpreter','tex');
        xlabel('x (mm)');
        ylabel('y (mm)');
    end

    title(tl, sprintf('DMD Spatial Modes 1 & 2: %s', caseLabels{i}), ...
          'Interpreter','tex','FontWeight','bold');

    fname = fullfile(outDir, sprintf('%02d_%s_DMD_SpatialModes_1and2.png', ...
                                     i, safen(caseLabels{i})));
    exportgraphics(f3, fname, 'Resolution', 300);
    close(f3);
end

%% ========== 2b) GLOBAL DMD AMPLITUDE SPECTRUM (ALL CASES) ===============
fprintf('Plotting global DMD amplitude spectrum (all cases)...\n');

% Legend labels for DMD amplitude spectrum based on S_a
Sa_labels = { ...
    'S_a = 5~\mum',  ...   % smooth
    'S_a = 12~\mum', ...
    'S_a = 20~\mum', ...
    'S_a = 30~\mum', ...
    'S_a = 53~\mum', ...
    'S_a = 80~\mum'  ...
};
if numel(Sa_labels) ~= nC
    warning('Sa_labels length (%d) != nC (%d). Update Sa_labels.', ...
        numel(Sa_labels), nC);
end

fAmp = figure('Position',[100 100 900 550],'Color','w');
axAmp = axes(fAmp); hold(axAmp,'on');

specColors = lines(nC);
colororder(axAmp, specColors);

% determine global max frequency to display
fMaxGlobal = 0;
for i = 1:nC
    if ~isfield(DMD(i),'freq_pos') || isempty(DMD(i).freq_pos), continue; end
    fMaxGlobal = max(fMaxGlobal, max(DMD(i).freq_pos));
end
if fMaxGlobal > 0
    fMaxPlot = fMaxGlobal * 1.05;
else
    warning('No valid DMD spectra found for global amplitude plot.');
    fMaxPlot = 0;
end

% --- SCATTER (no connecting lines) for each case ---
for i = 1:nC
    if ~isfield(DMD(i),'freq_pos') || isempty(DMD(i).freq_pos), continue; end

    f_pos   = DMD(i).freq_pos;
    amp_pos = DMD(i).amp_pos;

    if fMaxPlot > 0
        idx = f_pos >= 0 & f_pos <= fMaxPlot;
        f_plot = f_pos(idx);
        a_plot = amp_pos(idx);
    else
        f_plot = f_pos;
        a_plot = amp_pos;
    end

    if isempty(f_plot), continue; end

    scatter(axAmp, f_plot, a_plot, 28, specColors(i,:), ...
            'filled', 'MarkerFaceAlpha', 0.75);
end

grid(axAmp,'on'); grid(axAmp,'minor'); box(axAmp,'on');
set(axAmp, 'FontName','Times New Roman', ...
           'FontSize',12, ...
           'LineWidth',1.2, ...
           'XMinorTick','on', ...
           'YMinorTick','on');

xlabel(axAmp, '$f_k~\mathrm{(Hz)}$', 'Interpreter','latex');
ylabel(axAmp, '$|b_k|$',           'Interpreter','latex');
title(axAmp, 'DMD amplitude spectrum: comparison across roughness', ...
      'Interpreter','tex');

if fMaxPlot > 0
    xlim(axAmp, [0, fMaxPlot]);
end

if numel(Sa_labels) == nC
    lgAmp = legend(axAmp, Sa_labels, ...
                   'Location','northeastoutside', ...
                   'Interpreter','tex', ...
                   'Box','off');
    lgAmp.ItemTokenSize = [18 9];
else
    legend(axAmp, caseLabels, 'Location','northeastoutside', ...
           'Interpreter','tex', 'Box','off');
end

fname = fullfile(outDir, 'DMD_AmplitudeSpectrum_AllCases.png');
exportgraphics(fAmp, fname, 'Resolution', 300);
close(fAmp);

%% ========== 4) TIME EVOLUTION OF DOMINANT MODE (ALL CASES) ==============
fprintf('Plotting time evolution of dominant DMD mode amplitudes...\n');

f4 = figure('Position',[100 100 900 550],'Color','w');
ax4 = axes(f4); hold(ax4,'on');

for i = 1:nC
    if isempty(DMD(i).t) || isempty(DMD(i).a_dom), continue; end
    plot(ax4, DMD(i).t, abs(DMD(i).a_dom), ...
         'LineWidth',1.6, 'Color', modeColors(i,:));
end

grid(ax4,'on'); grid(ax4,'minor'); box(ax4,'on');
set(ax4,'FontName','Times New Roman','FontSize',12, ...
        'XMinorTick','on','YMinorTick','on');

xlabel(ax4, '$t$~(s)', 'Interpreter','latex');
ylabel(ax4, '$|a_{\mathrm{dom}}(t)|$', 'Interpreter','latex');
title(ax4, 'Time Evolution of Dominant DMD Mode Amplitude', 'Interpreter','tex');
legend(ax4, caseLabels, 'Location','best', 'Interpreter','tex', 'Box','off');

fname = fullfile(outDir, 'DMD_DominantMode_TimeEvolution_AllCases.png');
exportgraphics(f4, fname, 'Resolution', 300);
close(f4);

%% ========== 5) ENERGY OF FIRST nModesEnergy MODES (ALL CASES) ===========
fprintf('Plotting DMD energy content of dominant modes...\n');

f5 = figure('Position',[100 100 900 550],'Color','w');
ax5 = axes(f5); hold(ax5,'on');

for i = 1:nC
    if isempty(DMD(i).energy), continue; end
    [energy_sorted, ~] = sort(DMD(i).energy, 'descend');
    nShow = min(nModesEnergy, numel(energy_sorted));
    plot(ax5, 1:nShow, energy_sorted(1:nShow), ...
         '-o', 'LineWidth',1.6, 'MarkerSize',5, ...
         'Color', modeColors(i,:));
end

grid(ax5,'on'); box(ax5,'on');
set(ax5,'FontName','Times New Roman','FontSize',12, ...
        'XMinorTick','on','YMinorTick','on');

xlabel(ax5, 'Mode rank (sorted by energy)', 'Interpreter','tex');
ylabel(ax5, 'Energy fraction',            'Interpreter','tex');
title(ax5, sprintf('Energy of first %d DMD modes', nModesEnergy), ...
      'Interpreter','tex');
legend(ax5, caseLabels, 'Location','northeast', 'Interpreter','tex', 'Box','off');

fname = fullfile(outDir, 'DMD_ModeEnergy_First20_AllCases.png');
exportgraphics(f5, fname, 'Resolution', 300);
close(f5);

fprintf('=== DMD analysis complete. ===\n');

%% ===================== HELPER: REDBLUE COLORMAP =========================
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];
end




%% ========================== HELPERS =====================================
function mustHave(mf, names, file)
for k = 1:numel(names)
    if ~ismember(names{k}, who(mf))
        error('Variable %s missing in %s', names{k}, file);
    end
end
end

function r = resolveThroat(throatFile, x_mm, y_mm)
% Resolve throat location to physical (x_mm, y_mm)
S = load(throatFile);
r.x_mm = NaN; r.y_mm = NaN;

if isfield(S,'throatx_mm'), r.x_mm = S.throatx_mm; end
if isfield(S,'throaty_mm'), r.y_mm = S.throaty_mm; end

if isnan(r.x_mm) && isfield(S,'throatx_pix')
    % assume throatx_pix is 1-based img coordinate mapping to x_mm grid
    idx = min(max(round(S.throatx_pix),1), numel(x_mm));
    r.x_mm = x_mm(idx);
end
if isnan(r.y_mm) && isfield(S,'throaty_pix')
    idy = min(max(round(S.throaty_pix),1), numel(y_mm));
    r.y_mm = y_mm(idy);
end

if isnan(r.x_mm) && isfield(S,'throatx_idx')
    idx = min(max(S.throatx_idx,1), numel(x_mm));
    r.x_mm = x_mm(idx);
end
if isnan(r.y_mm) && isfield(S,'throaty_idx')
    idy = min(max(S.throaty_idx,1), numel(y_mm));
    r.y_mm = y_mm(idy);
end

% Fall back to something sane if missing
if isnan(r.x_mm)
    r.x_mm = x_mm(round(numel(x_mm)/2));
end
if isnan(r.y_mm)
    r.y_mm = y_mm(round(numel(y_mm)/2));
end
end

function s = safen(str)
% Make a label file-name-safe
s = regexprep(str,'[^\w\-]+','_');
end
