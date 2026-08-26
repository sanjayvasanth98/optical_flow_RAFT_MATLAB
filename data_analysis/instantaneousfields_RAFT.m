%% ========================================================================
%  make_instantaneous_vel_movies.m
%  Cluster-safe postprocessor for instantaneous velocity movies
%
%  For each velocity .mat (v7.3) with:
%    u_all, v_all, velMag_all       [Ny x Nx x Nt]
%    velMean_phys                   [Ny x Nx]  (for color scaling)
%    x_mm, y_mm                     (axes in mm)
%    fps                            (scalar)
%
%  This script:
%    - Finds global max(|V|) from velMean_phys across all cases
%    - For each case:
%         * extracts first 300 frames (or fewer if Nt < 300)
%         * writes a Motion JPEG .avi movie with:
%               imagesc(|V|) + turbo colormap + consistent caxis
%               quiver vectors overlaid (instantaneous u,v)
%         * saves a .mat file with u, v, |V| for those frames
%
%  Run on cluster, e.g.:
%    matlab -nodisplay -nosplash -nodesktop -batch "run('make_instantaneous_vel_movies.m')"
%
%  Author: ChatGPT helper for Sanjay (Nov 2025)
%% ========================================================================
clear; clc; close all;

% --- Load your personal toolbox path safely (optional) ---
userPathFile = fullfile(getenv('HOME'), 'matlab', 'pathdef.m');
if isfile(userPathFile)
    addpath(genpath(fileparts(userPathFile)));
    run(userPathFile);
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

% Case labels for plots (LaTeX/TeX-friendly)
caseLabels = {
    'Smooth'
    'SA 26.7~\mum'
    'SA 48.3~\mum'
    'SA 55~\mum'
    'SA 65~\mum'
    'SA 80~\mum'
};

% Number of instantaneous frames to export (max)
maxFramesToUse = 300;

% Output root folder
outBase = "/home/kbsanjayvasanth/VelocitydataRAFT/Processed_results/instantaneous_movies";
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outBase, ['InstMovies_' timestamp]);
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('Output directory: %s\n', outDir);

% Subfolder for movies + mats
instOutDir = fullfile(outDir, 'Instantaneous_first300');
if ~exist(instOutDir, 'dir'); mkdir(instOutDir); end

% Figure export resolution (if you later add PNG exports)
dpi = 300; %#ok<NASGU>

%% ========================== GLOBAL STYLE ================================
set(groot, 'defaultFigureVisible','off', ...
           'defaultAxesFontName','Times New Roman', ...
           'defaultTextFontName','Times New Roman', ...
           'defaultAxesFontSize',12, ...
           'defaultLineLineWidth',1.5, ...
           'defaultAxesLineWidth',1.0);

%% ====================== BASIC CHECKS ====================================
nC = numel(velMatFiles);
assert(nC == numel(caseLabels), 'caseLabels must match number of velMatFiles.');

%% ====================== PASS 1: LOAD AXES + MEAN FIELDS =================
fprintf('Pass 1: loading axes, fps, and global mean |V| for color scaling...\n');

X      = cell(nC,1);
Y      = cell(nC,1);
fps    = zeros(nC,1);
globMaxVel = -inf;

for i = 1:nC
    mf = matfile(velMatFiles{i});
    mustHave(mf, {'u_all','v_all','velMag_all','velMean_phys','x_mm','y_mm','fps'}, velMatFiles{i});

    x_mm = mf.x_mm;
    y_mm = mf.y_mm;
    [X{i}, Y{i}] = meshgrid(x_mm, y_mm);

    fps(i) = mf.fps;

    Qi_mean = mf.velMean_phys;
    globMaxVel = max(globMaxVel, max(Qi_mean(:),[],'omitnan'));

    fprintf('  Case %d: %s  (fps = %.3f, mean |V| max = %.3f m/s)\n', ...
            i, caseLabels{i}, fps(i), max(Qi_mean(:),[],'omitnan'));
end

if ~isfinite(globMaxVel) || globMaxVel <= 0
    error('Global max(|V|) from mean fields is non-positive. Check velMean_phys fields.');
end

fprintf('Global max of mean |V| across all cases: %.4f m/s\n', globMaxVel);

%% ====================== PASS 2: MAKE MOVIES + MAT SUBSETS ===============
fprintf('Pass 2: exporting first %d instantaneous frames to AVI + MAT...\n', maxFramesToUse);

for i = 1:nC
    fprintf('  Case %d/%d: %s\n', i, nC, caseLabels{i});

    % Load instantaneous data lazily via matfile (cluster-safe)
    mf = matfile(velMatFiles{i});
    mustHave(mf, {'u_all','v_all','velMag_all'}, velMatFiles{i});

    % Get size of u_all: [Ny x Nx x Nt]
    sz_u   = size(mf, 'u_all');
    Ny     = sz_u(1);
    Nx     = sz_u(2);
    if numel(sz_u) < 3
        error('u_all in %s must be Ny x Nx x Nt.', velMatFiles{i});
    end
    Nt_tot = sz_u(3);

    Nt_use = min(maxFramesToUse, Nt_tot);
    if Nt_use < 1
        warning('    No frames found in %s. Skipping.', velMatFiles{i});
        continue;
    end

    % Axes for plotting (from Pass 1)
    xi = X{i}(1,:);   % x in mm
    yi = Y{i}(:,1);   % y in mm

    % Coarsened grid for vectors (same for all frames in this case)
    stepX = max(1, round(Nx/25));
    stepY = max(1, round(Ny/25));
    xs = 1:stepX:Nx;
    ys = 1:stepY:Ny;

    % -------------------- Set up VideoWriter -----------------------------
    vidName = fullfile(instOutDir, ...
        sprintf('%02d_%s_InstVel_first%03d.avi', ...
                i, safen(caseLabels{i}), Nt_use));

    % Use Motion JPEG AVI for good compatibility
    vW = VideoWriter(vidName, 'Motion JPEG AVI');
    vW.FrameRate = fps(i);  % use actual sampling rate
    open(vW);

    % -------------------- Set up Figure & Graphics -----------------------
    f = figure('Position',[100 100 900 700], 'Color','w', 'Visible','off');
    ax = axes('Parent', f);
    hold(ax, 'on');

    % First frame: create graphics objects (imagesc + quiver)
    q0 = mf.velMag_all(:,:,1);
    u0 = mf.u_all(:,:,1);
    v0 = mf.v_all(:,:,1);

    hImg = imagesc(ax, xi, yi, q0);
    set(ax, 'YDir','normal'); axis(ax,'image');
    colormap(ax, turbo);
    caxis(ax, [0, globMaxVel]);
    cb = colorbar(ax);
    cb.Label.String = '|V| (m/s)';

    % Quiver: instantaneous vectors on coarsened grid
    hQuiv = quiver(ax, xi(xs), yi(ys), ...
                       u0(ys,xs), v0(ys,xs), ...
                       1.5, 'k');  % scale factor 1.5 like in DMD plot

    box(ax, 'on');
    set(ax, 'FontName','Times New Roman', ...
            'FontSize',12, ...
            'LineWidth',1.0, ...
            'XMinorTick','on', ...
            'YMinorTick','on');
    xlabel(ax, 'x (mm)');
    ylabel(ax, 'y (mm)');
    title(ax, sprintf('%s: Instantaneous |V| (frame 1/%d)', ...
                      caseLabels{i}, Nt_use), ...
          'Interpreter','tex');

    drawnow;
    frame = getframe(f);
    writeVideo(vW, frame);

    % -------------------- Allocate arrays for MAT subset -----------------
    u_firstN      = zeros(Ny, Nx, Nt_use, 'like', u0);
    v_firstN      = zeros(Ny, Nx, Nt_use, 'like', v0);
    velMag_firstN = zeros(Ny, Nx, Nt_use, 'like', q0);

    u_firstN(:,:,1)      = u0;
    v_firstN(:,:,1)      = v0;
    velMag_firstN(:,:,1) = q0;

    % -------------------- Loop over remaining frames ---------------------
    for k = 2:Nt_use
        if mod(k,50) == 0
            fprintf('    Case %d: frame %d/%d\n', i, k, Nt_use);
        end

        uk = mf.u_all(:,:,k);
        vk = mf.v_all(:,:,k);
        qk = mf.velMag_all(:,:,k);

        u_firstN(:,:,k)      = uk;
        v_firstN(:,:,k)      = vk;
        velMag_firstN(:,:,k) = qk;

        % Update plot objects (no re-creation)
        set(hImg, 'CData', qk);
        set(hQuiv, 'UData', uk(ys,xs), 'VData', vk(ys,xs));

        title(ax, sprintf('%s: Instantaneous |V| (frame %d/%d)', ...
                          caseLabels{i}, k, Nt_use), ...
              'Interpreter','tex');

        drawnow;
        frame = getframe(f);
        writeVideo(vW, frame);
    end

    % Close video and figure
    close(vW);
    close(f);

    fprintf('    Saved AVI: %s\n', vidName);

    % ------------------------ Save MAT subset ----------------------------
    matName = fullfile(instOutDir, ...
        sprintf('%02d_%s_InstVel_first%03d.mat', ...
                i, safen(caseLabels{i}), Nt_use));

    x_mm = xi;        %#ok<NASGU>
    y_mm = yi;        %#ok<NASGU>
    fps_case = fps(i); %#ok<NASGU>

    save(matName, ...
        'u_firstN', 'v_firstN', 'velMag_firstN', ...
        'x_mm', 'y_mm', 'fps_case', ...
        '-v7.3');

    fprintf('    Saved MAT: %s\n', matName);
end

fprintf('All done. Instantaneous movies and MAT subsets in: %s\n', instOutDir);

%% ========================== HELPERS =====================================
function mustHave(mf, names, file)
for k = 1:numel(names)
    if ~ismember(names{k}, who(mf))
        error('Variable %s missing in %s', names{k}, file);
    end
end
end

function s = safen(str)
% Make a label file-name-safe
s = regexprep(str,'[^\w\-]+','_');
end
