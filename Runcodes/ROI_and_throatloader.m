%% ------------------------------------------------------------
% Interactive ROI and Throat Selection Tool (Video or Image)
% ------------------------------------------------------------
clear; clc; close all;
fprintf('--- ROI and Throat Selection Tool ---\n');

%% --- Step 1: Select file (video OR image) ---
[fileName, filePath] = uigetfile({ ...
    '*.avi;*.mp4;*.mov;*.mkv;*.png;*.jpg;*.jpeg;*.tif;*.bmp', ...
    'Video or Image Files (*.avi, *.mp4, *.mov, *.mkv, *.png, *.jpg, *.jpeg, *.tif, *.bmp)'}, ...
    'Select a Video or Image File');

if isequal(fileName, 0)
    error('No file selected. Exiting.');
end

fullPath = fullfile(filePath, fileName);
[~, baseName, ext] = fileparts(fileName);
fprintf('Selected file: %s\n', fullPath);

%% --- Step 2: Create save folder ---
saveFolder = fullfile(filePath, baseName);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
    fprintf('Created folder: %s\n', saveFolder);
else
    fprintf('Folder already exists: %s\n', saveFolder);
end

%% --- Step 3: Load first frame or image ---
isVideo = ismember(lower(ext), {'.avi','.mp4','.mov','.mkv'});

if isVideo
    fprintf('File detected as VIDEO.\n');
    v = VideoReader(fullPath);
    frame1 = im2gray(readFrame(v));
else
    fprintf('File detected as IMAGE.\n');
    img = imread(fullPath);
    if size(img,3) == 3
        frame1 = rgb2gray(img);
    else
        frame1 = img;
    end
end

%% Display image for ROI selection
figure('Name', 'ROI Selection', 'NumberTitle', 'off');
imshow(frame1, []);
title('Select ROI using polygon (double-click to close)');
fprintf('Draw your ROI (double-click to confirm)...\n');

%% --- Step 4: ROI selection ---
maskROI = roipoly;
if isempty(maskROI)
    warning('No ROI selected. Using full image.');
    maskROI = true(size(frame1));
else
    fprintf('ROI selected.\n');
end

roiFile = fullfile(saveFolder, [baseName '_ROI.mat']);
save(roiFile, 'maskROI');
fprintf('ROI mask saved as: %s\n', roiFile);

%% --- Step 5: Select throat location ---
figure('Name', 'Select Throat Location', 'NumberTitle', 'off');
imshow(frame1, []);
title('Click to mark throat location, then press Enter');
fprintf('Click once to mark throat location, then press Enter.\n');

[x_throat, y_throat] = ginput(1);
hold on;
plot(x_throat, y_throat, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
text(x_throat + 10, y_throat, 'Throat', 'Color', 'r', 'FontWeight', 'bold');

throatFile = fullfile(saveFolder, [baseName '_throat.mat']);
save(throatFile, 'x_throat', 'y_throat');
fprintf('Throat location saved as: %s\n', throatFile);

%% --- Step 6: Completion ---
fprintf('\n--- Selection Complete ---\n');
fprintf('Saved ROI: %s\n', roiFile);
fprintf('Saved Throat Location: %s\n', throatFile);
fprintf('Data saved in folder: %s\n', saveFolder);


