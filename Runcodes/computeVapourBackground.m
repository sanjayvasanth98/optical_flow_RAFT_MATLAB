function [background, backgroundFrameNumbers, effectiveBgFrames] = computeVapourBackground(video, nBgFrames)
%COMPUTEVAPOURBACKGROUND Median grayscale background sampled over the whole video.
% Intensities are single precision on the 0..255 scale. Optional outputs
% retain the sample frame numbers and effective count for saved metadata.
totalFrames = countVideoFrames(video);
firstRawFrame = read(video, 1);
imageSize = [size(firstRawFrame, 1), size(firstRawFrame, 2)];
effectiveBgFrames = min(nBgFrames, totalFrames);
if effectiveBgFrames == 1
    backgroundFrameNumbers = round((1 + totalFrames) / 2);
else
    backgroundFrameNumbers = round(linspace(1, totalFrames, effectiveBgFrames));
end
backgroundSamples = zeros([imageSize effectiveBgFrames], 'single');
for sampleIndex = 1:effectiveBgFrames
    backgroundSamples(:, :, sampleIndex) = single(255 * im2double(im2gray( ...
        read(video, backgroundFrameNumbers(sampleIndex)))));
end
background = median(backgroundSamples, 3);
end

function totalFrames = countVideoFrames(video)
% Use the video frame count, not a rounded Duration * FrameRate estimate.
totalFrames = [];
try
    totalFrames = double(video.NumFrames);
catch
    % Some readers cannot supply a frame count; count decoded frames.
end
if isempty(totalFrames) || ~isscalar(totalFrames) || ...
        ~isfinite(totalFrames) || totalFrames < 1 || totalFrames ~= fix(totalFrames)
    video.CurrentTime = 0;
    totalFrames = 0;
    while hasFrame(video)
        readFrame(video);
        totalFrames = totalFrames + 1;
    end
    video.CurrentTime = 0;
end
assert(totalFrames > 0, 'Video contains no readable frames.');
end
