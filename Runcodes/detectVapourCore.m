function coreMask = detectVapourCore(rawFrame, background, roiMask, threshold, minArea)
%DETECTVAPOURCORE Detect dark vapour, returning the undilated logical core.
% Signed floating-point subtraction preserves dark-vapour contrast.
grayFrame = single(255 * im2double(im2gray(rawFrame)));
difference = background - grayFrame;
coreMask = (difference > threshold) & roiMask;
coreMask = bwareaopen(coreMask, minArea, 8);
coreMask = imfill(coreMask, 4, 'holes') & roiMask;
end
