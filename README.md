# Optical Flow RAFT — ARC Optimized (MATLAB)

ARC/HPC-friendly MATLAB pipeline to compute **RAFT optical flow** on high-speed videos and produce **calibrated velocity fields (m/s)**, **time-averaged plots**, and **throat-referenced profiles** with robust memory/I/O handling.

## What it does
- Runs **opticalFlowRAFT** (GPU-first, CPU fallback)
- Calibrates velocities using `mm_per_pixel` and `fps`
- Applies **ROI masking** (and ROI-respecting 5-point neighbor averaging)
- Streams results to disk with **chunked MAT writes** (default: every 100 frames)
- Generates publication-ready time-averaged plots:
  - Plot A: `|V|` + throat line
  - Plot B: `|V|` + **ROI-only arrowed streamlines** + dark-blue overlay outside ROI
- Saves ROI diagnostics as binary images

## Requirements
- MATLAB R2023b+ (recommended)
- Computer Vision Toolbox (for `opticalFlowRAFT`)
- Image Processing Toolbox
- Parallel Computing Toolbox (optional, for GPU)

## Inputs
- `videoPath` (AVI)
- ROI mask: `mat files/<video>_ROI.mat` containing `maskROI` (1 = ROI, 0 = outside)
- Throat file: `*_throat.mat` containing `x_throat`
- Calibration: `mm_per_pixel`, `fps`

## Outputs
- `plots/`
  - `<video>_TimeAvgVelMag_ThroatLine.png`
  - `<video>_TimeAvgVelMag_Streamlines_ROI.png`
  - `<video>_ROI_mask.png`, `<video>_ROI_outside.png`
- `mat files/`
  - `<video>_incremental.mat` (chunked frame outputs)
  - `<video>_means.mat` (time-averaged fields + profiles)
  - `<video>_ROI.mat` (persisted ROI)

## How to run
1. Edit `videoPath`, `mm_per_pixel`, `fps` in the script.
2. Ensure `mat files/<video>_ROI.mat` and `*_throat.mat` exist (ROI auto-created full-frame if missing).
3. Run:
   ```matlab
   raftmatlabsideview_ARC
