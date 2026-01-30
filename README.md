Optical Flow RAFT (ARC-Optimized MATLAB)

MATLAB pipeline to compute RAFT optical flow on high-speed videos with GPU-first execution, ROI-aware masking, and HPC-safe memory/I/O. Produces calibrated velocity fields (m/s), time-averaged magnitude plots (with/without streamlines), throat-referenced streamwise profiles, and chunked .mat outputs for long runs.

Features

Permanent fix for variable decoded frame sizes (forces all frames to reference H0×W0)

ROI-respecting 5-point neighbor averaging (no bleed across ROI boundary)

Chunked MAT writes every writeChunk frames (default 100)

Two plots: (A) |V| + throat line, (B) |V| + ROI-only arrowed streamlines + dark-blue outside-ROI overlay

Batch-safe (no figure popups), progress/ETA printing, optional GIF

Outputs

plots/: PNGs + ROI mask images

mat files/: *_incremental.mat, *_means.mat, *_ROI.mat

Run

Set videoPath, calibration (mm_per_pixel, fps)

Ensure mat files/<video>_ROI.mat and *_throat.mat exist (ROI auto-created as full-frame if missing)

Run: raftmatlabsideview_ARC.m
