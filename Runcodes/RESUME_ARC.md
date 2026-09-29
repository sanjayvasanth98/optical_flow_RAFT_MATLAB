# Resume an interrupted ARC run

Use `Runcodes/raftmatlabsideview_ARC.m` and keep the existing
`<video folder>/mat files/<video name>_velocity.mat` in place.

1. Set `contrastMode = "none"` to recover the original run. Set `videoPath` to the same video. Keep the original calibration, ROI, and RAFT settings.
2. For an old run without a completion counter, set `legacyProgressLog` near the top
   of the script to the interrupted job's `RAFT_SideView_<job id>.out` file.
   Alternatively, set `legacyCompletedFlowFrames` to the exact **N** in the last
   `Processed flow frame N/total (...)` line. Do not calculate N from the rounded 97%.
3. Copy the updated script to the path launched by `Runcodes/raftscript_ARC.sh`
   (or update that launcher to point to the updated script), then resubmit the job.

The script preserves completed velocity slices, reconstructs averages from those
slices, and decodes earlier video frames without running RAFT to restore the exact
frame position. It then processes the remaining pairs and regenerates the mean
fields, profiles, and plots. These preparation steps take some disk-reading time.
If all pairs were saved, only the output-generation steps run.

Future restarts use `completedFlowFrames` in the velocity MAT file automatically.
The counter advances after both components of each flow frame have been written;
an uncommitted pair is recomputed. The existing file is never reinitialized.
For a new video, leave both legacy recovery settings empty. To deliberately start
the same video afresh, move its old velocity file out of the output location first.

The old script did not initialize RAFT with video frame 1, so its first stored
flow slice is zero. Recovery preserves that existing slice. New runs initialize
RAFT with frame 1, and resumed runs initialize it with the last processed image.

Recovery requires a readable velocity MAT file and a log/count from the matching
run. The old file has no video identity or RAFT settings record, so those must be
checked by the operator on its first recovery. New runs record settings and video
size/modification time and reject mismatches. Run only one job per output file.

## CLAHE runs

`contrastMode = "clahe"` is now the default, matching the local runner's
`im2gray` followed by `adapthisteq` with default settings. Both the reference
frame used to initialize RAFT and each subsequent input receive this preprocessing.
All runs save to `mat files/` and `plots/` inside the video folder, regardless
of contrast mode. For a fresh run of the same video, move the old velocity MAT
file out of `mat files/` first. Derived MAT files and plots with matching names
are replaced when the new run generates its outputs.

A new CLAHE run processes the entire video; it cannot reuse the old 97% flow
results. Legacy log/count settings are ignored when creating a CLAHE run.
Interrupted CLAHE runs resume from their own completion counter and saved
contrast metadata. The existing ROI and throat files are still used.
