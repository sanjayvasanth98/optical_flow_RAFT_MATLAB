# ROI-only instantaneous velocity storage

New `<video>_velocity.mat` files from `raftmatlabsideview_ARC.m` store
`u_all` and `v_all` as **`nnz(maskROI) × numFlowFrames` single arrays**.
Each column is one flow frame; rows follow MATLAB's `find(maskROI)` order.
There are no stored velocity values for pixels outside the ROI.

The file retains `maskROI` and identifies this layout with
`instantaneousStorage = 'roi_pixels_by_frame_v1'`. Values are in m/s,
in image coordinates, just as in the previous full-frame format.
Only frames through `completedFlowFrames` are committed; preallocated
columns after that counter must not be used.

Reconstruct a frame without loading the entire recording:

```matlab
M = matfile('my_video_velocity.mat');
mask = logical(M.maskROI);
k = 1;
assert(k >= 1 && k <= M.completedFlowFrames);
u = nan(size(mask), 'single');
v = nan(size(mask), 'single');
u(mask) = M.u_all(:,k);
v(mask) = M.v_all(:,k);

% Optional: physical orientation, with the origin at the lower left.
u_phys = flipud(u);
v_phys = -flipud(v);
speed = hypot(u_phys, v_phys);
```

Use `zeros` instead of `nan` to reconstruct the previous zero-filled
exterior. Time-averaged fields, profiles, and plots retain their existing
format. `data_analysis/RAFT_Data_analysis.m` reads both storage layouts.
Other custom readers using `u_all(:,:,k)` must reconstruct packed frames
as above. The older `instantaneousfields_RAFT.m` also requires a separately
stored `velMag_all` and is not a reader for this ARC output format.

Existing full-frame velocity files resume in their original layout and
are not converted automatically. To use compact storage for a fresh run,
move the existing velocity file aside first. The uncompressed instantaneous
payload becomes `nnz(maskROI)/numel(maskROI)` of the full-grid payload;
actual MAT-file size savings depend on compression.
