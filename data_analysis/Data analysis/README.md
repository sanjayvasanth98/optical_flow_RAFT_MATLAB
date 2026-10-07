# Chunked RAFT analysis

Edit the configuration at the top of `RAFT_chunk_analysis.m`, then run it in
MATLAB. Each `matPaths` entry is a case; labels, range rows, and optional
per-case flow rates follow that original order.

## Choose phases and ranges

`analysisPhases` selects the phases for this run, in execution order:

```matlab
analysisPhases = ["pre-inception", "inception", "desinence", "fullcase"];
% Or run only one, even when all other ranges are configured:
analysisPhases = "desinence";
```

Set a separate range matrix for each event phase. Each matrix needs one
`[firstFrame lastFrame]` row per original MAT file. For example, with two cases:

```matlab
frameRangesByPhase.pre_inception = [1 1500; 1 1700];
frameRangesByPhase.inception = [1501 3000; 1701 3200];
frameRangesByPhase.desinence = [3001 4500; 3201 4700];
```

These numbers are examples. The script also includes six-row dummy range
matrices for pre-inception (`1:1999`) and desinence (`5001:6000`); replace them
with your actual ranges before enabling those phases. The existing inception
ranges are retained. Only ranges for cases
used by enabled analyses in requested phases are validated; unused rows can
remain `NaN`. Endpoints are included, and flow frame `k` describes the motion
between video frames `k` and `k+1`.

`fullcase` automatically uses `1:completedFlowFrames` for each analyzed case.
It excludes unwritten, preallocated frames and requires no manual range.
For older files without a completion counter, the existing
`assumeCompleteWhenNoCounter` option controls whether all stored frames may
be treated as complete.

## Choose analyses and cases

The four `analysisN_...` true/false switches still enable or disable an analysis
globally. Each enabled analysis runs only where its `analysisN_phases` list
intersects `analysisPhases`:

```matlab
% Run all phases, with profiles only for the three event phases.
analysisPhases = ["pre-inception", "inception", "desinence", "fullcase"];
analysis1_verticalProfiles = true;
analysis1_phases = ["pre-inception", "inception", "desinence"];
analysis4_frameSummary = true;
analysis4_phases = ["pre-inception", "inception", "desinence", "fullcase"];
```

The default profile phase list excludes `fullcase`; add it to
`analysis1_phases` if you want full-case profiles. The other analysis phase
lists include all four phases. A phase with no enabled analyses is skipped.
An empty phase list (`[]`) disables that analysis for every phase. If no
analysis can run at all, the script reports the configuration problem.

Each analysis also has an independent case filter:

```matlab
analysis1_caseIndices = [1 3 6]; % Profiles for these original cases only.
analysis4_caseIndices = [];    % Frame summaries for every case.
```

Filters apply to every phase selected for that analysis. A case unused by
every active analysis is never opened. Original case numbers are retained in
filenames, summary rows, metadata, and profile styles, even for reordered
subsets. For instantaneous exports, `instantaneousFramesToSave` frames are
taken from the beginning of each phase's selected range; the shared color
scale is calculated independently for each phase.

## Results and plotting

Each run has a unique folder, organized by phase and then analysis:

```text
results/<date>/run_<timestamp>/
  run_metadata.mat
  inception/
    run_metadata.mat
    vertical profiles/
      plot mat files/vertical_profile_data.mat
      station01_xplus_1.50mm/...
    instantaneous frames/
      plot mat files/case01_label.mat
      case01_label_frames/...
    profile station visualization/
      plot mat files/mean_speed_map_data.mat
      case01_label_mean_speed_map_....png
    frame summary/
      plot mat files/frame_summary_data.mat
  pre-inception/...
  desinence/...
  fullcase/...
```

Only requested phases with active analyses and their enabled analysis folders
are created. Run metadata records the phase lists, analysis switches, case
filters, and configured/resolved ranges. Saved analysis MAT files include
their phase, selected case labels, original case indices, and actual ranges.

Use the scripts in `plotting codes` to refine saved results without reopening
the velocity files. Set `options.analysisPhase = 'desinence'` in
`raft_plot_options.m` or an individual plotting script. When this option is
empty, the requested plot automatically uses the phase only if exactly one
phase contains that analysis; otherwise it asks you to select a phase.
See [the plotting instructions](plotting%20codes/README.md) for other controls.

## Verification

From this directory in MATLAB:

```matlab
addpath('tests');
test_phase_analysis;
```

The integration check uses temporary synthetic v7.3 files, including both
full-image and packed-ROI storage. It checks phase and case routing, saved
numerical results, completion limits, missing/unused cases, validation, and
new and legacy plotting sources. It does not read your experiment files.
