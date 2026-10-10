# Replot saved RAFT analysis

Open and run any `plot_*.m` script in MATLAB. No analysis workspace or original
velocity MAT files are needed. Edit `raft_plot_options.m` for shared settings,
or add overrides in an individual script's refinement section.

| Main analysis section | Plotting script | Saved input |
| --- | --- | --- |
| 1: all vertical profiles | `plot_01_vertical_profiles.m` | `vertical_profile_data.mat` |
| 1: mean streamwise velocity | `plot_01a_mean_streamwise_velocity.m` | Same profile MAT |
| 1: Reynolds shear | `plot_01b_reynolds_shear.m` | Same profile MAT |
| 1: streamwise fluctuations | `plot_01c_streamwise_fluctuations.m` | Same profile MAT |
| 1: wall-normal fluctuations | `plot_01d_wall_normal_fluctuations.m` | Same profile MAT |
| 1: 2-component fluctuation energy K_2C | `plot_01e_two_component_fluctuation_energy.m` | Same profile MAT |
| 2: instantaneous U and streamlines | `plot_02_instantaneous_velocity.m` | `caseNN_label.mat` subsets |
| 3: mean-speed maps and station lines | `plot_03_mean_speed_maps.m` | `mean_speed_map_data.mat` |
| 4: optional ROI frame summary | `plot_04_frame_summary.m` | `frame_summary_data.mat` |

Plot script filenames and headings match the analysis numbers in
`RAFT_chunk_analysis.m`. The frame-summary plot renders time traces from
the table saved by analysis 4.

The fifth profile metric is 2-component fluctuation energy,
`K_2C = 0.5 * (sigmaU^2 + sigmaV^2)`, plotted as `K_2C / Ub^2`.
Replotting also accepts the older `tkeInPlane_over_Ub2` saved field.

## Selecting data

By default, the scripts choose the newest dated results folder with saved
analysis data, then its newest `run_YYYYMMDD_HHMMSS` folder. New runs save
data at `<run>/<phase>/<analysis>/plot mat files`. Set
`options.analysisPhase` to `'pre-inception'`, `'inception'`, `'desinence'`, or
`'fullcase'` to choose a phase. An empty phase option auto-selects only when
exactly one phase contains the requested analysis; multiple matching phases
produce a clear error. Explicit date folders also select their newest run.
Existing results
directly under a date folder are also supported: the newest matching section MAT is used;
instantaneous cases come from one export folder. Legacy dates do not identify
a shared run across sections, so use an exact file when provenance matters.

Missing data in the chosen run produces a clear error; the scripts never
silently switch to an older date. For an older run or an exact saved file:

```matlab
options.runDir = fullfile(options.resultsRoot,'2026-10-01');
options.analysisPhase = 'inception';
% Or choose one exact input (for instantaneous plots, this selects one case):
options.dataFile = 'D:\path\to\run_YYYYMMDD_HHMMSS\inception\vertical profiles\plot mat files\vertical_profile_data.mat';
```

`runDir` also accepts an exact phase folder or analysis folder. An exact
`dataFile` takes precedence over phase and run selection. Legacy single-phase
runs with a shared `plot mat files` folder remain supported. A skipped
analysis in the chosen phase produces an error; it never uses another phase's
data or an older run automatically.

## Refining plots

Use `caseIndices`, `stationIndices`, and instantaneous `frameIndices` to select
data. `caseIndices` are positions among saved cases, which may be a subset of
the original MAT files. `frameIndices` are positions in the saved subset (for example, `1:3`),
not the original source frame numbers. Set `xLimits`, `yLimits`, `colorLimits`,
font sizes, figure size, legend location, or `showGrid` to adjust presentation.
Profile styles retain the original case order even when selecting cases.
`useSavedMarkerOptions = false` lets edited `markerOptions` override settings
saved by newer analysis runs. `useSavedColorLimits = false` lets instantaneous
plots calculate a scale using `negativeColorMin_mps`; an explicit
`colorLimits` always takes priority. Set `showStreamlines = false` to show
instantaneous velocity without streamlines.

Plots open for interactive editing and save PNG and editable MATLAB FIG files.
Set `savePDF = true` for PDF exports, or `closeAfterSave = true` for batch use.
Defaults match the original profile styles and map/instantaneous color schemes.
Profile scales stay consistent across saved stations and cases. Instantaneous
subsets are read one frame at a time and retain their saved physical coordinates.
Negative Reynolds stresses remain undefined in the square-root profile.

New outputs go to a unique `replots/<section>_<timestamp>` folder within the
chosen phase's analysis folder (or within the run for legacy inputs);
originals are preserved. Each vertical-profile output contains
station folders such as `station01_xplus_1.50mm`, with all requested metrics
for that station inside. The main analysis uses the same station folder names
under `<run>/<phase>/vertical profiles`.

`plot_settings.mat` records settings and source paths in a `plot mat files`
subfolder of each replot output. An explicit `outputDir` writes to that folder
and can overwrite existing plots with the same names.

## Saving optional sections in the main analysis

The optional sections remain disabled by default. To generate their plot data,
edit `RAFT_chunk_analysis.m` and run the analysis:

- `analysis4_frameSummary = true` saves the frame-summary data.
- `analysis3_profileStationVisualization = true` computes and saves maps.
- `analysis2_instantaneousFrameSaving = true` and
  `instantaneousFramesToSave = 9` save nine frames per case and phase.

`analysisPhases` chooses which phases run. Each `analysisN_phases` list and
`analysisN_caseIndices` filter controls that analysis independently. See
[the analysis configuration guide](../README.md) for ranges and examples.

Each future run saves provenance in `run_metadata.mat`, profiles with their
marker settings, maps with coordinates and actual station positions, and
instantaneous subsets with their original shared color scale.
Each phase has its own analysis folders and each analysis has its own
`plot mat files` folder. Instantaneous images go into per-case folders under
`<run>/<phase>/instantaneous frames`.
