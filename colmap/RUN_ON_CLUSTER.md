# Brief: run the COLMAP RGB mosaic pipeline on the first 10 transects

You are working on a compute cluster with a copy of the gref4hsi repository and a copy
of the Otter/Stavøya mission folder (`gref4hsi_full`, the 24 Oct 2024 survey). The
goal is to run the RGB orthomosaic pipeline in `colmap/` of the repository on the
first 10 transects, one transect at a time, and report what worked and what did not.
Do not modify the H5 files; everything the pipeline writes goes into a per-transect
`colmap/` folder that you create.

## What exists

* Scripts: `<repo>/colmap/01_extract_rgb.py` ... `05_export_tracks.py`. Read
  `<repo>/colmap/README.md` and `README_otter.md` first; they describe the method and
  the assumptions. Scripts 06 and 07 (feeding poses back into gref4hsi) are not part
  of this task.
* Mission data (adapt the root path to the cluster):
  * `gref4hsi_full/Input/H5/<transect>_<n>.h5` - all transects' chunks in ONE folder.
    Transect names look like `uhi_20241024_134634`; chunks are `_1.h5`, `_2.h5`, ...
    The first 10 transects by name are the first 10 distinct prefixes in sorted order.
  * `gref4hsi_full/Intermediate/pose.csv` - navigation for the whole mission
    (ECEF, 25 Hz). One file for all transects; the scripts interpolate by timestamp.
  * `gref4hsi_full/Input/GIS/<transect>/dem.tif` - altimeter DEM per transect,
    EPSG:25832, ellipsoidal heights. Required by stage 4.
* A reference configuration from the pilot transect: copy it from
  `<repo>/colmap/pipeline_reference.json` (same content as the pilot's
  `config/pipeline.json`). It already holds the camera intrinsics, focal override,
  offsets settings and mosaic settings that were tuned on transect
  `uhi_20241024_140601`.

## Requirements

* COLMAP on PATH (3.7 or newer). If the node has CUDA and COLMAP was built with it,
  set `"use_gpu": true` in the config; otherwise leave `false`.
* Python with numpy, scipy, pyproj, rasterio, opencv-python, h5py. Check with
  `python -c "import numpy, scipy, pyproj, rasterio, cv2, h5py"`.
* Per transect: ~0.5 GB disk, 4 to 8 cores, 8 GB RAM for the 2 mm mosaic (less at
  1 cm). The pilot transect (217 frames) took about 1 min extraction, 5 min
  features + matching, 2 to 4 min mapping, 10 s georeferencing, 1 min mosaic at 1 cm
  and 14 min at 2 mm, all on 4 CPU threads.

## Per-transect setup

For each transect `T` (first 10 by sorted name), create
`gref4hsi_full/colmap_<T>/config/pipeline.json` from the reference config and edit:

```json
"paths": {
  "h5_dir": "../Input/H5",
  "h5_pattern": "<T>_*.h5",          // only this transect's chunks
  ...unchanged...
},
"georeference": {
  "pose_csv": "../Intermediate/pose.csv",
  "epsg": 25832,
  "dem_path": "../Input/GIS/<T>/dem.tif",
  "window": 6,
  "nadir_constraint_length_m": 1.0,
  "altimeter_height_m": 1.21,
  "estimate_offsets": true          // estimate lever arm + clock per transect
},
"colmap": { ..., "num_threads": <cores you have>, "use_gpu": <true|false> },
"mosaic": { ..., "resolution_m": 0.01 }   // start at 1 cm; 2 mm only if 1 cm looks right
```

All relative paths resolve against `colmap_<T>/` (the parent of `config/`).

## Run, per transect, in this order

```bash
PY=python   # the environment with the packages above
CFG=gref4hsi_full/colmap_<T>/config/pipeline.json
$PY <repo>/colmap/01_extract_rgb.py --config $CFG
$PY <repo>/colmap/02_run_sparse.py all --config $CFG
$PY <repo>/colmap/03_georeference.py --config $CFG --estimate-offsets
$PY <repo>/colmap/04_create_orthomosaic.py --config $CFG
$PY <repo>/colmap/05_export_tracks.py --config $CFG
```

Run transects sequentially or in parallel across nodes; never two runs in the same
`colmap_<T>` folder. Each script prints a short summary; logs go to `colmap_<T>/logs/`.

## Checks to make and report for every transect

1. **Extraction**: number of frames, time span, whether the chunks are in
   chronological order (the script validates and raises otherwise).
2. **Matching**: run `02_run_sparse.py status --config $CFG` and record images,
   two_view_geometries. Adjacent frames should share hundreds of inliers.
3. **Mapper**: how many sparse models were produced in `workspace/sparse_fixed/`, and
   `colmap model_analyzer --path workspace/sparse_fixed/0` registered images and mean
   reprojection error. Expected: one model with (nearly) all images, 0.3 to 0.5 px.
   More than one model means the chain broke; note the frame ranges of each model
   (`colmap model_converter --output_type TXT` then read `images.txt`) and continue,
   the georeferencing handles several models.
4. **Georeferencing** (printed by script 03): the sliding-window xy residual RMS
   (expect 0.05 to 0.15 m), median tilt from nadir (expect < 6 deg), and the
   estimated offsets `forward_m`, `right_m`, `clock_s`. These offsets are physical
   constants of the vehicle, so they should be similar across transects. Report the
   10 values side by side; a consistent set is the most valuable result of this run.
   If they scatter, also run `--profile-clock` on two transects and save the tables.
5. **Mosaic** (printed by script 04 and in `output/report.json`): covered area,
   `overlap_luminance_std` (expect 4 to 6 DN; higher means misregistration) and the
   grid size. Open the GeoTIFF once (any viewer) to check it is not blank.
6. Anything that crashed: the traceback and the stage.

## Known pitfalls

* Transects with very few frames (a chunk of a few seconds) will fail in the mapper;
  report and skip.
* If the mapper produces many small models, do not change the intrinsics; note it.
  The intrinsics in the reference config are the only working set found so far.
* The scripts write `metadata/local_origin.json` on first run and reuse it; if you
  restart a transect from scratch, delete the whole `colmap_<T>` folder.
* The DEM must cover the camera footprint; if script 04 reports a low "DEM valid
  fraction" the transect's DEM is missing or offset. Report it.
* Do not run `02_run_sparse.py features` twice on the same database; it appends.
  Delete `workspace/database.db` to redo features.

## Deliverable

A short table with one row per transect: frames, registered images / models,
reprojection error, georef residual RMS, tilt, forward_m, right_m, clock_s, covered
area, overlap std, status. Plus the paths of the 10 `output/rgb_orthomosaic.tif` files
and the per-transect `logs/` folders. Keep all outputs; nothing needs to be cleaned up.
