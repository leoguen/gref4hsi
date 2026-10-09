# COLMAP RGB orthomosaic scripts (Otter / UHI RGB frames)

Scripts to turn the RGB frames stored in gref4hsi HDF5 files into a georeferenced
RGB orthomosaic, using COLMAP (3.7, CPU) for relative camera poses and the gref4hsi
navigation and DEM for georeferencing.

The scripts live here; **all data, config, logs and results live in a `colmap/`
folder inside the mission** (e.g. `<mission>/gref4hsi_test/colmap/`). Every script
takes `--config <mission>/colmap/config/pipeline.json` and resolves all paths in
that file relative to the `colmap/` folder. The mission folder's README documents
the findings for that dataset.

| Stage | Script | Writes (inside `<mission>/colmap/`) |
|---|---|---|
| 1 | `01_extract_rgb.py` | `images/extracted/*.jpg`, `metadata/frames.csv` |
| 2 | `02_run_sparse.py features|matches|mapper|all|status` | `workspace/database.db`, `workspace/sparse_fixed/` |
| 3 | `03_georeference.py [--estimate-offsets] [--profile-clock]` | `metadata/navigation.csv`, `metadata/camera_offsets.json`, `metadata/aligned_poses.json` |
| 4 | `04_create_orthomosaic.py [--resolution m] [--blend nadir|feather] [--focal px]` | `output/rgb_orthomosaic.tif` (+ `_count.tif`, `report.json`) |

```bash
PY=/home/leo/Documents/NTNU/PhD/UHI/gref4hsi/gref4hsi_venv/bin/python
CFG=/media/leo/NESP_1/.../gref4hsi_test/colmap/config/pipeline.json
$PY 01_extract_rgb.py --config $CFG
$PY 02_run_sparse.py all --config $CFG
$PY 03_georeference.py --config $CFG --estimate-offsets
$PY 04_create_orthomosaic.py --config $CFG
```

## Method in short

* **SfM with fixed intrinsics.** Self-calibration diverges on a flat seabed under
  forward motion, so the camera model is held fixed (`colmap.refine_intrinsics = false`).
* **Sliding-window georeferencing.** The incremental model drifts in scale and bends
  along the chain, so each image gets its own Sim3 fitted on its +-`window` neighbours
  against the navigation, with a near-nadir constraint (the navigation has no height and
  a straight track leaves the roll about the track undetermined).
* **Camera offsets.** `--estimate-offsets` searches the horizontal lever arm between the
  navigation reference (antenna) and the camera, plus a clock offset, by minimising the
  sliding-fit misfit. Because the vehicle yaw is periodic, the clock offset has periodic
  near-minima; `--profile-clock` prints misfit vs. clock offset so you can check that the
  chosen minimum is the deepest. Results are stored in `metadata/camera_offsets.json`
  and reused on later runs without `--estimate-offsets`.
* **Orthomosaic.** Images are projected onto the gref4hsi DEM (`georeference.dem_path`)
  through the OPENCV camera model, after one global flat-field correction. `nadir` blend
  lets the most central image win per cell (sharp, seams visible); `feather` averages
  (smooth, blurrier) and reports `overlap_luminance_std`, a registration quality number.
  `mosaic.focal_override_px` lets the projection use a different focal length than SfM.

Dependencies: COLMAP on PATH, numpy, scipy, pyproj, rasterio, opencv, h5py (all in the
gref4hsi venv).
