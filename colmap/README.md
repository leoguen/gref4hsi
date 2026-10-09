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
| 5 | `05_export_tracks.py` | `output/tracks/*.geojson` (antenna, corrected camera and SfM camera tracks for QGIS overlay) |
| 6 | `06_export_poses_to_h5.py --h5-dir <mission>/Input/H5` | `raw/nav_rgbsfm/{position_ecef,euler_angles,timestamp}` in each H5, `metadata/poses_body_ned.csv` |
| 7 | `07_run_gref4hsi.py --config <mission>/configuration.ini` | gref4hsi pose/georeference/orthorectify outputs of that mission |

```bash
PY=/home/leo/Documents/NTNU/PhD/UHI/gref4hsi/gref4hsi_venv/bin/python
CFG=/media/leo/NESP_1/.../gref4hsi_test/colmap/config/pipeline.json
$PY 01_extract_rgb.py --config $CFG
$PY 02_run_sparse.py all --config $CFG
$PY 03_georeference.py --config $CFG --estimate-offsets
$PY 04_create_orthomosaic.py --config $CFG
$PY 05_export_tracks.py --config $CFG
```

## Feeding the poses back into gref4hsi

To georeference the hyperspectral lines with the RGB-derived poses (the Løvås et al.
2022 loop: photogrammetry pose -> fixed HSI/RGB transform -> ray casting):

```bash
# 1. work on a copy of the mission so the navigation-based results stay
# 2. write the poses into the H5 files as a body pose series
$PY 06_export_poses_to_h5.py --config $CFG --h5-dir <mission_copy>/Input/H5
# 3. in <mission_copy>/configuration.ini point [HDF.raw_nav] eul_zyx/position/timestamp
#    at raw/nav_rgbsfm/..., and set the HSI lever arm in Input/Calib/HSI_2b.xml
#    (tx = -0.05: the HSI sits 5 cm behind the RGB camera along image-up)
# 4. run the stages without the Otter runner (it would reset [HDF.raw_nav])
$PY 07_run_gref4hsi.py --config <mission_copy>/configuration.ini
```

The body frame written is x = image-up, y = image-right, z = optical axis, with
roll/pitch/yaw relative to true-north NED (grid convergence applied). Timestamps are
the RGB frame times, which share the UHI clock with the HSI, so no clock offset is
involved. The existing HSI boresight (rz = -90 deg) is reused as a first guess; a
luminance-correlation calibration (paper, method 2) would refine it.

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
* **Orthomosaic (tiled).** Built in tiles so millimetre grids fit in memory; a coarse pass at `mosaic.stats_resolution_m` (1 cm) fixes the colour stretch and measures overlap consistency. Images are projected onto the gref4hsi DEM (`georeference.dem_path`)
  through the OPENCV camera model, after one global flat-field correction. `nadir` blend
  lets the most central image win per cell (sharp, seams visible); `feather` averages
  (smooth, blurrier) and reports `overlap_luminance_std`, a registration quality number.
  `mosaic.focal_override_px` lets the projection use a different focal length than SfM.

Dependencies: COLMAP on PATH, numpy, scipy, pyproj, rasterio, opencv, h5py (all in the
gref4hsi venv).
