# Ecotone HDF5 inspection helpers

Open `inspect_ecotone_h5.ipynb` with the `gref4hsi_venv` kernel. Set `SOURCE`
to either one HDF5 file or a directory containing sequential HDF5 chunks.

The notebook inventories the HDF5 tree, summarizes acquisition chunks, extracts
IMU and altimeter telemetry, reads calibration data, previews RGB and calibrated
hyperspectral imagery, and samples mean spectra without loading an entire cube.
Optional exports are written beneath `output/`, not beside the source data.
It can also create and display a compact, timestamped GIF sampled from the RGB
frames across the complete transect. RGB previews use a display-only percentile
stretch and gamma correction by default because the raw underwater frames are dark.

The default source is:

```text
/media/leo/NESP_1/NTNU/Fieldwork/UHI_Data/Slettvik_29092023/Transect3
```
