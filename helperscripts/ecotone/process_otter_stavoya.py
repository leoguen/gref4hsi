"""Prepare and process the 2024-10-24 Otter/Stavøya UHI test transect.

Run from the repository root with the project virtual environment::

    source gref4hsi_venv/bin/activate
    python helperscripts/ecotone/process_otter_stavoya.py prepare

The preparation stage modifies only the HDF5 working copies below, never the
files in the acquisition ``raw`` directory.
"""

from __future__ import annotations

import argparse
import configparser
from collections import namedtuple
from pathlib import Path

import numpy as np

from gref4hsi.scripts import georeference, orthorectification
from gref4hsi.utils import parsing_utils, uhi_parsing_utils
from gref4hsi.utils.config_utils import customize_config, prepend_data_dir_to_relative_paths


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_MISSION = Path(
    "/media/leo/NESP_1/NTNU/UHI/Malin/"
    "2024_24_Oct_Otter_Stavøya/24_10_24/gref4hsi_test"
)
CONFIG_TEMPLATE = REPOSITORY / "data/config_examples/configuration_uhi.ini"


SettingsPreprocess = namedtuple(
    "SettingsPreprocessing",
    [
        "dtype_datacube",
        "rotation_matrix_hsi_to_body",
        "translation_body_to_hsi",
        "rotation_matrix_alt_to_body",
        "translation_alt_to_body",
        "config_file_name",
        "time_offset_sec",
        "lon_lat_alt_origin",
        "resolution_dem",
        "agisoft_process",
        "smooth_navigation",
        "position_smoothing_window",
        "heading_smoothing_window",
        "attitude_smoothing_window",
        "smoothing_polyorder",
        "median_filter_window",
    ],
)


# Confirmed to use the same mounting convention as the UHI-DBE setup.
DBE_MOUNTING_ROTATION = np.array(
    [
        [0, 1, 0],
        [-1, 0, 0],
        [0, 0, 1],
    ],
    dtype=float,
)

def otter_settings(mission: Path) -> SettingsPreprocess:
    """Build preprocessing settings using the mission's median navigation origin."""
    import h5py

    latitudes = []
    longitudes = []
    h5_files = sorted((mission / "Input/H5").glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No HDF5 files found in {mission / 'Input/H5'}")
    for path in h5_files:
        with h5py.File(path, "r") as handle:
            navigation = handle["rawdata/navigation/external"]
            latitudes.append(navigation["Latitude"][::20])
            longitudes.append(navigation["Longitude"][::20])

    origin = np.array(
        [
            np.nanmedian(np.concatenate(longitudes)),
            np.nanmedian(np.concatenate(latitudes)),
            0.0,
        ]
    )
    return SettingsPreprocess(
        dtype_datacube=np.float32,
        rotation_matrix_hsi_to_body=DBE_MOUNTING_ROTATION,
        translation_body_to_hsi=np.zeros(3),
        rotation_matrix_alt_to_body=DBE_MOUNTING_ROTATION,
        translation_alt_to_body=np.zeros(3),
        config_file_name="configuration.ini",
        time_offset_sec=0.0,
        lon_lat_alt_origin=origin,
        resolution_dem=0.05,
        agisoft_process=False,
        smooth_navigation=True,
        position_smoothing_window=15,
        heading_smoothing_window=11,
        attitude_smoothing_window=7,
        smoothing_polyorder=2,
        median_filter_window=5,
    )


def create_configuration(mission: Path, overwrite: bool = False) -> configparser.ConfigParser:
    """Create the mission folder structure and its pipeline configuration."""
    config_file = mission / "configuration.ini"
    if overwrite or not config_file.exists():
        prepend_data_dir_to_relative_paths(
            config_path=str(CONFIG_TEMPLATE),
            DATA_DIR=str(mission),
        )

    custom_config = {
        "General": {
            "mission_dir": str(mission),
            "model_export_type": "dem_file",
            "max_ray_length": 20,
            "dem_per_transect": True,
        },
        "Coordinate Reference Systems": {
            "proj_epsg": 25832,
            "dem_epsg": 25832,
            "geocsc_epsg_export": 4978,
            "pos_epsg_orig": 4978,
            "dem_ref": "ellipsoid",
        },
        "Absolute Paths": {
            "geoid_path": str(REPOSITORY / "data/world/geoids/egm08_25.gtx"),
        },
        "HDF.raw_nav": {
            "eul_zyx": "raw/nav/euler_angles",
            "position": "raw/nav/position_ecef",
            "timestamp": "raw/nav/timestamp",
            "altitude": "raw/nav/altitude",
            "rotation_reference_type": "eul_ZYX",
            "is_global_rot": False,
            "eul_is_degrees": True,
        },
        "HDF.hyperspectral": {
            "datacube": "rawdata/hyperspectral/dataCube",
            "exposuretime": "rawdata/hyperspectral/exposureTime",
            "timestamp": "rawdata/hyperspectral/timestamp",
            "is_calibrated": False,
        },
        "HDF.calibration": {
            "band2wavelength": "rawdata/hyperspectral/calibration/spectral/band2Wavelength",
            "darkframe": "rawdata/hyperspectral/calibration/radiometric/darkFrame",
            "radiometricframe": "rawdata/hyperspectral/calibration/radiometric/radiometricFrame",
            "fov": "rawdata/hyperspectral/calibration/geometric/fieldOfView",
        },
        "HDF.rgb": {
            "rgb_frames": "rawdata/rgb/rgbFrames",
            "rgb_frames_timestamp": "rawdata/rgb/timestamp",
        },
        "Ancillary": {
            "pixel_nr_grid": "processed/georef/pixel_nr_grid",
            "unix_time_grid": "processed/georef/unix_time_grid",
        },
        "Orthorectification": {
            "resample_rgb_only": True,
            "resample_ancillary": False,
            "resolutionhyperspectralmosaic": 0.005,
            "raster_transform_method": "north_east",
            "mask_pixel_by_footprint": True,
        },
    }
    customize_config(str(config_file), custom_config)

    config = configparser.ConfigParser()
    config.read(config_file)
    return config


def run_stage(stage: str, mission: Path, reset_config: bool = False) -> None:
    config_file = mission / "configuration.ini"
    config = create_configuration(mission, overwrite=reset_config)
    settings = otter_settings(mission)
    print(f"Mission: {mission}")
    print(f"Navigation origin [lon, lat, h]: {settings.lon_lat_alt_origin}")

    if stage == "prepare":
        uhi_parsing_utils.uhi_otter(config=config, config_uhi=settings)
    elif stage == "pose":
        parsing_utils.export_pose(str(config_file))
    elif stage == "model":
        parsing_utils.export_model(str(config_file))
    elif stage == "georeference":
        georeference.main(str(config_file))
    elif stage == "orthorectify":
        orthorectification.main(str(config_file))
    elif stage == "all":
        uhi_parsing_utils.uhi_otter(config=config, config_uhi=settings)
        parsing_utils.export_pose(str(config_file))
        parsing_utils.export_model(str(config_file))
        georeference.main(str(config_file))
        orthorectification.main(str(config_file))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=("prepare", "pose", "model", "georeference", "orthorectify", "all"),
        help="Pipeline stage to run; start with 'prepare'.",
    )
    parser.add_argument(
        "--mission",
        type=Path,
        default=DEFAULT_MISSION,
        help=f"Working mission directory (default: {DEFAULT_MISSION}).",
    )
    parser.add_argument(
        "--reset-config",
        action="store_true",
        help="Recreate configuration.ini from the repository template.",
    )
    args = parser.parse_args()
    run_stage(args.stage, mission=args.mission.resolve(), reset_config=args.reset_config)


if __name__ == "__main__":
    main()
