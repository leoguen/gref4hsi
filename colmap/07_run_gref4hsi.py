#!/usr/bin/env python3
"""Run gref4hsi stages on a mission whose configuration points at RGB-derived poses.

The Otter runner (helperscripts/ecotone/process_otter_stavoya.py) rewrites the
[HDF.raw_nav] section on every call, which would undo the switch to the
``raw/nav_rgbsfm`` group written by 06_export_poses_to_h5.py. This runner calls the
gref4hsi stages directly with the configuration as it is on disk.

Usage: 07_run_gref4hsi.py --config <mission>/configuration.ini [--stages pose georeference orthorectify]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from gref4hsi.scripts import georeference, orthorectification  # noqa: E402
from gref4hsi.utils import parsing_utils  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True, help="gref4hsi configuration.ini of the mission")
    parser.add_argument("--stages", nargs="+", default=["pose", "georeference", "orthorectify"],
                        choices=["pose", "model", "georeference", "orthorectify"])
    args = parser.parse_args()
    cfg = str(args.config.expanduser().resolve())
    for stage in args.stages:
        t0 = time.time()
        print(f"=== {stage} ===", flush=True)
        if stage == "pose":
            parsing_utils.export_pose(cfg)
        elif stage == "model":
            parsing_utils.export_model(cfg)
        elif stage == "georeference":
            georeference.main(cfg)
        elif stage == "orthorectify":
            orthorectification.main(cfg)
        print(f"=== {stage} done in {time.time() - t0:.0f} s ===", flush=True)


if __name__ == "__main__":
    main()
