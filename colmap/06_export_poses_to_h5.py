#!/usr/bin/env python3
"""Write the RGB photogrammetry poses into gref4hsi H5 files as a navigation group.

gref4hsi (`pose_export_type = h5_embedded`) reads a body pose series from the H5:
ECEF position, roll/pitch/yaw of the body frame relative to NED (degrees) and unix
timestamps, then interpolates it onto the HSI line timestamps and applies the HSI
boresight and lever arm from the calibration XML.

Here the body frame is defined from the RGB camera: x = image-up (forward),
y = image-right, z = optical axis (down). This matches the vehicle body frame used in
the Otter preprocessing (x forward, z down) to within the unknown RGB mounting
rotation, so the existing HSI boresight (rz = -90 deg) can be reused as a first guess.
Timestamps are the RGB frame timestamps, which share the UHI clock with the HSI lines,
so no clock offset is needed.

The pose series is padded by half a second at both ends by linear extrapolation so
that HSI lines just outside the RGB span are covered.

Usage:
  06_export_poses_to_h5.py --config <mission>/colmap/config/pipeline.json \
      --h5-dir <other_mission>/Input/H5 [--group raw/nav_rgbsfm]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from pyproj import CRS, Proj, Transformer
from scipy.spatial.transform import Rotation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True, help="colmap pipeline.json of the mosaic workspace")
    parser.add_argument("--h5-dir", type=Path, required=True, help="Folder with the H5 files to write into")
    parser.add_argument("--group", default="raw/nav_rgbsfm", help="H5 group to create (default raw/nav_rgbsfm)")
    parser.add_argument("--pad-s", type=float, default=0.5, help="Extrapolated padding at both ends [s]")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    root = config_path.parents[1]
    meta = root / "metadata"
    origin_info = json.loads((meta / "local_origin.json").read_text())
    origin, epsg = np.array(origin_info["origin"]), int(origin_info["epsg"])
    poses = json.loads((meta / "aligned_poses.json").read_text())["poses"]
    frames = {Path(r["image"]).name: float(r["timestamp_unix"]) for r in csv.DictReader((root / "metadata/frames.csv").open())}
    poses.sort(key=lambda p: frames[p["image"]])
    t = np.array([frames[p["image"]] for p in poses])

    # --- positions: local grid -> lon/lat/h -> ECEF ---------------------------
    local = np.array([p["center"] for p in poses])
    E, N, h = local[:, 0] + origin[0], local[:, 1] + origin[1], local[:, 2] + origin[2]
    to_geo = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    lon, lat = to_geo.transform(E, N)
    to_ecef = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
    ecef = np.column_stack(to_ecef.transform(lon, lat, h))

    # --- orientations: camera (world->cam, grid ENU) -> body -> NED ----------
    # grid north differs from true north by the meridian convergence of the projection
    proj = Proj(CRS.from_epsg(epsg))
    conv = np.array([proj.get_factors(lo, la).meridian_convergence for lo, la in zip(lon, lat)])  # degrees
    body_from_cam = np.array([[0, -1, 0],   # body x = -cam y (image up)
                              [1, 0, 0],    # body y = +cam x (image right)
                              [0, 0, 1]])   # body z = +cam z (optical axis, down)
    enu_to_ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    rpy = np.zeros((len(poses), 3))
    for i, p in enumerate(poses):
        R_w2c = np.array(p["R"])                       # grid-ENU world -> camera
        R_body_to_enu = R_w2c.T @ body_from_cam.T      # body axes expressed in grid ENU
        R_body_to_ned = enu_to_ned @ R_body_to_enu     # grid NED
        # rotate from grid north to true north: true heading = grid heading + convergence
        R_grid_to_true = Rotation.from_euler("z", conv[i], degrees=True).as_matrix()
        R_body_to_ned = R_grid_to_true @ R_body_to_ned
        yaw, pitch, roll = Rotation.from_matrix(R_body_to_ned).as_euler("ZYX", degrees=True)
        rpy[i] = [roll, pitch, yaw]

    # --- pad both ends by linear extrapolation of position, constant attitude -
    def pad(tt, pos, att, dt):
        v0 = (pos[1] - pos[0]) / (tt[1] - tt[0])
        v1 = (pos[-1] - pos[-2]) / (tt[-1] - tt[-2])
        tt2 = np.concatenate([[tt[0] - dt], tt, [tt[-1] + dt]])
        pos2 = np.vstack([pos[0] - v0 * dt, pos, pos[-1] + v1 * dt])
        att2 = np.vstack([att[0], att, att[-1]])
        return tt2, pos2, att2

    t2, ecef2, rpy2 = pad(t, ecef, rpy, args.pad_s)

    out_csv = meta / "poses_body_ned.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["timestamp_unix", "ecef_x", "ecef_y", "ecef_z", "roll_deg", "pitch_deg", "yaw_deg"])
        for i in range(len(t2)):
            w.writerow([f"{t2[i]:.6f}", *[f"{v:.4f}" for v in ecef2[i]], *[f"{v:.3f}" for v in rpy2[i]]])
    print(f"{len(t)} poses (+2 padding): yaw {rpy[:, 2].min():.1f}..{rpy[:, 2].max():.1f} deg, "
          f"roll median {np.median(rpy[:, 0]):.1f}, pitch median {np.median(rpy[:, 1]):.1f} deg, "
          f"meridian convergence {conv.mean():.2f} deg -> {out_csv}")

    # --- write into every H5 (full series in each chunk, like the Otter preprocessing)
    for h5_path in sorted(args.h5_dir.glob("*.h5")):
        with h5py.File(h5_path, "r+") as f:
            if args.group in f:
                del f[args.group]
            g = f.create_group(args.group)
            g.create_dataset("position_ecef", data=ecef2)
            g.create_dataset("euler_angles", data=rpy2)
            g.create_dataset("timestamp", data=t2)
            g.attrs["source"] = "RGB COLMAP sliding-window poses (colmap/03_georeference.py)"
            g.attrs["body_frame"] = "x = RGB image-up, y = image-right, z = optical axis; angles roll,pitch,yaw vs true-north NED, degrees"
            th = f["rawdata/hyperspectral/timestamp"][:]
            inside = (th >= t2[0]) & (th <= t2[-1])
            # compare with the original navigation heading where available
            msg = ""
            if "raw/nav/euler_angles" in f:
                tn = f["raw/nav/timestamp"][:]
                yaw_nav = np.interp(t, tn, np.unwrap(np.radians(f["raw/nav/euler_angles"][:, 2])))
                d = np.degrees(np.angle(np.exp(1j * (np.radians(rpy[:, 2]) - yaw_nav))))
                msg = f"; yaw minus original nav yaw: median {np.median(d):+.1f} deg, IQR {np.percentile(d, 25):+.1f}..{np.percentile(d, 75):+.1f}"
            print(f"{h5_path.name}: wrote {args.group}; HSI lines covered {inside.sum()}/{len(th)}{msg}")


if __name__ == "__main__":
    main()
