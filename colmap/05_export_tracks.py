#!/usr/bin/env python3
"""Export the navigation and camera tracks as GeoJSON for overlay in QGIS.

Writes to ``<mission>/colmap/output/tracks/``:
  antenna_points.geojson  / antenna_line.geojson   navigation reference at each RGB frame time
  camera_nav_points.geojson / camera_nav_line.geojson  navigation shifted by the estimated lever arm + clock offset
  camera_sfm_points.geojson / camera_sfm_line.geojson  final per-image camera centres from the sliding-window alignment
Point attributes: frame index, image name, time, heading (deg, image-up direction) and fit residual.
GeoJSON is written in WGS84 (EPSG:4326) as the format requires; QGIS reprojects it onto the EPSG:25832 mosaic.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from pyproj import Transformer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True,
                        help="Mission pipeline.json, e.g. <mission>/colmap/config/pipeline.json")
    return parser.parse_args()


def write_geojson(path: Path, features: list) -> None:
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    root = config_path.parents[1]
    meta = root / "metadata"
    out = root / "output" / "tracks"
    out.mkdir(parents=True, exist_ok=True)
    origin_info = json.loads((meta / "local_origin.json").read_text())
    origin, epsg = np.array(origin_info["origin"]), int(origin_info["epsg"])
    to_wgs = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)

    nav = {r["image"]: r for r in csv.DictReader((meta / "navigation.csv").open())}
    cam_nav = {}
    if (meta / "navigation_camera.csv").exists():
        cam_nav = {r["image"]: r for r in csv.DictReader((meta / "navigation_camera.csv").open())}
    poses = json.loads((meta / "aligned_poses.json").read_text())
    offsets = poses.get("camera_offsets", {})
    sfm = {p["image"]: p for p in poses["poses"]}
    names = sorted(nav)

    def heading_deg(p) -> float:
        R = np.array(p["R"])
        f = -R.T[:, 1]  # image-up direction in world
        return float(np.degrees(np.arctan2(f[0], f[1])) % 360)

    layers = {"antenna": [], "camera_nav": [], "camera_sfm": []}
    for i, n in enumerate(names):
        r = nav[n]
        base = {"frame": i, "image": n, "time": float(r["timestamp_unix"])}
        if n in sfm:
            base["heading_deg"] = heading_deg(sfm[n])
        layers["antenna"].append((float(r["east"]), float(r["north"]), dict(base)))
        if n in cam_nav:
            c = cam_nav[n]
            layers["camera_nav"].append((float(c["cam_local_x"]) + origin[0], float(c["cam_local_y"]) + origin[1],
                                         dict(base, clock_s=offsets.get("clock_s"), forward_m=offsets.get("forward_m"),
                                              right_m=offsets.get("right_m"))))
        if n in sfm:
            p = sfm[n]
            layers["camera_sfm"].append((p["center"][0] + origin[0], p["center"][1] + origin[1],
                                         dict(base, residual_xy_m=p["residual_xy_m"], window_scale=p["window_scale"])))

    for name, pts in layers.items():
        if not pts:
            continue
        lon, lat = to_wgs.transform([p[0] for p in pts], [p[1] for p in pts])
        feats = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [float(lo), float(la)]},
                  "properties": dict(p[2], east=round(p[0], 3), north=round(p[1], 3))}
                 for lo, la, p in zip(lon, lat, pts)]
        write_geojson(out / f"{name}_points.geojson", feats)
        line = {"type": "Feature", "geometry": {"type": "LineString",
                                                "coordinates": [[float(lo), float(la)] for lo, la in zip(lon, lat)]},
                "properties": {"layer": name, "epsg_source": epsg}}
        write_geojson(out / f"{name}_line.geojson", [line])
        print(f"{name}: {len(pts)} points -> {out / (name + '_points.geojson')}")
    print("Load the *_line.geojson and *_points.geojson files in QGIS on top of output/rgb_orthomosaic.tif.")


if __name__ == "__main__":
    main()
