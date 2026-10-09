#!/usr/bin/env python3
"""Project aligned COLMAP images onto the gref4hsi DEM and blend an RGB orthomosaic.

For every registered image the DEM cells around the camera nadir are lifted to 3D
(E, N, h in the local frame), transformed into the camera, projected with the OPENCV
camera model and sampled from the flat-field corrected frame. Contributions are
blended with a border-feathering weight (``--blend feather``, weighted average) or
the most central image wins (``--blend nadir``). Output is an 8-bit RGBA GeoTIFF in
the configured EPSG, directly loadable in QGIS.

The mosaic is built tile by tile so that millimetre grids fit in memory. Two passes
are made: a coarse pass (``mosaic.stats_resolution_m``, default 1 cm) that fixes the
global colour stretch and measures the overlap consistency, and the full-resolution
pass that writes the GeoTIFF.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import cv2
import numpy as np
import rasterio
from rasterio.enums import ColorInterp
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module  # noqa: E402

geo = import_module("03_georeference")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True,
                        help="Mission pipeline.json, e.g. <mission>/colmap/config/pipeline.json; all paths resolve relative to its parent folder")
    parser.add_argument("--resolution", type=float, default=None, help="Mosaic cell size [m]")
    parser.add_argument("--blend", choices=("feather", "nadir"), default=None)
    parser.add_argument("--models", nargs="*", default=None, help="Aligned model names to use (default: all)")
    parser.add_argument("--no-flatfield", action="store_true")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--focal", type=float, default=None, help="Override focal length [px] for projection")
    parser.add_argument("--tile", type=int, default=None, help="Tile size in cells (default from config, 2048)")
    parser.add_argument("--source", choices=("raw", "normalized"), default=None,
                        help="Image set to mosaic: raw frames with flat field (default) or the grey-world normalised frames")
    return parser.parse_args()


def flat_field(image_paths: list[Path], sigma: float, vignette_threshold: float = 0.3,
               vignette_erode_px: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Smooth mean frame (unit mean per channel) and a per-pixel weight map.

    The weight map is the distance to the edge of the usable area, normalised to 1 at
    the frame centre. The usable area is where the mean frame is brighter than
    ``vignette_threshold`` times its central brightness (this removes the dark
    housing rim in the corners), eroded by a few pixels; outside it the weight is 0 so
    those pixels never contribute to the mosaic.
    """
    acc = None
    for p in image_paths:
        img = cv2.imread(str(p)).astype(np.float32)
        acc = img if acc is None else acc + img
    mean = acc / len(image_paths)
    k = int(6 * sigma + 1) | 1
    field = cv2.GaussianBlur(mean, (k, k), sigma)
    field /= field.reshape(-1, 3).mean(axis=0)  # unit-mean per channel
    lum = cv2.GaussianBlur(mean.mean(axis=2), (21, 21), 5)
    H, W = lum.shape
    centre = lum[H // 2 - H // 8:H // 2 + H // 8, W // 2 - W // 8:W // 2 + W // 8].mean()
    usable = (lum >= vignette_threshold * centre).astype(np.uint8)
    if vignette_erode_px > 0:
        usable = cv2.erode(usable, np.ones((2 * vignette_erode_px + 1,) * 2, np.uint8))
    dist = cv2.distanceTransform(usable, cv2.DIST_L2, 5)
    weight = dist / max(dist.max(), 1.0)
    print(f"usable image area {usable.mean() * 100:.1f}% (vignette threshold {vignette_threshold} x centre brightness)")
    return np.clip(field, 0.05, None), weight.astype(np.float32)


def sample_image(img: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Bilinear sampling of ``img`` at float pixel positions, in blocks (cv2.remap limit 32767)."""
    n = len(u)
    ncols = 4096
    nrows = -(-n // ncols)
    pad = nrows * ncols - n
    um = np.concatenate([u, np.zeros(pad, np.float32)]).reshape(nrows, ncols)
    vm = np.concatenate([v, np.zeros(pad, np.float32)]).reshape(nrows, ncols)
    return cv2.remap(img, um, vm, cv2.INTER_LINEAR).reshape(-1, 3)[:n]


class Mosaic:
    """Grid definition plus the per-tile accumulation."""

    def __init__(self, bounds_local, res, origin, epsg, dem_path, views, image_dir, field, weight_map, blend,
                 feather_power):
        self.xmin, self.ymin, self.xmax, self.ymax = bounds_local
        self.res = res
        self.origin, self.epsg = origin, epsg
        self.width = int(np.ceil((self.xmax - self.xmin) / res))
        self.height = int(np.ceil((self.ymax - self.ymin) / res))
        self.transform = from_origin(self.xmin + origin[0], self.ymax + origin[1], res, res)
        self.dem_path = dem_path
        self.views = views
        self.image_dir, self.field, self.weight_map = image_dir, field, weight_map
        self.blend, self.feather_power = blend, feather_power
        # per-image candidate window in grid cells (camera nadir +- search radius)
        self.boxes = []
        for _, im, _, _, radius in views:
            c = im["center"]
            c0 = int(np.clip((c[0] - radius - self.xmin) / res, 0, self.width))
            c1 = int(np.clip((c[0] + radius - self.xmin) / res, 0, self.width))
            r0 = int(np.clip((self.ymax - (c[1] + radius)) / res, 0, self.height))
            r1 = int(np.clip((self.ymax - (c[1] - radius)) / res, 0, self.height))
            self.boxes.append((r0, r1, c0, c1))

    def dem_tile(self, r0, r1, c0, c1) -> np.ndarray:
        h, w = r1 - r0, c1 - c0
        tf = from_origin(self.xmin + c0 * self.res + self.origin[0], self.ymax - r0 * self.res + self.origin[1],
                         self.res, self.res)
        dem = np.full((h, w), np.nan, dtype=np.float32)
        with rasterio.open(self.dem_path) as src:
            reproject(rasterio.band(src, 1), dem, dst_transform=tf, dst_crs=f"EPSG:{self.epsg}",
                      src_nodata=src.nodata, dst_nodata=np.nan, resampling=Resampling.bilinear)
        return dem - self.origin[2]

    def accumulate(self, r0, r1, c0, c1, image_cache: dict):
        """Blend all images touching the tile. Returns (rgb float32 BGR, count, overlap_var, covered)."""
        h, w = r1 - r0, c1 - c0
        acc = np.zeros((h, w, 3), np.float32)
        wsum = np.zeros((h, w), np.float32)
        acc_mean = np.zeros((h, w, 3), np.float32) if self.blend == "nadir" else None
        wsum_mean = np.zeros((h, w), np.float32) if self.blend == "nadir" else None
        acc2 = np.zeros((h, w), np.float32)
        count = np.zeros((h, w), np.uint16)
        dem = None
        xs = self.xmin + (np.arange(c0, c1) + 0.5) * self.res
        ys = self.ymax - (np.arange(r0, r1) + 0.5) * self.res
        for idx, (img_name, im, cam, _, _) in enumerate(self.views):
            br0, br1, bc0, bc1 = self.boxes[idx]
            ir0, ir1, ic0, ic1 = max(br0, r0), min(br1, r1), max(bc0, c0), min(bc1, c1)
            if ir1 <= ir0 or ic1 <= ic0:
                continue
            if dem is None:
                dem = self.dem_tile(r0, r1, c0, c1)
            gz = dem[ir0 - r0:ir1 - r0, ic0 - c0:ic1 - c0]
            valid = np.isfinite(gz)
            if not valid.any():
                continue
            gx, gy = np.meshgrid(xs[ic0 - c0:ic1 - c0], ys[ir0 - r0:ir1 - r0])
            pts = np.column_stack([gx[valid], gy[valid], gz[valid]])
            pc = (im["R"] @ pts.T).T + im["t"]
            front = pc[:, 2] > 0.05
            if not front.any():
                continue
            fx, fy, cx, cy, k1, k2, p1, p2 = cam["params"]
            K = np.array([[fx, 0, cx], [0, fy, cy], [0, 0, 1]])
            uv, _ = cv2.projectPoints(pc[front].reshape(-1, 1, 3), np.zeros(3), np.zeros(3), K,
                                      np.array([k1, k2, p1, p2]))
            uv = uv.reshape(-1, 2)
            W, H = cam["width"], cam["height"]
            inside = (uv[:, 0] >= 0) & (uv[:, 0] <= W - 1) & (uv[:, 1] >= 0) & (uv[:, 1] <= H - 1)
            if not inside.any():
                continue
            if img_name not in image_cache:
                if len(image_cache) > 40:
                    image_cache.clear()
                image_cache[img_name] = cv2.imread(str(self.image_dir / img_name)).astype(np.float32) / self.field
            img = image_cache[img_name]
            u, v = uv[inside, 0].astype(np.float32), uv[inside, 1].astype(np.float32)
            col = sample_image(img, u, v)
            # weight: distance to the edge of the usable (non-vignetted) area, 0 outside it
            wmap = self.weight_map[np.clip(np.rint(v).astype(int), 0, H - 1), np.clip(np.rint(u).astype(int), 0, W - 1)]
            keep = wmap > 0
            if not keep.any():
                continue
            col, u, v, wmap = col[keep], u[keep], v[keep], wmap[keep]
            wgt = np.clip(wmap, 1e-3, 1) ** self.feather_power
            rows, cols = np.nonzero(valid)
            rows, cols = rows[front][inside][keep] + (ir0 - r0), cols[front][inside][keep] + (ic0 - c0)
            lum = col.mean(axis=1)
            np.add.at(acc2, (rows, cols), wgt * lum ** 2)
            np.add.at(count, (rows, cols), 1)
            if self.blend == "feather":
                np.add.at(acc, (rows, cols), col * wgt[:, None])
                np.add.at(wsum, (rows, cols), wgt)
            else:
                np.add.at(acc_mean, (rows, cols), col * wgt[:, None])
                np.add.at(wsum_mean, (rows, cols), wgt)
                better = wgt > wsum[rows, cols]
                acc[rows[better], cols[better]] = col[better]
                wsum[rows[better], cols[better]] = wgt[better]
        covered = wsum > 0
        if self.blend == "feather":
            acc[covered] /= wsum[covered][:, None]
            mean_rgb, wm = acc, wsum
        else:
            acc_mean[covered] /= wsum_mean[covered][:, None]
            mean_rgb, wm = acc_mean, wsum_mean
        var = np.zeros((h, w), np.float32)
        multi = count >= 2
        var[multi] = acc2[multi] / wm[multi] - mean_rgb[multi].mean(axis=1) ** 2
        return acc, count, np.clip(var, 0, None), covered

    def tiles(self, size: int):
        for r0 in range(0, self.height, size):
            for c0 in range(0, self.width, size):
                yield r0, min(r0 + size, self.height), c0, min(c0 + size, self.width)


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parents[1]
    paths, mosaic_cfg = config["paths"], config.get("mosaic", {})
    res = args.resolution or mosaic_cfg.get("resolution_m", 0.01)
    stats_res = max(res, float(mosaic_cfg.get("stats_resolution_m", 0.01)))
    blend = args.blend or mosaic_cfg.get("blend", "nadir")
    feather_power = float(mosaic_cfg.get("feather_power", 2.0))
    search_radius = float(mosaic_cfg.get("search_radius_m", 2.5))
    tile = args.tile or int(mosaic_cfg.get("tile_cells", 2048))
    source = args.source or mosaic_cfg.get("image_source", "raw")
    if source == "normalized":
        image_dir = geo.resolve(root, paths.get("colmap_image_dir", "images/normalized"))
        args.no_flatfield = True  # the normalisation already removed the illumination pattern
    else:
        image_dir = geo.resolve(root, paths["image_dir"])
    dem_path = geo.resolve(root, config["georeference"]["dem_path"])
    origin_info = json.loads((root / "metadata/local_origin.json").read_text())
    origin, epsg = np.array(origin_info["origin"]), int(origin_info["epsg"])
    out_dir = root / "output"
    out_dir.mkdir(exist_ok=True)
    output = args.output or out_dir / mosaic_cfg.get("filename", "rgb_orthomosaic.tif")

    # --- per-image poses from the sliding-window alignment -------------------
    poses_file = root / "metadata/aligned_poses.json"
    if not poses_file.exists():
        raise SystemExit("metadata/aligned_poses.json missing; run 03_georeference.py first")
    aligned = json.loads(poses_file.read_text())
    views = []
    for p in aligned["poses"]:
        if args.models and p["model"] not in args.models:
            continue
        cam = dict(p["camera"])
        if cam["model"] != "OPENCV":
            raise SystemExit(f"Unsupported camera model {cam['model']}")
        focal_override = args.focal or mosaic_cfg.get("focal_override_px")
        if focal_override:
            ratio = float(focal_override) / cam["params"][0]
            cam["params"] = [cam["params"][0] * ratio, cam["params"][1] * ratio] + list(cam["params"][2:])
        im = {"R": np.array(p["R"]), "t": np.array(p["t"]), "center": np.array(p["center"])}
        views.append((p["image"], im, cam, p["model"], search_radius))
    use = sorted({v[3] for v in views})
    if not views:
        raise SystemExit("No aligned images found")
    views.sort(key=lambda v: v[0])
    print(f"{len(views)} registered images from models {use}; blend {blend}; focal {views[0][2]['params'][0]:.1f} px; "
          f"source {source} ({image_dir.name})")

    # --- flat field ---------------------------------------------------------
    all_images = sorted(image_dir.glob("*.jpg"))
    field, weight_map = flat_field(all_images, float(mosaic_cfg.get("flatfield_sigma_px", 60)),
                                   float(mosaic_cfg.get("vignette_threshold", 0.3)),
                                   int(mosaic_cfg.get("vignette_erode_px", 8)))
    if args.no_flatfield:
        field = np.ones((1, 1, 3), np.float32)
    else:
        cv2.imwrite(str(out_dir / "flatfield.png"), np.clip(field / field.max() * 255, 0, 255).astype(np.uint8))
    cv2.imwrite(str(out_dir / "pixel_weight.png"), (weight_map * 255).astype(np.uint8))

    centers = np.array([v[1]["center"] for v in views])
    bounds = (centers[:, 0].min() - search_radius, centers[:, 1].min() - search_radius,
              centers[:, 0].max() + search_radius, centers[:, 1].max() + search_radius)
    common = dict(origin=origin, epsg=epsg, dem_path=dem_path, views=views, image_dir=image_dir, field=field,
                  weight_map=weight_map, blend=blend, feather_power=feather_power)

    # --- pass 1: coarse statistics (colour stretch, overlap consistency) -----
    coarse = Mosaic(bounds, stats_res, **common)
    cache: dict = {}
    samples, var_sum, var_n = [], 0.0, 0
    for r0, r1, c0, c1 in coarse.tiles(tile):
        rgb, count, var, covered = coarse.accumulate(r0, r1, c0, c1, cache)
        samples.append(rgb[covered])
        multi = count >= 2
        var_sum += float(np.sqrt(var[multi]).sum())
        var_n += int(multi.sum())
    samples = np.concatenate(samples)
    overlap_std = var_sum / var_n if var_n else float("nan")
    print(f"stats pass at {stats_res} m: {len(samples)} covered cells, overlap luminance std "
          f"{overlap_std:.2f} DN over {var_n} cells (lower = better registration)")
    gamma = float(mosaic_cfg.get("gamma", 1.0))
    lo_pct, hi_pct = mosaic_cfg.get("stretch_low_pct", 0.5), mosaic_cfg.get("stretch_high_pct", 99.5)
    if mosaic_cfg.get("stretch", "common") == "per_channel":
        lo, hi = np.percentile(samples, lo_pct, axis=0), np.percentile(samples, hi_pct, axis=0)
    else:
        lo = np.repeat(np.percentile(samples.mean(axis=1), lo_pct), 3)
        hi = np.repeat(np.percentile(samples.max(axis=1), hi_pct), 3)
    del samples

    # --- pass 2: full resolution, written tile by tile ----------------------
    fine = Mosaic(bounds, res, **common)
    print(f"grid {fine.width} x {fine.height} cells at {res} m, tiles of {tile} cells")
    profile = dict(driver="GTiff", width=fine.width, height=fine.height, count=4, dtype="uint8",
                   crs=f"EPSG:{epsg}", transform=fine.transform, compress="deflate", tiled=True,
                   blockxsize=256, blockysize=256, photometric="RGB", alpha="unspecified", BIGTIFF="IF_SAFER")
    count_path = output.with_name(output.stem + "_count.tif")
    covered_cells = 0
    with rasterio.open(output, "w", **profile) as dst, \
            rasterio.open(count_path, "w", driver="GTiff", width=fine.width, height=fine.height, count=1,
                          dtype="uint16", crs=f"EPSG:{epsg}", transform=fine.transform, compress="deflate",
                          tiled=True, nodata=0, BIGTIFF="IF_SAFER") as dst_count:
        dst.colorinterp = [ColorInterp.red, ColorInterp.green, ColorInterp.blue, ColorInterp.alpha]
        n_tiles = sum(1 for _ in fine.tiles(tile))
        for k, (r0, r1, c0, c1) in enumerate(fine.tiles(tile)):
            rgb, count, _, covered = fine.accumulate(r0, r1, c0, c1, cache)
            win = Window(c0, r0, c1 - c0, r1 - r0)
            for b in range(3):  # source is BGR
                src_b = 2 - b
                band = np.clip((rgb[..., src_b] - lo[src_b]) / max(hi[src_b] - lo[src_b], 1e-3), 0, 1)
                if gamma != 1.0:
                    band **= 1 / gamma
                dst.write((band * 255).astype(np.uint8), b + 1, window=win)
            dst.write((covered * 255).astype(np.uint8), 4, window=win)
            dst_count.write(count, 1, window=win)
            covered_cells += int(covered.sum())
            if k % 10 == 0 or k == n_tiles - 1:
                print(f"  tile {k + 1}/{n_tiles}", flush=True)
        dst.update_tags(blend=blend, resolution_m=res, images=len(views), models=",".join(use),
                        local_origin=json.dumps(origin.tolist()))
    summary = {"output": str(output), "images": len(views), "models": use, "resolution_m": res, "blend": blend,
               "image_source": source,
               "covered_cells": covered_cells, "covered_area_m2": float(covered_cells * res * res),
               "grid": [fine.width, fine.height], "stretch_low": lo.tolist(), "stretch_high": hi.tolist(),
               "overlap_luminance_std": overlap_std, "stats_resolution_m": stats_res,
               "focal_px": float(views[0][2]["params"][0]),
               "camera_offsets": aligned.get("camera_offsets")}
    (out_dir / "report.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
