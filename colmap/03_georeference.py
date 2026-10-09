#!/usr/bin/env python3
"""Georeference COLMAP sparse models with the gref4hsi navigation solution.

Steps
1. Interpolate the ECEF reference positions in ``Intermediate/pose.csv`` onto the
   RGB frame timestamps from ``metadata/frames.csv``.
2. Convert the positions to EPSG:25832 easting/northing plus ellipsoidal height and
   subtract a local origin so COLMAP works with small numbers.
3. Run ``colmap model_aligner`` (Sim3, robust) for every sparse model that shares
   enough images with the navigation, writing aligned models to ``workspace/aligned``.
4. Report the residual between aligned camera centres and navigation positions.

The RGB camera lever arm relative to the navigation reference is not known, so the
alignment absorbs it as a constant offset in the Sim3 fit. Orientation is taken
entirely from the photogrammetry.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
from pyproj import Transformer
from scipy.spatial.transform import Rotation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True,
                        help="Mission pipeline.json, e.g. <mission>/colmap/config/pipeline.json; all paths resolve relative to its parent folder")
    parser.add_argument("--sparse-dir", type=Path, default=None, help="Override the sparse model directory")
    parser.add_argument("--max-error", type=float, default=0.5, help="Robust alignment inlier threshold [m]")
    parser.add_argument("--estimate-offsets", dest="estimate_offsets", action="store_true", default=None,
                        help="Estimate camera-antenna lever arm and clock offset (else reuse metadata/camera_offsets.json)")
    parser.add_argument("--profile-clock", action="store_true", help="With --estimate-offsets: print rms vs fixed clock offset and exit")
    parser.add_argument("--window", type=int, default=None,
                        help="Half-width (frames) of the sliding window for per-image alignment")
    return parser.parse_args()


def resolve(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (root / path).resolve()


def read_model_txt(model_dir: Path) -> tuple[dict, dict]:
    """Return cameras and images from a COLMAP TXT model."""
    cameras = {}
    for line in (model_dir / "cameras.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        cameras[int(parts[0])] = {"model": parts[1], "width": int(parts[2]), "height": int(parts[3]),
                                  "params": [float(p) for p in parts[4:]]}
    images = {}
    lines = [l for l in (model_dir / "images.txt").read_text().splitlines() if not l.startswith("#")]
    for line, pts_line in zip(lines[0::2], lines[1::2]):
        if not line.strip():
            continue
        p = line.split()
        qw, qx, qy, qz = map(float, p[1:5])
        t = np.array(list(map(float, p[5:8])))
        R = Rotation.from_quat([qx, qy, qz, qw]).as_matrix()  # world -> camera
        pts = pts_line.split()
        point_ids = np.array([int(v) for v in pts[2::3]]) if pts else np.array([], int)
        images[p[9]] = {"image_id": int(p[0]), "R": R, "t": t, "camera_id": int(p[8]), "center": -R.T @ t,
                        "point3D_ids": point_ids[point_ids >= 0]}
    return cameras, images


def read_points3d_txt(model_dir: Path) -> dict[int, np.ndarray]:
    points = {}
    for line in (model_dir / "points3D.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        points[int(p[0])] = np.array([float(p[1]), float(p[2]), float(p[3])])
    return points


def umeyama(src: np.ndarray, dst: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Weighted similarity transform dst = s * R @ src + t (Umeyama 1991)."""
    w = weights / weights.sum()
    mu_s, mu_d = (w[:, None] * src).sum(0), (w[:, None] * dst).sum(0)
    xs, xd = src - mu_s, dst - mu_d
    cov = (w[:, None] * xd).T @ xs
    U, S, Vt = np.linalg.svd(cov)
    D = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        D[2, 2] = -1
    R = U @ D @ Vt
    var_s = (w * (xs ** 2).sum(1)).sum()
    s = np.trace(np.diag(S) @ D) / var_s
    t = mu_d - s * R @ mu_s
    return float(s), R, t


def fit_sim3_with_nadir(src_c: np.ndarray, src_axis: np.ndarray, dst_c: np.ndarray, nadir_len: float,
                        max_error: float, iterations: int = 6) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Robust Sim3 from camera centres plus virtual points along the optical axis
    that must map to straight-down (0, 0, -nadir_len) offsets in the reference frame."""
    keep = np.ones(len(src_c), bool)
    s, R, t = umeyama(src_c, dst_c, keep.astype(float))
    s = max(s, 1e-6)
    for _ in range(iterations):
        src = np.vstack([src_c[keep], src_c[keep] + (nadir_len / s) * src_axis[keep]])
        dst = np.vstack([dst_c[keep], dst_c[keep] + np.array([0.0, 0.0, -nadir_len])])
        w = np.concatenate([np.ones(keep.sum()), 0.5 * np.ones(keep.sum())])
        s, R, t = umeyama(src, dst, w)
        s = max(s, 1e-6)  # degenerate reference (e.g. navigation clamped outside its time span)
        resid = np.linalg.norm((s * (R @ src_c.T)).T + t - dst_c, axis=1)
        thresh = max(max_error, 3.0 * np.median(resid[keep]))
        new_keep = resid <= thresh
        if new_keep.sum() < 5 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return s, R, t, keep


def export_txt(model_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["colmap", "model_converter", "--input_path", str(model_dir), "--output_path", str(out_dir),
                    "--output_type", "TXT"], check=True, capture_output=True)


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    config = json.loads(config_path.read_text())
    root = config_path.parents[1]
    paths = config["paths"]
    geo = config.setdefault("georeference", {})
    pose_csv = resolve(root, geo.get("pose_csv", "../Intermediate/pose.csv"))
    epsg = int(geo.get("epsg", 25832))
    sparse_dir = args.sparse_dir or resolve(root, paths.get("sparse_fixed_dir", paths["sparse_dir"]))
    aligned_dir = resolve(root, paths.get("aligned_dir", "workspace/aligned"))
    metadata_dir = root / "metadata"
    logs = root / "logs"
    logs.mkdir(exist_ok=True)

    # --- 1. navigation onto RGB timestamps ----------------------------------
    frames = list(csv.DictReader((root / paths["frame_manifest"]).open()))
    rgb_t = np.array([float(f["timestamp_unix"]) for f in frames])
    rgb_names = [Path(f["image"]).name for f in frames]

    nav = np.genfromtxt(pose_csv, delimiter=",", skip_header=1)
    nav_t, nav_xyz = nav[:, 0], nav[:, 1:4]
    order = np.argsort(nav_t)
    nav_t, nav_xyz = nav_t[order], nav_xyz[order]
    inside = (rgb_t >= nav_t[0]) & (rgb_t <= nav_t[-1])
    print(f"RGB frames: {len(rgb_t)}, inside navigation span: {inside.sum()}, "
          f"nav span {nav_t[0]:.3f}..{nav_t[-1]:.3f}, rgb span {rgb_t[0]:.3f}..{rgb_t[-1]:.3f}")
    to_llh = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)
    to_proj = Transformer.from_crs("EPSG:4979", f"EPSG:{epsg}", always_xy=True)

    def nav_enh_at(times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """ECEF antenna position interpolated at ``times`` -> (ecef, [E, N, h_ellipsoid])."""
        ecef = np.column_stack([np.interp(times, nav_t, nav_xyz[:, i]) for i in range(3)])
        lon, lat, h = to_llh.transform(ecef[:, 0], ecef[:, 1], ecef[:, 2])
        east, north = to_proj.transform(lon, lat)[:2]
        return ecef, np.column_stack([east, north, h])

    ecef, enh = nav_enh_at(rgb_t)
    origin_file = metadata_dir / "local_origin.json"
    if origin_file.exists():
        origin = np.array(json.loads(origin_file.read_text())["origin"])
    else:
        origin = np.round(np.median(enh, axis=0), 0)
        origin_file.write_text(json.dumps({"epsg": epsg, "origin": origin.tolist(),
                                           "note": "local = [E, N, h_ellipsoid] - origin"}, indent=2))
    local = enh - origin

    def nav_local_at(times: np.ndarray) -> np.ndarray:
        return nav_enh_at(times)[1] - origin

    with (metadata_dir / "navigation.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["image", "timestamp_unix", "ecef_x", "ecef_y", "ecef_z", "east", "north", "h_ellipsoid",
                    "local_x", "local_y", "local_z", "inside_nav_span"])
        for i, name in enumerate(rgb_names):
            w.writerow([name, f"{rgb_t[i]:.6f}", *[f"{v:.4f}" for v in ecef[i]], *[f"{v:.4f}" for v in enh[i]],
                        *[f"{v:.4f}" for v in local[i]], int(inside[i])])
    ref_path = metadata_dir / "ref_images_local.txt"
    with ref_path.open("w") as fh:
        for i, name in enumerate(rgb_names):
            if inside[i]:
                fh.write(f"{name} {local[i, 0]:.4f} {local[i, 1]:.4f} {local[i, 2]:.4f}\n")
    print(f"Local origin (E, N, h): {origin.tolist()} -> {ref_path}")

    # --- 2. align every sparse model ----------------------------------------
    # The navigation has no usable height (vehicle fixed at h=0) and the track is
    # nearly straight, so camera centres alone leave the roll about the track axis
    # undetermined. We therefore fit the Sim3 with camera centres AND a near-nadir
    # viewing-direction constraint (vehicle pitch/roll are within ~2 deg).
    models = sorted(p for p in sparse_dir.iterdir() if p.is_dir() and (p / "images.bin").exists())
    if aligned_dir.exists():
        shutil.rmtree(aligned_dir)
    aligned_dir.mkdir(parents=True)
    summary = []
    nav_lookup = {name: local[i] for i, name in enumerate(rgb_names) if inside[i]}
    nadir_len = float(geo.get("nadir_constraint_length_m", 1.0))
    min_images = int(geo.get("min_images", 5))
    for model in models:
        export_txt(model, model / "txt")
        cameras, images = read_model_txt(model / "txt")
        names = [n for n in images if n in nav_lookup]
        if len(names) < min_images:
            print(f"model {model.name}: only {len(images)} images, skipped")
            summary.append({"model": model.name, "aligned": False, "images": len(images)})
            continue
        src_c = np.array([images[n]["center"] for n in names])
        src_axis = np.array([images[n]["R"].T @ np.array([0.0, 0.0, 1.0]) for n in names])
        dst_c = np.array([nav_lookup[n] for n in names])
        s, R, t, keep = fit_sim3_with_nadir(src_c, src_axis, dst_c, nadir_len, args.max_error)
        T = np.eye(4)
        T[:3, :3] = s * R
        T[:3, 3] = t
        out = aligned_dir / model.name
        out.mkdir()
        np.savetxt(out / "sim3_transform.txt", T, fmt="%.12g")
        cmd = ["colmap", "model_transformer", "--input_path", str(model), "--output_path", str(out),
               "--transform_path", str(out / "sim3_transform.txt")]
        res = subprocess.run(cmd, capture_output=True, text=True)
        (logs / f"05_align_model_{model.name}.log").write_text(" ".join(cmd) + "\n" + res.stdout + res.stderr)
        if res.returncode or not (out / "images.bin").exists():
            raise RuntimeError(f"model_transformer failed for model {model.name}")
        export_txt(out, out / "txt")
        _, aligned_images = read_model_txt(out / "txt")
        resid = np.array([aligned_images[n]["center"] - nav_lookup[n] for n in names])
        d = np.linalg.norm(resid[:, :2], axis=1)
        axes = np.array([aligned_images[n]["R"].T @ np.array([0.0, 0.0, 1.0]) for n in names])
        tilt = np.degrees(np.arccos(np.clip(-axes[:, 2], -1, 1)))
        entry = {"model": model.name, "aligned": True, "images": len(images), "used_in_fit": int(keep.sum()),
                 "scale": float(s), "residual_xy_rms_m": float(np.sqrt(np.mean(d ** 2))),
                 "residual_xy_median_m": float(np.median(d)), "residual_xy_max_m": float(d.max()),
                 "residual_z_mean_m": float(resid[:, 2].mean()), "camera_tilt_from_nadir_deg_median": float(np.median(tilt)),
                 "camera_tilt_from_nadir_deg_max": float(tilt.max())}
        summary.append(entry)
        print(f"model {model.name}: {len(images)} images ({keep.sum()} inliers), scale {s:.4f}, "
              f"xy residual rms {entry['residual_xy_rms_m']:.3f} m, median {entry['residual_xy_median_m']:.3f} m, "
              f"max {entry['residual_xy_max_m']:.3f} m; tilt from nadir median {entry['camera_tilt_from_nadir_deg_median']:.1f} deg")
    # --- 3. sliding-window alignment -----------------------------------------
    # Incremental SfM drifts strongly in scale and orientation along this
    # forward-motion sequence, so a single Sim3 cannot place the whole chain. Each
    # image is instead aligned with a Sim3 fitted on its +-window neighbours. The
    # relative orientation between neighbours comes from photogrammetry, the
    # absolute position from navigation and the tilt from the nadir constraint.
    window = args.window if args.window is not None else int(geo.get("window", 6))
    altimeter_m = float(geo.get("altimeter_height_m", 1.21))
    model_data = []
    for model in models:
        if not (model / "txt" / "images.txt").exists():
            continue
        cameras, images = read_model_txt(model / "txt")
        points = read_points3d_txt(model / "txt")
        seq = sorted(n for n in images if n in nav_lookup)
        if len(seq) >= 2 * window + 1:
            model_data.append((model.name, cameras, images, points, seq))

    def run_sliding(nav_pos: dict, full: bool = False) -> tuple[list, float]:
        """Per-image Sim3 alignment against ``nav_pos`` (name -> local xyz). Returns poses and xy rms."""
        poses = []
        for model_name, cameras, images, points, seq in model_data:
            C = np.array([images[n]["center"] for n in seq])
            A = np.array([images[n]["R"].T @ np.array([0.0, 0.0, 1.0]) for n in seq])
            N = np.array([nav_pos[n] for n in seq])
            for i, name in enumerate(seq):
                lo, hi = max(0, i - window), min(len(seq), i + window + 1)
                s, R, t, keep = fit_sim3_with_nadir(C[lo:hi], A[lo:hi], N[lo:hi], nadir_len, args.max_error)
                im = images[name]
                R_new = im["R"] @ R.T               # world' -> cam, with world' = s R world + t
                center_new = s * R @ im["center"] + t
                entry = {"image": name, "model": model_name, "R": R_new, "center": center_new,
                         "window_scale": float(s), "residual_xy_m": float(np.linalg.norm(center_new[:2] - N[i, :2])),
                         "inliers": int(keep.sum())}
                if full:
                    ids = [k for k in im["point3D_ids"] if k in points]
                    depth = (np.median(((im["R"] @ np.array([points[k] for k in ids]).T).T + im["t"])[:, 2]) * s
                             if ids else float("nan"))
                    cam = cameras[im["camera_id"]]
                    entry.update({"t": (-R_new @ center_new), "height_above_points_m": float(depth),
                                  "camera": {"model": cam["model"], "width": cam["width"], "height": cam["height"],
                                             "params": cam["params"]}})
                poses.append(entry)
        r = np.array([p["residual_xy_m"] for p in poses]) if poses else np.array([np.nan])
        return poses, float(np.sqrt(np.mean(r ** 2)))

    # 3a. first pass: cameras placed on the antenna track
    poses, rms_antenna = run_sliding(nav_lookup)
    print(f"sliding window +-{window}, antenna positions: xy residual rms {rms_antenna:.3f} m")

    # 3b. estimate the horizontal camera-antenna offset (in the camera's own frame) and a
    #     clock offset. The vehicle yaws +-30 deg every few seconds, so a lever arm makes the
    #     camera path differ in shape from the antenna path; we search the offset that lets
    #     the sliding fit explain the SfM camera path best. Headings come from the first pass
    #     (rotations are insensitive to the position error).
    offsets = {"forward_m": 0.0, "right_m": 0.0, "clock_s": 0.0, "estimated": False}
    estimate = args.estimate_offsets if args.estimate_offsets is not None else bool(geo.get("estimate_offsets", False))
    offsets_file = metadata_dir / "camera_offsets.json"
    if not estimate and offsets_file.exists():
        offsets = json.loads(offsets_file.read_text())
    if estimate or offsets.get("forward_m") or offsets.get("right_m") or offsets.get("clock_s"):
        t_of = {n: rgb_t[i] for i, n in enumerate(rgb_names)}
        fwd = {}
        for p in poses:
            f = -np.array(p["R"]).T[:, 1]          # camera "up" (-y axis) in world = vehicle forward
            f = f[:2] / max(np.linalg.norm(f[:2]), 1e-9)
            fwd[p["image"]] = f
        names_all = [p["image"] for p in poses]
        times_all = np.array([t_of[n] for n in names_all])
        F = np.array([fwd[n] for n in names_all])
        RIGHT = np.column_stack([F[:, 1], -F[:, 0]])

        def corrected_nav(forward_m: float, right_m: float, clock_s: float) -> dict:
            P = nav_local_at(times_all + clock_s)
            P[:, :2] += forward_m * F + right_m * RIGHT
            return {n: P[i] for i, n in enumerate(names_all)}

        if estimate:
            from scipy.optimize import minimize
            def cost(theta):
                return run_sliding(corrected_nav(*theta))[1]
            if args.profile_clock:
                print("clock_s  forward_m  right_m  rms_m   (lever arm re-optimised for each fixed clock offset)")
                for c in np.arange(float(geo.get("profile_clock_min_s", -6.0)), float(geo.get("profile_clock_max_s", 6.0)) + 0.01, 0.5):
                    r2 = minimize(lambda th: cost([th[0], th[1], c]), [0.0, 0.0], method="Nelder-Mead",
                                  options={"initial_simplex": [[0, 0], [0.4, 0], [0, 0.4]], "xatol": 0.01, "fatol": 1e-4})
                    print(f"{c:+6.2f}  {r2.x[0]:+8.3f}  {r2.x[1]:+7.3f}  {r2.fun:6.3f}")
                return
            x0 = np.array([offsets["forward_m"], offsets["right_m"], offsets["clock_s"]])
            simplex = np.array([x0, x0 + [0.4, 0, 0], x0 + [0, 0.4, 0], x0 + [0, 0, 0.2]])
            res = minimize(cost, x0, method="Nelder-Mead",
                           options={"initial_simplex": simplex, "xatol": 0.005, "fatol": 1e-4, "maxfev": 300})
            offsets = {"forward_m": float(res.x[0]), "right_m": float(res.x[1]), "clock_s": float(res.x[2]),
                       "estimated": True, "rms_before_m": rms_antenna, "rms_after_m": float(res.fun),
                       "evaluations": int(res.nfev),
                       "note": "camera = antenna(t + clock_s) + forward_m * image-up direction + right_m * image-right"}
            offsets_file.write_text(json.dumps(offsets, indent=2))
            print(f"estimated offsets: forward {offsets['forward_m']:+.3f} m, right {offsets['right_m']:+.3f} m, "
                  f"clock {offsets['clock_s']:+.3f} s; rms {rms_antenna:.3f} -> {res.fun:.3f} m ({res.nfev} evals)")
        nav_cam = corrected_nav(offsets["forward_m"], offsets["right_m"], offsets["clock_s"])
        with (metadata_dir / "navigation_camera.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["image", "cam_local_x", "cam_local_y", "cam_local_z"])
            for n in names_all:
                w.writerow([n, *[f"{v:.4f}" for v in nav_cam[n]]])
    else:
        nav_cam = nav_lookup

    # 3c. final pass with the corrected camera positions
    poses, rms = run_sliding(nav_cam, full=True)
    for p in poses:
        p["R"], p["t"], p["center"] = p["R"].tolist(), p["t"].tolist(), p["center"].tolist()
    if poses:
        r = np.array([p["residual_xy_m"] for p in poses])
        hgt = np.array([p["height_above_points_m"] for p in poses])
        tilt = np.degrees(np.arccos(np.clip(-np.array([np.array(p["R"]).T[:, 2][2] for p in poses]), -1, 1)))
        focal_ratio = np.nanmedian(hgt) / altimeter_m
        f_cur = poses[0]["camera"]["params"][0]
        print(f"sliding window +-{window}: {len(poses)} images, xy residual rms {rms:.3f} m, "
              f"max {r.max():.3f} m; tilt from nadir median {np.median(tilt):.1f} deg; camera height above "
              f"points median {np.nanmedian(hgt):.2f} m vs altimeter {altimeter_m} m -> focal looks like "
              f"{f_cur / focal_ratio:.0f} px (current {f_cur:.0f})")
        (metadata_dir / "aligned_poses.json").write_text(json.dumps(
            {"epsg": epsg, "origin": origin.tolist(), "window": window, "nadir_len_m": nadir_len,
             "camera_offsets": offsets, "focal_suggested_px": float(f_cur / focal_ratio), "poses": poses}))
        for e in summary:
            if e.get("aligned"):
                e["sliding_window"] = {"window": window, "residual_xy_rms_m": rms,
                                       "residual_xy_rms_antenna_m": rms_antenna, "camera_offsets": offsets,
                                       "tilt_median_deg": float(np.median(tilt)),
                                       "height_above_points_median_m": float(np.nanmedian(hgt)),
                                       "focal_suggested_px": float(f_cur / focal_ratio)}
    (metadata_dir / "alignment_report.json").write_text(json.dumps(
        {"sparse_dir": str(sparse_dir), "aligned_dir": str(aligned_dir), "epsg": epsg, "origin": origin.tolist(),
         "max_error_m": args.max_error, "models": summary}, indent=2))
    print("Aligned models written to", aligned_dir)


if __name__ == "__main__":
    main()
