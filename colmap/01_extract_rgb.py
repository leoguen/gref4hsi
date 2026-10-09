#!/usr/bin/env python3
"""Extract a chronological RGB image sequence from the mission HDF5 chunks.

Optionally (``extraction.normalize = true``) also writes a second image set for COLMAP:
a per-pixel grey-world normalisation (Løvås et al. 2022, eq. 9: every pixel position
and channel is standardised with its mean and standard deviation over the whole
transect, then mapped to a common target mean/std) plus a feature mask that blanks
the dark housing rim. The raw frames stay in ``images/extracted`` for the mosaic.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import h5py
import numpy as np
from PIL import Image


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def grayworld_stretch(frame: np.ndarray, low_pct: float = 1.0, high_pct: float = 99.5, max_gain: float = 16.0) -> np.ndarray:
    """Gray-world white balance followed by a percentile contrast stretch (dark underwater frames).

    Statistics ignore the black vignette/border pixels. The stretch gain is capped to keep sensor noise in check."""
    img = frame.astype(np.float32)
    valid = img.max(axis=2) > 3
    if valid.sum() < 100:
        return frame
    means = img[valid].mean(axis=0)
    img *= (means.mean() / np.maximum(means, 1e-3))[None, None, :]
    lum = img.mean(axis=2)[valid]
    lo, hi = np.percentile(lum, [low_pct, high_pct])
    gain = min(255.0 / max(hi - lo, 1.0), max_gain)
    img = (img - lo) * gain
    img[~valid] = 0
    return np.clip(img, 0, 255).astype(np.uint8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True,
                        help="Mission pipeline.json, e.g. <mission>/colmap/config/pipeline.json; all paths resolve relative to its parent folder")
    parser.add_argument("--dry-run", action="store_true", help="Inspect inputs without writing images")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing extracted images")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    config = load_config(config_path)
    root = config_path.parents[1]
    h5_dir = resolve(root, config["paths"]["h5_dir"])
    image_dir = resolve(root, config["paths"]["image_dir"])
    manifest_path = resolve(root, config["paths"]["frame_manifest"])
    frame_dataset = config["hdf5"]["rgb_frames"]
    timestamp_dataset = config["hdf5"]["rgb_timestamps"]
    quality = int(config["extraction"]["jpeg_quality"])
    enhance = config["extraction"].get("enhance")   # optional: "grayworld"

    patterns = config["paths"].get("h5_pattern", "*.h5")  # e.g. "uhi_20241024_140601_*.h5" for one transect, or a list for several
    patterns = [patterns] if isinstance(patterns, str) else patterns
    files = sorted({path for pattern in patterns for path in h5_dir.glob(pattern)})
    if not files:
        raise FileNotFoundError(f"No HDF5 files found in {h5_dir}")

    inventory: list[tuple[Path, int, tuple[int, ...], float, float]] = []
    previous_end = -np.inf
    for h5_path in files:
        with h5py.File(h5_path, "r") as handle:
            if frame_dataset not in handle or timestamp_dataset not in handle:
                raise KeyError(f"Missing RGB datasets in {h5_path.name}")
            frames = handle[frame_dataset]
            timestamps = handle[timestamp_dataset]
            if frames.ndim != 4 or frames.shape[-1] != 3 or frames.dtype != np.uint8:
                raise ValueError(f"Unexpected RGB array in {h5_path.name}: {frames.shape}, {frames.dtype}")
            if len(frames) != len(timestamps) or not len(frames):
                raise ValueError(f"Frame/timestamp mismatch in {h5_path.name}")
            times = timestamps[()]
            if np.any(np.diff(times) <= 0):
                raise ValueError(f"Timestamps are not strictly increasing in {h5_path.name}")
            start, end = float(times[0]), float(times[-1])
            if start <= previous_end:
                raise ValueError(f"Timestamp overlap before {h5_path.name}")
            previous_end = end
            inventory.append((h5_path, len(frames), tuple(frames.shape[1:]), start, end))

    total = sum(item[1] for item in inventory)
    print(f"Input directory: {h5_dir}")
    for h5_path, count, shape, start, end in inventory:
        print(f"  {h5_path.name}: {count} frames, {shape}, {start:.6f} .. {end:.6f}")
    print(f"Total: {total} frames")
    if args.dry_run:
        return

    image_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_manifest = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
    fieldnames = ["sequence", "image", "timestamp_unix", "source_h5", "source_frame"]

    written = skipped = 0
    sequence = 0
    with temporary_manifest.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for h5_path, _, _, _, _ in inventory:
            with h5py.File(h5_path, "r") as handle:
                frames = handle[frame_dataset]
                timestamps = handle[timestamp_dataset]
                for source_frame in range(len(frames)):
                    timestamp = float(timestamps[source_frame])
                    filename = f"frame_{sequence:06d}_{timestamp:.6f}.jpg"
                    destination = image_dir / filename
                    if destination.exists() and not args.overwrite:
                        skipped += 1
                    else:
                        temporary_image = destination.with_suffix(".jpg.tmp")
                        Image.fromarray(grayworld_stretch(frames[source_frame]) if enhance == "grayworld" else frames[source_frame], mode="RGB").save(
                            temporary_image, format="JPEG", quality=quality, subsampling=0
                        )
                        os.replace(temporary_image, destination)
                        written += 1
                    writer.writerow(
                        {
                            "sequence": sequence,
                            "image": destination.relative_to(root),
                            "timestamp_unix": f"{timestamp:.9f}",
                            "source_h5": h5_path.name,
                            "source_frame": source_frame,
                        }
                    )
                    sequence += 1

    os.replace(temporary_manifest, manifest_path)
    print(f"Images: {image_dir}")
    print(f"Manifest: {manifest_path}")
    print(f"Written: {written}; already present: {skipped}")

    if config["extraction"].get("normalize", False):
        normalize_images(config, root, image_dir, quality)


def normalize_images(config: dict, root: Path, image_dir: Path, quality: int) -> None:
    """Per-pixel grey-world normalisation and vignette mask for feature matching."""
    import cv2

    ex = config["extraction"]
    out_dir = resolve(root, config["paths"].get("colmap_image_dir", "images/normalized"))
    mask_dir = resolve(root, config["paths"].get("mask_dir", "images/masks"))
    out_dir.mkdir(parents=True, exist_ok=True)
    mask_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(image_dir.glob("*.jpg"))
    stack = np.stack([cv2.imread(str(f)).astype(np.float32) for f in files])  # N,H,W,3 (BGR)
    mu = stack.mean(axis=0)
    sd = stack.std(axis=0)
    sigma_px = float(ex.get("normalize_stats_sigma_px", 5))
    if sigma_px > 0:  # smooth the statistics so fixed-pattern noise is not amplified
        k = int(6 * sigma_px + 1) | 1
        mu = cv2.GaussianBlur(mu, (k, k), sigma_px)
        sd = cv2.GaussianBlur(sd, (k, k), sigma_px)
    sd = np.maximum(sd, float(ex.get("normalize_min_std", 2.0)))
    target_mean = float(ex.get("normalize_target_mean", 110.0))
    target_std = float(ex.get("normalize_target_std", 40.0))

    # usable area: mean luminance above a fraction of the central brightness, eroded
    lum = cv2.GaussianBlur(mu.mean(axis=2), (21, 21), 5)
    H, W = lum.shape
    centre = lum[H // 2 - H // 8:H // 2 + H // 8, W // 2 - W // 8:W // 2 + W // 8].mean()
    usable = (lum >= float(ex.get("vignette_threshold", 0.3)) * centre).astype(np.uint8)
    erode = int(ex.get("vignette_erode_px", 8))
    if erode > 0:
        usable = cv2.erode(usable, np.ones((2 * erode + 1,) * 2, np.uint8))
    mask = usable * 255

    for f, img in zip(files, stack):
        out = (img - mu) / sd * target_std + target_mean
        out[usable == 0] = 0
        cv2.imwrite(str(out_dir / f.name), np.clip(out, 0, 255).astype(np.uint8),
                    [cv2.IMWRITE_JPEG_QUALITY, quality])
        cv2.imwrite(str(mask_dir / (f.name + ".png")), mask)  # COLMAP: <image name>.png, 0 = ignore
    cv2.imwrite(str(mask_dir.parent / "normalize_mean.png"), np.clip(mu, 0, 255).astype(np.uint8))
    print(f"Normalised images: {out_dir} ({len(files)} files, target mean {target_mean}, std {target_std}, "
          f"stats smoothed sigma {sigma_px} px)")
    print(f"Feature masks: {mask_dir} (usable area {usable.mean() * 100:.1f}%)")


if __name__ == "__main__":
    main()
