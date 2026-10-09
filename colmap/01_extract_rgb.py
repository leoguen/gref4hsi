#!/usr/bin/env python3
"""Extract a chronological RGB image sequence from the mission HDF5 chunks."""

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

    files = sorted(h5_dir.glob(config["paths"].get("h5_pattern", "*.h5")))  # e.g. "uhi_20241024_140601_*.h5" for one transect
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
                        Image.fromarray(frames[source_frame], mode="RGB").save(
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


if __name__ == "__main__":
    main()
