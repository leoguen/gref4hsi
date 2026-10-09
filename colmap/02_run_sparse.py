#!/usr/bin/env python3
"""Run inspectable COLMAP sparse-reconstruction stages."""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=("features", "matches", "mapper", "all", "status"),
        help="Run one stage, all stages, or report database/model status",
    )
    parser.add_argument("--config", type=Path, required=True,
                        help="Mission pipeline.json, e.g. <mission>/colmap/config/pipeline.json; all paths resolve relative to its parent folder")
    return parser.parse_args()


def resolve(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (root / path).resolve()


def database_status(database: Path) -> dict[str, int]:
    if not database.exists():
        return {"images": 0, "keypoints": 0, "descriptors": 0, "matches": 0, "two_view_geometries": 0}
    result: dict[str, int] = {}
    with sqlite3.connect(database) as connection:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        for table in ("images", "keypoints", "descriptors", "matches", "two_view_geometries"):
            result[table] = (
                connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                if table in tables else 0
            )
    return result


def report(database: Path, sparse_dir: Path) -> None:
    status = database_status(database)
    print("Database:", database)
    for key, value in status.items():
        print(f"  {key}: {value}")
    models = sorted(path for path in sparse_dir.glob("*") if path.is_dir()) if sparse_dir.exists() else []
    print(f"Sparse models: {len(models)}")
    for model in models:
        print(f"  {model}")


def run(command: list[str], log_path: Path) -> None:
    print("Running:", " ".join(command), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        stamp = datetime.now(timezone.utc).isoformat()
        log.write(f"\n[{stamp}] {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log.write(line)
        returncode = process.wait()
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    root = config_path.parents[1]
    image_dir = resolve(root, config["paths"].get("colmap_image_dir", config["paths"]["image_dir"]))
    mask_dir = resolve(root, config["paths"]["mask_dir"]) if config["paths"].get("mask_dir") else None
    database = resolve(root, config["paths"]["database"])
    options = config["colmap"]
    refine = bool(options.get("refine_intrinsics", True))
    sparse_key = "sparse_dir" if refine else "sparse_fixed_dir"
    sparse_dir = resolve(root, config["paths"].get(sparse_key, config["paths"]["sparse_dir"]))
    logs = root / "logs"

    colmap = shutil.which("colmap")
    if colmap is None:
        raise FileNotFoundError("COLMAP is not available on PATH")
    # COLMAP 3.13 renamed SiftExtraction/SiftMatching use_gpu and num_threads to FeatureExtraction/FeatureMatching
    help_text = subprocess.run([colmap, "feature_extractor", "-h"], capture_output=True, text=True)
    new_names = "FeatureExtraction.use_gpu" in help_text.stdout + help_text.stderr
    extract_prefix, match_prefix = ("FeatureExtraction", "FeatureMatching") if new_names else ("SiftExtraction", "SiftMatching")
    image_count = len(list(image_dir.glob("*.jpg")))
    if image_count == 0:
        raise FileNotFoundError(f"No JPEG images found in {image_dir}")
    database.parent.mkdir(parents=True, exist_ok=True)
    sparse_dir.mkdir(parents=True, exist_ok=True)

    if args.stage == "status":
        report(database, sparse_dir)
        return

    stages = ("features", "matches", "mapper") if args.stage == "all" else (args.stage,)
    if "features" in stages:
        camera_params = ",".join(str(value) for value in options["camera_params"])
        run(
            [
                colmap, "feature_extractor",
                "--database_path", str(database),
                "--image_path", str(image_dir),
                "--ImageReader.camera_model", str(options["camera_model"]),
                "--ImageReader.camera_params", camera_params,
                "--ImageReader.single_camera", "1" if options["single_camera"] else "0",
                f"--{extract_prefix}.use_gpu", "1" if options["use_gpu"] else "0",
                f"--{extract_prefix}.num_threads", str(options.get("num_threads", 4)),
            ],
            logs / "02_features.log",
        )
    if "matches" in stages:
        if database_status(database)["keypoints"] != image_count:
            raise RuntimeError("Feature extraction is incomplete; run the features stage first")
        run(
            [
                colmap, "sequential_matcher",
                "--database_path", str(database),
                f"--{match_prefix}.use_gpu", "1" if options["use_gpu"] else "0",
                f"--{match_prefix}.num_threads", str(options.get("num_threads", 4)),
                "--SequentialMatching.overlap", str(options["sequential_overlap"]),
            ],
            logs / "03_matches.log",
        )
    if "mapper" in stages:
        if database_status(database)["two_view_geometries"] == 0:
            raise RuntimeError("No verified image pairs; run the matches stage first")
        run(
            [
                colmap, "mapper",
                "--database_path", str(database),
                "--image_path", str(image_dir),
                "--output_path", str(sparse_dir),
                "--Mapper.num_threads", str(options.get("num_threads", 4)),
                "--Mapper.ba_refine_focal_length", "1" if refine else "0",
                "--Mapper.ba_refine_principal_point", "0",
                "--Mapper.ba_refine_extra_params", "1" if refine else "0",
            ] + [str(v) for kv in options.get("mapper_extra_args", {}).items() for v in (f"--Mapper.{kv[0]}", kv[1])],
            logs / ("04_mapper.log" if refine else "04_mapper_fixed_intrinsics.log"),
        )
    report(database, sparse_dir)


if __name__ == "__main__":
    main()
