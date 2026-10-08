"""Memory-conscious inspection helpers for Ecotone/Coralis UHI HDF5 files."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw


HSI = "/rawdata/hyperspectral"
RGB = "/rawdata/rgb"
IMU = "/rawdata/navigation/imu"
ALTITUDE = "/rawdata/navigation/altitude"
EXTERNAL = "/rawdata/navigation/external"


def _natural_key(path: Path) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", path.name)]


def find_h5(source: str | Path) -> list[Path]:
    """Return naturally sorted HDF5 files from a file or directory."""
    source = Path(source).expanduser()
    if source.is_file():
        if source.suffix.lower() not in {".h5", ".hdf5"}:
            raise ValueError(f"Not an HDF5 file: {source}")
        return [source]
    if not source.is_dir():
        raise FileNotFoundError(source)
    files = sorted((*source.glob("*.h5"), *source.glob("*.hdf5")), key=_natural_key)
    if not files:
        raise FileNotFoundError(f"No .h5 or .hdf5 files found in {source}")
    return files


def _plain(value):
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def file_inventory(path: str | Path) -> pd.DataFrame:
    """List groups/datasets without reading dataset contents."""
    rows: list[dict] = []
    with h5py.File(path, "r") as handle:
        rows.append({"path": "/", "kind": "group", "shape": None, "dtype": None,
                     "chunks": None, "compression": None, "attributes": dict(handle.attrs)})

        def visitor(name, obj):
            is_dataset = isinstance(obj, h5py.Dataset)
            rows.append({
                "path": f"/{name}",
                "kind": "dataset" if is_dataset else "group",
                "shape": obj.shape if is_dataset else None,
                "dtype": str(obj.dtype) if is_dataset else None,
                "chunks": obj.chunks if is_dataset else None,
                "compression": obj.compression if is_dataset else None,
                "attributes": {key: _plain(value) for key, value in obj.attrs.items()},
            })

        handle.visititems(visitor)
    return pd.DataFrame(rows)


def summarize_files(files: Iterable[str | Path]) -> pd.DataFrame:
    """Build a compact per-chunk acquisition summary."""
    rows = []
    for item in files:
        path = Path(item)
        with h5py.File(path, "r") as handle:
            cube = handle[f"{HSI}/dataCube"]
            hsi_time = handle[f"{HSI}/timestamp"]
            rgb_frames = handle.get(f"{RGB}/rgbFrames")
            imu_time = handle.get(f"{IMU}/TimestampMeasured")
            external_time = handle.get(f"{EXTERNAL}/TimestampMeasured")
            altitude_time = handle.get(f"{ALTITUDE}/TimestampMeasured")
            start = float(hsi_time[0]) if hsi_time.size else np.nan
            end = float(hsi_time[-1]) if hsi_time.size else np.nan
            rows.append({
                "file": path.name,
                "size_GB": path.stat().st_size / 1e9,
                "hsi_lines": cube.shape[0],
                "spatial_pixels": cube.shape[1],
                "spectral_bands": cube.shape[2],
                "hsi_start_unix": start,
                "hsi_end_unix": end,
                "duration_s": end - start,
                "rgb_frames": 0 if rgb_frames is None else rgb_frames.shape[0],
                "imu_samples": 0 if imu_time is None else imu_time.shape[0],
                "external_nav_samples": 0 if external_time is None else external_time.shape[0],
                "altitude_samples": 0 if altitude_time is None else altitude_time.shape[0],
            })
    result = pd.DataFrame(rows)
    result["hsi_start_utc"] = pd.to_datetime(result["hsi_start_unix"], unit="s", utc=True)
    result["hsi_end_utc"] = pd.to_datetime(result["hsi_end_unix"], unit="s", utc=True)
    return result


def read_calibration(path: str | Path) -> dict[str, object]:
    """Read the small calibration arrays and instrument metadata."""
    with h5py.File(path, "r") as handle:
        calibration = handle[f"{HSI}/calibration"]
        return {
            "file_attributes": {key: _plain(value) for key, value in handle.attrs.items()},
            "instrument_attributes": {key: _plain(value) for key, value in calibration.attrs.items()},
            "wavelength_nm": handle[f"{HSI}/calibration/spectral/band2Wavelength"][()],
            "field_of_view": handle[f"{HSI}/calibration/geometric/fieldOfView"][()],
            "dark_frame": handle[f"{HSI}/calibration/radiometric/darkFrame"][()],
            "radiometric_frame": handle[f"{HSI}/calibration/radiometric/radiometricFrame"][()],
        }


def read_telemetry(files: Iterable[str | Path], group: str) -> pd.DataFrame:
    """Concatenate all one-dimensional datasets in a telemetry group."""
    chunks = []
    for item in files:
        path = Path(item)
        with h5py.File(path, "r") as handle:
            if group not in handle:
                continue
            h5_group = handle[group]
            values = {
                name: dataset[()]
                for name, dataset in h5_group.items()
                if isinstance(dataset, h5py.Dataset) and dataset.ndim == 1
            }
            if not values:
                continue
            frame = pd.DataFrame(values)
            frame.insert(0, "source_file", path.name)
            chunks.append(frame)
    if not chunks:
        return pd.DataFrame()
    result = pd.concat(chunks, ignore_index=True)
    if "TimestampMeasured" in result:
        result["datetime_utc"] = pd.to_datetime(result["TimestampMeasured"], unit="s", utc=True)
        result = result.sort_values("TimestampMeasured", kind="stable").reset_index(drop=True)
    return result


def nearest_band_indices(wavelengths: np.ndarray, targets_nm=(640.0, 550.0, 460.0)) -> np.ndarray:
    wavelengths = np.asarray(wavelengths)
    return np.asarray([np.abs(wavelengths - target).argmin() for target in targets_nm])


def percentile_stretch(image: np.ndarray, low=2.0, high=98.0) -> np.ndarray:
    """Stretch each channel independently to [0, 1]."""
    image = np.asarray(image, dtype=np.float32)
    output = np.empty_like(image)
    channels = 1 if image.ndim == 2 else image.shape[-1]
    for channel in range(channels):
        plane = image if image.ndim == 2 else image[..., channel]
        lo, hi = np.nanpercentile(plane, (low, high))
        stretched = np.clip((plane - lo) / max(hi - lo, np.finfo(np.float32).eps), 0, 1)
        if image.ndim == 2:
            output = stretched
        else:
            output[..., channel] = stretched
    return output


def enhance_rgb(image: np.ndarray, low=1.0, high=99.5, gamma=0.65) -> np.ndarray:
    """Return a display-enhanced RGB image using stretch, white balance, and gamma.

    This is intended only for visualization. Per-channel percentile limits make
    dark, strongly tinted underwater frames easier to inspect.
    """
    stretched = percentile_stretch(image, low=low, high=high)
    corrected = np.power(stretched, gamma)
    return np.round(255 * corrected).clip(0, 255).astype(np.uint8)


def hsi_rgb_preview(
    path: str | Path,
    targets_nm=(640.0, 550.0, 460.0),
    max_lines=500,
    calibrated=True,
) -> tuple[np.ndarray, np.ndarray]:
    """Read only three bands and return a stretched along-track HSI preview."""
    with h5py.File(path, "r") as handle:
        wavelengths = handle[f"{HSI}/calibration/spectral/band2Wavelength"][()]
        indices = nearest_band_indices(wavelengths, targets_nm)
        cube = handle[f"{HSI}/dataCube"]
        step = max(1, int(np.ceil(cube.shape[0] / max_lines)))
        # h5py requires increasing fancy indices; restore requested RGB order afterwards.
        order = np.argsort(indices)
        inverse = np.argsort(order)
        data = cube[::step, :, indices[order]].astype(np.float32)[..., inverse]
        if calibrated:
            dark = handle[f"{HSI}/calibration/radiometric/darkFrame"][:, indices[order]][()][..., inverse]
            gain = handle[f"{HSI}/calibration/radiometric/radiometricFrame"][:, indices[order]][()][..., inverse]
            exposure = handle[f"{HSI}/exposureTime"][::step].astype(np.float32)
            data = (data - dark[None, :, :]) * gain[None, :, :] / (exposure[:, None, None] / 1000.0)
        return percentile_stretch(data), wavelengths[indices]


def rgb_frame(path: str | Path, index=None) -> tuple[np.ndarray, float]:
    """Read one RGB frame and its Unix timestamp."""
    with h5py.File(path, "r") as handle:
        frames = handle[f"{RGB}/rgbFrames"]
        index = frames.shape[0] // 2 if index is None else index
        return frames[index], float(handle[f"{RGB}/timestamp"][index])


def create_rgb_gif(
    files: Iterable[str | Path],
    output_path: str | Path,
    frame_step=3,
    max_frames=180,
    width=640,
    fps=8,
    enhance=True,
    low=1.0,
    high=99.5,
    gamma=0.65,
) -> Path:
    """Create a compact, timestamped GIF sampled across HDF5 RGB frames."""
    files = [Path(item) for item in files]
    candidates: list[tuple[Path, int, float]] = []
    for path in files:
        with h5py.File(path, "r") as handle:
            if f"{RGB}/rgbFrames" not in handle:
                continue
            timestamps = handle[f"{RGB}/timestamp"][()]
            candidates.extend((path, index, float(timestamps[index])) for index in range(0, len(timestamps), frame_step))
    if not candidates:
        raise ValueError("No RGB frames were found")
    if len(candidates) > max_frames:
        selected = np.linspace(0, len(candidates) - 1, max_frames, dtype=int)
        candidates = [candidates[index] for index in selected]

    images = []
    open_path = None
    handle = None
    try:
        for path, index, timestamp in candidates:
            if path != open_path:
                if handle is not None:
                    handle.close()
                handle = h5py.File(path, "r")
                open_path = path
            frame = handle[f"{RGB}/rgbFrames"][index]
            if enhance:
                frame = enhance_rgb(frame, low=low, high=high, gamma=gamma)
            image = Image.fromarray(frame).convert("RGB")
            height = round(image.height * width / image.width)
            image = image.resize((width, height), Image.Resampling.LANCZOS)
            draw = ImageDraw.Draw(image)
            label = pd.to_datetime(timestamp, unit="s", utc=True).strftime("%Y-%m-%d %H:%M:%S.%f UTC")[:-3]
            text_box = draw.textbbox((0, 0), label)
            draw.rectangle((8, 8, text_box[2] + 16, text_box[3] + 16), fill=(0, 0, 0))
            draw.text((12, 10), label, fill=(255, 255, 255))
            images.append(image)
    finally:
        if handle is not None:
            handle.close()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=round(1000 / fps),
        loop=0,
        optimize=False,
    )
    return output_path


def mean_spectrum(path: str | Path, line_step=25, sample_step=16) -> tuple[np.ndarray, np.ndarray]:
    """Compute a sampled, radiometrically calibrated mean spectrum."""
    with h5py.File(path, "r") as handle:
        cube = handle[f"{HSI}/dataCube"]
        raw = cube[::line_step, ::sample_step, :].astype(np.float32)
        dark = handle[f"{HSI}/calibration/radiometric/darkFrame"][::sample_step, :]
        gain = handle[f"{HSI}/calibration/radiometric/radiometricFrame"][::sample_step, :]
        exposure = handle[f"{HSI}/exposureTime"][::line_step].astype(np.float32)
        radiance = (raw - dark[None, :, :]) * gain[None, :, :] / (exposure[:, None, None] / 1000.0)
        wavelengths = handle[f"{HSI}/calibration/spectral/band2Wavelength"][()]
        return wavelengths, np.nanmean(radiance, axis=(0, 1))
