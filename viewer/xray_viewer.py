from __future__ import annotations

import argparse
from collections import deque
import json
import math
import os
import subprocess
import struct
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageTk

from xray_physics import ConeBeamCamera, close_planar_boundaries, path_lengths
from wire_raster import wire_transmission


APP_TITLE = "Guidewire Projection Viewer"
BG = (139, 143, 145)
PANEL_BG = "#151a1d"
TEXT = "#d9e0e3"
MUTED = "#829096"
ACCENT = "#d7e9ee"
DISPLAY_ORIGIN = np.zeros(3, dtype=np.float64)
IODINE_DELTA_MU_PER_MM = 0.08


@dataclass(frozen=True)
class ViewerConfig:
    display_scale: float = 0.5
    source_object_mm: float = 750.0
    source_detector_mm: float = 1000.0
    detector_width_mm: float = 150.0
    detector_height_mm: float = 150.0
    detector_pixel_pitch_mm: float = 0.154
    matrix: int = 974
    wire_display_diameter_mm: float = 0.30
    wire_attenuation_coefficient_per_mm: float = 4.0
    wire_depth_attenuation_per_mm: float = 0.0
    vessel_depth_attenuation_per_mm: float = 0.0
    view_mode: str = "overview"
    overview_fov_object_mm: float = 112.5
    tip_tracking_fov_object_mm: float = 100.0
    frame_rate_fps: float = 7.5
    display_delay_ms: float = 0.0
    show_reaction_force: bool = True
    show_contact_points: bool = False
    force_bar_max_n: float = 10.0


def load_viewer_config(path: Path | None = None) -> ViewerConfig:
    path = path or Path(os.environ.get("GUIDEWIRE_VIEWER_CONFIG",
                                      Path(__file__).resolve().with_name("viewer_config.json")))
    defaults = ViewerConfig()
    try:
        values = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(values, dict):
            raise ValueError("configuration root must be an object")
    except FileNotFoundError:
        values = {}
    except (OSError, json.JSONDecodeError, ValueError) as error:
        print(f"Viewer config ignored ({error}); using defaults", flush=True)
        values = {}

    xray = values.get("xray", {})
    wire = values.get("wire", {})
    vessel = values.get("vessel", {})
    view = values.get("view", {})
    ui = values.get("ui", {})
    if not isinstance(xray, dict):
        xray = {}
    if not isinstance(wire, dict):
        wire = {}
    if not isinstance(vessel, dict):
        vessel = {}
    if not isinstance(view, dict):
        view = {}
    if not isinstance(ui, dict):
        ui = {}

    def positive_float(section: dict, key: str, default: float) -> float:
        try:
            value = float(section.get(key, default))
            return value if math.isfinite(value) and value > 0 else default
        except (TypeError, ValueError):
            return default

    def positive_root_float(key: str, default: float) -> float:
        try:
            value = float(values.get(key, default))
            return value if math.isfinite(value) and value > 0 else default
        except (TypeError, ValueError):
            return default

    def nonnegative_float(section: dict, key: str, default: float) -> float:
        try:
            value = float(section.get(key, default))
            return value if math.isfinite(value) and value >= 0 else default
        except (TypeError, ValueError):
            return default

    source_object_mm = positive_float(xray, "source_object_mm", defaults.source_object_mm)
    source_detector_mm = positive_float(xray, "source_detector_mm", defaults.source_detector_mm)
    detector_pixel_pitch_mm = positive_float(
        xray, "detector_pixel_pitch_mm", defaults.detector_pixel_pitch_mm)
    overview_fov_object_mm = positive_float(
        view, "overview_fov_object_mm", defaults.overview_fov_object_mm)
    tip_tracking_fov_object_mm = positive_float(
        view, "tip_tracking_fov_object_mm", defaults.tip_tracking_fov_object_mm)
    detector_width_mm = overview_fov_object_mm * source_detector_mm / source_object_mm
    detector_height_mm = detector_width_mm
    matrix = max(1, int(round(detector_width_mm / detector_pixel_pitch_mm)))
    view_mode = view.get("mode", defaults.view_mode)
    if not isinstance(view_mode, str) or view_mode not in {"overview", "tip_tracking"}:
        view_mode = defaults.view_mode

    return ViewerConfig(
        display_scale=positive_root_float("display_scale", defaults.display_scale),
        source_object_mm=source_object_mm,
        source_detector_mm=source_detector_mm,
        detector_width_mm=detector_width_mm,
        detector_height_mm=detector_height_mm,
        detector_pixel_pitch_mm=detector_pixel_pitch_mm,
        matrix=matrix,
        wire_display_diameter_mm=positive_float(wire, "display_diameter_mm", defaults.wire_display_diameter_mm),
        wire_attenuation_coefficient_per_mm=positive_float(
            wire, "attenuation_coefficient_per_mm", defaults.wire_attenuation_coefficient_per_mm),
        wire_depth_attenuation_per_mm=nonnegative_float(
            wire, "depth_attenuation_per_mm", defaults.wire_depth_attenuation_per_mm),
        vessel_depth_attenuation_per_mm=nonnegative_float(
            vessel, "depth_attenuation_per_mm", defaults.vessel_depth_attenuation_per_mm),
        view_mode=view_mode,
        overview_fov_object_mm=overview_fov_object_mm,
        tip_tracking_fov_object_mm=tip_tracking_fov_object_mm,
        frame_rate_fps=positive_float(ui, "frame_rate_fps", defaults.frame_rate_fps),
        display_delay_ms=nonnegative_float(ui, "display_delay_ms", defaults.display_delay_ms),
        show_reaction_force=(ui.get("show_reaction_force", defaults.show_reaction_force)
                             is True),
        show_contact_points=(ui.get("show_contact_points", defaults.show_contact_points)
                             is True),
        force_bar_max_n=positive_float(ui, "force_bar_max_n", defaults.force_bar_max_n),
    )


VIEWER_CONFIG = load_viewer_config()
GUIDEWIRE_DIAMETER_MM = VIEWER_CONFIG.wire_display_diameter_mm
WIRE_DELTA_MU_PER_MM = VIEWER_CONFIG.wire_attenuation_coefficient_per_mm
DISPLAY_SCALE = VIEWER_CONFIG.display_scale
SOURCE_OBJECT_MM = VIEWER_CONFIG.source_object_mm
SOURCE_DETECTOR_MM = VIEWER_CONFIG.source_detector_mm
XRAY_MATRIX = VIEWER_CONFIG.matrix
DETECTOR_WIDTH_MM = VIEWER_CONFIG.detector_width_mm
DETECTOR_HEIGHT_MM = VIEWER_CONFIG.detector_height_mm
DETECTOR_PIXEL_PITCH_MM = VIEWER_CONFIG.detector_pixel_pitch_mm


def to_view_coordinates(points: np.ndarray) -> np.ndarray:
    """Map SOFA world coordinates into viewer-only display coordinates."""
    points = np.asarray(points, dtype=np.float64)
    return DISPLAY_ORIGIN + DISPLAY_SCALE * (points - DISPLAY_ORIGIN)


def _font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    names = [
        "C:/Windows/Fonts/seguisb.ttf" if bold else "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/meiryo.ttc",
    ]
    for name in names:
        try:
            return ImageFont.truetype(name, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def load_stl(path: str | Path) -> np.ndarray:
    """Return STL triangles as (n, 3, 3) float array."""
    path = Path(path)
    raw = path.read_bytes()
    if len(raw) >= 84:
        count = struct.unpack_from("<I", raw, 80)[0]
        expected = 84 + count * 50
        if expected == len(raw):
            dtype = np.dtype(
                [
                    ("normal", "<f4", (3,)),
                    ("vertices", "<f4", (3, 3)),
                    ("attr", "<u2"),
                ]
            )
            records = np.frombuffer(raw, dtype=dtype, count=count, offset=84)
            return np.asarray(records["vertices"], dtype=np.float64)

    vertices: list[list[float]] = []
    for line in raw.decode("utf-8", errors="ignore").splitlines():
        fields = line.strip().split()
        if len(fields) == 4 and fields[0].lower() == "vertex":
            vertices.append([float(fields[1]), float(fields[2]), float(fields[3])])
    if not vertices or len(vertices) % 3:
        raise ValueError(f"STLを読み込めません: {path}")
    return np.asarray(vertices, dtype=np.float64).reshape(-1, 3, 3)


def load_wire_log(path: str | Path, max_frames: int = 5000) -> tuple[np.ndarray, list[np.ndarray]]:
    """Read BeamAdapter T=/X= position logs."""
    path = Path(path)
    times: list[float] = []
    frames: list[np.ndarray] = []
    current_time = 0.0
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            text = line.strip()
            if text.startswith("T="):
                try:
                    current_time = float(text.split("=", 1)[1])
                except ValueError:
                    pass
            elif text.startswith("X="):
                values = np.fromstring(text.split("=", 1)[1], sep=" ")
                if values.size >= 6 and values.size % 3 == 0:
                    times.append(current_time)
                    frames.append(values.reshape(-1, 3))
    if not frames:
        raise ValueError(f"有効なX=フレームがありません: {path}")
    if len(frames) > max_frames:
        indices = np.linspace(0, len(frames) - 1, max_frames).round().astype(int)
        frames = [frames[i] for i in indices]
        times = [times[i] for i in indices]
    return np.asarray(times, dtype=np.float64), frames


def read_latest_wire_frame(path: str | Path) -> tuple[float, np.ndarray] | None:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return None
    with path.open("rb") as stream:
        size = stream.seek(0, os.SEEK_END)
        stream.seek(max(0, size - 262144))
        block = stream.read().decode("utf-8", errors="ignore")
    lines = block.splitlines()
    current_time = 0.0
    latest: np.ndarray | None = None
    latest_time = 0.0
    for line in lines:
        text = line.strip()
        if text.startswith("T="):
            try:
                current_time = float(text.split("=", 1)[1])
            except ValueError:
                pass
        elif text.startswith("X="):
            values = np.fromstring(text.split("=", 1)[1], sep=" ")
            if values.size >= 6 and values.size % 3 == 0:
                latest = values.reshape(-1, 3)
                latest_time = current_time
    if latest is None:
        return None
    return latest_time, latest


def clean_wire(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=np.float64)
    points = points[np.isfinite(points).all(axis=1)]
    if len(points) <= 1:
        return points
    keep = np.ones(len(points), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(points, axis=0), axis=1) > 1e-7
    result = points[keep]
    return result


def read_live_snapshot(path: str | Path) -> tuple[float, float, float, np.ndarray] | None:
    """Read the compact snapshot produced by interactive_xray_scene.py."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return None
    if len(lines) < 2:
        return None
    header = lines[0].split()
    values = np.fromstring(lines[1], sep=" ")
    if len(header) != 3 or values.size < 6 or values.size % 3:
        return None
    try:
        return float(header[0]), float(header[1]), float(header[2]), values.reshape(-1, 3)
    except ValueError:
        return None


@dataclass
class Projection:
    center: np.ndarray
    horizontal: np.ndarray
    vertical: np.ndarray
    depth: np.ndarray
    angle_deg: float = 0.0
    flipped: bool = False

    @classmethod
    def from_vertices(cls, vertices: np.ndarray) -> "Projection":
        pts = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
        center = (pts.min(axis=0) + pts.max(axis=0)) * 0.5
        centered = pts - pts.mean(axis=0)
        covariance = centered.T @ centered / max(1, len(centered) - 1)
        values, vectors = np.linalg.eigh(covariance)
        order = np.argsort(values)[::-1]
        vertical = vectors[:, order[0]]
        horizontal = vectors[:, order[1]]
        depth = vectors[:, order[2]]
        if np.dot(np.cross(horizontal, vertical), depth) < 0:
            depth = -depth
        return cls(center=center, horizontal=horizontal, vertical=vertical, depth=depth)

    def basis(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        angle = math.radians(self.angle_deg)
        horizontal = math.cos(angle) * self.horizontal + math.sin(angle) * self.depth
        depth = -math.sin(angle) * self.horizontal + math.cos(angle) * self.depth
        vertical = -self.vertical if self.flipped else self.vertical
        return horizontal, vertical, depth

    def project(self, points: np.ndarray) -> np.ndarray:
        horizontal, vertical, depth = self.basis()
        relative = np.asarray(points, dtype=np.float64) - self.center
        return np.column_stack((relative @ horizontal, relative @ vertical, relative @ depth))


class XRayRenderer:
    def __init__(self, triangles: np.ndarray):
        self.triangles = to_view_coordinates(triangles)
        self.closed_triangles = close_planar_boundaries(self.triangles)
        self.projection = Projection.from_vertices(self.triangles)
        self.overview_center = self.projection.center.copy()
        self.view_mode = VIEWER_CONFIG.view_mode
        self.object_fov_mm = self._fov_for_mode(self.view_mode)
        self.detector_fov_mm = self.object_fov_mm * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM
        self.base_detector_fov_mm = (VIEWER_CONFIG.overview_fov_object_mm
                                     * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM)
        self.last_tip_position: np.ndarray | None = None
        self._cached_key: tuple | None = None
        self._cached_vessel: np.ndarray | None = None
        self._cached_camera: ConeBeamCamera | None = None
        self._log_view_state()

    @staticmethod
    def _fov_for_mode(mode: str) -> float:
        return (VIEWER_CONFIG.overview_fov_object_mm if mode == "overview"
                else VIEWER_CONFIG.tip_tracking_fov_object_mm)

    def set_view_mode(self, mode: str) -> None:
        if mode not in {"overview", "tip_tracking"}:
            return
        previous_base_fov = self.base_detector_fov_mm
        self.view_mode = mode
        self.object_fov_mm = self._fov_for_mode(mode)
        self.detector_fov_mm = self.object_fov_mm * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM
        self.base_detector_fov_mm = (VIEWER_CONFIG.overview_fov_object_mm
                                     * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM)
        if not math.isclose(previous_base_fov, self.base_detector_fov_mm):
            self.invalidate()
        if mode == "overview":
            self.projection.center = self.overview_center.copy()
        self._log_view_state()

    def update_tracking_tip(self, points: np.ndarray) -> None:
        if not len(points):
            return
        tip = np.asarray(points[-1], dtype=np.float64)
        if not np.isfinite(tip).all():
            return
        self.last_tip_position = tip.copy()

    def _log_view_state(self) -> None:
        overview_detector_fov = (VIEWER_CONFIG.overview_fov_object_mm
                                 * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM)
        tracking_detector_fov = (VIEWER_CONFIG.tip_tracking_fov_object_mm
                                 * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM)
        overview_matrix = max(1, int(round(
            overview_detector_fov / DETECTOR_PIXEL_PITCH_MM)))
        tracking_crop_pixels = max(1, int(round(
            tracking_detector_fov / DETECTOR_PIXEL_PITCH_MM)))
        print(
            f"Viewer: DISPLAY_SCALE={DISPLAY_SCALE:g}, mode={self.view_mode}, "
            f"object FOV={self.object_fov_mm:.3f} mm, "
            f"detector FOV={self.detector_fov_mm:.3f} mm, "
            f"detector pixel pitch={DETECTOR_PIXEL_PITCH_MM:.3f} mm/pixel, "
            f"overview object FOV={VIEWER_CONFIG.overview_fov_object_mm:.3f} mm, "
            f"overview detector-side FOV={overview_detector_fov:.3f} mm, "
            f"overview matrix={overview_matrix} x {overview_matrix}, "
            f"tip tracking object FOV={VIEWER_CONFIG.tip_tracking_fov_object_mm:.3f} mm, "
            f"tip tracking detector-side FOV={tracking_detector_fov:.3f} mm, "
            f"tip tracking crop pixel count={tracking_crop_pixels}, "
            f"SOD={SOURCE_OBJECT_MM:.3f} mm, SID={SOURCE_DETECTOR_MM:.3f} mm, "
            f"projection center={np.array2string(self.projection.center, precision=3)}, "
            f"guidewire tip={('unavailable' if self.last_tip_position is None else np.array2string(self.last_tip_position, precision=3))}",
            flush=True,
        )

    def invalidate(self) -> None:
        self._cached_key = None
        self._cached_vessel = None
        self._cached_camera = None

    def _camera(self, center: np.ndarray | None = None) -> ConeBeamCamera:
        horizontal, vertical, forward = self.projection.basis()
        camera_center = self.projection.center if center is None else center
        return ConeBeamCamera(
            source=camera_center - SOURCE_OBJECT_MM * forward,
            horizontal=horizontal,
            vertical=vertical,
            forward=forward,
            source_object_mm=SOURCE_OBJECT_MM,
            source_detector_mm=SOURCE_DETECTOR_MM,
            detector_width_mm=self.base_detector_fov_mm,
            detector_height_mm=self.base_detector_fov_mm,
            width=XRAY_MATRIX,
            height=XRAY_MATRIX,
        )

    def _vessel_projection(self) -> tuple[np.ndarray, ConeBeamCamera]:
        # Roadmap geometry is fixed; tip tracking only crops its composited image.
        key = (round(self.projection.angle_deg, 3), self.projection.flipped,
               tuple(np.round(self.overview_center, 6)),
               round(self.base_detector_fov_mm, 6), SOURCE_OBJECT_MM,
               SOURCE_DETECTOR_MM, VIEWER_CONFIG.vessel_depth_attenuation_per_mm)
        if (self._cached_key == key and self._cached_vessel is not None
                and self._cached_camera is not None):
            return self._cached_vessel.copy(), self._cached_camera
        camera = self._camera(center=self.overview_center)
        depth_mu = VIEWER_CONFIG.vessel_depth_attenuation_per_mm
        if depth_mu > 0:
            vessel_mm, first_hit = path_lengths(self.closed_triangles, camera,
                                                 return_first_hit=True)
            depth_tau = depth_mu * np.where(np.isfinite(first_hit), first_hit, 0.0)
        else:
            vessel_mm = path_lengths(self.closed_triangles, camera)
            depth_tau = 0.0
        transmission = np.exp(-IODINE_DELTA_MU_PER_MM * vessel_mm - depth_tau)
        # Inverted contrast-only DSA map: neutral grey background, white lumen.
        vessel = np.clip(150.0 + 105.0 * (1.0 - transmission), 0, 255).astype(np.uint8)
        self._cached_key = key
        self._cached_vessel = vessel.copy()
        self._cached_camera = camera
        return vessel, camera

    @staticmethod
    def _detector_pixels(points: np.ndarray, camera: ConeBeamCamera,
                         offset_x: float, offset_y: float, side: float,
                         crop: tuple[int, int, int] | None = None) -> np.ndarray:
        uv = camera.project_mm(points)[:, :2]
        pixels = np.empty((len(uv), 2), dtype=np.float64)
        detector_x = (uv[:, 0] / camera.detector_width_mm + 0.5) * camera.width
        detector_y = (0.5 - uv[:, 1] / camera.detector_height_mm) * camera.height
        if crop is None:
            crop_left, crop_top, crop_size = 0, 0, camera.width
        else:
            crop_left, crop_top, crop_size = crop
        pixels[:, 0] = offset_x + (detector_x - crop_left) / crop_size * side
        pixels[:, 1] = offset_y + (detector_y - crop_top) / crop_size * side
        return pixels

    def _tracking_crop(self, image: Image.Image, camera: ConeBeamCamera,
                       tip: np.ndarray | None) -> tuple[Image.Image, tuple[int, int, int]]:
        center = tip if tip is not None else self.overview_center
        detector_uv = camera.project_mm(np.asarray(center, dtype=np.float64)[None, :])[0, :2]
        center_x = (detector_uv[0] / camera.detector_width_mm + 0.5) * camera.width
        center_y = (0.5 - detector_uv[1] / camera.detector_height_mm) * camera.height
        crop_size = max(1, int(round(self.detector_fov_mm
                                     / self.base_detector_fov_mm * camera.width)))
        left = int(round(center_x - crop_size * 0.5))
        top = int(round(center_y - crop_size * 0.5))
        crop = Image.new("RGB", (crop_size, crop_size), (150, 150, 150))
        src_left, src_top = max(0, left), max(0, top)
        src_right = min(camera.width, left + crop_size)
        src_bottom = min(camera.height, top + crop_size)
        if src_left < src_right and src_top < src_bottom:
            region = image.crop((src_left, src_top, src_right, src_bottom))
            crop.paste(region, (src_left - left, src_top - top))
        return crop, (left, top, crop_size)

    @staticmethod
    def _wire_absorption(points: np.ndarray, camera: ConeBeamCamera) -> np.ndarray:
        """Integrate the union of wire intervals at four detector sub-samples."""
        return wire_transmission(points, camera, XRAY_MATRIX, GUIDEWIRE_DIAMETER_MM,
                                 WIRE_DELTA_MU_PER_MM,
                                 VIEWER_CONFIG.wire_depth_attenuation_per_mm)

    def render(
        self,
        width: int,
        height: int,
        wire: np.ndarray | None,
        frame_time: float = 0.0,
        frame_index: int = 0,
        frame_count: int = 0,
        show_nodes: bool = False,
        show_roadmap: bool = True,
    ) -> Image.Image:
        ss = 1
        points = (to_view_coordinates(clean_wire(wire)) if wire is not None
                  else np.empty((0, 3)))
        self.update_tracking_tip(points)
        vessel, camera = self._vessel_projection()
        detector = vessel.astype(np.float32) if show_roadmap else np.full_like(vessel, 150, dtype=np.float32)
        if len(points) >= 2:
            detector *= self._wire_absorption(points, camera)
        detector = np.clip(detector, 0, 255).astype(np.uint8)
        detector_image = Image.fromarray(detector, mode="L").convert("RGB")
        crop = None
        if self.view_mode == "tip_tracking":
            detector_image, crop = self._tracking_crop(
                detector_image, camera, self.last_tip_position)

        canvas = Image.new("RGB", (width * ss, height * ss), BG)
        side = min(width, height) * ss
        offset_x = (width * ss - side) * 0.5
        offset_y = (height * ss - side) * 0.5
        # Reuse this exact camera/crop for display-only diagnostic overlays.
        self.last_display_geometry = (camera, offset_x, offset_y, side, crop)
        # Monitor resampling only; keep detector transmission samples unchanged.
        display_filter = (Image.Resampling.BOX if side < detector_image.width
                          else Image.Resampling.BILINEAR)
        panel = detector_image.resize((side, side), display_filter)
        canvas.paste(panel, (round(offset_x), round(offset_y)))
        canvas = canvas.convert("RGBA")
        draw = ImageDraw.Draw(canvas, "RGBA")

        if show_nodes and len(points) >= 2:
            pixels = self._detector_pixels(points, camera, offset_x, offset_y, side, crop)
            radius = 1.8 * ss
            for x, y in pixels:
                draw.ellipse((x-radius, y-radius, x+radius, y+radius),
                             outline=(19, 66, 96, 220), width=max(1, ss))

        small_font = _font(10 * ss)
        object_field_mm = self.object_fov_mm
        scale = side / object_field_mm
        target_px = 90 * ss
        target_units = target_px / max(scale, 1e-9)
        magnitude = 10 ** math.floor(math.log10(max(target_units, 1e-9)))
        candidates = np.asarray([1, 2, 5, 10], dtype=float) * magnitude
        bar_units = float(candidates[np.argmin(np.abs(candidates - target_units))])
        bar_px = bar_units * scale
        x2 = offset_x + side - 22 * ss
        x1 = x2 - bar_px
        y = height * ss - 27 * ss
        draw.line((x1, y, x2, y), fill=(30, 33, 35, 220), width=2 * ss)
        draw.line((x1, y-4*ss, x1, y+4*ss), fill=(30, 33, 35, 220), width=2 * ss)
        draw.line((x2, y-4*ss, x2, y+4*ss), fill=(30, 33, 35, 220), width=2 * ss)
        label = f"{bar_units:g} mm"
        bbox = draw.textbbox((0, 0), label, font=small_font)
        draw.text(((x1+x2-(bbox[2]-bbox[0]))/2, y-18*ss), label, font=small_font, fill=(37, 40, 42, 210))

        result = canvas.convert("RGB")
        return result

    def project_diagnostic_points(self, sofa_points: np.ndarray) -> tuple[np.ndarray, tuple]:
        camera, offset_x, offset_y, side, crop = self.last_display_geometry
        pixels = self._detector_pixels(to_view_coordinates(sofa_points), camera,
                                       offset_x, offset_y, side, crop)
        return pixels, (offset_x, offset_y, offset_x + side, offset_y + side)


class ViewerApp:
    def __init__(self, root, stl_path: Path, wire_path: Path | None):
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk

        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.root = root
        self.root.title(APP_TITLE)
        self.root.configure(bg="#0b0f11")
        self.root.geometry("1180x820")
        self.root.minsize(900, 650)

        self.base_dir = Path(__file__).resolve().parent
        self.runtime_dir = Path(os.environ.get("GUIDEWIRE_RUNTIME_DIR", self.base_dir / "runtime"))
        self.runtime_dir.mkdir(parents=True, exist_ok=True)
        self.command_path = self.runtime_dir / "control.json"
        self.force_state_path = self.runtime_dir / "live_force_state.json"
        self.force_state_mtime = -1
        self.latest_force_state = None
        self.force_state_buffer = deque()
        self.last_force_received = 0.0
        self.force_session_started_at = 0.0
        # Each SOFA launch gets a new snapshot path, so a previous run can
        # never appear during this run's startup.
        self.snapshot_path = self.runtime_dir / f"live_wire_{time.time_ns()}.txt"
        self.scene_path = self.base_dir / "interactive_xray_scene.py"
        self.engine_path = self.base_dir / "native_engine.py"
        self.sim_process: subprocess.Popen | None = None
        self.sofa_log_stream = None
        self.pressed_actions: set[str] = set()
        self.live_xtip = 0.0
        self.live_rotation = 0.0
        self.last_command_write = 0.0
        self.stl_path = stl_path
        self.wire_path = wire_path
        self.renderer = XRayRenderer(load_stl(stl_path))
        self.view_mode_var = tk.StringVar(value=self.renderer.view_mode)
        self.frame_rate_fps = VIEWER_CONFIG.frame_rate_fps
        self.display_delay_ms = VIEWER_CONFIG.display_delay_ms
        self.config_path = Path(os.environ.get("GUIDEWIRE_VIEWER_CONFIG", self.base_dir / "viewer_config.json"))
        self.settings_window = None
        self.times = np.asarray([0.0])
        self.frames: list[np.ndarray] = []
        self.frame_index = 0
        self.playing = False
        self.live_mtime = -1.0
        self.live_state_buffer: deque[tuple[float, float, np.ndarray]] = deque()
        self.last_tick = time.monotonic()
        self.poll_after_id = None
        self.render_after_id = None

        self.main = tk.Frame(root, bg="#0b0f11")
        self.main.pack(fill="both", expand=True)
        self.image_label = tk.Label(self.main, bg="#0a0e10", bd=0, highlightthickness=0)
        self.image_label.pack(fill="both", expand=True)
        self.reaction_var = tk.BooleanVar(value=VIEWER_CONFIG.show_reaction_force)
        self.contacts_var = tk.BooleanVar(value=VIEWER_CONFIG.show_contact_points)
        self.force_bar_max_n = VIEWER_CONFIG.force_bar_max_n
        # 投影への重ね描きから、OSテーマを使う下部の計測欄へ移す。
        self.force_panel = ttk.Frame(root, padding=(10, 5), relief="sunken", borderwidth=1)
        ttk.Label(self.force_panel, text="入口支持力").pack(side="left", padx=(0, 12))
        self.force_bar = ttk.Progressbar(self.force_panel, orient="horizontal",
                                         mode="determinate", length=180, maximum=1.)
        self.force_bar.pack(side="left", padx=(0, 12))
        self.force_value_label = ttk.Label(self.force_panel, text="データ待ち", width=28, anchor="w")
        self.force_value_label.pack(side="left")
        self.force_range_label = ttk.Label(self.force_panel, anchor="e")
        self.force_range_label.pack(side="right")
        self._update_force_bar()

        # Keep callback targets alive but unmanaged; their controls live only in Settings.
        self.sim_button = tk.Button(root, text="Start simulation")
        self.sim_status = tk.Label(root, text="Stopped")
        self.play_button = tk.Button(root, text="Play")
        self.slider = ttk.Scale(root, from_=0, to=0, orient="horizontal", command=self.on_slider)
        self.time_label = tk.Label(root, text="No wire data")
        self.nodes_var = tk.BooleanVar(value=False)
        self.roadmap_var = tk.BooleanVar(value=True)
        self.file_label = None

        self.live_var = tk.BooleanVar(value=False)
        self.create_menu_bar()
        self.root.bind("<space>", lambda _e: self.toggle_play())
        self.root.bind("j", lambda _e: self.step_frame(-1))
        self.root.bind("l", lambda _e: self.step_frame(1))
        for key, action in (("w", "insert"), ("s", "retract"), ("a", "rotate_left"), ("d", "rotate_right"),
                            ("Up", "insert"), ("Down", "retract"), ("Left", "rotate_left"), ("Right", "rotate_right")):
            self.root.bind(f"<KeyPress-{key}>", lambda _e, a=action: self.press_action(a))
            self.root.bind(f"<KeyRelease-{key}>", lambda _e, a=action: self.release_action(a))
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.image_label.bind("<Configure>", lambda _e: self.redraw())

        if wire_path and wire_path.exists():
            self.load_wire(wire_path)
        self.update_file_label()
        print(
            f"Frame rate: {self.frame_rate_fps:g} fps | "
            f"Display delay: {self.display_delay_ms:g} ms", flush=True)
        self.poll_after_id = self.root.after(10, self.tick)
        self.render_after_id = self.root.after(self._tick_interval_ms(), self.render_tick)
        self.root.after(350, self.start_simulation)

    def create_menu_bar(self) -> None:
        # rootのmenuに接続すると、Windows標準のウィンドウメニューになる。
        # 投影画像内へボタンや疑似ツールバーを重ねない。
        tk = self.tk
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff=False)
        self.file_menu.add_command(label="血管STLを開く…", command=self.open_stl)
        self.file_menu.add_command(label="ワイヤーログを開く…", command=self.open_wire)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="投影画像を保存…", command=self.save_image)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="終了", command=self.on_close)
        self.menu_bar.add_cascade(label="ファイル", menu=self.file_menu)

        self.view_menu = tk.Menu(self.menu_bar, tearoff=False)
        for label, mode in (("全体表示", "overview"), ("先端追従", "tip_tracking")):
            self.view_menu.add_radiobutton(label=label, variable=self.view_mode_var,
                                           value=mode, command=self.on_view_mode_change)
        self.view_menu.add_separator()
        self.view_menu.add_checkbutton(label="血管ロードマップ", variable=self.roadmap_var,
                                       command=self.redraw)
        self.view_menu.add_checkbutton(label="入口支持力", variable=self.reaction_var,
                                       command=self.redraw)
        self.view_menu.add_checkbutton(label="接触点", variable=self.contacts_var,
                                       command=self.redraw)
        self.view_menu.add_command(label="投影方向を反転", command=self.flip_view)
        self.menu_bar.add_cascade(label="表示", menu=self.view_menu)

        self.simulation_menu = tk.Menu(self.menu_bar, tearoff=False)
        self.simulation_menu.add_command(label="開始", command=self.start_simulation)
        self.simulation_menu.add_command(label="停止", command=self.stop_simulation)
        self.simulation_menu.add_separator()
        self.simulation_menu.add_command(label="ログを再生／一時停止", command=self.toggle_play,
                                         accelerator="Space")
        self.simulation_menu.add_command(label="ログの前フレーム", command=lambda: self.step_frame(-1),
                                         accelerator="J")
        self.simulation_menu.add_command(label="ログの次フレーム", command=lambda: self.step_frame(1),
                                         accelerator="L")
        self.menu_bar.add_cascade(label="シミュレーション", menu=self.simulation_menu)

        self.settings_menu = tk.Menu(self.menu_bar, tearoff=False)
        self.settings_menu.add_command(label="Viewer設定…", command=self.open_settings)
        self.menu_bar.add_cascade(label="設定", menu=self.settings_menu)
        self.help_menu = tk.Menu(self.menu_bar, tearoff=False)
        self.help_menu.add_command(label="キーボード操作", command=self.show_keyboard_help)
        self.menu_bar.add_cascade(label="ヘルプ", menu=self.help_menu)
        for menu in (self.file_menu, self.view_menu, self.simulation_menu,
                     self.settings_menu, self.help_menu):
            menu.configure(postcommand=self.release_menu_input)
        self.root.configure(menu=self.menu_bar)

    def release_menu_input(self) -> None:
        # メニューにフォーカスが移るとKeyReleaseが届かない場合がある。
        # 入力を保持したまま設定操作中に挿入が続くことを防ぐ。
        if self.pressed_actions:
            self.pressed_actions.clear()
            self.write_command()

    def show_keyboard_help(self) -> None:
        self.messagebox.showinfo(
            "キーボード操作",
            "投影画面をクリックしてから操作してください。\n\n"
            "W ／ ↑：挿入\nS ／ ↓：引き戻し\n"
            "A ／ ←：左回転\nD ／ →：右回転\n\n"
            "挿入・回転はキーを押している間だけ実行します。\n"
            "ログ再生時：Spaceで再生／一時停止、J・Lで前後のフレーム。",
            parent=self.root,
        )

    def open_settings(self) -> None:
        if self.settings_window is not None and self.settings_window.winfo_exists():
            self.settings_window.lift()
            return

        tk = self.tk
        window = tk.Toplevel(self.root)
        self.settings_window = window
        window.title("Viewer設定")
        window.transient(self.root)
        window.resizable(False, False)
        window.configure(bg=PANEL_BG)
        body = tk.Frame(window, bg=PANEL_BG, padx=18, pady=16)
        body.pack(fill="both", expand=True)

        def section(title, row):
            tk.Label(body, text=title, bg=PANEL_BG, fg="#aab5b9",
                     font=("Segoe UI Semibold", 9)).grid(
                         row=row, column=0, columnspan=3, sticky="w", pady=(10, 4))

        section("表示設定", 0)
        mode_var = tk.StringVar(value=self.renderer.view_mode)
        choices = tk.Frame(body, bg=PANEL_BG)
        choices.grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 3))
        for label, mode in (("全体表示", "overview"), ("先端追従", "tip_tracking")):
            tk.Radiobutton(
                choices, text=label, variable=mode_var, value=mode,
                bg=PANEL_BG, fg=TEXT, selectcolor="#20282c",
                activebackground=PANEL_BG, activeforeground="white",
                font=("Segoe UI", 9),
            ).pack(side="left", padx=(0, 12))

        entries = {}
        settings_rows = (
            ("Frame rate", f"{self.frame_rate_fps:g}", "fps"),
            ("Display delay", f"{self.display_delay_ms:g}", "ms"),
            ("Overview FOV", f"{VIEWER_CONFIG.overview_fov_object_mm:g}", "mm"),
            ("Tip tracking FOV", f"{VIEWER_CONFIG.tip_tracking_fov_object_mm:g}", "mm"),
            ("Wire attenuation coefficient", f"{WIRE_DELTA_MU_PER_MM:g}", "mm⁻¹"),
            ("Force bar max", f"{self.force_bar_max_n:g}", "N*"),
        )
        japanese_labels = {
            "Frame rate": "描画フレームレート", "Display delay": "表示遅延",
            "Overview FOV": "全体表示の視野幅", "Tip tracking FOV": "先端追従の視野幅",
            "Wire attenuation coefficient": "ワイヤーの減弱係数", "Force bar max": "支持力の表示レンジ上限",
        }
        for row, (label, value, unit) in enumerate(settings_rows, start=2):
            tk.Label(body, text=japanese_labels[label], bg=PANEL_BG, fg=TEXT,
                     font=("Segoe UI", 9)).grid(row=row, column=0, sticky="w", pady=5)
            entry = tk.Entry(body, width=11, bg="#20282c", fg=TEXT,
                             insertbackground=TEXT, relief="flat", justify="right")
            entry.insert(0, value)
            entry.grid(row=row, column=1, sticky="e", padx=(12, 6), ipady=4)
            tk.Label(body, text=unit, bg=PANEL_BG, fg=MUTED,
                     font=("Segoe UI", 9), width=4, anchor="w").grid(row=row, column=2, sticky="w")
            entries[label] = entry

        section("力学設定", 8)
        damping = self._current_damping()
        tk.Label(body, text="剛性減衰 β", bg=PANEL_BG, fg=TEXT,
                 font=("Segoe UI", 9)).grid(row=9, column=0, sticky="w", pady=5)
        damping_entry = tk.Entry(body, width=11, bg="#20282c", fg=TEXT,
                                 insertbackground=TEXT, relief="flat", justify="right")
        damping_entry.insert(0, f"{damping:g}" if damping is not None else "—")
        damping_entry.grid(row=9, column=1, sticky="e", padx=(12, 6), ipady=4)
        tk.Label(body, text="s", bg=PANEL_BG, fg=MUTED,
                 font=("Segoe UI", 9)).grid(row=9, column=2, sticky="w")
        if damping is not None:
            entries["Rayleigh stiffness"] = damping_entry
            damping_note = "減衰を変更して適用すると、挿入状態をリセットして計算を再開始します。"
        else:
            damping_entry.configure(state="disabled")
            damping_note = "減衰の変更には、材料設定付きの起動バットを使用してください。"
        tk.Label(body, text=damping_note, bg=PANEL_BG, fg=MUTED,
                 wraplength=330, justify="left", font=("Segoe UI", 8)).grid(
                     row=10, column=0, columnspan=3, sticky="w", pady=(2, 4))

        section("データ", 12)
        data_buttons = tk.Frame(body, bg=PANEL_BG)
        data_buttons.grid(row=13, column=0, columnspan=3, sticky="ew")
        self._settings_button(data_buttons, "血管STLを開く…", self.open_stl).pack(fill="x", pady=(0, 5))
        self._settings_button(data_buttons, "ワイヤーログを開く…", self.open_wire).pack(fill="x")
        self.file_label = tk.Label(body, text="", bg=PANEL_BG, fg=MUTED,
                                   justify="left", anchor="w", wraplength=300, font=("Segoe UI", 8))
        self.file_label.grid(row=14, column=0, columnspan=3, sticky="ew", pady=(6, 0))
        self.update_file_label()

        section("表示項目", 15)
        tk.Checkbutton(
            body, text="BeamAdapterのサンプル点", variable=self.nodes_var,
            command=self.redraw, bg=PANEL_BG, fg=TEXT, selectcolor="#20282c",
            activebackground=PANEL_BG, activeforeground="white", font=("Segoe UI", 9),
        ).grid(row=16, column=0, columnspan=3, sticky="w")
        for row, label, variable in (
                (17, "入口支持力を表示", self.reaction_var),
                (18, "接触点を表示", self.contacts_var)):
            tk.Checkbutton(
                body, text=label, variable=variable, command=self.redraw,
                bg=PANEL_BG, fg=TEXT, selectcolor="#20282c",
                activebackground=PANEL_BG, activeforeground="white", font=("Segoe UI", 9),
            ).grid(row=row, column=0, columnspan=3, sticky="w")
        tk.Label(body, text="* 単位は起動プロファイルに依存します。実測校正は未実施です。",
                 bg=PANEL_BG, fg=MUTED, font=("Segoe UI", 8)).grid(
                     row=19, column=0, columnspan=3, sticky="w")

        section("出力", 20)
        self._settings_button(body, "投影画像を保存…", self.save_image).grid(
            row=21, column=0, columnspan=3, sticky="ew")

        buttons = tk.Frame(body, bg=PANEL_BG)
        buttons.grid(row=22, column=0, columnspan=3, sticky="e", pady=(16, 0))
        tk.Button(
            buttons, text="適用", command=lambda: self.apply_settings(
                mode_var.get(), entries, window), bg="#222a2e", fg=TEXT,
            activebackground="#303a3f", activeforeground="white", relief="flat",
            bd=0, padx=14, pady=6, font=("Segoe UI", 9),
        ).pack(side="left", padx=(0, 7))
        tk.Button(
            buttons, text="閉じる", command=window.destroy, bg="#222a2e", fg=TEXT,
            activebackground="#303a3f", activeforeground="white", relief="flat",
            bd=0, padx=14, pady=6, font=("Segoe UI", 9),
        ).pack(side="left")
        window.protocol("WM_DELETE_WINDOW", window.destroy)

    def _settings_button(self, parent, text: str, command):
        return self.tk.Button(
            parent, text=text, command=command, bg="#222a2e", fg=TEXT,
            activebackground="#303a3f", activeforeground="white",
            relief="flat", bd=0, anchor="w", padx=10, pady=6,
            font=("Segoe UI", 9),
        )

    def _current_damping(self):
        path = os.environ.get("GUIDEWIRE_PARAMETERS")
        if not path:
            return None
        try:
            profile = json.loads(Path(path).read_text(encoding="utf-8"))
            value = float(profile["mechanics_code"]["rayleigh_stiffness"])
            return value if math.isfinite(value) and value >= 0 else None
        except (OSError, KeyError, TypeError, ValueError):
            return None

    def _save_damping(self, value: float) -> bool:
        if not math.isfinite(value) or value < 0:
            raise ValueError("減衰は0以上の有限値を指定してください。")
        path = Path(os.environ["GUIDEWIRE_PARAMETERS"])
        profile = json.loads(path.read_text(encoding="utf-8"))
        current = float(profile["mechanics_code"]["rayleigh_stiffness"])
        if value == current:
            return False
        profile["mechanics_code"]["rayleigh_stiffness"] = value
        if profile.get("name", "").startswith("低減衰試用（β="):
            profile["name"] = f"低減衰試用（β={value:g}）：" + profile["name"].split("：", 1)[-1]
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(profile, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, path)
        return True

    def apply_settings(self, mode: str, entries: dict, window) -> None:
        global VIEWER_CONFIG, XRAY_MATRIX, DETECTOR_WIDTH_MM, DETECTOR_HEIGHT_MM
        global WIRE_DELTA_MU_PER_MM
        try:
            frame_rate = float(entries["Frame rate"].get())
            display_delay_ms = float(entries["Display delay"].get())
            overview_fov = float(entries["Overview FOV"].get())
            tracking_fov = float(entries["Tip tracking FOV"].get())
            wire_attenuation = float(entries["Wire attenuation coefficient"].get())
            force_bar_max = float(entries["Force bar max"].get())
            damping = (float(entries["Rayleigh stiffness"].get())
                       if "Rayleigh stiffness" in entries else None)
            if damping is not None and (not math.isfinite(damping) or damping < 0):
                raise ValueError
            if (mode not in {"overview", "tip_tracking"}
                    or not all(math.isfinite(value) and value > 0
                               for value in (frame_rate, overview_fov, tracking_fov,
                                             wire_attenuation, force_bar_max))
                    or not math.isfinite(display_delay_ms)
                    or display_delay_ms < 0):
                raise ValueError
        except (KeyError, TypeError, ValueError):
            self.messagebox.showerror(
                "設定値を確認してください",
                "フレームレート・視野幅・減弱係数・支持力上限は正の有限値、表示遅延・減衰は0以上の有限値を指定してください。",
                parent=window,
            )
            return

        damping_changed = False
        if damping is not None:
            try:
                damping_changed = self._save_damping(damping)
            except (OSError, KeyError, TypeError, ValueError) as error:
                self.messagebox.showerror("減衰設定の保存", str(error), parent=window)
                return
        restart_simulation = (damping_changed and self.sim_process is not None
                              and self.sim_process.poll() is None)

        detector_fov = overview_fov * SOURCE_DETECTOR_MM / SOURCE_OBJECT_MM
        matrix = max(1, int(round(detector_fov / VIEWER_CONFIG.detector_pixel_pitch_mm)))
        VIEWER_CONFIG = replace(
            VIEWER_CONFIG,
            view_mode=mode,
            overview_fov_object_mm=overview_fov,
            tip_tracking_fov_object_mm=tracking_fov,
            frame_rate_fps=frame_rate,
            display_delay_ms=display_delay_ms,
            wire_attenuation_coefficient_per_mm=wire_attenuation,
            show_reaction_force=self.reaction_var.get(),
            show_contact_points=self.contacts_var.get(),
            force_bar_max_n=force_bar_max,
            detector_width_mm=detector_fov,
            detector_height_mm=detector_fov,
            matrix=matrix,
        )
        XRAY_MATRIX = matrix
        DETECTOR_WIDTH_MM = detector_fov
        DETECTOR_HEIGHT_MM = detector_fov
        self.frame_rate_fps = frame_rate
        self.display_delay_ms = display_delay_ms
        WIRE_DELTA_MU_PER_MM = wire_attenuation
        self.force_bar_max_n = force_bar_max
        self.last_tick = time.monotonic()
        self.renderer.set_view_mode(mode)
        self.view_mode_var.set(mode)
        print(
            f"Frame rate: {self.frame_rate_fps:g} fps | "
            f"Display delay: {self.display_delay_ms:g} ms", flush=True)
        if not self._save_settings_config(mode, overview_fov, tracking_fov,
                                          frame_rate, display_delay_ms,
                                          wire_attenuation):
            self.messagebox.showwarning(
                "設定の保存", "設定は適用されましたが、viewer_config.jsonへの保存に失敗しました。",
                parent=window,
            )
        self.redraw()
        if self.render_after_id is not None:
            try:
                self.root.after_cancel(self.render_after_id)
            except self.tk.TclError:
                pass
        self.render_after_id = self.root.after(0, self.render_tick)
        if restart_simulation:
            self.stop_simulation()
            self.start_simulation()

    def _save_settings_config(self, mode: str, overview_fov: float,
                              tracking_fov: float, frame_rate: float,
                              display_delay_ms: float,
                              wire_attenuation: float) -> bool:
        try:
            if self.config_path.exists():
                values = json.loads(self.config_path.read_text(encoding="utf-8"))
                if not isinstance(values, dict):
                    return False
            else:
                values = {}
            view = values.setdefault("view", {})
            ui = values.setdefault("ui", {})
            if not isinstance(view, dict) or not isinstance(ui, dict):
                return False
            view.update({
                "mode": mode,
                "overview_fov_object_mm": overview_fov,
                "tip_tracking_fov_object_mm": tracking_fov,
            })
            ui["frame_rate_fps"] = frame_rate
            ui["display_delay_ms"] = display_delay_ms
            ui["show_reaction_force"] = self.reaction_var.get()
            ui["show_contact_points"] = self.contacts_var.get()
            ui["force_bar_max_n"] = self.force_bar_max_n
            wire = values.setdefault("wire", {})
            if not isinstance(wire, dict):
                return False
            wire["attenuation_coefficient_per_mm"] = wire_attenuation
            xray = values.setdefault("xray", {})
            if not isinstance(xray, dict):
                return False
            xray["detector_pixel_pitch_mm"] = VIEWER_CONFIG.detector_pixel_pitch_mm
            temporary = self.config_path.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(values, indent=2) + "\n", encoding="utf-8")
            os.replace(temporary, self.config_path)
            return True
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            return False

    def _tick_interval_ms(self) -> int:
        return max(1, round(1000.0 / self.frame_rate_fps))

    def _hold_button(self, parent, text: str, action: str):
        button = self.tk.Button(
            parent, text=text, bg="#222a2e", fg=TEXT,
            activebackground="#435158", activeforeground="white",
            relief="flat", bd=0, padx=5, pady=8, font=("Segoe UI Semibold", 8),
        )
        button.bind("<ButtonPress-1>", lambda _e: self.press_action(action))
        button.bind("<ButtonRelease-1>", lambda _e: self.release_action(action))
        button.bind("<Leave>", lambda _e: self.release_action(action))
        return button

    def find_runsofa(self) -> Path | None:
        candidates: list[Path] = []
        env_root = os.environ.get("SOFA_ROOT")
        if env_root:
            candidates.append(Path(env_root) / "bin" / "runSofa.exe")
        for parent in self.base_dir.parents:
            candidates.extend(parent.glob("v*/bin/runSofa.exe"))
        return next((path for path in candidates if path.exists()), None)

    def write_command(self) -> None:
        insert = int("insert" in self.pressed_actions) - int("retract" in self.pressed_actions)
        rotate = int("rotate_right" in self.pressed_actions) - int("rotate_left" in self.pressed_actions)
        payload = {"insert": insert, "rotate": rotate, "updated": time.time()}
        temporary = self.command_path.with_suffix(".tmp")
        try:
            temporary.write_text(json.dumps(payload), encoding="utf-8")
            os.replace(temporary, self.command_path)
            self.last_command_write = time.monotonic()
        except OSError:
            pass

    def press_action(self, action: str) -> None:
        if action not in self.pressed_actions:
            self.pressed_actions.add(action)
            self.write_command()

    def release_action(self, action: str) -> None:
        if action in self.pressed_actions:
            self.pressed_actions.discard(action)
            self.write_command()

    def toggle_simulation(self) -> None:
        if self.sim_process and self.sim_process.poll() is None:
            self.stop_simulation()
        else:
            self.start_simulation()

    def start_simulation(self) -> None:
        if self.sim_process and self.sim_process.poll() is None:
            return
        if not getattr(sys, "frozen", False) and not self.engine_path.exists():
            self.sim_status.configure(text="Native engine file is missing", fg="#e7b38f")
            return
        self.pressed_actions.clear()
        self.write_command()
        # Never leave a recorded/sample wire on screen while waiting for the
        # first live SOFA frame.  At xtip=0 BeamAdapter has no visible segment,
        # so the correct display is an empty vessel rather than stale data.
        self.times = np.asarray([0.0])
        self.frames = []
        self.frame_index = 0
        self.live_state_buffer.clear()
        self.latest_force_state = None
        self.force_state_buffer.clear()
        self.force_state_mtime = -1
        self.force_session_started_at = time.time()
        self.slider.configure(to=0)
        self.slider.set(0)
        self.redraw()
        log_path = self.runtime_dir / "sofa.log"
        self.sofa_log_stream = log_path.open("w", encoding="utf-8")
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        sofa_env = os.environ.copy()
        sofa_env["GUIDEWIRE_VIEWER_SNAPSHOT"] = str(self.snapshot_path)
        sofa_env["GUIDEWIRE_VIEWER_FORCE_STATE"] = str(self.force_state_path)
        sofa_env["GUIDEWIRE_VESSEL_STL"] = str(self.stl_path)
        self.sim_process = subprocess.Popen(
            ([sys.executable, "--worker"] if getattr(sys, "frozen", False)
             else [sys.executable, str(self.engine_path)]),
            cwd=str(self.base_dir), stdout=self.sofa_log_stream, stderr=subprocess.STDOUT,
            creationflags=creationflags, env=sofa_env,
        )
        self.live_var.set(True)
        self.wire_path = self.snapshot_path
        self.live_mtime = -1.0
        self.sim_button.configure(text="Stop simulation")
        self.sim_status.configure(text="Starting extracted core… waiting for its first wire frame", fg="#9fc9b0")
        self.update_file_label()

    def stop_simulation(self) -> None:
        self.pressed_actions.clear()
        self.write_command()
        if self.sim_process and self.sim_process.poll() is None:
            self.sim_process.terminate()
            try:
                self.sim_process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                self.sim_process.kill()
        self.sim_process = None
        self.latest_force_state = None
        self.force_state_buffer.clear()
        self._update_force_bar()
        if self.sofa_log_stream:
            self.sofa_log_stream.close()
            self.sofa_log_stream = None
        self.sim_button.configure(text="Start simulation")
        self.sim_status.configure(text="Stopped", fg=MUTED)

    def on_close(self) -> None:
        self.stop_simulation()
        self.root.destroy()

    def update_file_label(self) -> None:
        if self.file_label is None:
            return
        try:
            if not self.file_label.winfo_exists():
                self.file_label = None
                return
        except self.tk.TclError:
            self.file_label = None
            return
        vessel = self.stl_path.name
        wire = self.wire_path.name if self.wire_path else "未読込み"
        self.file_label.configure(text=f"血管：{vessel}\nワイヤー：{wire}")

    def load_wire(self, path: Path) -> None:
        self.times, self.frames = load_wire_log(path)
        self.wire_path = path
        self.frame_index = len(self.frames) - 1
        self.slider.configure(to=max(0, len(self.frames) - 1))
        self.slider.set(self.frame_index)
        self.update_file_label()
        self.redraw()

    def open_stl(self) -> None:
        name = self.filedialog.askopenfilename(title="血管STLを開く", filetypes=[("STLメッシュ", "*.stl"), ("すべてのファイル", "*.*")])
        if not name:
            return
        self.stl_path = Path(name)
        self.renderer = XRayRenderer(load_stl(self.stl_path))
        self.update_file_label()
        self.redraw()

    def open_wire(self) -> None:
        name = self.filedialog.askopenfilename(title="ワイヤーログを開く", filetypes=[("ワイヤーログ", "*.txt"), ("すべてのファイル", "*.*")])
        if name:
            self.load_wire(Path(name))

    def set_angle(self, angle: float) -> None:
        self.renderer.projection.angle_deg = angle
        self.renderer.invalidate()
        self.redraw()

    def on_view_mode_change(self) -> None:
        _frame_time, wire = self.current_wire()
        if wire is not None:
            self.renderer.update_tracking_tip(to_view_coordinates(clean_wire(wire)))
        self.renderer.set_view_mode(self.view_mode_var.get())
        self.redraw()

    def flip_view(self) -> None:
        self.renderer.projection.flipped = not self.renderer.projection.flipped
        self.renderer.invalidate()
        self.redraw()

    def on_slider(self, value: str) -> None:
        if not self.frames or self.live_var.get():
            return
        self.frame_index = max(0, min(len(self.frames) - 1, int(float(value))))
        self.redraw()

    def toggle_play(self) -> None:
        if not self.frames or self.live_var.get():
            return
        self.playing = not self.playing
        self.play_button.configure(text="Pause" if self.playing else "Play")
        self.last_tick = time.monotonic()

    def step_frame(self, amount: int) -> None:
        if not self.frames or self.live_var.get():
            return
        self.playing = False
        self.play_button.configure(text="Play")
        self.frame_index = max(0, min(len(self.frames) - 1, self.frame_index + amount))
        self.slider.set(self.frame_index)
        self.redraw()

    def on_live_toggle(self) -> None:
        self.playing = False
        self.play_button.configure(text="Play")
        self.slider.configure(state="disabled" if self.live_var.get() else "normal")
        self.live_mtime = -1.0
        self.redraw()

    def current_wire(self) -> tuple[float, np.ndarray | None]:
        if not self.frames:
            return 0.0, None
        index = max(0, min(self.frame_index, len(self.frames) - 1))
        return float(self.times[index]), self.frames[index]

    def redraw(self) -> None:
        width = max(500, self.image_label.winfo_width())
        height = max(500, self.image_label.winfo_height())
        frame_time, wire = self.current_wire()
        image = self.renderer.render(
            width,
            height,
            wire,
            frame_time=frame_time,
            frame_index=self.frame_index,
            frame_count=len(self.frames),
            show_nodes=self.nodes_var.get(),
            show_roadmap=self.roadmap_var.get(),
        )
        # Diagnostics are composited ONLY on the GUI copy, never on the X-ray
        # renderer output, roadmap cache, attenuation calculation or saved image.
        if self.contacts_var.get() and self.live_var.get():
            state = self._force_state_for_frame(frame_time)
            contacts = state.get("contacts", []) if state else []
            if contacts:
                positions = np.asarray([contact["position"] for contact in contacts], dtype=float)
                try:
                    pixels, bounds = self.renderer.project_diagnostic_points(positions)
                    image = image.copy()
                    draw = ImageDraw.Draw(image)
                    for x, y in pixels:
                        if (np.isfinite(x) and np.isfinite(y)
                                and bounds[0] + 3 <= x <= bounds[2] - 3
                                and bounds[1] + 3 <= y <= bounds[3] - 3):
                            draw.ellipse((x-3, y-3, x+3, y+3),
                                         outline=(230, 35, 35), width=1)
                except ValueError:
                    pass
        self.photo = ImageTk.PhotoImage(image)
        self.image_label.configure(image=self.photo)
        self._update_force_bar()
        if self.frames:
            self.time_label.configure(text=f"t = {frame_time:8.2f} s")
        else:
            self.time_label.configure(text="No wire data")

    def save_image(self) -> None:
        name = self.filedialog.asksaveasfilename(
            title="投影画像を保存",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png")],
        )
        if not name:
            return
        frame_time, wire = self.current_wire()
        self.renderer.render(
            1400,
            1000,
            wire,
            frame_time=frame_time,
            frame_index=self.frame_index,
            frame_count=len(self.frames),
            show_nodes=self.nodes_var.get(),
            show_roadmap=self.roadmap_var.get(),
        ).save(name)

    def tick(self) -> None:
        now = time.monotonic()
        self._poll_force_state(now)
        if self.sim_process:
            return_code = self.sim_process.poll()
            if return_code is not None:
                self.sim_status.configure(text=f"Simulation stopped (code {return_code}).\nSee runtime/sofa.log", fg="#e7b38f")
                self.sim_button.configure(text="Start simulation")
                self.sim_process = None
                if self.sofa_log_stream:
                    self.sofa_log_stream.close()
                    self.sofa_log_stream = None
            elif self.snapshot_path.exists():
                self.sim_status.configure(
                    text=f"Running   insertion {self.live_xtip * float(os.environ.get('GUIDEWIRE_INSERTION_DISPLAY_SCALE', '1')):6.1f} mm\nrotation {math.degrees(self.live_rotation):7.1f}°",
                    fg="#9fc9b0",
                )
            if self.sim_process and now - self.last_command_write >= 0.1:
                self.write_command()
        if self.live_var.get() and self.wire_path and self.wire_path.exists():
            mtime = self.wire_path.stat().st_mtime
            if mtime != self.live_mtime:
                if self.wire_path == self.snapshot_path:
                    snapshot = read_live_snapshot(self.wire_path)
                    latest = None if snapshot is None else (snapshot[0], snapshot[3])
                    if snapshot is not None:
                        self.live_xtip, self.live_rotation = snapshot[1], snapshot[2]
                else:
                    latest = read_latest_wire_frame(self.wire_path)
                self.live_mtime = mtime
                if latest is not None:
                    frame_time, frame = latest
                    self.live_state_buffer.append((now, frame_time, frame.copy()))
        self.poll_after_id = self.root.after(10, self.tick)

    def _poll_force_state(self, now: float) -> None:
        if self.live_var.get():
            try:
                mtime = self.force_state_path.stat().st_mtime_ns
                if mtime != self.force_state_mtime:
                    state = json.loads(self.force_state_path.read_text(encoding="utf-8"))
                    if (state.get("session") == str(self.snapshot_path)
                            and float(state["updated"]) >= self.force_session_started_at
                            and 0 <= time.time() - float(state["updated"]) < 1.0):
                        self.latest_force_state = state
                        self.last_force_received = now
                        self.force_state_buffer.append(state)
                        cutoff = float(state["time"]) - max(2.0, self.display_delay_ms / 1000.0 + 1.0)
                        while len(self.force_state_buffer) > 2 and self.force_state_buffer[1]["time"] < cutoff:
                            self.force_state_buffer.popleft()
                    self.force_state_mtime = mtime
            except (OSError, ValueError, KeyError, TypeError):
                pass
        if not self.live_var.get() or now - self.last_force_received > 1.0:
            self.latest_force_state = None
            self.force_state_buffer.clear()
        self._update_force_bar()

    def _force_state_for_frame(self, frame_time: float):
        for state in reversed(self.force_state_buffer):
            if state["time"] <= frame_time + 1e-7:
                return state
        return None

    def _update_force_bar(self) -> None:
        if not self.reaction_var.get():
            self.force_panel.pack_forget()
            return
        if not self.force_panel.winfo_manager():
            self.force_panel.pack(side="bottom", fill="x", before=self.main)
        state = self.latest_force_state
        support = state.get("insertion_support_force") if state else None
        magnitude = None
        if support:
            vector = np.asarray(support.get("force", []), dtype=float)
            if vector.shape == (3,) and np.all(np.isfinite(vector)):
                magnitude = float(np.linalg.norm(vector))
        fraction = min(1.0, magnitude / self.force_bar_max_n) if magnitude is not None else 0.0
        self.force_bar.configure(value=fraction)
        unit = "N" if state and str(state.get("force_unit", "")).startswith("N (SI") else "シーン単位"
        self.force_value_label.configure(
            text=f"{magnitude:.3g} {unit}" if magnitude is not None else "データ待ち")
        timing = "直近の計算速度：データ待ち"
        if state and "wall_time" in state:
            wall_time = float(state["wall_time"])
            simulation_time = float(state["time"])
            if np.isfinite(wall_time) and np.isfinite(simulation_time):
                previous = None
                # 実時間で約1秒前の状態と比較し、累積遅れを含めない。
                for sample in reversed(self.force_state_buffer):
                    if "wall_time" not in sample:
                        continue
                    sample_wall = float(sample["wall_time"])
                    sample_time = float(sample["time"])
                    if not np.isfinite(sample_wall) or not np.isfinite(sample_time):
                        continue
                    if sample_wall < wall_time and sample_time <= simulation_time:
                        previous = sample
                        if wall_time - sample_wall >= 1.0:
                            break
                if previous is not None:
                    elapsed = wall_time - float(previous["wall_time"])
                    if elapsed >= 0.25:
                        speed = (simulation_time - float(previous["time"])) / elapsed
                        pace_ms = (1.0 - speed) * 1000.0
                        direction = "遅れ" if pace_ms >= 0.0 else "先行"
                        timing = (f"直近{elapsed:.2f}秒：計算速度 {speed * 100:.0f}％"
                                  f" / {direction} {abs(pace_ms):.0f} ms/秒")
        self.force_range_label.configure(text=timing)

    def render_tick(self) -> None:
        now = time.monotonic()
        if self.live_var.get():
            selected = self._delayed_live_state(now)
            if selected is not None:
                frame_time, frame = selected
                self.times = np.asarray([frame_time])
                self.frames = [frame]
                self.frame_index = 0
                self.redraw()
        elif self.playing and self.frames:
            elapsed = now - self.last_tick
            if elapsed >= 1 / self.frame_rate_fps:
                increment = max(1, int(elapsed * self.frame_rate_fps))
                self.frame_index = (self.frame_index + increment) % len(self.frames)
                self.slider.set(self.frame_index)
                self.last_tick = now
                self.redraw()
        self.render_after_id = self.root.after(self._tick_interval_ms(), self.render_tick)

    def _delayed_live_state(self, now: float) -> tuple[float, np.ndarray] | None:
        if not self.live_state_buffer:
            return None
        target_time = now - self.display_delay_ms / 1000.0
        if target_time < self.live_state_buffer[0][0]:
            return None
        selected = self.live_state_buffer[0]
        for state in self.live_state_buffer:
            if state[0] > target_time:
                break
            selected = state

        # Keep the selected predecessor plus enough recent history for the delay.
        cutoff = target_time - 1.0
        while len(self.live_state_buffer) > 2 and self.live_state_buffer[1][0] < cutoff:
            self.live_state_buffer.popleft()
        return selected[1], selected[2]


def make_sample(source: Path, destination: Path, count: int) -> None:
    times, frames = load_wire_log(source, max_frames=max(count, 2))
    if len(frames) > count:
        indices = np.linspace(0, len(frames) - 1, count).round().astype(int)
        times = times[indices]
        frames = [frames[i] for i in indices]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        for frame_time, frame in zip(times, frames):
            stream.write(f"T= {frame_time:.8g}\n")
            values = " ".join(f"{value:.9g}" for value in frame.reshape(-1))
            stream.write(f"  X= {values}\n")


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Native guidewire projection viewer")
    parser.add_argument("--stl", type=Path, default=here / "assets" / "vessel.stl")
    parser.add_argument("--wire", type=Path, default=None)
    parser.add_argument("--render-preview", type=Path)
    parser.add_argument("--make-sample", nargs=2, metavar=("SOURCE", "DESTINATION"))
    parser.add_argument("--sample-count", type=int, default=240)
    args = parser.parse_args()

    if args.make_sample:
        make_sample(Path(args.make_sample[0]), Path(args.make_sample[1]), args.sample_count)
        return 0

    if not args.stl.exists():
        raise FileNotFoundError(f"Vessel STL not found: {args.stl}")
    triangles = load_stl(args.stl)
    renderer = XRayRenderer(triangles)
    times = np.asarray([0.0])
    frames: list[np.ndarray] = []
    preview_wire = args.wire
    if args.render_preview and preview_wire is None:
        preview_wire = here / "assets" / "sample_wire_positions.txt"
    if preview_wire and preview_wire.exists():
        times, frames = load_wire_log(preview_wire)

    if args.render_preview:
        wire = frames[-1] if frames else None
        renderer.render(
            1080,
            760,
            wire,
            frame_time=float(times[-1]),
            frame_index=max(0, len(frames) - 1),
            frame_count=len(frames),
        ).save(args.render_preview)
        return 0

    import tkinter as tk

    root = tk.Tk()
    ViewerApp(root, args.stl, args.wire if args.wire and args.wire.exists() else None)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
