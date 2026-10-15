"""Cone-beam path-length projection used by the native guidewire viewer.

Geometry is expressed in millimetres.  This module contains no UI and does not
modify the SOFA/BeamAdapter mechanical model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ConeBeamCamera:
    source: np.ndarray
    horizontal: np.ndarray
    vertical: np.ndarray
    forward: np.ndarray
    source_object_mm: float
    source_detector_mm: float
    detector_width_mm: float
    detector_height_mm: float
    width: int
    height: int

    def project_mm(self, points: np.ndarray) -> np.ndarray:
        delta = np.asarray(points, dtype=np.float64) - self.source
        depth = delta @ self.forward
        if np.any(depth <= 0):
            raise ValueError("Projection geometry lies behind the X-ray source")
        scale = self.source_detector_mm / depth
        return np.column_stack((scale * (delta @ self.horizontal),
                                scale * (delta @ self.vertical), depth))

    def ray_directions(self) -> np.ndarray:
        xs = (np.arange(self.width) + 0.5 - self.width / 2) * (
            self.detector_width_mm / self.width)
        ys = (self.height / 2 - np.arange(self.height) - 0.5) * (
            self.detector_height_mm / self.height)
        u, v = np.meshgrid(xs, ys)
        rays = (self.source_detector_mm * self.forward[None, None, :]
                + u[..., None] * self.horizontal[None, None, :]
                + v[..., None] * self.vertical[None, None, :])
        return rays / np.linalg.norm(rays, axis=2, keepdims=True)


def close_planar_boundaries(triangles: np.ndarray) -> np.ndarray:
    """Close simple planar STL openings so inside path lengths are defined."""
    triangles = np.asarray(triangles, dtype=np.float64)
    vertices: list[np.ndarray] = []
    ids: dict[tuple[float, float, float], int] = {}
    edge_uses: dict[tuple[int, int], int] = {}
    for triangle in triangles:
        face: list[int] = []
        for point in triangle:
            key = tuple(float(value) for value in point)
            if key not in ids:
                ids[key] = len(vertices)
                vertices.append(np.asarray(point, dtype=np.float64))
            face.append(ids[key])
        for a, b in zip(face, face[1:] + face[:1]):
            edge = (min(a, b), max(a, b))
            edge_uses[edge] = edge_uses.get(edge, 0) + 1
    if any(count > 2 for count in edge_uses.values()):
        raise ValueError("Vessel STL has non-manifold edges")
    boundary = {edge for edge, count in edge_uses.items() if count == 1}
    if not boundary:
        return triangles
    neighbors: dict[int, list[int]] = {}
    for a, b in boundary:
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)
    if any(len(values) != 2 for values in neighbors.values()):
        raise ValueError("Vessel STL has a broken boundary")
    caps: list[np.ndarray] = []
    while boundary:
        start, second = next(iter(boundary))
        loop = [start, second]
        boundary.remove((min(start, second), max(start, second)))
        previous, current = start, second
        while current != start:
            nxt = next(value for value in neighbors[current] if value != previous)
            if nxt == start:
                break
            edge = (min(current, nxt), max(current, nxt))
            if edge not in boundary:
                raise ValueError("Vessel STL boundary is not a simple loop")
            boundary.remove(edge)
            loop.append(nxt)
            previous, current = current, nxt
        closing = (min(current, start), max(current, start))
        if closing not in boundary:
            raise ValueError("Vessel STL boundary does not close")
        boundary.remove(closing)
        points = np.stack([vertices[index] for index in loop])
        center = points.mean(axis=0)
        _, singular, _ = np.linalg.svd(points - center)
        if singular[-1] > 1e-3 * max(singular[0], 1.0):
            raise ValueError("Vessel STL opening is not planar")
        for a, b in zip(points, np.roll(points, -1, axis=0)):
            caps.append(np.stack((a, b, center)))
    return np.concatenate((triangles, np.asarray(caps)), axis=0)


def wire_tube(points: np.ndarray, diameter_mm: float, sides: int = 6) -> np.ndarray:
    """Convert BeamAdapter VisualPoints to a closed display-only wire tube."""
    points = np.asarray(points, dtype=np.float64)
    keep = np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1) > 1e-7]
    points = points[keep]
    if len(points) < 2:
        return np.empty((0, 3, 3), dtype=np.float64)
    radius = diameter_mm * 0.5
    rings: list[np.ndarray] = []
    for index, point in enumerate(points):
        tangent = points[min(index + 1, len(points) - 1)] - points[max(index - 1, 0)]
        tangent /= np.linalg.norm(tangent)
        reference = np.array([0.0, 1.0, 0.0])
        if abs(float(reference @ tangent)) > 0.9:
            reference = np.array([1.0, 0.0, 0.0])
        first = np.cross(tangent, reference)
        first /= np.linalg.norm(first)
        second = np.cross(tangent, first)
        angles = 2 * np.pi * np.arange(sides) / sides
        rings.append(point + radius * (np.cos(angles)[:, None] * first
                                       + np.sin(angles)[:, None] * second))
    faces: list[np.ndarray] = []
    for first, second in zip(rings[:-1], rings[1:]):
        for index in range(sides):
            nxt = (index + 1) % sides
            faces.extend((np.stack((first[index], first[nxt], second[nxt])),
                          np.stack((first[index], second[nxt], second[index]))))
    for ring, center in ((rings[0], points[0]), (rings[-1], points[-1])):
        for index in range(sides):
            faces.append(np.stack((center, ring[index], ring[(index + 1) % sides])))
    return np.asarray(faces)


def path_lengths(triangles: np.ndarray, camera: ConeBeamCamera,
                 max_hits: int = 32, return_first_hit: bool = False):
    """Return material path length in mm for every finite-source detector ray."""
    height, width = camera.height, camera.width
    directions = camera.ray_directions()
    hits = np.full((max_hits, height, width), np.inf, dtype=np.float32)
    counts = np.zeros((height, width), dtype=np.uint8)
    flat_hits = hits.reshape(max_hits, -1)
    flat_counts = counts.ravel()
    dx = camera.detector_width_mm / width
    dy = camera.detector_height_mm / height
    projected = camera.project_mm(np.asarray(triangles).reshape(-1, 3))[:, :2].reshape(-1, 3, 2)
    for triangle, uv in zip(triangles, projected):
        u0, u1 = uv[:, 0].min(), uv[:, 0].max()
        v0, v1 = uv[:, 1].min(), uv[:, 1].max()
        x0 = max(0, int(math.floor(u0 / dx + width / 2 - 1)))
        x1 = min(width, int(math.ceil(u1 / dx + width / 2 + 1)))
        y0 = max(0, int(math.floor(height / 2 - v1 / dy - 1)))
        y1 = min(height, int(math.ceil(height / 2 - v0 / dy + 1)))
        if x0 >= x1 or y0 >= y1:
            continue
        local_rays = directions[y0:y1, x0:x1].reshape(-1, 3)
        a, b, c = triangle
        edge1, edge2 = b - a, c - a
        pvec = np.cross(local_rays, edge2)
        determinant = pvec @ edge1
        stable = np.abs(determinant) > 1e-9
        inverse = np.zeros_like(determinant)
        inverse[stable] = 1.0 / determinant[stable]
        tvec = camera.source - a
        bary_u = (pvec @ tvec) * inverse
        qvec = np.cross(tvec, edge1)
        bary_v = (local_rays @ qvec) * inverse
        distance = float(edge2 @ qvec) * inverse
        valid = (stable & (bary_u >= -1e-8) & (bary_v >= -1e-8)
                 & (bary_u + bary_v <= 1 + 1e-8) & (distance > 0))
        if not np.any(valid):
            continue
        local_y, local_x = np.divmod(np.flatnonzero(valid), x1 - x0)
        indices = (local_y + y0) * width + local_x + x0
        levels = flat_counts[indices].astype(np.intp)
        usable = levels < max_hits
        indices, levels = indices[usable], levels[usable]
        distances = distance[valid][usable]
        flat_hits[levels, indices] = distances
        flat_counts[indices] += 1
    hits.sort(axis=0)
    unique = np.full_like(hits, np.inf)
    unique_counts = np.zeros((height, width), dtype=np.uint8)
    flat_unique = unique.reshape(max_hits, -1)
    flat_unique_counts = unique_counts.ravel()
    last = np.full(height * width, -np.inf, dtype=np.float32)
    for layer in hits.reshape(max_hits, -1):
        distinct = np.isfinite(layer) & (layer - last > 1e-3)
        indices = np.flatnonzero(distinct)
        flat_unique[flat_unique_counts[indices], indices] = layer[indices]
        flat_unique_counts[indices] += 1
        last[indices] = layer[indices]
    lengths = np.zeros((height, width), dtype=np.float32)
    for level in range(0, max_hits - 1, 2):
        valid = unique_counts >= level + 2
        lengths[valid] += unique[level + 1, valid] - unique[level, valid]
    if return_first_hit:
        return lengths, np.where(lengths > 0, unique[0], np.inf)
    return lengths


def capsule_path_lengths(points: np.ndarray, radius_mm: float,
                         camera: ConeBeamCamera, u: np.ndarray,
                         v: np.ndarray) -> np.ndarray:
    """Exact ray path lengths through the union of centerline capsules.

    ``u`` and ``v`` are detector-plane sample coordinates in millimetres.
    Each consecutive VisualPoints pair forms a finite cylinder with spherical
    end caps; overlapping intervals are merged so joints are not double-counted.
    """
    points = np.asarray(points, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64).ravel()
    v = np.asarray(v, dtype=np.float64).ravel()
    rays = (camera.source_detector_mm * camera.forward[None, :]
            + u[:, None] * camera.horizontal[None, :]
            + v[:, None] * camera.vertical[None, :])
    rays /= np.linalg.norm(rays, axis=1, keepdims=True)
    n = len(rays)
    detector_tmax = camera.source_detector_mm / np.maximum(
        rays @ camera.forward, 1e-12)
    starts: list[np.ndarray] = []
    ends: list[np.ndarray] = []

    def add_sphere(center: np.ndarray) -> None:
        m = camera.source - center
        b = rays @ m
        c = float(m @ m) - radius_mm * radius_mm
        disc = b * b - c
        root = np.sqrt(np.maximum(disc, 0.0))
        lo, hi = -b - root, -b + root
        valid = disc >= 0.0
        starts.append(np.where(valid, np.maximum(lo, 0.0), np.inf))
        ends.append(np.where(valid, np.minimum(hi, detector_tmax), -np.inf))

    if len(points) == 1:
        add_sphere(points[0])
    else:
        for a, bpoint in zip(points[:-1], points[1:]):
            axis = bpoint - a
            length = float(np.linalg.norm(axis))
            if length <= 1e-10:
                continue
            e = axis / length
            m = camera.source - a
            de = rays @ e
            me = float(m @ e)
            dp = rays - de[:, None] * e
            mp = m - me * e
            qa = np.einsum("ij,ij->i", dp, dp)
            qb = 2.0 * (dp @ mp)
            qc = float(mp @ mp) - radius_mm * radius_mm
            disc = qb * qb - 4.0 * qa * qc
            valid = (disc >= 0.0) & (qa > 1e-14)
            root = np.sqrt(np.maximum(disc, 0.0))
            lo = (-qb - root) / np.maximum(2.0 * qa, 1e-14)
            hi = (-qb + root) / np.maximum(2.0 * qa, 1e-14)
            axial_lo = np.full(n, -np.inf)
            axial_hi = np.full(n, np.inf)
            moving = np.abs(de) > 1e-14
            t0 = -me / np.where(moving, de, 1.0)
            t1 = (length - me) / np.where(moving, de, 1.0)
            axial_lo[moving] = np.minimum(t0[moving], t1[moving])
            axial_hi[moving] = np.maximum(t0[moving], t1[moving])
            parallel_inside = (~moving) & (me >= 0.0) & (me <= length)
            valid &= moving | parallel_inside
            lo = np.maximum.reduce((lo, axial_lo, np.zeros(n)))
            hi = np.minimum(np.minimum(hi, axial_hi), detector_tmax)
            valid &= hi > lo
            starts.append(np.where(valid, lo, np.inf))
            ends.append(np.where(valid, hi, -np.inf))
            add_sphere(a)
            add_sphere(bpoint)

    if not starts:
        return np.zeros(n, dtype=np.float32)
    lo = np.stack(starts)
    hi = np.stack(ends)
    order = np.argsort(lo, axis=0)
    lo = np.take_along_axis(lo, order, axis=0)
    hi = np.take_along_axis(hi, order, axis=0)
    total = np.zeros(n, dtype=np.float64)
    current_lo = lo[0].copy()
    current_hi = hi[0].copy()
    for index in range(1, len(starts)):
        next_lo, next_hi = lo[index], hi[index]
        disjoint = next_lo > current_hi
        total += np.where(disjoint, np.maximum(current_hi - current_lo, 0.0), 0.0)
        current_lo = np.where(disjoint, next_lo, current_lo)
        current_hi = np.where(disjoint, next_hi, np.maximum(current_hi, next_hi))
    total += np.maximum(current_hi - current_lo, 0.0)
    return total.astype(np.float32)
