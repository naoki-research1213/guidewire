# Extracted from current interactive_xray_scene.py; see reference_manifest.json.
from __future__ import annotations
import os,json,time,csv
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
RUNTIME=Path(os.environ.get("GUIDEWIRE_RUNTIME_DIR",HERE.parent/"runtime/native_viewer"))
COMMAND=RUNTIME/"control.json"
SNAPSHOT=Path(os.environ.get("GUIDEWIRE_VIEWER_SNAPSHOT",RUNTIME/"live_wire.txt"))
FORCE_STATE=Path(os.environ.get("GUIDEWIRE_VIEWER_FORCE_STATE",RUNTIME/"live_force_state.json"))
class DiagnosticsBase:
    def __init__(self,*args,**kwargs):pass
    def getContext(self):return self.context
def scalar(value, default=0.0):
    try:
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            return float(value[0]) if len(value) else float(default)
        return float(value)
    except Exception:
        return float(default)


def support_force_components(position, rest_position, velocity, stiffness, beta):
    """Opposite of the translational spring + its implicit Rayleigh damping.

    EulerImplicitSolver uses -beta*K_s*v. RestShapeSprings' rest target is
    stationary during that solve; use the solver velocity, NOT an insertion
    speed or a finite difference across BeamAdapter's node redistribution.
    This is a support force, not the generalized force conjugate to xtip.
    """
    elastic = stiffness * (np.asarray(position) - np.asarray(rest_position))
    damping = beta * stiffness * np.asarray(velocity)
    return elastic, damping, elastic + damping


class ExternalDriveController(DiagnosticsBase):
    """File-based bridge between the native viewer and BeamAdapter."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.deploy = kwargs["deploy"]
        self.collision_dofs = kwargs["collision_dofs"]
        self.beam_dofs = kwargs["beam_dofs"]
        self.visual_dofs = kwargs["visual_dofs"]
        self.constraint_solver = kwargs["constraint_solver"]
        self.rest_springs = kwargs["rest_springs"]
        self.ode_solver = kwargs["ode_solver"]
        self.contact_visual = kwargs["contact_visual"]
        self.force_scale = float(kwargs.get("force_scale", 1.0))
        self.force_unit = kwargs.get("force_unit", "SOFA scene unit; nominal N display, not SI-calibrated")
        self.position_unit = kwargs.get("position_unit", "scene mm")
        self.insert_rate = float(kwargs.get("insert_rate", 10.0))
        self.rotation_rate = float(kwargs.get("rotation_rate", 0.5))
        self.insert_axis = 0
        self.rotate_axis = 0
        self.last_command_mtime = -1
        self.last_snapshot = 0.0
        self.force_warning_reported = False
        self.last_force_snapshot = 0.0
        self.contact_force_source = "unavailable"
        RUNTIME.mkdir(parents=True, exist_ok=True)
        self.contact_log = (RUNTIME / "contact_force_log.csv").open(
            "w", newline="", encoding="utf-8")
        self.contact_writer = csv.writer(self.contact_log)
        self.contact_log.write(
            f"# Units: time=s, position={self.position_unit}, force={self.force_unit}\n")
        self.contact_writer.writerow((
            "time", "contact_index", "x", "y", "z", "fx", "fy", "fz", "force_magnitude"))
        self.contact_log.flush()
        self.reaction_log = (RUNTIME / "insertion_support_force_log.csv").open(
            "w", newline="", encoding="utf-8")
        self.reaction_writer = csv.writer(self.reaction_log)
        self.reaction_log.write(
            "# Translational support spring + its solver Rayleigh damping; NOT total xtip actuator reaction.\n"
            f"# Units: time=s, position={self.position_unit}, force={self.force_unit}\n")
        self.reaction_writer.writerow((
            "time", "fx", "fy", "fz", "force_magnitude", "axial_force",
            "elastic_fx", "elastic_fy", "elastic_fz",
            "damping_fx", "damping_fy", "damping_fz",
            "support_index", "stiffness", "rayleigh_beta", "vx", "vy", "vz"))
        self.reaction_log.flush()

    def read_command(self):
        try:
            mtime = COMMAND.stat().st_mtime_ns
            if mtime != self.last_command_mtime:
                command = json.loads(COMMAND.read_text(encoding="utf-8"))
                updated = float(command.get("updated", 0.0))
                if time.time() - updated <= 0.5:
                    self.insert_axis = max(-1, min(1, int(command.get("insert", 0))))
                    self.rotate_axis = max(-1, min(1, int(command.get("rotate", 0))))
                else:
                    self.insert_axis = self.rotate_axis = 0
                self.last_command_mtime = mtime
        except Exception:
            self.insert_axis = self.rotate_axis = 0

    def onAnimateBeginEvent(self, _event):
        self.read_command()
        dt = float(self.getContext().dt.value)
        if self.insert_axis:
            xtip = scalar(self.deploy.xtip.value)
            xtip = max(0.0, min(1000.0, xtip + self.insert_axis * self.insert_rate * dt))
            self.deploy.xtip.value = [xtip]
        if self.rotate_axis:
            rotation = scalar(self.deploy.rotationInstrument.value)
            self.deploy.rotationInstrument.value = [rotation + self.rotate_axis * self.rotation_rate * dt]

    def onAnimateEndEvent(self, _event):
        simulation_time = float(self.getContext().getTime())
        self.record_force_diagnostics(simulation_time)
        now = time.monotonic()
        if now - self.last_snapshot < 1.0 / 30.0:
            return
        try:
            # This is BeamAdapter's own spline evaluation, generated by
            # AdaptiveBeamMapping.  The viewer must not invent a second spline.
            positions = self.visual_dofs.position.value
            values = []
            for point in positions:
                values.extend((float(point[0]), float(point[1]), float(point[2])))
            if len(values) < 6:
                return
            xtip = scalar(self.deploy.xtip.value)
            rotation = scalar(self.deploy.rotationInstrument.value)
            text = f"{simulation_time:.8g} {xtip:.8g} {rotation:.8g}\n"
            text += " ".join(f"{value:.9g}" for value in values) + "\n"
            temporary = SNAPSHOT.with_suffix(".tmp")
            temporary.write_text(text, encoding="utf-8")
            os.replace(temporary, SNAPSHOT)
            self.last_snapshot = now
        except Exception:
            pass

    @staticmethod
    def _flat_values(value):
        try:
            return np.asarray(value, dtype=np.float64).reshape(-1)
        except (TypeError, ValueError):
            return np.fromstring(str(value), sep=" ", dtype=np.float64)

    def _world_contact_forces(self):
        """Get world contact forces on collision DOFs, without adding constraints."""
        positions = np.asarray(self.collision_dofs.position.value, dtype=np.float64).reshape(-1, 3)
        forces = np.zeros_like(positions)
        dt = float(self.getContext().dt.value)
        if dt <= 0.0 or not np.isfinite(dt):
            return positions, forces
        # Prefer the mapped MechanicalObject reaction vector when this solver
        # actually exposes it. The current LCP solver does not allocate lambda.
        # EulerImplicit's correction RHS is an impulse (velocity factor=1,
        # position factor=dt); convert impulse to time-step-averaged force.
        direct = self.collision_dofs.getData("lambda")
        if direct is not None:
            values = np.asarray(direct.value, dtype=np.float64)
            if values.shape == positions.shape and np.all(np.isfinite(values)):
                self.contact_force_source = "CollisionDOFs.lambda / dt"
                return positions, values / dt
        # Compressed constraint data may be exposed as scipy.sparse by newer
        # bindings; its text accessor is sufficient for the official compact
        # constraint-map format and avoids requiring SciPy in the SOFA Python.
        raw_constraints = self.collision_dofs.getData("constraint").getValueString()
        constraints = self._flat_values(raw_constraints)
        impulses = self._flat_values(self.constraint_solver.constraintForces.value)
        self.contact_force_source = "CollisionDOFs.constraint^T * ConstraintSolver.constraintForces / dt"
        if not constraints.size:
            return positions, forces
        # SOFA's MatrixDeriv stream starts with the number of stored rows;
        # then each row is: global_constraint_id, count, (DOF_id, x,y,z)*count.
        # Global IDs index solver.constraintForces, not the compact row number.
        row_count = int(round(constraints[0]))
        cursor = 1
        for _row in range(row_count):
            if cursor + 2 > constraints.size:
                raise ValueError("Truncated CollisionDOFs.constraint row")
            constraint_id = int(round(constraints[cursor]))
            dof_count = int(round(constraints[cursor + 1]))
            record_end = cursor + 2 + 4 * dof_count
            if dof_count < 0 or record_end > constraints.size:
                raise ValueError("Invalid CollisionDOFs.constraint row")
            if 0 <= constraint_id < impulses.size:
                force_scalar = impulses[constraint_id] / dt
                offset = cursor + 2
                for _ in range(dof_count):
                    point_id = int(round(constraints[offset]))
                    direction = constraints[offset + 1:offset + 4]
                    if 0 <= point_id < len(forces):
                        forces[point_id] += direction * force_scalar
                    offset += 4
            cursor = record_end
        if cursor != constraints.size:
            raise ValueError("Unexpected CollisionDOFs.constraint stream format")
        return positions, forces

    @staticmethod
    def _insert_direction(pose):
        """Rotate the beam's local +X insertion axis by its SOFA quaternion."""
        qx, qy, qz, qw = (float(value) for value in pose[3:7])
        norm = (qx*qx + qy*qy + qz*qz + qw*qw) ** 0.5
        if norm <= 1e-12:
            return None
        qx, qy, qz, qw = qx/norm, qy/norm, qz/norm, qw/norm
        direction = np.asarray((
            1.0 - 2.0 * (qy*qy + qz*qz),
            2.0 * (qx*qy + qw*qz),
            2.0 * (qx*qz - qw*qy),
        ), dtype=np.float64)
        direction_norm = np.linalg.norm(direction)
        return direction / direction_norm if direction_norm > 1e-12 else None

    def record_force_diagnostics(self, simulation_time):
        try:
            positions, forces = self._world_contact_forces()
            forces = forces * self.force_scale
            active = []
            contacts = []
            for contact_index, (position, force) in enumerate(zip(positions, forces)):
                magnitude = float(np.linalg.norm(force))
                if not np.isfinite(magnitude) or magnitude <= 1e-6:
                    continue
                active.append(position.tolist())
                contacts.append({"index": contact_index, "position": position.tolist(),
                                 "force": force.tolist(), "magnitude": magnitude})
                self.contact_writer.writerow((
                    f"{simulation_time:.9g}", contact_index,
                    *(f"{float(value):.9g}" for value in position),
                    *(f"{float(value):.9g}" for value in force),
                    f"{magnitude:.9g}"))
            self.contact_log.flush()
            # A non-colliding MechanicalObject keeps this display-only: it adds no
            # collision geometry or force to the simulated system.
            self.contact_visual.position.value = active

            index = int(round(scalar(self.deploy.indexFirstNode.value, -1.0)))
            current = np.asarray(self.beam_dofs.position.value, dtype=np.float64).reshape(-1, 7)
            rest = np.asarray(self.beam_dofs.rest_position.value, dtype=np.float64).reshape(-1, 7)
            stiffness = scalar(self.rest_springs.stiffness.value, 1e8)
            velocity = np.asarray(self.beam_dofs.velocity.value, dtype=np.float64).reshape(-1, 6)
            beta = scalar(self.ode_solver.rayleighStiffness.value)
            # ForceField-local Rayleigh stiffness contributes in addition to
            # EulerImplicitSolver's global coefficient (MechanicalParams).
            local_beta = self.rest_springs.getData("rayleighStiffness")
            if local_beta is not None:
                beta += scalar(local_beta.value)
            support = None
            if 0 <= index < min(len(current), len(rest), len(velocity)):
                # IRC updates indexFirstNode and fixes preceding virtual nodes.
                # Measure this actual spring boundary, not a hard-coded node 0.
                elastic, damping, reaction = support_force_components(
                    current[index, :3], rest[index, :3], velocity[index, :3], stiffness, beta)
                elastic, damping, reaction = (value * self.force_scale for value in (elastic, damping, reaction))
                direction = self._insert_direction(current[index])
                magnitude = float(np.linalg.norm(reaction))
                axial = float(np.dot(reaction, direction)) if direction is not None else float("nan")
                self.reaction_writer.writerow((
                    f"{simulation_time:.9g}", *(f"{float(v):.9g}" for v in reaction),
                    f"{magnitude:.9g}", f"{axial:.9g}",
                    *(f"{float(v):.9g}" for v in elastic),
                    *(f"{float(v):.9g}" for v in damping),
                    index, f"{stiffness:.9g}", f"{beta:.9g}",
                    *(f"{float(v):.9g}" for v in velocity[index, :3])))
                if np.all(np.isfinite(reaction)):
                    support = {"index": index, "position": current[index, :3].tolist(),
                               "force": reaction.tolist(), "magnitude": magnitude,
                               "elastic_force": elastic.tolist(),
                               "damping_force": damping.tolist(), "rayleigh_beta": beta,
                               "axial_force": axial if np.isfinite(axial) else None,
                               "source": "support spring + implicit Rayleigh damping (not total xtip actuator reaction)"}
            else:
                self.reaction_writer.writerow((f"{simulation_time:.9g}", *([""] * 17)))
            self.reaction_log.flush()
            now = time.monotonic()
            if now - self.last_force_snapshot >= 1.0 / 30.0:
                payload = {"time": simulation_time, "updated": time.time(),
                           "session": str(SNAPSHOT),
                           "force_unit": self.force_unit,
                           "position_unit": self.position_unit,
                           "contact_force_source": self.contact_force_source,
                           "contacts": contacts, "insertion_support_force": support}
                wall_time = getattr(self.getContext(), "wall_time", None)
                if wall_time is not None:
                    payload["wall_time"] = float(wall_time)
                temporary = FORCE_STATE.with_suffix(".tmp")
                temporary.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
                try:
                    os.replace(temporary, FORCE_STATE)
                except PermissionError:
                    # Windows can briefly deny replacement while the viewer
                    # reads the previous file. Keep that complete state and
                    # retry on the next simulation step, without blocking it.
                    return
                self.last_force_snapshot = now
        except Exception as error:
            if not self.force_warning_reported:
                print(f"Force diagnostics unavailable: {error}", flush=True)
                self.force_warning_reported = True

