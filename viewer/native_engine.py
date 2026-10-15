"""現行のファイル入出力・診断処理と抽出計算核を接続するワーカー。"""
from pathlib import Path
import os,time
from native_core import NativeCore,Model,Data
from reference_diagnostics import ExternalDriveController

class Context:
    dt=Data(.005)
    time=0.
    def getTime(self):return self.time

def run():
    mesh=Path(os.environ.get('GUIDEWIRE_VESSEL_STL',Path(__file__).parent/'assets/vessel.stl'))
    core=NativeCore(mesh,os.environ.get('GUIDEWIRE_PARAMETERS'));context=Context()
    settings=core.profile
    dt=settings['numerics']['dt_s'] if settings else .005
    stiffness=settings['mechanics_code']['support_stiffness'] if settings else 1e8
    beta=settings['mechanics_code']['rayleigh_stiffness'] if settings else .2
    force_scale=settings.get('force_scale_n_per_code',1.) if settings else 1.
    context.dt=Data(dt)
    deploy=Model(xtip=[0.],rotationInstrument=[0.],indexFirstNode=[60])
    beam=Model(position=core.read(0,7),rest_position=core.read(1,7),velocity=core.read(2,6))
    collision=Model(position=core.read(3,3),constraint=core.read(7))
    visual=Model(position=core.read(4,3));solver=Model(constraintForces=[])
    diagnostics=ExternalDriveController(deploy=deploy,beam_dofs=beam,collision_dofs=collision,visual_dofs=visual,
        constraint_solver=solver,rest_springs=Model(stiffness=[stiffness],rayleighStiffness=[0.]),
        ode_solver=Model(rayleighStiffness=[beta]),contact_visual=Model(position=[]),
        force_scale=force_scale,force_unit='N (SI unit conversion; model not experimentally validated)' if settings and 'material_si' in settings else 'SOFA scene unit; nominal N display, not SI-calibrated',
        position_unit=f"raw scene coordinate (1 unit = {settings['internal_length_m']} m)" if settings and 'material_si' in settings else 'scene mm')
    diagnostics.context=context
    print('Extracted SOFA / BeamAdapter core running; no SOFA DLL is loaded',flush=True)
    if settings:
        import json
        from reference_diagnostics import RUNTIME
        RUNTIME.mkdir(parents=True,exist_ok=True)
        (RUNTIME/'effective_parameters.json').write_text(json.dumps(settings,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f"Parameter profile: {settings['name']}; dt={dt}s",flush=True)
    simulation_started = time.perf_counter()
    context.wall_time = 0.0
    try:
        while True:
            # WindowsのこのPython環境ではmonotonicの分解能が15.625 ms。
            # 2 ms刻みの残り時間は高精度なperf_counterで測る。
            start=time.perf_counter();diagnostics.read_command();core.step(diagnostics.insert_axis,diagnostics.rotate_axis)
            state=core.read(6);context.time=float(state[0]);deploy.xtip.value=[state[1]];deploy.rotationInstrument.value=[state[2]];deploy.indexFirstNode.value=[state[3]]
            beam.position.value=core.read(0,7);beam.rest_position.value=core.read(1,7);beam.velocity.value=core.read(2,6)
            collision.position.value=core.read(3,3);collision.constraint.value=core.read(7);visual.position.value=core.read(4,3);solver.constraintForces.value=core.read(5)
            context.wall_time = time.perf_counter() - simulation_started
            diagnostics.onAnimateEndEvent(None)
            remaining=dt-(time.perf_counter()-start)
            if remaining>0:time.sleep(remaining)
    finally:
        diagnostics.contact_log.close();diagnostics.reaction_log.close();core.close()
if __name__=='__main__':run()
