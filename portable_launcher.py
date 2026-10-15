"""Portable entry point; numerical methods are in the unchanged native core."""
from pathlib import Path
import json
import os
import sys
import traceback
import shutil
import tempfile

def root_path():
    return Path(sys.executable).resolve().parent if getattr(sys, 'frozen', False) else Path(__file__).resolve().parent

def data_path():
    override = os.environ.get('GUIDEWIRE_PORTABLE_DATA')
    base = Path(override) if override else Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir()))/'GuidewireSimulator'
    base.mkdir(parents=True, exist_ok=True)
    return base

def main():
    root = root_path()
    data = data_path()
    sys.path.insert(0, str(root / 'viewer'))
    os.environ['GUIDEWIRE_CORE_LIBRARY'] = str(root / 'runtime/guidewire_core.dll')
    if '--worker' in sys.argv:
        from native_engine import run
        run()
        return 0
    # Explicit paths prevent accidental use of another local installation.
    runtime = data / 'sessions' / str(os.getpid())
    runtime.mkdir(parents=True, exist_ok=True)
    os.environ['GUIDEWIRE_RUNTIME_DIR'] = str(runtime)
    config = data / 'config'
    config.mkdir(parents=True, exist_ok=True)
    for name in ('parameters.json', 'viewer.json'):
        if not (config/name).exists():
            shutil.copyfile(root/'config'/name, config/name)
    os.environ['GUIDEWIRE_PARAMETERS'] = str(config / 'parameters.json')
    os.environ['GUIDEWIRE_VIEWER_CONFIG'] = str(config / 'viewer.json')
    os.environ['GUIDEWIRE_VESSEL_STL'] = str(root / 'assets/public_sample.stl')
    profile = json.loads((config/'parameters.json').read_text(encoding='utf-8'))
    os.environ['GUIDEWIRE_INSERTION_DISPLAY_SCALE'] = str(profile['internal_length_m']*1000)
    os.environ.pop('GUIDEWIRE_VIEWER_SNAPSHOT',None)
    os.environ.pop('GUIDEWIRE_VIEWER_FORCE_STATE',None)
    if '--self-test' in sys.argv:
        import numpy as np
        import PIL
        import tkinter as tk
        from native_core import NativeCore
        core = NativeCore(root/'assets/public_sample.stl', config/'parameters.json')
        try:
            for _ in range(3): core.step(1.,0.)
            state=core.read(6)
            assert np.all(np.isfinite(core.read(0))) and state[0]>0 and state[1]>0
            result={'ok': True, 'python':sys.version.split()[0], 'numpy':np.__version__,
                    'pillow':PIL.__version__, 'simulation_time':float(state[0]),
                    'tip_coordinate':float(state[1]), 'beta':profile['mechanics_code']['rayleigh_stiffness']}
            screen=tk.Tk();screen.withdraw();result['tk']=screen.tk.call('info','patchlevel');screen.destroy()
            (data/'self-test-result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
        finally: core.close()
        return 0
    if '--smoke-test' in sys.argv:
        import tkinter as tk
        from xray_viewer import ViewerApp
        screen=tk.Tk()
        app=ViewerApp(screen,root/'assets/public_sample.stl',None)
        # Exercise the actual executable's worker and live state exchange.
        attempts=[0]
        def check():
            attempts[0]+=1
            state=app.latest_force_state
            if state and state.get('time',0)>0 and 'wall_time' in state:
                app.open_settings()
                screen.update_idletasks()
                screen.grab_set() if False else None
                from PIL import ImageGrab
                x,y=screen.winfo_rootx(),screen.winfo_rooty()
                ImageGrab.grab(bbox=(x,y,x+screen.winfo_width(),y+screen.winfo_height())).save(data/'smoke-test-viewer.png')
                (data/'smoke-test-result.json').write_text(json.dumps({'ok':True,'time':state['time'],'wall_time':state['wall_time'],'settings_opened':True}),encoding='utf-8')
                app.on_close()
            elif attempts[0]>90:
                app.on_close()
                raise RuntimeError('No live state from portable worker within startup timeout')
            else: screen.after(500,check)
        screen.after(1000,check);screen.mainloop()
        if not (data/'smoke-test-result.json').exists():raise RuntimeError('Portable GUI smoke test failed')
        return 0
    sys.argv=['xray_viewer.py','--stl',str(root/'assets/public_sample.stl')]
    from xray_viewer import main as viewer_main
    return viewer_main()

if __name__=='__main__':
    root=root_path()
    log=data_path()/('startup-'+str(os.getpid())+'.log')
    # PyInstaller windowed executables do not provide console streams.
    if sys.stdout is None:sys.stdout=log.open('a',encoding='utf-8',buffering=1)
    if sys.stderr is None:sys.stderr=sys.stdout
    try:raise SystemExit(main())
    except Exception:
        with log.open('a',encoding='utf-8') as stream:traceback.print_exc(file=stream)
        if '--worker' not in sys.argv:
            import ctypes
            ctypes.windll.user32.MessageBoxW(None,'起動に失敗しました。ログ: '+str(log),'Guidewire Simulator',16)
        raise
