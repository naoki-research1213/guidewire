"""原実装を抽出したC++計算核へのデータ受け渡し。力学式はここに書かない。"""
from pathlib import Path
import ctypes as ct
import os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
class NativeCore:
    def __init__(self,mesh,parameter_path=None):
        path=Path(os.environ.get('GUIDEWIRE_CORE_LIBRARY', ROOT/'runtime/guidewire_core.dll'))
        if not path.exists():raise FileNotFoundError('Run build_native.ps1 -FetchToolchain before launching')
        self.dll=ct.CDLL(str(path))
        self.dll.gw_create.argtypes=[ct.c_char_p];self.dll.gw_create.restype=ct.c_void_p
        self.dll.gw_destroy.argtypes=[ct.c_void_p]
        self.dll.gw_step.argtypes=[ct.c_void_p,ct.c_double,ct.c_double];self.dll.gw_step.restype=ct.c_int
        self.dll.gw_read.argtypes=[ct.c_void_p,ct.c_int,ct.POINTER(ct.c_double),ct.c_int];self.dll.gw_read.restype=ct.c_int
        self.dll.gw_error.restype=ct.c_char_p
        self.profile=None
        mesh_bytes=str(Path(mesh).resolve()).encode('mbcs' if os.name=='nt' else 'utf-8')
        if parameter_path:
            from scene_parameters import load_parameters
            self.profile,values=load_parameters(parameter_path)
            self.dll.gw_create_parameters.argtypes=[ct.c_char_p,ct.POINTER(ct.c_double),ct.c_int]
            self.dll.gw_create_parameters.restype=ct.c_void_p
            array=(ct.c_double*len(values))(*values)
            self.handle=self.dll.gw_create_parameters(mesh_bytes,array,len(values))
        else:self.handle=self.dll.gw_create(mesh_bytes)
        if not self.handle:raise RuntimeError(self.dll.gw_error().decode('utf-8',errors='replace'))
    def read(self,kind,columns=None):
        n=self.dll.gw_read(self.handle,kind,None,0)
        result=np.empty(n,dtype=np.float64)
        self.dll.gw_read(self.handle,kind,result.ctypes.data_as(ct.POINTER(ct.c_double)),n)
        return result.reshape(-1,columns) if columns else result
    def step(self,insert,rotate):
        if self.dll.gw_step(self.handle,float(insert),float(rotate)):
            raise RuntimeError(self.dll.gw_error().decode('utf-8',errors='replace'))
    def close(self):
        if self.handle:self.dll.gw_destroy(self.handle);self.handle=None

class Data:
    def __init__(self,value):self.value=value
    def getValueString(self):return ' '.join(format(float(x),'.17g') for x in np.asarray(self.value).reshape(-1))
class Model:
    def __init__(self,**values):
        for name,value in values.items():setattr(self,name,Data(value))
    def getData(self,name):return getattr(self,name,None)
