from source_files import source_file, copy_source
from extract_native import ROOT
from extract_solver import COMMIT
import re,shutil,json,hashlib

def extract(src,pattern):
    masked=re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]),src)
    m=re.search(pattern,masked);assert m,pattern
    opening=masked.index('{',m.end());end=opening+1;level=1
    while level:
        level+=(masked[end]=='{')-(masked[end]=='}');end+=1
    return src[m.start():opening],src[opening:end],src[:m.start()].count('\n')+1

def main():
    helper=ROOT/'vendor/sofa/Sofa/framework/Helper/src/sofa/helper'
    dest=ROOT/'upstream/include/sofa/helper';dest.mkdir(parents=True,exist_ok=True)
    copy_source(helper/'LCPcalc.h',dest/'LCPcalc.h')
    shim=ROOT/'native/shim/sofa/helper/system/thread';shim.mkdir(parents=True,exist_ok=True)
    (shim/'CTime.h').write_text('#pragma once\n#define SOFA_HELPER_API\n')
    src=source_file((helper/'LCPcalc.cpp')).read_text(encoding='utf-8');out=[src.split('#include')[0],'\nnamespace sofa::helper {\n']
    manifest=[]
    for name in ('compute','New_GS_State'):
        sig,body,line=extract(src,r'void\s+LocalBlock33::'+name+r'\s*\(')
        out.extend([f'\n// upstream:{line}; {COMMIT}\n',sig,body,'\n'])
        manifest.append(dict(source='LCPcalc.cpp',function='LocalBlock33::'+name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    out.append('\n}\n')
    (ROOT/'native/generated/local_block_functions.inc').write_text(''.join(out),encoding='utf-8')
    lcpath=ROOT/'vendor/sofa/Sofa/Component/Constraint/Lagrangian/Solver/src/sofa/component/constraint/lagrangian/solver/LCPConstraintSolver.cpp'
    lsrc=source_file(lcpath).read_text(encoding='utf-8')
    sig,body,line=extract(lsrc,r'int\s+LCPConstraintSolver::nlcp_gaussseidel_unbuilt\s*\(')
    (ROOT/'native/generated/contact_solver_functions.inc').write_text(lsrc.split('#include')[0]+f'\n// upstream:{line}; {COMMIT}\n'+sig.replace('LCPConstraintSolver','ContactSolver')+body+'\n',encoding='utf-8')
    manifest.append(dict(source='LCPConstraintSolver.cpp',function='nlcp_gaussseidel_unbuilt',line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    d=ROOT/'upstream/sofa/constraint';d.mkdir(parents=True,exist_ok=True)
    copy_source(lcpath,d/lcpath.name);copy_source(helper/'LCPcalc.cpp',d/'LCPcalc.cpp')
    (ROOT/'native/generated/contact_solver_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
