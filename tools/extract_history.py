from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
import json,hashlib
def main():
    path=ROOT/'vendor/sofa/Sofa/Component/Constraint/Lagrangian/Solver/src/sofa/component/constraint/lagrangian/solver/LCPConstraintSolver.cpp'
    src=source_file(path).read_text(encoding='utf-8');out=[src.split('#include')[0]];manifest=[]
    for name in ('computeInitialGuess','keepContactForcesValue'):
        sig,body,line=extract(src,r'void\s+LCPConstraintSolver::'+name+r'\(')
        out.extend([f'\n// LCPConstraintSolver.cpp:{line}\n',sig.replace('LCPConstraintSolver::','ContactSolver::'),body,'\n'])
        manifest.append(dict(function=name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    (ROOT/'native/generated/history_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/history_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
