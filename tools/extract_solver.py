from source_files import source_file, copy_source
"""SOFA BTDLinearSolverとBTDMatrix::Blockを原文から抽出する。"""
from extract_native import ROOT,functions
import re,hashlib,json,shutil

SOURCE=ROOT/'vendor/sofa/Sofa/Component/LinearSolver/Direct/src/sofa/component/linearsolver/direct/BTDLinearSolver.inl'
COMMIT='7c18e95d5c5f2839079892c69e7d89a313c79603'

def main():
    src=source_file(SOURCE).read_text(encoding='utf-8')
    masked=re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]),src)
    out=[src.split('#pragma once')[0]];manifest=[]
    names='my_identity|invert|solve|computeMinvBlock|getMinvElement|init_partial_inverse|init_partial_solve|bwdAccumulateRHinBloc|bwdAccumulateLHGlobal|fwdAccumulateRHGlobal|fwdComputeLHinBloc|partial_solve'
    for m in re.finditer(r'(?:void|double)\s+BTDLinearSolver<Matrix,Vector>::('+names+r')\s*\(',masked):
        opening=masked.index('{',m.end());end=opening+1;level=1
        while level:
            level+=(masked[end]=='{')-(masked[end]=='}');end+=1
        sig=src[m.start():opening];body=src[opening:end]
        adapted=sig.replace('BTDLinearSolver<Matrix,Vector>','BTD')
        line=src[:m.start()].count('\n')+1
        out.extend([f'\n// upstream BTDLinearSolver.inl:{line}, {COMMIT}\n',adapted,body,'\n'])
        manifest.append(dict(function=m[1],line=line,body_sha256=hashlib.sha256(body.encode()).hexdigest(),body_changed=False,commit=COMMIT))
    (ROOT/'native/generated/btd_functions.inc').write_text(''.join(out),encoding='utf-8')
    blockpath=ROOT/'vendor/sofa/Sofa/framework/LinearAlgebra/src/sofa/linearalgebra/BTDMatrix.h'
    blocksrc=source_file(blockpath).read_text(encoding='utf-8')
    start=blocksrc.index('    class TransposedBlock');end=blocksrc.index('    typedef Block SubMatrixType;')
    (ROOT/'native/generated/btd_block.inc').write_text(blocksrc.split('#pragma once')[0]+'\nconstexpr size_t BSIZE=6;\nusing T=double;\n'+blocksrc[start:end],encoding='utf-8')
    dest=ROOT/'upstream/sofa/linearalgebra';dest.mkdir(parents=True,exist_ok=True)
    copy_source(blockpath,dest/blockpath.name);copy_source(SOURCE,dest/SOURCE.name)
    vectorpath=blockpath.parent/'BlockVector.h';vsrc=source_file(vectorpath).read_text(encoding='utf-8')
    start=vsrc.index('    class Block :');end=vsrc.index('    typedef Block SubVectorType;')
    (ROOT/'native/generated/btd_vector.inc').write_text(vsrc.split('#pragma once')[0]+'\nconstexpr size_t N=6;\n'+vsrc[start:end].replace('class Block :','class SubVector :'),encoding='utf-8')
    copy_source(vectorpath,dest/vectorpath.name)
    (ROOT/'native/generated/solver_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    spring=ROOT/'vendor/sofa/Sofa/Component/SolidMechanics/Spring/src/sofa/component/solidmechanics/spring/RestShapeSpringsForceField.inl'
    sig,body,line=functions(spring,'RestShapeSpringsForceField','addForce')[0]
    (ROOT/'native/generated/spring_functions.inc').write_text(source_file(spring).read_text(encoding='utf-8').split('#pragma once')[0]+f'\n// upstream:{line}; {COMMIT}\ntemplate<class DataTypes>\n'+sig.replace('RestShapeSpringsForceField','Spring')+body+'\n',encoding='utf-8')
    copy_source(spring,dest/spring.name)
    posepath=ROOT/'vendor/sofa/Sofa/framework/DefaultType/src/sofa/defaulttype/RigidCoord.h'
    psrc=source_file(posepath).read_text(encoding='utf-8');start=psrc.index('    void operator +=(const Deriv& dg)')
    opening=psrc.index('{',start);end=opening+1;level=1
    while level:
        level+=(psrc[end]=='{')-(psrc[end]=='}');end+=1
    body=psrc[opening:end]
    (ROOT/'native/generated/rigid_integrate.inc').write_text(psrc.split('#pragma once')[0]+'\n// Original RigidCoord<3,real>::operator+=(Deriv), unchanged body.\nvoid integratePose(Pose& pose,const RigidDerivativeView& dg)\n{\nusing real=double;auto& center=pose.p;auto& orientation=pose.q;\n'+body+'\n}\n',encoding='utf-8')
    copy_source(posepath,dest/posepath.name)
    ccpath=ROOT/'vendor/sofa/Sofa/Component/Constraint/Lagrangian/Correction/src/sofa/component/constraint/lagrangian/correction/LinearSolverConstraintCorrection.inl'
    ccsrc=source_file(ccpath).read_text(encoding='utf-8');out=[ccsrc.split('#pragma once')[0]]
    for name in ('resetForUnbuiltResolution','addConstraintDisplacement','setConstraintDForce','getBlockDiagonalCompliance','applyContactForce'):
        sig,body,line=functions(ccpath,'LinearSolverConstraintCorrection',name)[0]
        out.extend([f'\n// upstream:{line}; {COMMIT}\ntemplate<class DataTypes>\n',sig.replace('LinearSolverConstraintCorrection','Correction'),body,'\n'])
    (ROOT/'native/generated/correction_functions.inc').write_text(''.join(out),encoding='utf-8')
    copy_source(ccpath,dest/ccpath.name)
    # BTD's templated sparse scalar J * Minv * Jt multiplication.
    start=src.index('bool BTDLinearSolver<Matrix,Vector>::addJMInvJt(RMatrix&')
    opening=src.index('{',start);end=opening+1;level=1
    while level:
        level+=(src[end]=='{')-(src[end]=='}');end+=1
    (ROOT/'native/generated/compliance_functions.inc').write_text(src.split('#pragma once')[0]+'\ntemplate<class RMatrix,class JMatrix>\n'+src[start:end].replace('BTDLinearSolver<Matrix,Vector>','BTD')+'\n',encoding='utf-8')
if __name__=='__main__':main()
