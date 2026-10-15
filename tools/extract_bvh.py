from source_files import source_file, copy_source
from extract_native import ROOT,functions
from extract_contact_solver import extract
from extract_solver import COMMIT
import shutil,json,hashlib

def main():
    base=ROOT/'vendor/sofa/Sofa/Component/Collision'
    paths=[base/'Geometry/src/sofa/component/collision/geometry/CubeCollisionModel.cpp',
           base/'Detection/Algorithm/src/sofa/component/collision/detection/algorithm/BVHNarrowPhase.cpp']
    out=[];manifest=[];dest=ROOT/'upstream/sofa/collision'
    for path in paths:
        src=source_file(path).read_text(encoding='utf-8');cls=path.stem;out.append(src.split('#include')[0])
        names={'CubeCollisionModel':['resize','setParentOf','addCube','updateCube','updateCubes','getInternalChildren','getExternalChildren','computeBoundingTree'],
               'BVHNarrowPhase':['isSelfCollision','addCollisionPair','initializeExternalCells','processExternalCell','processInternalCell','visitCollisionElements','visitExternalChildren','finalCollisionPairs','getCollisionModelsFromTestPair','isRangeEmpty']}[cls]
        # A function may have multiple overloads (setParentOf).
        import re
        masked=re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',lambda m:re.sub(r'[^\n]',' ',m[0]),src)
        pattern=r'(?:void|bool|Index|std::pair<core::CollisionElementIterator,core::CollisionElementIterator>|std::pair<core::CollisionModel\*, core::CollisionModel\*>)\s+'+cls+r'::(?:'+'|'.join(names)+r')\s*\('
        for m in re.finditer(pattern,masked):
            opening=masked.index('{',m.end());end=opening+1;level=1
            while level:
                level+=(masked[end]=='{')-(masked[end]=='}');end+=1
            sig=src[m.start():opening];body=src[opening:end];line=src[:m.start()].count('\n')+1
            out.extend([f'\n// upstream {path.name}:{line}; {COMMIT}\n',sig,body,'\n'])
            manifest.append(dict(source=path.name,signature=sig,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
        copy_source(path,dest/path.name)
    prox=base/'Detection/Intersection/src/sofa/component/collision/detection/intersection/BaseProximityIntersection.cpp'
    src=source_file(prox).read_text(encoding='utf-8');sig,body,line=extract(src,r'bool\s+BaseProximityIntersection::testIntersection\(Cube&')
    out.extend([src.split('#include')[0],sig.replace('BaseProximityIntersection::testIntersection','LMDIntersector::testCubes'),body,'\n'])
    copy_source(prox,dest/prox.name)
    for cls in ('PointCollisionModel','LineCollisionModel','TriangleCollisionModel'):
        p=base/('Geometry/src/sofa/component/collision/geometry/'+cls+'.inl');sig,body,line=functions(p,cls,'computeBoundingTree')[0]
        out.extend([source_file(p).read_text(encoding='utf-8').split('#pragma once')[0],f'\n// upstream {p.name}:{line}\ntemplate<class DataTypes>\n',sig,body,'\n'])
        copy_source(p,dest/p.name)
    (ROOT/'native/generated/bvh_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/bvh_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
