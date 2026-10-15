from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
from extract_solver import COMMIT
import shutil,json,hashlib

def main():
    path=ROOT/'vendor/sofa/Sofa/Component/Collision/Detection/Intersection/src/sofa/component/collision/detection/intersection/LocalMinDistance.cpp'
    src=source_file(path).read_text(encoding='utf-8');out=[src.split('#include')[0]];manifest=[]
    patterns=[r'int\s+LocalMinDistance::computeIntersection\(Triangle& e2, Point& e1,',
              r'int\s+LocalMinDistance::computeIntersection\(Line& e2, Point& e1,',
              r'int\s+LocalMinDistance::computeIntersection\(Line& e1, Line& e2,']
    patterns += [r'bool\s+LocalMinDistance::testValidity\('+typ+r'\s*&' for typ in ('Point','Line','Triangle')]
    for pattern in patterns:
        sig,body,line=extract(src,pattern);out.extend([f'\n// upstream:{line}; {COMMIT}\n',sig,body,'\n'])
        manifest.append(dict(signature=sig,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest(),commit=COMMIT))
    (ROOT/'native/generated/collision_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/collision_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    dest=ROOT/'upstream/sofa/collision';dest.mkdir(parents=True,exist_ok=True);copy_source(path,dest/path.name)
if __name__=='__main__':main()
