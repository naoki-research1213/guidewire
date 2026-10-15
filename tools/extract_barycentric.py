from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
import shutil,json,hashlib,re
def main():
    path=ROOT/'vendor/sofa/Sofa/Component/Mapping/Linear/src/sofa/component/mapping/linear/BarycentricMappers/BarycentricMapperMeshTopology.inl'
    src=source_file(path).read_text(encoding='utf-8');out=[src.split('#include')[0]];manifest=[]
    for name in ('addPointInLine','addPointInTriangle','createPointInLine','createPointInTriangle','apply','applyJT'):
        ret=r'typename\s+BarycentricMapperMeshTopology<In, Out>::Index\s*' if name not in ('apply','applyJT') else r'void\s+'
        pattern=ret+r'BarycentricMapperMeshTopology<In,Out>::'+name+r'\s*\('
        if name=='applyJT':pattern+=r' typename In::MatrixDeriv'
        sig,body,line=extract(src,pattern)
        sig=re.sub(r'typename\s+BarycentricMapperMeshTopology<In, Out>::Index', 'unsigned',sig).replace('BarycentricMapperMeshTopology<In,Out>::','BaryMapper::')
        out.extend([f'\n// BarycentricMapperMeshTopology.inl:{line}\n',sig,body,'\n'])
        manifest.append(dict(function=name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    (ROOT/'native/generated/barycentric_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/barycentric_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    copy_source(path,ROOT/'upstream/sofa/collision'/path.name)
if __name__=='__main__':main()
