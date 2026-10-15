from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
from extract_solver import COMMIT
import shutil,json,hashlib

def main():
    base=ROOT/'vendor/sofa/Sofa';out=[];manifest=[];dest=ROOT/'upstream/sofa/collision'
    def record(path,line,body):manifest.append(dict(source=path.name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    stl=base/'Component/IO/Mesh/src/sofa/component/io/mesh/MeshSTLLoader.cpp'
    src=source_file(stl).read_text(encoding='utf-8');out.append(src.split('#include')[0])
    for pattern in (r'bool\s+isBinarySTLValid\(',r'bool\s+MeshSTLLoader::readBinarySTL\(',r'bool\s+MeshSTLLoader::readSTL\('):
        sig,body,line=extract(src,pattern);out.extend([f'\n// upstream {stl.name}:{line}; {COMMIT}\n',sig,body,'\n']);record(stl,line,body)
    copy_source(stl,dest/stl.name)
    loader=base/'framework/Core/src/sofa/core/loader/MeshLoader.cpp';src=source_file(loader).read_text(encoding='utf-8')
    sig,body,line=extract(src,r'void\s+MeshLoader::addTriangle\(type::vector<Triangle\s*>& pTriangles, const Triangle& p\)')
    record(loader,line,body)
    out.extend([src.split('#include')[0],sig.replace('MeshLoader::','MeshSTLLoader::'),body,'\n']);copy_source(loader,dest/loader.name)
    mesh=base/'Component/Topology/Container/Constant/src/sofa/component/topology/container/constant/MeshTopology.cpp';src=source_file(mesh).read_text(encoding='utf-8');out.append(src.split('#include')[0])
    for name in ('EdgeUpdate::updateFromSurface','createEdgesAroundVertexArray','createEdgesInTriangleArray','createTrianglesAroundVertexArray','createTrianglesAroundEdgeArray'):
        sig,body,line=extract(src,r'void\s+MeshTopology::'+name+r'\s*\(');out.extend([f'\n// upstream {mesh.name}:{line}; {COMMIT}\n',sig,body,'\n']);record(mesh,line,body)
    copy_source(mesh,dest/mesh.name)
    (ROOT/'native/generated/mesh_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/mesh_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
