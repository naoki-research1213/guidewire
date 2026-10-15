from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
import shutil,json,hashlib
def main():
    path=ROOT/'vendor/BeamAdapter/src/BeamAdapter/component/mapping/AdaptiveBeamMapping.inl'
    src=source_file(path).read_text(encoding='utf-8');out=[src.split('#include')[0]];manifest=[]
    for name in ('computeIdxAndBaryCoordsForAbs','computeDistribution','applyJTonPoint','applyJonPoint','apply','applyJ','applyJT'):
        pattern=r'void\s+AdaptiveBeamMapping< TIn, TOut>::'+name+r'\('
        if name=='applyJT':pattern+=r'const core::ConstraintParams\*'
        sig,body,line=extract(src,pattern)
        sig=sig.replace('AdaptiveBeamMapping< TIn, TOut>::','Mapping::').replace('const core::MechanicalParams*','const MappingParams*').replace('const MechanicalParams*','const MappingParams*')
        out.extend([f'\n// AdaptiveBeamMapping.inl:{line}\n',sig,body,'\n'])
        manifest.append(dict(function=name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    (ROOT/'native/generated/mapping_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/mapping_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    dest=ROOT/'upstream/BeamAdapter/src/BeamAdapter/component/mapping';dest.mkdir(parents=True,exist_ok=True);copy_source(path,dest/path.name)
if __name__=='__main__':main()
