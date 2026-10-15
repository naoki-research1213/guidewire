from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
import shutil,json,hashlib
def main():
    path=ROOT/'vendor/sofa/Sofa/Component/Mapping/Linear/src/sofa/component/mapping/linear/SubsetMapping.inl'
    src=source_file(path).read_text(encoding='utf-8');out=[src.split('#include')[0]];manifest=[]
    for name,param in (('apply','MechanicalParams'),('applyJT','ConstraintParams')):
        sig,body,line=extract(src,r'void\s+SubsetMapping<TIn, TOut>::'+name+r'\s*\( const core::'+param+r'\*?\s*')
        sig=sig.replace('SubsetMapping<TIn, TOut>::','SubsetMapper::').replace('core::MechanicalParams','MappingParams')
        out.extend([f'\n// SubsetMapping.inl:{line}\n',sig,body,'\n']);manifest.append(dict(function=name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    (ROOT/'native/generated/subset_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/subset_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8');copy_source(path,ROOT/'upstream/sofa/collision'/path.name)
if __name__=='__main__':main()
