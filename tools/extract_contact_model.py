from source_files import source_file, copy_source
from extract_native import ROOT
from extract_contact_solver import extract
import shutil,json,hashlib
def main():
    base=ROOT/'vendor/sofa/Sofa/Component';out=[];manifest=[]
    path=base/'Constraint/Lagrangian/Model/src/sofa/component/constraint/lagrangian/model/BaseContactLagrangianConstraint.inl'
    src=source_file(path).read_text(encoding='utf-8');out.append(src.split('#include')[0])
    patterns=[r'void\s+BaseContactLagrangianConstraint<DataTypes, ContactParams>::addContact\(const ContactParams& parameters, Deriv norm, '+tail for tail in (r'Coord P, Coord Q, Real contactDistance, int m1, int m2, long',r'Real contactDistance',r'Coord P, Coord Q, Real contactDistance, int m1, int m2, Coord')]
    patterns += [r'void\s+BaseContactLagrangianConstraint<DataTypes, ContactParams>::'+name+r'\(' for name in ('buildConstraintMatrix','getPositionViolation','getConstraintInfo')]
    for pattern in patterns:
        sig,body,line=extract(src,pattern);out.extend([f'\n// {path.name}:{line}\n',sig.replace('BaseContactLagrangianConstraint<DataTypes, ContactParams>::','ContactConstraint::'),body,'\n'])
        manifest.append(dict(source=path.name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    copy_source(path,ROOT/'upstream/sofa/collision'/path.name)
    path=base/'Collision/Response/Contact/src/sofa/component/collision/response/contact/BaseUnilateralContactResponse.inl';src=source_file(path).read_text(encoding='utf-8');out.append(src.split('#include')[0])
    for name in ('setDetectionOutputs','activateMappers','createResponse'):
        sig,body,line=extract(src,r'void\s+BaseUnilateralContactResponse<TCollisionModel1, TCollisionModel2, ConstraintParameters, ResponseDataTypes>::'+name+r'\(')
        out.extend([f'\n// {path.name}:{line}\n',sig.replace('BaseUnilateralContactResponse<TCollisionModel1, TCollisionModel2, ConstraintParameters, ResponseDataTypes>::','ContactResponse::'),body,'\n'])
        manifest.append(dict(source=path.name,line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest()))
    copy_source(path,ROOT/'upstream/sofa/collision'/path.name)
    (ROOT/'native/generated/contact_model_functions.inc').write_text(''.join(out),encoding='utf-8')
    (ROOT/'native/generated/contact_model_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
