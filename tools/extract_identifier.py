from pathlib import Path
import json,hashlib
from extract_native import ROOT
from extract_contact_solver import extract
def main():
    src=ROOT/'upstream/sofa/collision';dest=ROOT/'native/generated'
    header=(src/'ContactIdentifier.h').read_text(encoding='utf-8')
    adapted=header.replace('#include <sofa/component/collision/response/contact/config.h>','').replace('#include <sofa/core/collision/DetectionOutput.h>','').replace('SOFA_COMPONENT_COLLISION_RESPONSE_CONTACT_API','')
    (dest/'contact_identifier.h').write_text(adapted,encoding='utf-8')
    text=(src/'ContactIdentifier.cpp').read_text(encoding='utf-8')
    (dest/'contact_identifier.inc').write_text(text.replace('#include <sofa/component/collision/response/contact/ContactIdentifier.h>','#include "contact_identifier.h"'),encoding='utf-8')
    sig,body,line=extract(text,r'ContactIdentifier::ContactIdentifier\(')
    (dest/'identifier_manifest.json').write_text(json.dumps([dict(source='ContactIdentifier.cpp',line=line,body_changed=False,body_sha256=hashlib.sha256(body.encode()).hexdigest())],indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
