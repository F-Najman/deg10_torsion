#!/usr/bin/env python3
"""Complete scoped C65 ordinary branch certificate; not a full C65 exclusion."""
from pathlib import Path
import json,hashlib,subprocess,sys,time
C=Path(__file__).resolve().parent;R=C.parents[1];H=C/'composite65';t=time.monotonic();m=json.loads((H/'ordinary_manifest.json').read_text());assert m['degree']==9 and m['group']==[1,65]and m['status']=='COMPLETE EXCLUSION OF ORDINARY DEGREE7_8_9 REDUCTION BRANCHES'
for a in m['files']:assert hashlib.sha256((R/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
z=subprocess.run([sys.executable,str(H/'verify_all_components_readonly.py')]);assert z.returncode==0
print('NONIC65_ORDINARY_FROZEN_CERTIFICATE_PASS',json.dumps({'scope':'ordinary degree7/8/9 only; CM6 and allbad separate','files':len(m['files']),'selected_native_exits_zero':13,'binary_cases':336,'binary_rank':121,'manifest_sha256':hashlib.sha256((H/'ordinary_manifest.json').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t},sort_keys=True))
