"""Read-only certificate for C66 exclusion via complete nonic X1(33) classification."""
from pathlib import Path
import json,hashlib,subprocess,sys,shutil
p=Path(__file__).resolve().parent;sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();J=lambda r:json.loads((p/r).read_text())
m=J('manifest66.json')
for r,v in m['files'].items():
 q=(p/r).resolve();assert q.is_file()and q.stat().st_size==v['bytes']and sha(q)==v['sha256'],r
for f in ['rebound/classification33','frobenius2_verify33','short_classes33','replay33','full_torsion33','F5_rr/basis_jets33','F5_rr/verify33','F5_rr/replay33']:
 st=J(f+'.exit.json');assert st['actual_exit_code']==0 and st['completed_marker'],f
rootmodel=J('model_retained/build33_pool_complete.json');assert rootmodel['all_complete']and len(rootmodel['jobs'])==1;st=rootmodel['jobs'][0];assert st['actual_exit_code']==0 and st['completion_marker'];assert st['source_sha256']==sha(p/'model_retained/build33.m')and st['log_sha256']==sha(p/'model_retained/build33.log')
ident=J('model_retained/source_package_identity33.json');assert ident['complete']and len(ident['files'])==195
for z in ident['files']:
 assert sha(p/z['retained_path'])==z['sha256']
assert J('independent_replay33.json')['finite_class_tests']==155770
assert J('F5_rr/independent_certificate33.json')['all_ten_h0_zero']
assert J('full_torsion33.json')['full_global_torsion_order']==33
sage=shutil.which('sage');assert sage,'Sage is required for independent arithmetic replays'
for script in ['replay33.py','F5_rr/replay33.py','full_torsion33.py']:
 r=subprocess.run([sage,'-python',str(p/script),'--check'],capture_output=True,text=True);assert r.returncode==0,(script,r.stdout,r.stderr);print(r.stdout.strip())
print('NONIC66_COMPLETE_SOURCE33_CERTIFICATE_PASS',len(m['files']),'files; C66 excluded, full nonic33 torsion classified')
