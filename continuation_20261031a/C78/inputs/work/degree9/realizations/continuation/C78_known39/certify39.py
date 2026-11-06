"""Read-only known nonic39 global/full-torsion and positive-reduction binding audit."""
from pathlib import Path
import json,hashlib,shutil,subprocess
p=Path(__file__).resolve().parent;J=lambda f:json.loads((p/f).read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();m=J('manifest39.json')
for r,v in m['files'].items():
 q=p/r;assert q.stat().st_size==v['bytes']and sha(q)==v['sha256'],r
for f in ['verify39','replay39']:
 d=J(f+'.exit.json');assert d['actual_exit_code']==0 and d['completed_marker']
b=J('positive_divisors39.json');assert b['global_closed_nonic_divisors']==4 and len(b['known_positive_divisors'])==4
assert J('independent39.json')['full_global_torsion_order']==39
sage=shutil.which('sage');assert sage;r=subprocess.run([sage,'-python',str(p/'replay39.py'),'--check'],capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr);print(r.stdout.strip())
print('NONIC78_KNOWN39_READONLY_CERTIFICATE_PASS',len(m['files']),'files; known4nonic39divisors only, no C78exclusion asserted')
