#!/usr/bin/env python3
"""Read-only audit of the complete independently frozen C37 exclusion."""
from pathlib import Path
import json,hashlib,subprocess,sys,re
H=Path(__file__).resolve().parent;B=H.parent.parent;P=H/'P37';M=P/'manifest.json';z=json.loads(M.read_text())
for f in z['files']:
 q=B/f['path'];assert q.is_file()and hashlib.sha256(q.read_bytes()).hexdigest()==f['sha256'],q
for f in z['files']:
 if f['path'].endswith('.exit.json'):assert json.loads((B/f['path']).read_text())['native_exit_code']==0,f
s=(P/'good_branches_64GB.log').read_text();actual={(int(d),int(i)):v for d,i,v in re.findall(r'NONIC_PRIME_TEST_RESULT N=37 degree=(\d+) index=(\d+) q=3 principal=(true|false)',s)}
assert actual=={(6,1):'true',(6,2):'false',(7,1):'false',(9,1):'false',(9,2):'false',(9,3):'true',(9,4):'true'}
assert 'NONIC_PRIME_GOOD_BRANCHES_COMPLETED N=37 summary=[ <6, 21, 2, 1>, <7, 18, 1, 0>, <9, 39, 4, 2> ]'in s
s=(P/'finite_Hecke_controls_v2.log').read_text();assert s.count('agrees=true')==7 and 'NONIC_FINITE_AND_HECKE_CONTROLS_COMPLETE N=37 inputs=7'in s
assert 'NONIC37_NC9_SURVIVORS_ALL_DELTA4_FIXED'in(P/'finite_nc9_survivors.log').read_text()
assert 'NONIC37_SEXTIC_LIFTS_COMPLETE all_three_residue_classes_lift=true rational37isogeny_verified=true'in(P/'sextic_lifts_control_v2.log').read_text()
assert 'NONIC37_RATIONAL_ISOGENY_FIBERS_COMPLETE'in(P/'rational_isogeny_fibers.log').read_text()
assert 'NONIC37_SEXTIC_CUSP3_NATIVE_RR_CONTROL rank=9 complete=true'in(P/'sextic_cusp3_export.log').read_text()
for cmd in [[sys.executable,str(H/'verify_FI_prefixes.py'),str(P/'FI_rows.txt')],[sys.executable,str(P/'replay_rational_two.py')],[sys.executable,str(P/'replay_sextic_cusp3.py'),'--readonly']]:
 r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,(cmd,r.stdout,r.stderr);print(r.stdout.strip())
print('NONIC37_COMPLETE_READ_ONLY_CERTIFICATE_PASS',len(z['files']),'hashes; actual exits0; complete FI and dyadic replays; seven native Hecke controls; explicit sextic lifts; full native RR control')
