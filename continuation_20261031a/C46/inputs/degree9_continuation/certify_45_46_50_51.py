#!/usr/bin/env python3
"""Read-only replay of four completed nonic cyclic exclusions."""
from pathlib import Path
import json,hashlib,re,subprocess,sys,ast
H=Path(__file__).resolve().parent;B=H.parents[1]
M=json.loads((H/'packet_45_46_50_51.json').read_text());assert M['degree']==9 and M['complete'] and M['excluded_orders']==[45,46,50,51]
for e in M['files']:assert hashlib.sha256((B/e['path']).read_bytes()).hexdigest()==e['sha256'],e['path']
def native(tag,marker):
 p=H/(tag+'.log');d=json.loads((H/(tag+'_status.json')).read_text());s=p.read_text()
 assert d['native_exit_code']==0 and d['log_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
 assert marker in s and not re.search(r'Runtime error|User error|Assertion failed|syntax error|failed memory request',s,re.I)
 return s
native('independent_counts_46_50_51','DEGREE9_FAST_COUNTS_DONE complete=true')
s=native('native_controls_short_fields','NONIC_NATIVE_CONTROLS_COMPLETE levels=46,50,51 source_degrees=9,6,9')
assert all(f'NONIC_NATIVE_HECKE_EQUAL N={n}'in s for n in[46,50,51])
s=native('N51_good9_and_pure9','NONIC51_GOOD9_AND_PURE9_DONE all_places=80 representatives=5 complete=true')
assert 'NONIC51_PURE9_CUSP_CONTROL rank=9 covered=16 dimensionU=62'in s
assert [(int(i),q,f)for i,q,f in re.findall(r'NONIC51_GOOD9_CASE index=(\d+) q=(\d+) principal=(\w+)',s)]==[(i,'5','false')for i in range(1,6)]
s=native('N45_good9_and_pure9_retry','NONIC45_GOOD9_AND_PURE9_DONE all_places=60 representatives=5 complete=true')
assert 'NONIC45_PURE9_CUSP_CONTROL rank=9 covered=12 dimensionU=39'in s and 'NONIC45_FULL_NATIVE_EQUAL degree=72'in s
assert [(int(i),q,f)for i,q,f in re.findall(r'NONIC45_GOOD9_CASE index=(\d+) q=(\d+) principal=(\w+)',s)]==[(i,'7','false')for i in range(1,6)]
native('independent_counts45','DEGREE9_FAST_COUNTS_DONE complete=true')
s=native('N45_quotient_injectivity_control','INJECTIVE_REDUCTION_3_15_DONE')
assert 'EQUAL_TWO_PRIMARY_ORDERS 8192 8192'in s
for name in ['N45','N51']:
 snap=json.loads((H/(name+'_geometry_manifest.json')).read_text());items=snap['files'];items=items.items()if isinstance(items,dict)else[(e['path'],e['sha256'])for e in items]
 for p,d in items:assert hashlib.sha256((H/p).read_bytes()).hexdigest()==d,p
 data=H/(name+'_geometry/cyclic_sieves/'+name+('_completion_jets.m'if name=='N45'else'_frozen_jets.m'))
 cmd=[sys.executable,str(H/(name+'_geometry/cyclic_sieves/replay_cartier_projection.py')),str(data)]
 r=subprocess.run(cmd,capture_output=True,text=True,check=True);print(r.stdout.strip())
 r=subprocess.run([sys.executable,str(H/(name+'_replay_degree9.py')),'--readonly'],capture_output=True,text=True,check=True)
 assert name+'_DEGREE9_CUSP_REPLAY_DONE'in r.stdout;print(r.stdout.splitlines()[-1])
 j=json.loads((H/(name+'_degree9_FI_replay.json')).read_text())
 if name=='N51':assert j['tested_divisors']==3333248 and j['exceptions_count']==0
 else:
  assert j['tested_divisors']==377348 and j['exceptions_count']==8
  weights=ast.literal_eval(re.search(r'ds:=(.*?);',data.read_text())[1])
  assert all(e['rank']==8 and sum(a*w for a,w in zip(e['vector'],weights))==9 and all(w<=2 for a,w in zip(e['vector'],weights)if a)for e in j['exceptions'])
r=subprocess.run([sys.executable,str(H/'N45_quotient_cusp_labels.py')],capture_output=True,text=True,check=True);assert 'NONIC45_QUOTIENT_RATIONAL_CUSPS_IDENTIFIED_PASS'in r.stdout
r=subprocess.run([sys.executable,str(H/'audit_46_50_D9.py')],capture_output=True,text=True,check=True);assert 'NONIC46_50_READONLY_COVERAGE_PASS'in r.stdout;print(r.stdout.strip())
s=(H/'N45_geometry/cyclic_sieves/N45_p2_good_only.log').read_text();assert 'RESULT type=[ 8 ] total=5 survivors=0'in s and 'total_survivors=0'in s and 'CC complement degree <= 1:'in s
print('NONIC45_46_50_51_READONLY_CERTIFICATE_PASS complete=true orders=45,46,50,51')
