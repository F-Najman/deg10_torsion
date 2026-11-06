from pathlib import Path
import json,hashlib,ast
C=Path(__file__).resolve().parent;B=C.parents[1];m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,78],[3,39]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
r=json.loads((C/'independent_genuine_and_known_lifts.json').read_text());assert r['status']=='PASS' and r['actual_geometric5_kernels']==60 and r['known_decic_full_torsion_order']==39
v=ast.literal_eval((C/'native_replay.m').read_text().replace('[*','[').replace('*]',']'))
for ci,d,kp,kd,hec,h0 in v:
 x=next(x for x in r['source_cases']if x['case']==ci);assert x['degree']==d and kp==x['source_key']and kd==x['diamond5_key']and hec==x['target_keys']and h0==x['native_h0']
assert [x[5]for x in v]==[1,0,0,0,0,0,0,0,0,1]
assert [(x['degree'],x['principal_source_case'])for x in r['principal_orbit_known_lift_bindings']]==[(9,1),(10,10)]
for row in m['jobs']:
 st=json.loads((B/row['status_path']).read_text());assert st['state']=='EXITED' and st['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker']in(B/st['log']).read_text()
print('D31A_C78_M3N39_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' all_decic_X139_known=true fulltorsion39=true actual5kernels60=true')
