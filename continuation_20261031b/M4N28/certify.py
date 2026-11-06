from pathlib import Path
import hashlib,json,ast,re
C=Path(__file__).resolve().parent;H=C.parents[1]
def rd(p):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',p.read_text().replace('[*','[').replace('*]',']'))))
m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS'and m['parent_mathematical_review']=='PASS'and m['excluded_groups']==[[4,28]]
for p,s in m['files_sha256'].items():assert hashlib.sha256((H/p).read_bytes()).hexdigest()==s,p
for j in m['jobs']:
 s=json.loads((H/j['status_path']).read_text());assert s['state']=='EXITED'and s['actual_exit_code']==0;assert j['accepted_completion_marker']in(H/s['log']).read_text();assert hashlib.sha256((H/s['script']).read_bytes()).hexdigest()==s['script_sha256']
a=json.loads((C/'parent_independent_fourmap_v3.json').read_text());assert a['status']=='PASS'and a['source_models']==174 and a['genuine_marked_maps']==696 and a['remaining_cases']==19
u=json.loads((C/'parent_independent_union.json').read_text());assert u['status']=='PASS'and u['cases']==19
S=C.parent/'Y228_shared';a=rd(S/'parent_native_all45_and_neighbors_v2.m');assert len(a[0])==45 and len(a[1])==4 and a[2:]==[272,129]
left={r[0]for r in a[0]if any(r[4])};ids={v['union_case']for v in u['checks']};assert ids&left=={27,42,43,45};assert {r[0]for r in a[1]}==ids&left and all(r[1]==1 and r[5]==0 for r in a[1])
g=[json.loads((S/f'independent_genuine_v2_s{i}.json').read_text())for i in range(2)];assert all(x['status']=='PASS'for x in g);assert sum(x['component_tests']for x in g)==53;assert sum('neighbor_case'in v for x in g for v in x['checks'])==4
print('D31B_M428_FINAL_CERTIFICATE_COMPLETE',len(m['files_sha256']),'eligible174=true fourmaps696=true nativeclasses656=true genuine5kernels342=true final4negative=true')
