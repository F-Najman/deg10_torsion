from pathlib import Path
import json,hashlib,re
C=Path(__file__).resolve().parent;B=C.parents[1]
m=json.loads((C/'final_certificate_summary.json').read_text())
assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,56]]
assert m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():
 assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
for r in m['jobs']:
 s=json.loads((B/r['status_path']).read_text())
 assert s['state']=='EXITED' and s['actual_exit_code']==r['actual_exit_code']
 if r['accepted_completion_marker']:
  assert r['accepted_completion_marker'] in (B/s['log']).read_text()
for case in [65,271,327,491,578,818,837,978,1009]:
 q=11 if case==65 else 5
 v=json.loads((C/'kummer_rr'/f'direct_case{case}_A{q}_result.m').read_text().replace('[*','[').replace('*]',']'))
 assert v==([65,11,96,72,60,192,132,0] if case==65 else [case,5,60,36,24,120,60,0])
f=json.loads((C/'inputs/independent_fourmap.json').read_text())
assert f['status']=='PASS' and f['tested_support_patterns']==1229 and f['remaining_support_patterns']==135
g=json.loads((C/'sieve_independent_genuine_isogenies.json').read_text())
assert g['status']=='PASS' and g['source_contexts']==11 and g['total_geometric_kernels']==72
c=json.loads((C/'cusp_A11_degree2_audit.json').read_text())
assert c['status']=='PASS' and c['tested_cusps_degree_at_most2']==27 and c['all_A11_columns_zero_literally']
for n in m['independent_report_paths']:
 assert json.loads((B/n).read_text())['status']=='PASS',n
print('D31B_C56_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' source_patterns1229=true source_left135=true finalRR9=true genuine72kernels=true known27_nonsquare=true')
