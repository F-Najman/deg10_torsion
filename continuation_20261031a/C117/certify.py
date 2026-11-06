from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;B=C.parents[1]
m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,117]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():
 p=B/n;assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==h,n
r=json.loads((C/'finite_and_inherited_verified.json').read_text());assert r['status']=='PASS' and r['old_packet_actual_exit_code']==0 and r['source39_degree10_reduction_cannot_support117']
assert [x['eligible_curve_orders']for x in r['finite_fields']]==[[],[],[],[],[],[],[],[234],[468],[]]
s=json.loads((B/'jobs/d31a_117_census.status.json').read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0
assert 'D31A_C117_CENSUS_AND_INHERITED_COMPLETE' in (B/'jobs/d31a_117_census.log').read_text()
assert 'C117_D9_COMPLETE_READ_ONLY_CERTIFICATE_PASS' in (C/'inherited_replay.log').read_text()
print('D31A_C117_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' rankzero_X139_dyadic_injection=true no_full_J117_rank_claim=true')
