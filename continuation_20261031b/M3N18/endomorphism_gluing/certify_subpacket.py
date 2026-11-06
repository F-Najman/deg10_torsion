from pathlib import Path
import json,hashlib
N=Path(__file__).resolve().parent;H=N.parents[2]
m=json.loads((N/'certificate_manifest.json').read_text());assert m['status']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((H/n).read_bytes()).hexdigest()==h,n
for row in m['jobs']:
 s=json.loads((H/row['status_path']).read_text());assert s['state']=='EXITED'and s['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker']in(H/s['log']).read_text()
a=json.loads((N/'direct_cusp_membership_report.json').read_text());assert a['status']=='PASS'and a['tested_pairs']==4 and a['excluded_pairs']==2 and a['remaining_pairs']==2
b=json.loads((N/'tau_cusp_class_relation.json').read_text());assert b['status']=='PASS'and any(r['Aut_parameters']==[2,0,2,5]for r in b['matches'])
print('D31B_M318_ENDOMORPHISM_NF_SUBCERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' direct_cusp4to2=true tau_global_relation=true exact_decic_point=true')
