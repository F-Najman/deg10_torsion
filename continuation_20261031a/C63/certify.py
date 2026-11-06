from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;B=C.parents[1];m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,63]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
r=json.loads((C/'independent_final_controls.json').read_text());assert r['status']=='PASS' and r['native_RR_tests_bound_to_genuine_isogenies']==21 and r['total_pure_cusp_patterns']==239238 and r['actual_geometric5_kernels']==36 and r['manual_Velu_checks']==72
for row in m['jobs']:
 s=json.loads((B/row['status_path']).read_text());assert s['state']=='EXITED' and s['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker']in(B/s['log']).read_text()
assert m['failed_diagnostics_not_certificates']==['d31a_63_independent_bind_boundary']
print('D31A_C63_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' cusp_patterns239238=true independent_RR21=true genuine5kernels36=true')
