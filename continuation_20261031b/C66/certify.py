from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;B=C.parents[1];m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,66]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
for row in m['jobs']:
 s=json.loads((B/row['status_path']).read_text());assert s['state']=='EXITED' and s['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker'] in (B/s['log']).read_text()
a=json.loads((C/'F2_independent_verified.json').read_text());assert a['status']=='PASS' and a['native_points_bound_to_actual_33torsion']==250 and a['independent_rational_image_classes']==321210 and a['independent_effective_divisors']==321230
b=json.loads((C/'gonal_genuine_independent.json').read_text());assert b['status']=='PASS' and b['genuine_marked_2kernels']==23 and b['all_paired_fibre_degrees_at_most']==4
r=json.loads((C/'F5_common_basis58_independent.json').read_text());assert r['status']=='PASS' and r['test_count']==100 and r['independent_coefficient_rank']==38 and r['all_h0_zero']
print('D31B_C66_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' finite250=true classes321210=true signedRR100=true genuine2kernels23=true paired_fibre_degree4=true')
