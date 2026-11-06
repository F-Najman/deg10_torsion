from pathlib import Path
import json,hashlib,ast
C=Path(__file__).resolve().parent;B=C.parents[1];m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,74]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
for row in m['jobs']:
 s=json.loads((B/row['status_path']).read_text());assert s['state']=='EXITED' and s['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker']in(B/s['log']).read_text()
r=json.loads((C/'parent_prefix_ranks_verified.json').read_text());assert r['status']=='PASS' and r['rational_cusp_patterns']==8436285 and r['sextic_plus_cusp_patterns']==5985 and r['rank']==10
for n,h in r['source_hashes'].items():assert hashlib.sha256((C/n).read_bytes()).hexdigest()==h,n
c=json.loads((C/'parent_independent_census.json').read_text());assert c['status']=='PASS' and c['closed_marked_points']==99 and c['all99_native_prime_keys_match_independent_marked_curves']
n=ast.literal_eval((C/'parent_native_summary_v2.m').read_text().replace('[*','[').replace('*]',']'));assert n==[99,18,10,38]
assert json.loads((B/'continuation_20261031b/C66/final_certificate_summary.json').read_text())['parent_mathematical_review']=='PASS'
print('D31B_C74_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' cusp8436285=true mixed5985=true nativeNC10_99=true independent2046curves=true')
