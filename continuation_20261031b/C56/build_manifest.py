from pathlib import Path
import json,hashlib,argparse,re
P=argparse.ArgumentParser();P.add_argument('--report',action='append',default=[]);A=P.parse_args()
C=Path(__file__).resolve().parent;H=C.parents[1];R=C.parent;S=R/'Y228_shared'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
review=json.loads((C/'mathematical_review.json').read_text());assert review['status']=='PASS'
for n in A.report:assert json.loads((H/n).read_text())['status']=='PASS',n
files=set()
def add(p):
 assert p.is_file(),p
 files.add(p)
def tree(p):
 for q in p.rglob('*'):
  if q.is_file() and q.suffix not in ['.pdf','.aux','.toc','.out','.synctex'] and q.name!='final_certificate_summary.json':add(q)
# C56 primary, exact inherited snapshots, and independent controls.
tree(C);tree(H/'continuation_20261031a/C56')
# Stable accepted shared quotient model, ranks and genuine marked maps.
old=R/'M4N28/final_certificate_summary.json';od=json.loads(old.read_text());assert od['status']=='PASS';add(old)
for n,h in od['files_sha256'].items():assert sha(H/n)==h,n;add(H/n)
k=S/'known_decic_C56_obstruction_packet.json';kd=json.loads(k.read_text());assert kd['status']=='PASS' and kd['parent_review']=='PASS';add(k)
for n,h in kd['files_sha256'].items():assert sha(H/n)==h,n;add(H/n)
for n in ['fast_inverse_verified.json','fast_primitive_inverse.m']:add(S/n)
for pat in ['d31b_C56*','d31a_C56*']:
 for q in (H/'computations').glob(pat):
  if q.is_file():add(q)
add(H/'run_job.py')
primary={
 'd31b_C56_project_to_Y28','d31b_C56_final_source_subset','d31b_C56_quadratic_cover_profile',
 'd31b_C56_base_operator_s0','d31b_C56_base_operator_s1','d31b_C56_independent_class_ideals',
 'd31b_C56_kummer_direct_case491_v3','d31b_C56_kummer_cross271','d31b_C56_kummer_cross327','d31b_C56_kummer_cross65',
 *[f'd31b_C56_kummer_case{i}'for i in[578,818,837,978,1009]],
 'd31a_C56_fourmap_filter_v4','d31a_C56_independent_fourmap_v3'}
jobs=[]
for pat in ['d31b_C56*.status.json','d31a_C56*.status.json']:
 for q in sorted((H/'jobs').glob(pat)):
  s=json.loads(q.read_text());assert s['state']=='EXITED',q
  add(q);add(H/s['log']);add(H/s['script']);assert sha(H/s['script'])==s['script_sha256'],q
  stem=q.name.removesuffix('.status.json');text=(H/s['log']).read_text();markers=[z for z in text.splitlines()if re.match(r'^D31[AB]_C56_[A-Z0-9_]*COMPLETE(?:\s|$)',z)]
  accepted=s['actual_exit_code']==0 and markers and (stem in primary or 'sieve_'in stem or 'parent_basis_v3_'in stem)
  jobs.append({'status_path':str(q.relative_to(H)),'actual_exit_code':s['actual_exit_code'],'accepted_completion_marker':markers[-1].split()[0]if accepted else None,'scope':'completed final proof component'if accepted else'diagnostic, auxiliary exploration, or superseded computation; not a certificate'})
for n in A.report:add(H/n)
m={'status':'PASS','degree':10,'excluded_groups':[[1,56]],'parent_mathematical_review':'PASS','method':'Four X1(28) rational Picard images and genus25 Hecke sieve; exact quadratic X1(56) Kummer restrictions on nine divisors; known decic quotient and four nonsquare 2-descent classes','files_sha256':{str(q.relative_to(H)):sha(q)for q in sorted(files)},'jobs':jobs,'independent_report_paths':A.report,'final_source_case_partition':{'global_known_quotient_nonsquare':[492,737],'A5_exact_RR_zero':[271,327,491,578,818,837,978,1009],'A11_exact_RR_zero':[65]},'historical_stage_status':'The frozen 31a C56 notes were partial; the present proof supplies its remaining exclusions. Earlier incomplete or failed native jobs are retained only as diagnostics. All final RR tests have separate completed Kummer jobs.','diagnostics_are_not_certificates':True}
f=C/'final_certificate_summary.json';assert not f.exists();f.write_text(json.dumps(m,indent=2)+'\n');print('C56_MANIFEST_FROZEN',len(files),sha(f))
