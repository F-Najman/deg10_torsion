from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[3];W=R/'work/deg10';accept=set();diag=set()
def add(p):assert p.is_file();accept.add(p)
for p in (H/'inputs').rglob('*'):
 if p.is_file()and'__pycache__'not in p.parts:add(p)
for s in ['resolution.tex','new_finite_jets.m','new_finite_jets.raw.jsonseq','new_finite_jets.normalized.jsonl','prepare_new_sieve.py',
 'new_sieve_inputs.m','new_matrix_sieve.m','new_matrix_results.raw.jsonseq','q7_residual_jets.m','q7_finite_jets.raw.jsonseq',
 'q7_finite_jets.normalized.jsonl','q7_sieve_inputs.m','q7_matrix_sieve.m','q7_matrix_results.raw.jsonseq',
 'cusp_audit.py','cusp_audit.json','degree10_frobenius_audit.py','degree10_frobenius_audit.json','completion_summary.json','prepare_completion.py']:add(H/s)
required={'d31a_C44_new_finite_jets':'D31A_C44_NEW_FINITE_JETS_COMPLETE cases=956 distinct10=5854 distinct5=60',
 'd31a_C44_new_matrix_sieve':'D31A_C44_CANONICAL22_COMPLETE cases=956 survivors=[ 510, 511 ] all_rank105=false',
 'd31a_C44_q7_residual_jets':'D31A_C44_Q7_FINITE_JETS_COMPLETE cases=2',
 'd31a_C44_q7_matrix_sieve':'D31A_C44_Q7_CANONICAL22_COMPLETE cases=2 survivors=[] all_rank105=true'}
for stem,marker in required.items():
 for p in [W/'computations'/(stem+'.m'),W/'jobs'/(stem+'.log'),W/'jobs'/(stem+'.status.json')]:add(p)
 a=json.loads((W/'jobs'/(stem+'.status.json')).read_text());assert a['state']=='EXITED'and a['actual_exit_code']==0
 assert marker in(W/'jobs'/(stem+'.log')).read_text()
nested=H/'independent_certificate_manifest.json';x=json.loads(nested.read_text());assert x['status']=='PASS'and x['all_owned_workers_exited'];add(nested)
base=R/x['path_base']
for q,sha in x['accepted_files_sha256'].items():
 p=base/q;assert hashlib.sha256(p.read_bytes()).hexdigest()==sha;add(p)
for q,sha in x['diagnostic_files_sha256'].items():
 p=base/q;assert hashlib.sha256(p.read_bytes()).hexdigest()==sha;diag.add(p)
assert json.loads((H/'independent_final_verification.json').read_text())['status']=='PASS'
statuses=[]
for p in (W/'jobs').glob('d31a_C44*.status.json'):
 a=json.loads(p.read_text());assert a['state']=='EXITED';statuses.append(dict(path=str(p.relative_to(R)),actual_exit_code=a['actual_exit_code']))
life=dict(status='PASS',all_case_jobs_exited=True,all_case_statuses=statuses,required_primary_producers=list(required),
 independent_accepted_jobs=x['jobs'],retired_or_failed_jobs_are_not_certificates=True,unlaunched_prefix_search_helpers_are_diagnostics=True)
(H/'lifecycle.json').write_text(json.dumps(life,indent=2)+chr(10));add(H/'lifecycle.json')
for p in H.rglob('*'):
 if p.is_file()and p.name!='certificate_manifest.json'and p not in accept:diag.add(p)
for directory in['computations','jobs']:
 for p in(W/directory).glob('d31a_C44*'):
  if p.is_file()and p not in accept:diag.add(p)
diag-=accept
# Verify every retained source snapshot once more, preserving the old proof packets.
for name in['inputs/inherited_canonical_provenance.json','inputs/lower_snapshot_provenance.json']:
 for a in json.loads((H/name).read_text())['files']:
  assert hashlib.sha256((R/a['snapshot']).read_bytes()).hexdigest()==a['sha256']
  assert hashlib.sha256((R/a['source']).read_bytes()).hexdigest()==a['sha256']
def ent(p):return dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
a=dict(status='PASS',scope='Complete decic C44 exclusion: Varivoda gonality boundary and exact cusp Gal5=Gal3^8; all956 new support tests combined with inherited lower tests; full canonical/jet/genuine-kernel controls',files=[ent(p)for p in sorted(accept)],diagnostic_files=[ent(p)for p in sorted(diag)])
(H/'certificate_manifest.json').write_text(json.dumps(a,indent=2)+chr(10))
for z in a['files']+a['diagnostic_files']:assert hashlib.sha256((R/z['path']).read_bytes()).hexdigest()==z['sha256']
print('C44_FROZEN',len(a['files']),len(a['diagnostic_files']),hashlib.sha256((H/'certificate_manifest.json').read_bytes()).hexdigest())
