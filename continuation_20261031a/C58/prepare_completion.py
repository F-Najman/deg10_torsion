from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[3];W=R/'work/deg10'
required={'d31a_C58_class_filter':'D31A_DEGREE10_CLASS_FILTER_COMPLETE level=29 points=5971 tests=293 survivors=2',
 'd31a_C58_known_lifts_v3':'D31A_C58_KNOWN_LIFT_BINDING_COMPLETE',
 'd31a_C58_independent_replay':'D31A_C58_INDEPENDENT_REPLAY_COMPLETE',
 'd31a_C58_inherited_replay':'NONIC58_INDEPENDENT_REPLAY',
 'd31a_58_parent_class_replay_v2':'C58_PARENT_CLASS_REPLAY_COMPLETE 293 [ 462, 1748 ]'}
accept=set();diagnostics=set()
def add(p):assert p.is_file();accept.add(p)
for p in (H/'inputs').rglob('*'):
 if p.is_file():add(p)
for s in ['resolution.tex','class_common.m','degree10_filter.m','degree10_census.m','degree10_classes.m','degree10_complete.m',
 'independent_decic_replay.py','independent_decic_replay.json','known_lift_binding_job_v3.m','known_lift_binding.m','parent_class_replay.m',
 'inherited_provenance.json','completion_summary.json','prepare_completion.py']:add(H/s)
statuses=[]
for stem,marker in required.items():
 suffix='.py'if 'independent_replay'in stem or 'inherited_replay'in stem else'.m'
 add(W/'computations'/(stem+suffix));add(W/'jobs'/(stem+'.log'));add(W/'jobs'/(stem+'.status.json'))
 a=json.loads((W/'jobs'/(stem+'.status.json')).read_text());assert a['state']=='EXITED'and a['actual_exit_code']==0
 assert marker in (W/'jobs'/(stem+'.log')).read_text()
for p in list((W/'jobs').glob('d31a_C58*.status.json'))+list((W/'jobs').glob('d31a_58_parent*.status.json')):
 a=json.loads(p.read_text());assert a['state']=='EXITED';statuses.append(dict(path=str(p.relative_to(R)),actual_exit_code=a['actual_exit_code']))
life=dict(status='PASS',all_case_jobs_exited=True,accepted_producers=list(required),all_case_statuses=statuses,failed_preflights_are_not_certificates=True)
(H/'lifecycle.json').write_text(json.dumps(life,indent=2)+chr(10));add(H/'lifecycle.json')
for p in H.rglob('*'):
 if p.is_file()and p.name!='certificate_manifest.json'and p not in accept:diagnostics.add(p)
for directory in ['computations','jobs']:
 for p in list((W/directory).glob('d31a_C58*'))+list((W/directory).glob('d31a_58_parent*')):
  if p.is_file()and p not in accept:diagnostics.add(p)
def ent(p):return dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
a=dict(status='PASS',scope='Complete C58 decic exclusion via proven genus22 Picard image, exact new5971-point census,293 primary tests, both global decic rigidity lifts, independent Sage and native ideal controls; old2324 tests replayed with new complement argument',files=[ent(p)for p in sorted(accept)],diagnostic_files=[ent(p)for p in sorted(diagnostics)])
(H/'certificate_manifest.json').write_text(json.dumps(a,indent=2)+chr(10))
for e in a['files']+a['diagnostic_files']:assert hashlib.sha256((R/e['path']).read_bytes()).hexdigest()==e['sha256']
print('C58_FROZEN',len(a['files']),len(a['diagnostic_files']),hashlib.sha256((H/'certificate_manifest.json').read_bytes()).hexdigest())
