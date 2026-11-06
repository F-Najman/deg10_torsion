#!/usr/bin/env python3
"""Read-only integrity and scope audit for the isolated degree-ten project."""
from pathlib import Path
import json,hashlib,math
H=Path(__file__).resolve().parent;ROOT=H.parent.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
L=json.loads((H/'certified_results.json').read_text());D=json.loads((H/'candidates.json').read_text());evidence=set()
for r in L['realizations']+L['exclusions']+L['auxiliary_results']:
 assert r['verified'] is True
 for e in r['evidence']:
  assert sha(H/e['path'])==e['sha256'],e['path'];evidence.add(e['path'])
for f,base in [('inputs/global_archive_manifest.json','inputs/global_archive'),('inputs/C115_archive_manifest.json','inputs/C115_archive'),('inputs/M3N15_octic_archive_manifest.json','inputs/M3N15_octic_archive'),('shared/mdsage_manifest.json','shared/mdsage')]:
 m=json.loads((H/f).read_text())
 for p,s in m['files'].items():assert sha(H/base/p)==s,(f,p)
m=json.loads((H/'shared/mdmagma_manifest.json').read_text())
for e in m['files']:assert sha(H/e['path'])==e['sha256'],e['path']
for e in json.loads((H/'inputs/provenance.json').read_text()):assert sha(H/e['snapshot'])==e['sha256'],e['snapshot']
old=json.loads((H/'earlier_degrees_manifest.json').read_text())
for e in old['files']:assert sha(ROOT/e['path'])==e['sha256'],e['path']
known=set(map(tuple,L['inherited_occurrences']['groups']))|{tuple(r['group']) for r in L['realizations']};excluded={tuple(r['group']) for r in L['exclusions']}
assert sorted(known)==D['known_lower_bound'] or sorted(known)==list(map(tuple,D['known_lower_bound']))
assert sorted(excluded)==list(map(tuple,D['certified_exclusions']))
sub=lambda a,b:b[0]%a[0]==0 and b[1]%a[1]==0
assert not any(sub(e,g) for e in excluded for g in known)
frontier=set(map(tuple,D['minimal_outside_known']));covered={g for g in frontier if any(sub(e,g) for e in excluded)}
assert sorted(frontier-covered)==list(map(tuple,D['remaining_minimal_obstructions']))
complete=bool(L['classification_complete'])
assert len(known)==(69 if complete else 68) and len(frontier)==(110 if complete else 109)
assert len(excluded)==len(covered) and len(covered)+len(frontier-covered)==len(frontier)
complete=not(frontier-covered)
assert L['classification_complete'] is D['classification_complete'] is complete
assert D['classification']==([list(g) for g in sorted(known)] if complete else None)
layout=H/'paper/layout/provenance.json'
if layout.exists():
 for e in json.loads(layout.read_text()).values():
  assert sha(H/e['source'])==e['source_sha256']
  assert sha(H/e['typeset_copy'])==e['typeset_sha256']
bad=[];pending=[];statuses=[]
for p in sorted((H/'jobs').glob('*.status.json')):
 s=json.loads(p.read_text());statuses.append(s)
 if s['state']!='EXITED':pending.append(p.name)
 if 'script_sha256' in s:assert sha(H/s['script'])==s['script_sha256'],s['script']
 if s.get('actual_exit_code')!=0:bad.append({'job':s['purpose'],'actual_exit_code':s.get('actual_exit_code'),'reason':s.get('termination_reason')})
assert not pending,pending
log=(H/'paper/main.log').read_text();assert 'Undefined control sequence' not in log and 'LaTeX Error' not in log and 'Overfull' not in log
report={'status':'PASS; CLASSIFICATION COMPLETE' if complete else 'PASS; CLASSIFICATION INCOMPLETE','known_full_groups':len(known),'certified_excluded_subgroups':len(excluded),'minimal_obstructions':len(frontier),'remaining_minimal_obstructions':len(frontier-covered),'evidence_files_checked':len(evidence),'protected_earlier_files_unchanged':len(old['files']),'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'unsuccessful_or_terminated_probes_not_promoted':bad,'paper_compiles_without_errors_or_overfull_boxes':True}
print(json.dumps(report,indent=2))
