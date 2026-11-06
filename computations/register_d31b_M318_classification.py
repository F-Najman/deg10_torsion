"""Register the frozen occurrence and complete low-degree X1(3,18) classification."""
from pathlib import Path
import hashlib,json,subprocess
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;C=H/'continuation_20261031b/M3N18'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
mf=C/'final_certificate_summary.json';d=json.loads(mf.read_text())
assert d['status']=='PASS' and d['parent_mathematical_review']=='PASS'
assert d['complete_low_degree_classification'] is True
assert d['realized_groups']==[[3,18]] and sorted(d['excluded_groups'])==[[3,36],[6,18]]
ev={}
if 'files_sha256' in d:items=d['files_sha256'].items()
else:items=((e['path'],e['sha256'])for e in d['files'])
for rel,digest in items:
 p=ROOT/rel if rel.startswith('work/')else H/rel
 assert p.is_relative_to(H) and sha(p)==digest,p
 ev[str(p.relative_to(H))]=digest
review=C/'registration_review.json';assert not review.exists()
paths=[H/'certified_results.json',H/'candidates.json',H/'candidates.txt',H/'continuation_20261031b/round_status.json']
before={p:p.read_bytes() for p in paths}
oldledger=json.loads(before[H/'certified_results.json'])
assert not any(r['group']in[[3,18],[3,36],[6,18]]for r in oldledger['realizations']+oldledger['exclusions'])
try:
 review.write_text(json.dumps({'status':'PASS','degree':10,'parent_mathematical_review':'PASS','realized_groups':[[3,18]],'excluded_groups':[[3,36],[6,18]],'complete_low_degree_classification':True,'frozen_summary_sha256':sha(mf)},indent=2)+'\n')
 for p in[mf,review]:ev[str(p.relative_to(H))]=sha(p)
 evidence=[{'path':p,'sha256':s}for p,s in sorted(ev.items())]
 p=H/'certified_results.json';L=json.loads(p.read_text());proof=str((C/'resolution.tex').relative_to(H))
 assert not any(r['group']==[3,18]for r in L['realizations']+L['exclusions'])
 L['realizations'].append({'group':[3,18],'verified':True,'method':'Exact characteristic-zero degree-five point over Q(zeta3); exact independent marked-point and full-torsion verification','proof':proof,'evidence':evidence})
 for group in[[3,36],[6,18]]:
  assert not any(r['group']==group for r in L['exclusions'])
  L['exclusions'].append({'group':group,'verified':True,'method':'Complete classification of noncuspidal degree-five divisors on X1(3,18) over Q(zeta3); every resulting elliptic curve has exact full torsion C3 x C18','proof':proof,'evidence':evidence})
 L['auxiliary_results'].append({'name':'decic_M3N18_complete_degree5_classification','verified':True,'statement':'All noncuspidal degree-five divisors are accounted for by the two certified cuspidal classes and their marking translates. Full torsion equals C3 x C18.','evidence':evidence})
 p.write_text(json.dumps(L,indent=2)+'\n')
 subprocess.run(['python3',str(H/'candidates.py')],cwd=ROOT,check=True,stdout=subprocess.PIPE)
 D=json.loads((H/'candidates.json').read_text());assert len(D['known_lower_bound'])==69 and len(D['minimal_outside_known'])==110 and len(D['certified_exclusions'])==110 and not D['remaining_minimal_obstructions']
 L['classification_complete']=True;p.write_text(json.dumps(L,indent=2)+'\n')
 original=[[1,48],[1,54],[1,56],[1,60],[1,66],[1,74],[3,18],[3,27],[3,33],[4,28]]
 new=[r['group']for r in L['exclusions']if r.get('proof','').startswith('continuation_20261031b/')]
 assert set(map(tuple,original))<=set(map(tuple,new+[[3,18]]))
 (H/'continuation_20261031b/round_status.json').write_text(json.dumps({'original_task_groups':original,'target_resolved_groups':10,'registered_new_exclusions':sorted(new),'registered_new_occurrences':[[3,18]],'completed_count':10,'round_complete':True,'classification_complete':True,'known_full_groups':69,'frontier_size':110,'excluded':110,'remaining':0,'note':'The last original case occurs. Its complete low-degree classification excludes two newly exposed larger groups.'},indent=2)+'\n')
except BaseException:
 for p,data in before.items():p.write_bytes(data)
 if review.exists():review.unlink()
 raise
print('REGISTERED M318 occurrence and M3N36/M6N18 exclusions;',len(evidence),'evidence hashes; Phi(10) complete')
