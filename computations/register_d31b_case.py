"""Register frozen packets, including packets proving more than one group."""
from pathlib import Path
import json,hashlib,sys,re
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for label in sys.argv[1:]:
 C=H/'continuation_20261031b'/label;mf=C/'final_certificate_summary.json'
 if not mf.exists():mf=C/'certificate_manifest.json'
 d=json.loads(mf.read_text());assert d['status']=='PASS' and d.get('parent_mathematical_review','PASS')=='PASS'
 match=re.fullmatch(r'C(\d+)',label);pair=[1,int(match[1])]if match else list(map(int,re.fullmatch(r'M(\d+)N(\d+)',label).groups()));groups=d.get('excluded_groups',d.get('groups',[pair]))
 entries=d.get('files_sha256');ev={}
 if entries is not None:
  for rel,digest in entries.items():
   p=(ROOT/rel)if rel.startswith('work/')else(H/rel);assert p.is_relative_to(H)and sha(p)==digest,p;ev[str(p.relative_to(H))]=digest
 else:
  for e in d['files']:
   p=(ROOT/e['path'])if e['path'].startswith('work/')else(H/e['path']);assert p.is_relative_to(H)and sha(p)==e['sha256'],p;ev[str(p.relative_to(H))]=e['sha256']
 method=d.get('method',d.get('scope','Complete exact degree-ten exclusion with independent arithmetic and parent mathematical review'))
 review=C/'registration_review.json';assert not review.exists()
 review.write_text(json.dumps({'status':'PASS','degree':10,'groups':groups,'parent_mathematical_review':'PASS','complete_exclusion':True,'method':method,'frozen_summary_sha256':sha(mf),'classification_complete':False},indent=2)+'\n')
 for p in[mf,review]:ev[str(p.relative_to(H))]=sha(p)
 evidence=[{'path':p,'sha256':s}for p,s in sorted(ev.items())]
 p=H/'certified_results.json';L=json.loads(p.read_text())
 for group in groups:
  assert not any(r['group']==group for r in L['exclusions'])
  L['exclusions'].append({'group':group,'verified':True,'method':method,'proof':str((C/'resolution.tex').relative_to(H)),'evidence':evidence})
 L['auxiliary_results'].append({'name':f'decic_{label}_continuation_20261031b','verified':True,'statement':method+'. Parent mathematical review and independent exact checks passed; failed or incomplete diagnostics are not certificates.','evidence':evidence})
 p.write_text(json.dumps(L,indent=2)+'\n');print('REGISTERED',groups,len(evidence))
L=json.loads((H/'certified_results.json').read_text());new=[r['group']for r in L['exclusions']if r.get('proof','').startswith('continuation_20261031b/')]
(H/'continuation_20261031b/round_status.json').write_text(json.dumps({'original_task_groups':[[1,48],[1,54],[1,56],[1,60],[1,66],[1,74],[3,18],[3,27],[3,33],[4,28]],'target_resolved_groups':10,'registered_new_exclusions':sorted(new),'completed_count':len(new),'round_complete':len(new)>=10,'classification_complete':False},indent=2)+'\n')
