from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;C=H/'continuation_20261031a/C70';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();mf=C/'certificate_manifest.json';digest=sha(mf)
assert digest=='eadb712f7b96415986835801baefecbaea589cebd80ab9f1229917cb768a9e59';d=json.loads(mf.read_text());assert d['status']=='PASS'
ev={}
for e in d['files']:
 p=ROOT/e['path'];assert p.is_relative_to(H)and sha(p)==e['sha256'];ev[str(p.relative_to(H))]=e['sha256']
review=C/'registration_review_v2.json';assert not review.exists();review.write_text(json.dumps({'status':'PASS','degree':10,'group':[1,70],'parent_mathematical_review':'PASS','complete_exclusion':True,'frozen_summary_sha256':digest,'replaces':'registration_review.json','replacement_reason':'The original registration occurred during manifest finalization. Only the manifest hash changed; the mathematical proof and all previously accepted arithmetic files are unchanged. This review binds the explicitly confirmed final freeze, including its additional record.','classification_complete':False},indent=2)+'\n')
for p in [mf,review]:ev[str(p.relative_to(H))]=sha(p)
evidence=[{'path':p,'sha256':s}for p,s in sorted(ev.items())]
p=H/'certified_results.json';L=json.loads(p.read_text())
for r in L['exclusions']:
 if r['group']==[1,70]:r['evidence']=evidence
for r in L['auxiliary_results']:
 if r['name']=='decic_C70_continuation_20261031a':r['evidence']=evidence
p.write_text(json.dumps(L,indent=2)+'\n');print('C70_FINAL_HANDOFF_REPINNED',len(evidence))
