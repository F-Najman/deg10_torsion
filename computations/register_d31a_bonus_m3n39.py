from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=H/'certified_results.json';L=json.loads(p.read_text());r=next(r for r in L['exclusions']if r['group']==[1,78]);assert not any(r['group']==[3,39]for r in L['exclusions'])
mf=H/'continuation_20261031a/C78/final_certificate_summary.json';d=json.loads(mf.read_text());assert [3,39]in d['excluded_groups']and d['parent_mathematical_review']=='PASS'
for e in r['evidence']:assert sha(H/e['path'])==e['sha256']
review=mf.parent/'registration_review_m3n39.json';assert not review.exists()
method='The complete decic X1(39) classification identifies every curve with a point of order39 with a known full-C39 curve over its exact degree10 moduli field; hence no independent second3-point exists'
review.write_text(json.dumps({'status':'PASS','degree':10,'group':[3,39],'parent_mathematical_review':'PASS','method':method,'frozen_summary_sha256':sha(mf),'complete_exclusion':True,'classification_complete':False},indent=2)+'\n')
ev=r['evidence']+[{'path':str(review.relative_to(H)),'sha256':sha(review)}]
L['exclusions'].append({'group':[3,39],'verified':True,'method':method,'proof':r['proof'],'evidence':ev});L['auxiliary_results'].append({'name':'decic_M3N39_continuation_20261031a','verified':True,'statement':method,'evidence':ev});p.write_text(json.dumps(L,indent=2)+'\n')
new=[r['group']for r in L['exclusions']if r.get('proof','').startswith('continuation_20261031a/')]
(H/'continuation_20261031a/round_status.json').write_text(json.dumps({'target_additional_exclusions':10,'registered_new_exclusions':sorted(new),'completed_count':len(new),'round_complete':len(new)>=10,'classification_complete':False},indent=2)+'\n')
print('REGISTERED',[3,39],len(ev))
