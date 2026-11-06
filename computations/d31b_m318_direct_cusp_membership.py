from sage.all import *
from pathlib import Path
import json,time
N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');exec(compile((N/'framed_prefix.snapshot.py').read_text(),str(N/'framed_prefix.snapshot.py'),'exec'));N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');START=time.monotonic();sur=json.loads((P/'paired_final_survivors.json').read_text())['survivors'];H13,_=hnf(r13+c13);out=[];short=json.loads((N/'short_cusp_divisors.json').read_text())['data'];short={v['target_id']:v for v in short}
for rec in sur:
 for sid,h in rec['source_matches']:
  if sid not in[3514,3561]:continue
  desc,sv=source[sid-1];co=vector(ZZ,sv)*HC.inverse();assert all(z in ZZ for z in co);cc=list((vector(ZZ,co)*UHC)[len(r7):]);assert kr7(list(vector(ZZ,cc)*matrix(ZZ,c7)-vector(ZZ,sv)))==kr7([0]*n7);perm=per[h-1];hv=sum((cc[i]*vector(ZZ,c7[perm[i]-1])for i in range(24)),vector(ZZ,[0]*n7))+(5-sum(ds[i]*cc[i]for i in range(24)))*vector(ZZ,c7[perm[base]-1]);tv=vector(ZZ,rec['target_class']);tmem=all(v in ZZ for v in tv*H13.inverse());row={'target_id':rec['target_id'],'target_descriptor':rec['target_descriptor'],'source_pattern':sid,'source_descriptor':desc,'normalizing_Aut':h,'source_class_in_global_cusp_image':True,'source_cusp_coefficients':list(map(int,cc)),'normalized_source_class':list(map(int,hv)),'target_class':list(map(int,tv)),'target_class_in_global_cusp_image':bool(tmem)}
  if tmem:
   rep=short[rec['target_id']];vv=vector(ZZ,rep['actual_degree5_cusp_divisor_coefficients']);assert kr7(list(vv*matrix(ZZ,c7)-hv))==kr7([0]*n7);assert key13(list(vv*matrix(ZZ,c13)-tv))==key13([0]*n13);row['global_degree5_cusp_divisor_coefficients']=list(map(int,vv));row['actual_unmultiplied_classes_match_both_primes']=True
  out.append(row);print('M318_DIRECT_CUSP_MEMBERSHIP',rec['target_id'],sid,h,tmem,flush=True)
assert len(out)==4 and sum(v['target_class_in_global_cusp_image']for v in out)==2
rep={'status':'PASS','scope':'FourfinalcandidatepairsafterindependentquarticandHeckefilters; globaltorsioninject7 forcescuspidalityfromsourceHmembership; two13classesnotcuspidal areexcluded','tested_pairs':4,'excluded_pairs':2,'remaining_pairs':2,'global_cusp_group_order':20253807,'source_class_membership_uses_unmultiplied_classes':True,'normalization_uses_only_actual24cusppermutations_and_affine5basecorrection':True,'observations':out};(N/'direct_cusp_membership_report.json').write_text(json.dumps(rep,indent=2)+'\n');print('D31B_M318_DIRECT_CUSP_MEMBERSHIP_COMPLETE',4,2,2,time.monotonic()-START,flush=True)
