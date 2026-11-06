from sage.all import *
from pathlib import Path
import json,time
N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');exec(compile((N/'framed_prefix.snapshot.py').read_text(),str(N/'framed_prefix.snapshot.py'),'exec'));N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');START=time.monotonic();rows=json.loads((N/'short_cusp_divisors.json').read_text())['data'];A1=vector(ZZ,rows[0]['actual_degree5_cusp_divisor_coefficients']);A2=vector(ZZ,rows[1]['actual_degree5_cusp_divisor_coefficients']);tau=rd('tau_cusp_permutation_p7.m')[1];assert tau==rd('tau_cusp_permutation_p13.m')[1];AUT=rd('cusp_automorphisms.m');gp=AUT[2]
def mat(par):a,b,c,d=par;return matrix(ZZ,[[a,c],[6*b,d]])
def norm(M):
 a=int(M[0,0]%3);b=int((M[1,0]%18)//6);c=int(M[0,1]%3);d=int(M[1,1]%18)
 if d not in[1,5,7]:a=(-a)%3;b=(-b)%3;c=(-c)%3;d=(-d)%18
 assert d in[1,5,7]and a==inverse_mod(d,3)
 return[a,b,c,d]
D={tuple(range(1,25)):[1,0,0,1]};q=list(D)
while q:
 pp=q.pop(0)
 for gg,gpar in zip(gp,[[1,1,0,1],[1,0,1,1],[2,0,0,5]]):
  vv=tuple(gg[i-1]for i in pp);par=norm(mat(D[pp])*mat(gpar))
  if vv in D:assert D[vv]==par
  else:D[vv]=par;q.append(vv)
assert len(D)==27
matches=[]
for h,pp in enumerate(per,1):
 v7=sum((A1[i]*vector(ZZ,c7[tau[i]-1])-A2[i]*vector(ZZ,c7[pp[i]-1])for i in range(24)),vector(ZZ,[0]*n7));v13=sum((A1[i]*vector(ZZ,c13[tau[i]-1])-A2[i]*vector(ZZ,c13[pp[i]-1])for i in range(24)),vector(ZZ,[0]*n13));z7=kr7(list(v7))==kr7([0]*n7);z13=key13(list(v13))==key13([0]*n13);assert z7==z13
 if z7:matches.append({'Aut_index':h,'Aut_parameters':D[tuple(pp)],'cusp_permutation':pp,'difference_class7':list(map(int,v7)),'difference_class13':list(map(int,v13)),'both_finite_class_relations_zero':True})
assert matches
rep={'status':'PASS','statement':'tau(A1) islinearlyequivalent to g(A2)globally: bothcuspidal, ManinDrinfeldtorsion, injectivesplitodd7reduction; exactrelationat13alsoverified','matches':matches,'tau_cusp_permutation':tau,'source_A1':list(map(int,A1)),'source_A2':list(map(int,A2)),'global_second_section_follows_from_first':True};(N/'tau_cusp_class_relation.json').write_text(json.dumps(rep,indent=2)+'\n');print('D31B_M318_TAU_CUSP_CLASS_RELATION_COMPLETE',[(m['Aut_index'],m['Aut_parameters'])for m in matches],time.monotonic()-START,flush=True)
