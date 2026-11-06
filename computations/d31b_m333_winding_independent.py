from sage.all import *
from pathlib import Path
import ast,json,time
T=time.monotonic();R=Path('work/deg10/continuation_20261031b/M3N33');a=ast.literal_eval((R/'winding_antih10_N99_H34.m').read_text().replace('[*','[').replace('*]',']'));assert a[:4]==[99,34,81,40];out=[]
for index,order,d,p,v,ops in a[4]:
 assert order in[2,10]and d==8;K=GF(p);V=VectorSpace(K,d);W=V.subspace([V(v)]);As=[matrix(K,d,d,x)for x in ops]
 while True:
  old=W.dimension();W=V.subspace(list(W.basis())+[b*A for A in As for b in W.basis()])
  if W.dimension()==old:break
 assert W.dimension()==d;out.append({'character_order':order,'dimension_over_character_field':d,'finite_prime':p,'full_Hecke_orbit_rank':int(W.dimension())})
assert sorted(x['character_order']for x in out)==[2,10];result={'status':'PASS','quotient_source_genus':81,'h10_invariant_quotient_genus':41,'anti_invariant_dimension':40,'components':out,'twist_containment':'Quadratic character of conductor3 preserves nebentype and twists each level-dividing99 factor to level dividing99; Derickx-Sutherland Theorem4.1 proof (Atkin-Li). Hence all h10-minus factors and their -3 twists have analytic rank0.','seconds':time.monotonic()-T};(R/'winding_antih10_independent.json').write_text(json.dumps(result,indent=2)+'\n');print('D31B_M333_WINDING_INDEPENDENT_COMPLETE',40,time.monotonic()-T,flush=True)
