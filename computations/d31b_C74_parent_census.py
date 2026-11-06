from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time
C=Path('work/deg10/continuation_20261031b/C74');start=time.monotonic();Fp=GF(2);K=GF(1024,'a');tau=next(z for z in K if z.trace()==1);set_random_seed(313174)
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def xy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return[(s-r)/t,t/(s*s-s-r+1)]
def key2(vals):
 L=vals[0].parent();d=lcm([v.minpoly().degree()for v in vals]);ok=False
 for bits in[[1,0],[0,1],[1,1]]:
  g=bits[0]*vals[0]+bits[1]*vals[1];mp=g.minpoly()
  if mp.degree()==d:chosen=bits;ok=True;break
 assert ok
 vv=lambda v:vector(Fp,[L(v).polynomial()[i]for i in range(L.degree())]);M=matrix(Fp,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),chosen,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(v))))for v in vals]]
# Every supersingular trace over F_(2^10) is in this containing list.
ss_traces=[-64,-32,0,32,64];assert all((1025-t)%37 for t in ss_traces)
keys=Counter();orders=Counter();ncurves=0
for j in K:
 if j==0:continue
 for a2 in[0,tau]:
  E=EllipticCurve(K,[1,a2,0,0,1/j]);assert E.j_invariant()==j;n=E.cardinality();ncurves+=1
  if n%37:continue
  assert n in[962,1036]and (1025-n)%2 and 1023%37!=0;orders[int(n)]+=1
  while True:
   P=(n//37)*E.random_point()
   if P!=E(0):break
  assert 37*P==E(0)
  for u in range(1,19):
   kk=key2(xy(E,u*P));assert kk[0]==10;keys[json.dumps(kk)]+=1
assert ncurves==2046 and orders==Counter({962:5,1036:50})and len(keys)==99 and set(keys.values())=={10}
native=ast.literal_eval((C/'parent_native_point_keys.m').read_text().replace('[*','[').replace('*]',']'));assert len(native)==99 and set(map(json.dumps,native))==set(keys)
out={'status':'PASS','ordinary_isomorphism_classes_enumerated':ncurves,'supersingular_j0_excluded_by_complete_possible_trace_list':ss_traces,'eligible_orders':dict(orders),'closed_marked_points':99,'every_closed_point_degree10':True,'all99_native_prime_keys_match_independent_marked_curves':True,'geometric_marked_points':990,'seconds':time.monotonic()-start};(C/'parent_independent_census.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C74_PARENT_CENSUS_COMPLETE',99,dict(orders),time.monotonic()-start,flush=True)
