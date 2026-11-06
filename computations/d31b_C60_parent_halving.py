from sage.all import *
from pathlib import Path
import json,ast,time,re
C=Path('work/deg10/continuation_20261031b/C60');start=time.monotonic();Fp=GF(7);RP=PolynomialRing(Fp,'x')
def read(n):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',(C/n).read_text().replace('[*','[').replace('*]',']'))))
points,left=read('three_maps_remaining.m');oldobs,oldtests,oldleft=read('final_halving.m');assert len(left)==24 and len(oldobs)==72 and oldleft==[];obs={int(r[0]):r for r in oldobs};result={};rows=[]
for pid in sorted(obs):
 desc=points[pid-1];e,bits,pol,co=desc;assert len(pol)==e+1;K=GF(7**e,'a',modulus=RP(pol));a=K.gen()
 if len(co)==2 and isinstance(co[0],list):bb=RP(co[0])(a);cc=RP(co[1])(a)
 else:bb=RP(co[:e])(a);cc=RP(co[e:])(a)
 E=EllipticCurve(K,[0,cc,0,bb*(1-bb-cc),0]);P2=E(0,0);Q=E(bb,bb);assert 2*P2==E(0)and P2!=E(0)and P2!=15*Q and 30*Q==E(0)and 15*Q!=E(0)and 10*Q!=E(0)and 6*Q!=E(0)
 phi=E.isogeny(P2);A=phi.codomain();S=phi(Q);dual=phi.dual();assert dual(S)==2*Q;assert(Q[0]+(Q+P2)[0]-P2[0],Q[1]+(Q+P2)[1]-P2[1])==S.xy()
 R=PolynomialRing(K,'z');z=R.gen();b2,b4,b6,b8=A.b_invariants();f=z**4-b4*z*z-2*b6*z-b8-S[0]*(4*z**3+b2*z*z+2*b4*z+b6);fac=sorted(int(g.degree())for g,k in f.factor()for _ in range(k));halves=[T for sx in f.roots(multiplicities=False)for T in A.lift_x(sx,all=True)if 2*T==S];marked=[T for T in halves if dual(T)==Q];n=len(marked);assert n in[0,2]
 assert all(60*T==A(0)and 30*T!=A(0)and phi(dual(T))==2*T for T in marked)
 assert e==obs[pid][1]and fac==sorted(obs[pid][2])and n==obs[pid][3];result[pid]=n;rows.append({'pointID':pid,'degree':int(e),'half_quartic_factor_degrees':fac,'marked_lifts':n})
assert all(cusp==0 and all(mul==1 for pid,mul in parts)and any(result[int(pid)]==0 for pid,mul in parts)for cid,parts,cusp,flags in left)
out={'status':'PASS','active_source_points':72,'patterns_rejected':24,'original_marked_halving_lifts_verified':True,'manual_order2_Velu_coordinates_verified':True,'dual_composition_and_orders_verified':True,'all_multiplicity_one_obstructions_verified':True,'points':rows,'seconds':time.monotonic()-start};(C/'parent_independent_halving.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C60_PARENT_HALVING_COMPLETE',72,24,time.monotonic()-start,flush=True)
