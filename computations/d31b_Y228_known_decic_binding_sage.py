from sage.all import *
from pathlib import Path
import ast,json
C=Path(__file__).resolve().parents[1]/'continuation_20261031b/Y228_shared';s=(C/'union_cases.m').read_text();fd=ast.literal_eval(s.split('UnionFields:=')[1].split(';')[0].replace('[*','[').replace('*]',']'));cases=ast.literal_eval(s.split('UnionCases:=')[1].split(';')[0].replace('[*','[').replace('*]',']'));p=next(z for z in fd if z[0]==10)[1];F=GF(3);R=PolynomialRing(F,'X');K=GF(3**10,name='z',modulus=R(p));z=K.gen();roots=R([2,2,-1,-2,-1,1]).change_ring(K).roots(multiplicities=False);assert len(roots)==5;a=roots[0];b=sum(K(v)*a**i for i,v in enumerate([1032,-239,-222,-759,419]));c=sum(K(v)*a**i for i,v in enumerate([42,-11,-8,-31,17]));E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0);assert P.order()==28;T=14*P;ai=E.a_invariants();S=PolynomialRing(K,'t');t=S.gen();cub=4*t**3+(ai[0]**2+4*ai[1])*t**2+2*(ai[0]*ai[2]+2*ai[3])*t+(ai[2]**2+4*ai[4]);rs=[r for r in cub.roots(multiplicities=False)if r!=T[0]];assert len(rs)==2;P2=E(rs[0],-(ai[0]*rs[0]+ai[2])/2)
def bc(U,Q):
 r,s=U.xy();x,y=Q.xy();X=x-r;Y=y-s+ai[0]*X/2;return(X**3/Y**2,(3*r+ai[1]+ai[0]**2/4)*X**2/Y**2)
def key(vals):return min(tuple(int(w.polynomial()[i])for w in[v**(3**j)for v in vals]for i in range(10))for j in range(10))
want=tuple(cases[26][2][0][2]);matches=[]
for B in[0,1]:
 for Cc in[0,1]:
  for u in[1,3,5,9,11,13]:
   vals=bc(P2+B*T,Cc*P2+u*P)
   if key(vals)==want:matches.append([B,Cc,u])
assert matches;vals=bc(P2,P);assert lcm(v.minpoly().degree()for v in vals)==10
out={'status':'PASS','known_quintic28_good_reduction_at3':True,'full2_extension_residue_degree':10,'global_trace_reduces_to_union_case':27,'actual_automorphisms_matching_canonical_case':matches,'canonical_F3_orbit_key_matches':True};(C/'known_decic_binding_independent.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_Y228_KNOWN_DECIC_BINDING_INDEPENDENT_COMPLETE',matches,flush=True)
