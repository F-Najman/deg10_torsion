"""Independent exhaustive source census and all genuine marked five-isogenies."""
from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time,hashlib
ROOT=Path('work/deg10/continuation_20261030c/C99');START=time.monotonic()
def read(name):return ast.literal_eval((ROOT/name).read_text().replace('[*','[').replace('*]',']'))
F2=GF(2);R2=PolynomialRing(F2,'T')
source=read('source99.m');primary=read('quotient33_RR.m')
raw=ROOT/'independent_inputs/all_curves.m';inputcheck=json.loads((ROOT/'independent_inputs/independent_point_counts.json').read_text())
assert inputcheck['status']=='PASS' and hashlib.sha256(raw.read_bytes()).hexdigest()==inputcheck['input_sha256']
allcurves=ast.literal_eval(raw.read_text().replace('[*','[').replace('*]',']'))
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y
 n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return [-n2**3/b3**2,-(n1*n2-b3)/b3]
def unitxy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def vec(K,v):return vector(F2,[K(v).polynomial()[i] for i in range(K.degree())])
def key2(v):
 K=v[0].parent();d=lcm([a.minpoly().degree() for a in v]);found=False
 for bits in [[1,0],[0,1],[1,1]]:
  g=bits[0]*v[0]+bits[1]*v[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found
 M=matrix(F2,[vec(K,g**i) for i in range(d)]);assert M.rank()==d
 co=[list(map(int,M.transpose().solve_right(vec(K,a)))) for a in v]
 return [int(d),bits,list(map(int,mp.list())),co]
def frobkey(v,d):
 K=v[0].parent()
 return min(tuple(int(a.polynomial()[i]) for a in [v[0]**(2**j),v[1]**(2**j)] for i in range(d)) for j in range(d))

# Independently inspect all relevant curve groups and every rational point.
closed=[];rational=[];census=[];sourceorbits=set();allkeys=set()
for d,pol,curves in allcurves:
 K=GF(2**d,'a',modulus=R2(pol));a=K.gen();levels=0;eligible=0;pointorders=[]
 for jj,aa,card in curves:
  if card%99:continue
  E=EllipticCurve(K,[R2(v)(a) for v in aa]);assert E.cardinality()==card
  pts=[P for P in E.points() if P!=E(0) and 99*P==E(0) and 33*P!=E(0) and 9*P!=E(0)]
  if not pts:continue
  assert d==10 and card==990 and len(E.automorphisms())==2 and len(pts)==60
  eligible+=1;levels+=len(pts)//2
  # There is one cyclic subgroup of order99. Enumerate every marked generator.
  P=pts[0];ks={frobkey(pointbc(E,u*P),d) for u in range(1,50) if gcd(u,99)==1}
  assert len(ks)==30;allkeys.update(ks);sourceorbits.add(min(ks));pointorders.append([int(card),len(pts)])
 rational.append(levels);v=levels-sum(e*closed[e-1] for e in range(1,d) if d%e==0);assert v%d==0;closed.append(v//d)
 census.append({'degree':d,'eligible_curves':eligible,'marked_points_mod_sign':levels,'closed_places':v//d})
assert closed==[0]*9+[120] and len(sourceorbits)==4 and len(allkeys)==120
assert sorted(sourceorbits)==sorted(tuple(ky) for d,ky in source[2])
print('C99_INDEPENDENT_EXHAUSTIVE_SOURCE_PASS',closed,4,flush=True)

pol=next(row[1] for row in source[0] if row[0]==10)
K=GF(2**10,'a',modulus=R2(pol));a=K.gen()
out=[]
for ci,(d,ky) in enumerate(source[2],1):
 assert d==10
 b=R2(ky[:d])(a);c=R2(ky[d:])(a);E0=EllipticCurve(K,[1-c,-b,-b,0,0]);P0=E0(0,0)
 assert 99*P0==E0(0) and 33*P0!=E0(0) and 9*P0!=E0(0)
 # All five-torsion points are rational over F_(2^20).
 L=GF(2**20,'u');RP=PolynomialRing(L,'x');x=RP.gen();root=RP(pol).roots(multiplicities=False)[0]
 emb=lambda v:sum((L(v.polynomial()[i])*root**i for i in range(10)),L(0))
 E=EllipticCurve(L,list(map(emb,E0.a_invariants())));P=E(list(map(emb,P0.xy())));Q=3*P
 assert 99*P==E(0) and 33*Q==E(0) and 11*Q!=E(0) and 3*Q!=E(0)
 xs=E.division_polynomial(5).roots(multiplicities=False);assert len(xs)==12
 pts=[T for sx in xs for T in E.lift_x(sx,all=True)];assert len(pts)==24 and all(5*T==E(0) and T!=E(0) for T in pts)
 kernels={}
 for T in pts:
  ker=(x-T[0])*(x-(2*T)[0]);assert ker.degree()==2 and ker.is_squarefree()
  kernels[tuple(ker.list())]=T
 assert len(kernels)==6
 actual=Counter();observations=[]
 for coeffs,T in kernels.items():
  ker=RP(list(coeffs));phi=E.isogeny(ker);ET=phi.codomain();QT=phi(Q)
  assert phi.degree()==5 and 33*QT==ET(0) and 11*QT!=ET(0) and 3*QT!=ET(0)
  nonzero=[j*T for j in range(1,5)]
  vx=Q[0]+sum(((Q+S)[0]-S[0] for S in nonzero),L(0))
  vy=Q[1]+sum(((Q+S)[1]-S[1] for S in nonzero),L(0))
  assert (vx,vy)==QT.xy()
  key=key2(unitxy(ET,QT));actual[json.dumps(key)]+=1
 # Parse every positive primary target over its explicit absolute field.
 expected=Counter();pobs=[]
 for row in primary[2][ci-1][5]:
  orbit,norm,pd,co,mp=row;A=GF(2**(len(mp)-1),'v',modulus=R2(mp));v=A.gen()
  xy=[sum((A(c)*v**i for i,c in enumerate(z)),A(0)) for z in co]
  key=key2(xy);assert key[0]==pd[0] and 10*orbit==norm*pd[0]
  expected[json.dumps(key)]+=orbit;pobs.append([int(orbit),int(norm),int(pd[0]),key])
 assert actual==expected and sorted(v[0] for v in pobs)==[1,1,2,2]
 assert sorted(v[2] for v in pobs)==[10,10,20,20] and all(v[1]==1 for v in pobs)
 out.append({'case':ci,'source_key':key2(unitxy(E,Q)),'diamond5_key':key2(unitxy(E,5*Q)),'kernel_orbit_observations':pobs,'geometric_cyclic_kernels':6,'individual_nonzero_kernel_points':24,'manual_Velu_coordinate_checks':6})
 print('C99_INDEPENDENT_GENUINE_FIVE_KERNELS_PASS',ci,[v[0] for v in pobs],flush=True)
result={'status':'PASS','method':'Exhaustive rational-point source census and six actual quadratic5-kernels, Sage isogenies, four-point manual Velu sums, physical33-level and exact primary residue-key/norm binding','finite_NC_closed_counts':closed,'source_closed_places':120,'source_diamond_orbits':4,'source_degree10_census':census,'normalized_source_records_exact':True,'source_cases':out,'total_geometric_cyclic5_kernels':24,'all_primary_RR_h0':[v[6] for v in primary[2]],'elapsed_seconds':time.monotonic()-START}
(ROOT/'independent_genuine_isogenies.json').write_text(json.dumps(result,indent=2)+'\n')
print('C99_INDEPENDENT_SOURCE_AND_GENUINE_ISOGENY_COMPLETE',24,time.monotonic()-START,flush=True)
