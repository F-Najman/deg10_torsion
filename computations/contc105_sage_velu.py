from sage.all import *
from pathlib import Path
import ast,json,time
ROOT=Path('work/deg10/continuation_20261030c/C105');START=time.monotonic()
def read(name):return ast.literal_eval((ROOT/name).read_text().replace('[*','[').replace('*]',']'))
F2=GF(2);R2=PolynomialRing(F2,'T')
rows=read('kernel3_explicit.m');source=read('source105.m');native=read('parent_native_RR.m')
def unitxy(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y
 n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 b=-n2**3/b3**2;c=-(n1*n2-b3)/b3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def key2(v):
 K=v[0].parent();n=K.degree();d=lcm([a.minpoly().degree() for a in v]);found=False
 for bits in [[1,0],[0,1],[1,1]]:
  g=bits[0]*v[0]+bits[1]*v[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found
 vv=lambda a:vector(F2,[K(a).polynomial()[i] for i in range(n)])
 M=matrix(F2,[vv(g**i) for i in range(d)]);assert M.rank()==d
 co=[list(map(int,M.transpose().solve_right(vv(a)))) for a in v]
 return [int(d),bits,list(map(int,mp.list())),co]
checks=[];seen={};degrees={}
for row in rows:
 ci,e,poly,zbits,eco,p105co,qco,sco,tco,qtco,xyco,td,norm=row
 K=GF(2**(len(poly)-1),'a',modulus=R2(poly),impl='pari_ffelt');a=K.gen()
 dec=lambda v:sum((K(c)*a**i for i,c in enumerate(v)),K(0))
 E=EllipticCurve(K,list(map(dec,eco)));P=E(list(map(dec,p105co)));Q=E(list(map(dec,qco)));s=dec(sco)
 assert 105*P==E(0) and all((105//p)*P!=E(0) for p in [3,5,7]) and Q==3*P
 z=dec(zbits);d,ky=source[2][ci-1]
 b=sum((K(c)*z**i for i,c in enumerate(ky[:d])),K(0));c=sum((K(c)*z**i for i,c in enumerate(ky[d:])),K(0))
 assert list(E.a_invariants())==[1-c,-b,-b,K(0),K(0)] and P==E(0,0)
 RP=PolynomialRing(K,'x');x=RP.gen();assert E.division_polynomial(3)(s)==0
 phi=E.isogeny(x-s);ET=phi.codomain();QT=phi(Q)
 assert list(ET.a_invariants())==list(map(dec,tco)) and list(QT.xy())==list(map(dec,qtco))
 # Reconstruct Velu coordinates directly from both nonzero kernel points.
 ry=PolynomialRing(K,'v');v=ry.gen();aa=E.a_invariants();yp=v*v+(aa[0]*s+aa[2])*v-(s**3+aa[1]*s*s+aa[3]*s+aa[4]);ys=yp.roots(multiplicities=False)
 if ys:
  EE=E;QQ=Q;T=EE(s,ys[0]);targetxy=QT.xy()
 else:
  L=K.extension(yp,'u');EE=E.change_ring(L);QQ=EE(Q);T=EE(L(s),L.gen());targetxy=tuple(map(L,QT.xy()))
 assert 3*T==EE(0) and T!=EE(0)
 manualx=QQ[0]+(QQ+T)[0]+(QQ-T)[0]-2*T[0]
 manualy=QQ[1]+(QQ+T)[1]+(QQ-T)[1]-T[1]-(-T)[1]
 assert (manualx,manualy)==targetxy
 xy=unitxy(ET,QT);assert xy==list(map(dec,xyco));key=key2(xy)
 assert key[0]==td and td*norm==10*e
 nativecase=native[ci-1];assert nativecase[0]==ci
 assert key2(unitxy(E,Q))==nativecase[1] and key2(unitxy(E,3*Q))==nativecase[2]
 assert [e,norm,td,key] in nativecase[3]
 seen.setdefault(ci,[]).append([e,norm,td,key]);degrees[ci]=degrees.get(ci,0)+e
 checks.append({'case':ci,'kernel_orbit_degree':e,'target_closed_degree':td,'norm':norm,'exact_source_order':105,'Sage_isogeny_and_manual_Velu_coordinates_match':True,'native_positive_prime_intrinsic_key_matches':True})
 print('C105_SAGE_GENUINE_KERNEL_POINT_PASS',ci,e,td,norm,flush=True)
assert len(rows)==12 and len(seen)==6 and all(v==4 for v in degrees.values())
for ci,obs in seen.items():assert sorted(obs,key=str)==sorted(native[ci-1][3],key=str)
out={'status':'PASS','cases':6,'actual_kernel_orbits':12,'all_four_cyclic3_kernels_covered_per_case':True,'manual_Velu_coordinates':True,'exact_native_RR_prime_bindings':True,'checks':checks,'elapsed_seconds':time.monotonic()-START}
(ROOT/'independent_sage_velu.json').write_text(json.dumps(out,indent=2)+'\n')
print('C105_SAGE_VELU_ALL_SIX_PASS',time.monotonic()-START,flush=True)
