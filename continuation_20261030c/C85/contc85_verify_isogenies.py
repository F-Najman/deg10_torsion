from pathlib import Path
from collections import Counter
import json,time
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
from arithmetic10 import lit
P=Path(__file__).resolve().parent;t0=time.monotonic();s=(P/'jets10.m').read_text();ss=lit(s,'SourceRecords');ts=lit(s,'SourceHeckeRecords');K=GF(2**30,'a');a=K.gen();RK=PolynomialRing(K,'x');roots={}
def decode(rec):
 pol,aa,pt,zz,uv,jt=rec;key=tuple(pol)
 if key not in roots:roots[key]=RK(pol).roots(multiplicities=False)[0]
 r=roots[key]
 cv=lambda v:sum(K(b)*r**i for i,b in enumerate(v))
 E=EllipticCurve(K,[cv(v)for v in aa]);Q=E(*[cv(v)for v in pt]);assert 85*Q==E(0)and 5*Q!=E(0)and 17*Q!=E(0)
 return E,Q,[cv(v)for v in uv]
def unit(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=a3+a1*x;b2=x+a2;b4=x*x+a4-a1*y;sl=b4/b3;b2-=sl*sl+a1*sl;b=-b2**3/b3**2;c=1-a1*b2/b3;r=b/c;s=c*c/(b-c);t=r*s+1;return[(s-r)/t,t/(s*s-s-r+1)]
def bits(v):return sum(int(c)<<i for i,c in enumerate(v._vector_()))
def key(vs):
 d=int(lcm([v.minpoly().degree()for v in vs]));return d,min(tuple(bits(v**(2**i))for v in vs)for i in range(d))
checks=[]
for i,(sr,target)in enumerate(zip(ss,ts),1):
 E,P0,uv=decode(sr);assert unit(E,P0)==uv;assert key(uv)[0]==10
 native=Counter()
 for e,d,m,rec in target:
  EG,PG,units=decode(rec);assert unit(EG,PG)==units and key(units)[0]==d and e*10==d*m;native[key(units)]+=m
 xs=RK(E.division_polynomial(3)).roots(multiplicities=False);assert len(xs)==4;geo=Counter();RY=PolynomialRing(K,'y');Y=RY.gen();a1,a2,a3,a4,a6=E.a_invariants()
 for x in xs:
  ys=(Y*Y+(a1*x+a3)*Y-(x**3+a2*x*x+a4*x+a6)).roots(multiplicities=False);assert len(ys)==2;T=E(x,ys[0]);assert T!=E(0)and 3*T==E(0);T.set_order(3);phi=E.isogeny(T,check=True);assert phi.degree()==3;PG=phi(P0);assert 85*PG==phi.codomain()(0);geo[key(unit(phi.codomain(),PG))]+=1
 expected=Counter()
 for ky,n in geo.items():assert 10*n%ky[0]==0;expected[ky]=10*n//ky[0]
 assert native==expected and sum(ky[0]*m for ky,m in native.items())==40
 checks.append(dict(source=i,all_four_geometric_kernels=True,target_closed_degrees=sorted(ky[0]for ky in native),actual_point_images_verified=True,norms_verified=True))
out=dict(status='PASS',sources=4,geometric_kernels=16,physical_isogenies_independently_rebuilt=True,checks=checks,seconds=time.monotonic()-t0);(P/'isogenies10_independent.json').write_text(json.dumps(out,indent=2)+'\n');print('C85_DEG10_INDEPENDENT_ALL16_ACTUAL_3_ISOGENIES_PASS',time.monotonic()-t0,flush=True)
