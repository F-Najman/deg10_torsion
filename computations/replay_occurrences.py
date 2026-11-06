from pathlib import Path
import re,json,math
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,EllipticCurve
H=Path(__file__).resolve().parents[1];R=PolynomialRing(QQ,'t')
arr=lambda s:[QQ(x.strip()) for x in s.split(',')]
for n in [29,31,35,37,38,39,40,42,45,50]:
 stem='witness_C'+str(n)+('_v2' if n==42 else '');s=(H/f'computations/{stem}.m').read_text();log=(H/f'jobs/{stem}.log').read_text()
 f=R(arr(re.search(r'f:=R!\[(.*?)\]',s)[1]));bb=R(arr(re.search(r'bb:=R!\[(.*?)\]',s)[1]));cc=R(arr(re.search(r'cc:=R!\[(.*?)\]',s)[1]));assert f.degree()==10 and f.is_irreducible();K=NumberField(f,'a');a=K.gen();b,c=bb(a),cc(a);E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E([0,0]);assert n*P==E(0)
 for p in ZZ(n).prime_divisors():assert (n//p)*P!=E(0)
 bound=0
 for p,r,count in re.findall(r'N=\d+ p=(\d+) root=(\d+) residue_degree=1 count=(\d+)',log):
  p,r,count=map(int,(p,r,count));F=GF(p);ff=f.change_ring(F);assert ff(F(r))==0 and ff.derivative()(F(r))!=0
  b=bb.change_ring(F)(F(r));c=cc.change_ring(F)(F(r));Er=EllipticCurve(F,[1-c,-b,-b,0,0]);assert Er.cardinality()==count;bound=math.gcd(bound,count)
 assert bound==n;print('DECIC_CYCLIC_WITNESS_REPLAY_PASS',n,flush=True)
for name,n in [('28',28),('30a',30)]:
 stem='quintic_full2_'+name+'_v2';s=(H/f'computations/{stem}.m').read_text();log=(H/f'jobs/{stem}.log').read_text()
 f=R(arr(re.search(r'f:=R!\[(.*?)\]',s)[1]));bb=R(arr(re.search(r'b:=Evaluate\(R!\[(.*?)\]',s)[1]));cc=R(arr(re.search(r'c:=Evaluate\(R!\[(.*?)\]',s)[1]));assert f.is_irreducible() and f.degree()==5
 K=NumberField(f,'a');a=K.gen();b,c=bb(a),cc(a);E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E([0,0]);assert n*P==E(0)
 for p in ZZ(n).prime_divisors():assert (n//p)*P!=E(0)
 factors=E.division_polynomial(2).factor();hs=[h.monic() for h,e in factors if h.degree()==2];assert len(hs)==1 and hs[0].is_irreducible();h=hs[0]
 p,r,w,count=map(int,re.search(r'p=(\d+) quintic_root=(\d+) relative_root=(\d+) count=(\d+)',log).groups());F=GF(p);red=lambda z:sum((F(v)*F(r)**j for j,v in enumerate(K(z).list())),F(0))
 ff=f.change_ring(F);assert ff(F(r))==0 and ff.derivative()(F(r))!=0
 S=PolynomialRing(F,'x');hh=S([red(v) for v in h.list()]);assert hh(F(w))==0 and hh.derivative()(F(w))!=0
 Er=EllipticCurve(F,[red(z) for z in E.a_invariants()]);assert Er.cardinality()==count==2*n
 print('DECIC_FULL2_TOWER_REPLAY_PASS',n,flush=True)
print('DECIC_OCCURRENCES_REPLAY_COMPLETE count=12',flush=True)
