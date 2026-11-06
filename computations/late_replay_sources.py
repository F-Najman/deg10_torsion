from pathlib import Path
import ast,json,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
H=Path('work/deg10');P=H/'continuation_20261029b';start=time.monotonic()
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def key(pair,d):
 out=[]
 for _ in range(d):
  out.append(tuple(int(c) for z in pair for c in z._vector_()));pair=[z**3 for z in pair]
 return min(out)
def bc(E,T,Q):
 x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
reports=[]
for n,g in [(38,64),(40,57)]:
 C=P/f'C2{n}';fp,cd,nc=read(C/('frobenius.m' if n==38 else 'source_frobenius.txt'));assert len(fp)==2*g+1
 trace=[];closed=[]
 for d in range(1,11):
  tr=-d*fp[2*g-d]-sum(fp[2*g-j]*trace[d-j-1] for j in range(1,d));trace.append(tr)
  num=3**d+1-tr-sum(e*closed[e-1] for e in range(1,d) if d%e==0);assert num%d==0;closed.append(num//d)
 assert [closed[d-1]-Counter(cd)[d] for d in range(1,11)]==nc
 counts=[];orbits=[]
 R=PolynomialRing(GF(3),'a')
 for d in range(4,11):
  a=read(C/(f'source{d}.m' if n==38 else f'M2N40_D10_direct_full{d}_arithmetic.txt'));assert a[0]==d and a[2]
  K=GF(3**d,'a',modulus=R(a[1]));z=K.gen();rows=a[6] if n==38 else a[8];seen=set();total=Counter();own=0
  for deg,ky in rows:
   b=sum(K(v)*z**i for i,v in enumerate(ky[:d]));c=sum(K(v)*z**i for i,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
   assert 2*T==E(0) and n*Q==E(0) and (n//2)*Q!=E(0) and (n//(19 if n==38 else 5))*Q!=E(0) and T!=(n//2)*Q
   if n==40:assert 20*Q!=E(0)
   actual=min(int(e) for e in divisors(d) if b**(3**int(e))==b and c**(3**int(e))==c);assert actual==deg
   ims=set()
   if n==38:
    ts=[T,19*Q,T+19*Q];V=20*Q
    for i in range(3):
     for j in range(3):
      if i==j:continue
      for k in range(1,10):ims.add(key(bc(E,ts[i],ts[j]+k*V),d))
   else:
    for bb in range(2):
     for cc in range(2):
      for aa in [1,3,7,9,11,13,17,19]:ims.add(key(bc(E,T+20*bb*Q,cc*T+aa*Q),d))
   assert tuple(ky)==min(ims) and seen.isdisjoint(ims);seen.update(ims);total[deg]+=len(ims);own+=deg==d
  native=a[5] if n==38 else a[7]
  assert [total[int(e)] for e in divisors(d)]==native and all(total[int(e)]==nc[int(e)-1] for e in divisors(d))
  counts.append(total[d]);orbits.append(own)
  print('INDEPENDENT_SOURCE_REPLAY',n,d,total[d],own,time.monotonic()-start,flush=True)
 reports.append({'group':[2,n],'prime':3,'source_genus':g,'closed_non_cuspidal_counts_through_10':nc,'cusp_signature':dict(Counter(cd)),'full_automorphism_orbit_counts_4_through_10':orbits,'independently_verified_full_coverage':True})
(P/'independent_sources.json').write_text(json.dumps(reports,indent=2)+'\n')
print('BOTH_INDEPENDENT_SOURCE_REPLAYS_COMPLETE',time.monotonic()-start,flush=True)
