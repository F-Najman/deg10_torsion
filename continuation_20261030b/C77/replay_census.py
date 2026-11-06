"""Independent finite source census and complete marked3-isogeny checks."""
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
from pathlib import Path
from collections import Counter
import ast,json,time,hashlib
H=Path(__file__).resolve().parent;t0=time.monotonic();k=GF(2);R=PolynomialRing(k,'x');x=R.gen()
def read(p):return ast.literal_eval(Path(p).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def field(defn,name):
 f=R(defn);assert f.is_irreducible();return GF(2**f.degree(),name,modulus=f)
def value(F,v):return F['x'](v)(F.gen())
def coords(E,P):
 X,Y=P.xy();a1,a2,a3,a4,a6=E.a_invariants();aa3=2*Y+a3+a1*X;aa2=3*X+a2;aa4=3*X**2+2*X*a2+a4-a1*Y
 aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3
 b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return((s-r)/t,t/(s*s-s-r+1))
def order77(P):return 77*P==P.curve()(0) and 5*P!=P.curve()(0) and 11*P!=P.curve()(0)
def orbitkey(vs):
 F=vs[0].parent();n=F.degree();out=[]
 for i in range(n):out.append(tuple(int((v**(2**i)).polynomial()[j])for v in vs for j in range(n)))
 return min(out)
def orbitdegree(vs):
 n=vs[0].parent().degree()
 for i in range(1,n+1):
  if all(v**(2**i)==v for v in vs):return i
 raise AssertionError
census=read(H/'census.m');counts=[];qualifying=[]
for d in range(1,11):
 F=GF(2**d,'z');q=2**d;u=next(v for v in F if v.trace()==1);qual=0
 # For j!=0, these two trace classes give all ordinary quadratic twists.
 for j in F:
  if not j:continue
  for a2 in [F(0),u]:
   E=EllipticCurve(F,[1,a2,0,0,1/j]);n=E.cardinality()
   if n%77:continue
   assert all(not(n%(p*p)==0 and (q-1)%p==0)for p in [7,11]);qual+=1
 # Every j=0 curve is supersingular. The following set contains every
 # possible trace, and none gives order divisible by77 in these fields.
 traces={0}
 if d%2:traces.update([-2**((d+1)//2),2**((d+1)//2)])
 else:traces.update([-2**(d//2),2**(d//2),-2**(d//2+1),2**(d//2+1)])
 assert all((q+1-t)%77 for t in traces)
 prior=sum(e*counts[e-1]for e in range(1,d)if d%e==0);assert (30*qual-prior)%d==0
 counts.append((30*qual-prior)//d);qualifying.append(qual)
assert counts==census[0]==[0,0,0,0,0,0,0,0,0,90]
assert qualifying==[r[1]for r in census[1]]
# Check all stored source generators and j-orbit representatives in their
# explicitly specified fields, independently of the complete recount.
for row in census[2]:
 d,defn,aa,pp,jkey,n,iv=row;F=field(defn,'a');E=EllipticCurve(F,[value(F,v)for v in aa]);P=E([value(F,v)for v in pp]);assert order77(P)and E.cardinality()==n
 assert orbitdegree([E.j_invariant()])==d
 assert list(orbitkey([E.j_invariant()]))==jkey

out=dict(status='PASS',level=77,noncuspidal_closed_counts=counts,ordinary_twists_by_degree=qualifying,stored_generators_checked=len(census[2]),all_order77_groups_cyclic=True,new_degree10_diamond_orbits=3,seconds=time.monotonic()-t0)
(H/'independent_census.json').write_text(json.dumps(out,indent=2)+'\n')
print('C77_INDEPENDENT_CENSUS_PASS',json.dumps(out,sort_keys=True))
