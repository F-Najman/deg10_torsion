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
def order55(P):return 55*P==P.curve()(0) and 5*P!=P.curve()(0) and 11*P!=P.curve()(0)
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
   if n%55:continue
   assert all(not(n%(p*p)==0 and (q-1)%p==0)for p in [5,11]);qual+=1
 # Every j=0 curve is supersingular. The following set contains every
 # possible trace, and none gives order divisible by55 in these fields.
 traces={0}
 if d%2:traces.update([-2**((d+1)//2),2**((d+1)//2)])
 else:traces.update([-2**(d//2),2**(d//2),-2**(d//2+1),2**(d//2+1)])
 assert all((q+1-t)%55 for t in traces)
 prior=sum(e*counts[e-1]for e in range(1,d)if d%e==0);assert (20*qual-prior)%d==0
 counts.append((20*qual-prior)//d);qualifying.append(qual)
assert counts==census[0]==[0,0,0,0,0,0,20,0,40,100]
assert qualifying==[r[1]for r in census[1]]
# Check all stored source generators and j-orbit representatives in their
# explicitly specified fields, independently of the complete recount.
for row in census[2]:
 d,defn,aa,pp,jkey,n,iv=row;F=field(defn,'a');E=EllipticCurve(F,[value(F,v)for v in aa]);P=E([value(F,v)for v in pp]);assert order55(P)and E.cardinality()==n
 assert orbitdegree([E.j_invariant()])==d
 assert list(orbitkey([E.j_invariant()]))==jkey
outputs=[];geometric=0
for idx in range(1,6):
 packet=read(H/f'source_case{idx}_packet.m');row=packet[6];obs=packet[7]
 n=lcm([len(v[4][0])-1 for v in obs]);F=GF(2**n,'b');T=PolynomialRing(F,'t');t=T.gen()
 def embed(defn):
  f=T(defn);roots=f.roots(multiplicities=False);assert len(roots)==len(defn)-1;return roots[0]
 a=embed(row[1]);vv=lambda v:T(v)(a)
 E=EllipticCurve(F,[vv(v)for v in row[2]]);P=E([vv(v)for v in row[3]]);assert order55(P)
 kerpoly=E.division_polynomial(3);roots=kerpoly.roots(multiplicities=False);assert len(roots)==4
 expected=Counter()
 for r in roots:
  phi=E.isogeny(t-r);PP=phi(P);assert order55(PP);expected[orbitkey(coords(phi.codomain(),PP))]+=1
 actual=Counter();details=[]
 for number,w,d,m,data in obs:
  assert d==10*w and m==1
  beta=embed(data[0]);val=lambda v:T(v)(beta);vals=[val(v)for v in data[1]]
  EE=EllipticCurve(F,vals[:5]);PP=EE(vals[5:7]);assert order55(PP)
  xy=coords(EE,PP);assert list(xy)==vals[7:9]and orbitdegree(xy)==d
  actual[orbitkey(xy)]+=w;details.append([int(w),int(d),int(m)])
 assert expected==actual and sum(actual.values())==4
 # Recover the complete division polynomial from all four distinct kernels.
 product=T(1)
 for r in roots:product*=t-r
 assert T(kerpoly)==kerpoly.leading_coefficient()*product
 geometric+=4;outputs.append(dict(case=idx,geometric_kernels=4,orbits=details))
out=dict(status='PASS',independent_complete_finite_census=True,
 ordinary_twists_by_degree=qualifying,noncuspidal_closed_counts=counts,
 stored_source_generators_checked=len(census[2]),new_degree10_diamond_representatives=5,
 independent_marked_Velu3_kernels=geometric,
 complete_three_division_polynomials_recovered=5,all_norm_weights_and_marked_Tate_coordinates_checked=True,
 cases=outputs,seconds=time.monotonic()-t0,
 evidence=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in [H/'census.m']+[H/f'source_case{i}_packet.m'for i in range(1,6)]])
(H/'independent_finite_geometry.json').write_text(json.dumps(out,indent=2)+'\n')
print('C55_INDEPENDENT_CENSUS_AND_MARKED_ISOGENIES_PASS',json.dumps(out,sort_keys=True))
