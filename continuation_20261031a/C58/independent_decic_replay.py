from pathlib import Path
import ast,re,json,time,hashlib
from collections import Counter
from sage.all import GF,QQ,PolynomialRing,EllipticCurve,matrix,vector,ZZ,identity_matrix
H=Path(__file__).resolve().parent;D=H/'inputs/C58_X129';tt=time.monotonic();k=GF(3);R=PolynomialRing(k,'z')
def lit(path,name):
 s=path.read_text();v=re.search(r'\b'+re.escape(name)+r':=(.*?);',s,re.S).group(1)
 v=re.sub(r'\[\s*\*','[',v);v=re.sub(r'\*\s*\]',']',v);v=re.sub(r'\btrue\b','True',v);v=re.sub(r'\bfalse\b','False',v)
 return ast.literal_eval(v)
def counts(pc,q,n=10):
 aa=list(reversed(pc));power=[0];closed=[]
 for j in range(1,n+1):
  power.append(-j*aa[j]-sum(aa[i]*power[j-i]for i in range(1,j)))
  z=q**j+1-power[j]-sum(i*closed[i-1]for i in range(1,j)if j%i==0);assert z%j==0;closed.append(z//j)
 return closed
G=json.loads((D/'geometry58.json').read_text());C=json.loads((D/'characters58.json').read_text())
PK=lit(H/'degree10_census.m','PointKeys');nr=lit(H/'degree10_census.m','EligibleTwoCounts');reps=lit(H/'degree10_census.m','Representatives');orbits=lit(H/'degree10_census.m','Orbits');tests=lit(H/'degree10_complete.m','Tests')
assert len(PK)==counts(G['N29_frobenius_coefficients'],3)[9]==5971
flat=[i for o in orbits for i in o];assert sorted(flat)==list(range(1,5972))and len(orbits)==len(reps)==430 and all(reps[j]in o for j,o in enumerate(orbits))
fields={};Bs={};objects={};keys={}
def coeff(v):return [int(a)for a in v.polynomial().list()]+[0]*(10-len(v.polynomial().list()))
def tuplekey(v):return tuple(tuple(int(a)for a in row)for row in v)
def key(x,y):
 for h in [x,y,x+y,x-y]:
  if h.minpoly().degree()==10:break
 else:raise AssertionError('no primitive coordinate')
 A=matrix(k,[coeff(h**i)for i in range(10)]);AI=A.inverse()
 return (tuple(int(a)for a in h.minpoly().list()),tuple(int(a)for a in vector(k,coeff(x))*AI),tuple(int(a)for a in vector(k,coeff(y))*AI))
def xy(E,P):
 a1,a2,a3,a4,a6=E.ainvs();x,y=P.xy();v3=2*y+a3+a1*x;v2=3*x+a2;v4=3*x*x+2*x*a2+a4-a1*y
 z1=2*v4/v3+a1;z2=v2-(v4/v3)**2-a1*v4/v3;b=-z2**3/v3**2;c=-(z1*z2-v3)/v3
 r=b/c;s=c*c/(b-c);v=r*s-2*r+1
 return (s-r)/v,v/(s*s-s-r+1)
source_total=0
for i,entry in enumerate(PK,1):
 cp,xx,yy=entry;cp=tuple(cp)
 if cp not in fields:
  f=R(cp);assert f.degree()==10 and f.is_irreducible();fields[cp]=GF(3**10,'a',modulus=f)
 F=fields[cp];a=F.gen();x=R(xx)(a);y=R(yy)(a);ky=key(x,y);assert ky==tuplekey(entry)and ky not in keys;keys[ky]=i
 r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1)
 E=EllipticCurve(F,[1-c,-b,-b,0,0]);P=E(0,0);assert E.discriminant()!=0 and 29*P==E(0)and P!=E(0)
 assert sum(int(co)*x**ex[0]*y**ex[1]for ex,co in zip(G['N29_affine_exponents'],G['N29_affine_coefficients']))==0
 S=PolynomialRing(F,'u');u=S.gen();a1,a2,a3,a4,a6=E.ainvs();two=4*u**3+(a1*a1+4*a2)*u*u+2*(a1*a3+2*a4)*u+a3*a3+4*a6
 num=len(two.roots());assert num==nr[i-1]and num in [0,1,3];source_total+=num
 if i in set(reps):objects[i]=(E,P)
 if i%1000==0:print('D31A_C58_INDEPENDENT_POINTS',i,round(time.monotonic()-tt,2),flush=True)
assert source_total==counts(G['N58_frobenius_coefficients'],3)[9]
for j,i in enumerate(reps):
 E,P=objects[i];Q=P;seen=set()
 for a in range(1,15):
  xx,yy=xy(E,Q);seen.add(keys[key(xx,yy)]);Q+=P
 assert seen==set(orbits[j])and all(nr[z-1]==nr[i-1]for z in seen)
print('D31A_C58_INDEPENDENT_ORBITS',len(reps),flush=True)
# Exact lattice tests for the two complete primary images.
checks=[]
for mod in [8,9]:
 A=matrix(ZZ,C['matrix'+str(mod)]);B=A.matrix_from_columns([i-1 for i in C['cusp_FB_indices']]);HH=B.augment(mod*identity_matrix(ZZ,A.nrows())).transpose().hermite_form(include_zero_rows=False).transpose();checks.append((mod,A,mod*HH.inverse()))
for a in tests:
 i,ims,yes,pr,wit=a;extra,ff,inds,es=wit;assert i in reps and nr[i-1]>0
 for (mod,A,HH),im in zip(checks,ims):
  v=vector(ZZ,[int(-sum(e*A[j,z-1]for z,e in zip(inds,es)))%mod for j in range(A.nrows())]);assert list(v)==im
 membership=all(all(v in ZZ for v in HH*vector(QQ,im)/mod)for (mod,A,HH),im in zip(checks,ims));assert membership==yes
assert {a[0]for a in tests}=={i for i in reps if nr[i-1]>0}and len(tests)==293
assert [a[0]for a in tests if a[2]]==[462,1748]
# Independent checks of both exact known global models, reductions, and odd cardinalities.
known=[]
for label,rep,root29 in [('a',1748,11),('b',462,5)]:
 A=json.loads((H/'inputs'/('1.29.10.'+label+'.json')).read_text());QX=PolynomialRing(QQ,'x');ff=QX(list(map(QQ,A['field']['poly_coeffs'])));assert ff.is_irreducible();assert R(ff).is_irreducible()
 F=GF(3**10,'b',modulus=R(ff));alpha=F.gen()
 def evalrat(co,x):
  F=x.parent();z=F(0)
  for c in reversed(co):v=QQ(c);z=z*x+F(v.numerator())/F(v.denominator())
  return z
 b=evalrat(A['curve']['tate_b_coeffs'],alpha);c=evalrat(A['curve']['tate_c_coeffs'],alpha);E=EllipticCurve(F,[1-c,-b,-b,0,0]);P=E(0,0);assert 29*P==E(0)
 seen=set();Q=P
 for j in range(14):xx,yy=xy(E,Q);seen.add(keys[key(xx,yy)]);Q+=P
 assert seen==set(orbits[reps.index(rep)])
 L=GF(29);a=L(root29);assert ff.change_ring(L)(a)==0 and ff.derivative().change_ring(L)(a)!=0
 b=evalrat(A['curve']['tate_b_coeffs'],a);c=evalrat(A['curve']['tate_c_coeffs'],a);E=EllipticCurve(L,[1-c,-b,-b,0,0]);assert E.discriminant()!=0 and E.cardinality()==29
 known.append(dict(label=label,representative=rep,orbit_size=len(seen),good_prime=29,simple_root=root29,odd_cardinality=29))
out=dict(status='PASS',target_degree10_places=5971,source_degree10_places=source_total,diamond_orbits=430,eligible_class_tests=293,survivor_representatives=[462,1748],known_lifts=known,seconds=time.monotonic()-tt,
 complete_marked_point_and_2division_census=True,all_diamond_orbits_reconstructed=True,all_primary_sums_and_memberships_verified=True,
 native_ideal_identities_are_separately_replayed=True)
(H/'independent_decic_replay.json').write_text(json.dumps(out,indent=2)+chr(10));print('D31A_C58_INDEPENDENT_REPLAY_COMPLETE',out,flush=True)
