from pathlib import Path
import ast,re,json,time,hashlib
from collections import Counter
from sage.all import GF,QQ,PolynomialRing,EllipticCurve,matrix,vector,ZZ,identity_matrix
H=Path(__file__).resolve().parent;D=H/'inputs/C62_X131';tt=time.monotonic();k=GF(3);R=PolynomialRing(k,'z')
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
G=json.loads((D/'geometry62.json').read_text());C=json.loads((D/'characters62.json').read_text())
PK=lit(H/'degree10_resume_v2_census.m','PointKeys');nr=lit(H/'degree10_resume_v2_census.m','EligibleTwoCounts');reps=lit(H/'degree10_resume_v2_census.m','Representatives');orbits=lit(H/'degree10_resume_v2_census.m','Orbits');tests=lit(H/'degree10_resume_v2_complete.m','Tests')
assert len(PK)==counts(G['N31_frobenius_coefficients'],3)[9]==5685
flat=[i for o in orbits for i in o];assert sorted(flat)==list(range(1,5686))and len(orbits)==len(reps)==379 and all(reps[j]in o for j,o in enumerate(orbits))
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
 E=EllipticCurve(F,[1-c,-b,-b,0,0]);P=E(0,0);assert E.discriminant()!=0 and 31*P==E(0)and P!=E(0)
 assert sum(int(co)*x**ex[0]*y**ex[1]for ex,co in zip(G['N31_affine_exponents'],G['N31_affine_coefficients']))==0
 S=PolynomialRing(F,'u');u=S.gen();a1,a2,a3,a4,a6=E.ainvs();two=4*u**3+(a1*a1+4*a2)*u*u+2*(a1*a3+2*a4)*u+a3*a3+4*a6
 num=len(two.roots());assert num==nr[i-1]and num in [0,1,3];source_total+=num
 if i in set(reps):objects[i]=(E,P)
 if i%1000==0:print('D31A_C62_INDEPENDENT_POINTS',i,round(time.monotonic()-tt,2),flush=True)
assert source_total+45==counts(G['N62_frobenius_coefficients'],3)[9]
for j,i in enumerate(reps):
 E,P=objects[i];Q=P;seen=set()
 for a in range(1,16):
  xx,yy=xy(E,Q);seen.add(keys[key(xx,yy)]);Q+=P
 assert seen==set(orbits[j])and all(nr[z-1]==nr[i-1]for z in seen)
print('D31A_C62_INDEPENDENT_ORBITS',len(reps),flush=True)
# Exact lattice tests for the two complete primary images.
checks=[]
for mod in [16,31,61]:
 A=matrix(ZZ,C['matrix'+str(mod)]);B=A.matrix_from_columns([i-1 for i in C['cusp_FB_indices']]);HH=B.augment(mod*identity_matrix(ZZ,A.nrows())).transpose().hermite_form(include_zero_rows=False).transpose();checks.append((mod,A,mod*HH.inverse()))
for a in tests:
 i,ims,yes,pr,wit=a;extra,ff,inds,es=wit;assert i in reps and nr[i-1]>0
 for (mod,A,HH),im in zip(checks,ims):
  v=vector(ZZ,[int(-sum(e*A[j,z-1]for z,e in zip(inds,es)))%mod for j in range(A.nrows())]);assert list(v)==im
 membership=all(all(v in ZZ for v in HH*vector(QQ,im)/mod)for (mod,A,HH),im in zip(checks,ims));assert membership==yes
assert {a[0]for a in tests}=={i for i in reps if nr[i-1]>0}and len(tests)==253
assert [a[0]for a in tests if a[2]]==[]

# Reconstruct the exact known nonic31 trace and all five surviving lower-degree reductions.
old=json.loads((D/'classification62.json').read_text());surv=set(old['survivor_indices']);assert len(surv)==5
QQX=PolynomialRing(QQ,'w');ff=QQX([-1,-2,19,-28,14,-1,1,0,-2,1]);assert ff.is_irreducible()and R(ff).is_irreducible()
L=GF(3**9,'v',modulus=R(ff));v=L.gen();xx=v;yy=sum(L(QQ(a).numerator())/L(QQ(a).denominator())*v**i for i,a in enumerate([QQ(n)/119 for n in[150,386,-440,125,-11,46,16,-28,9]]))
r=(xx*xx*yy-xx*yy+yy-1)/(xx*xx*yy-xx);ss=(xx*yy-yy+1)/(xx*yy);b=r*ss*(r-1);c=ss*(r-1);E=EllipticCurve(L,[1-c,-b,-b,0,0]);P=E(0,0);assert 31*P==E(0)
def key9(x,y):
 for h in[x,y,x+y,x-y]:
  if h.minpoly().degree()==9:break
 else:raise AssertionError('nonicprimitive')
 def co(v):a=v.polynomial().list();return list(a)+[0]*(9-len(a))
 A=matrix(k,[co(h**i)for i in range(9)]);AI=A.inverse()
 return(tuple(map(int,h.minpoly().list())),tuple(map(int,vector(k,co(x))*AI)),tuple(map(int,vector(k,co(y))*AI)))
known=set();Q=P
for a in range(1,16):x,y=xy(E,Q);known.add(key9(x,y));Q+=P
assert len(known)==5
seen=set()
for a in map(json.loads,(D/'finite_points62.jsonseq').read_text().splitlines()):
 if a.get('type')!='point'or a['FB_index']not in surv:continue
 L=GF(3**9,'u',modulus=R(a['field_polynomial']));v=L.gen();b=R(a['b'])(v);c=R(a['c'])(v);r=b/c;ss=c*c/(b-c);z=r*ss-2*r+1;x=(ss-r)/z;y=z/(ss*ss-ss-r+1);seen.add(key9(x,y))
assert seen==known
out=dict(status='PASS',target_degree10_places=5685,source_degree10_places=source_total+45,diamond_orbits=379,eligible_class_tests=253,survivors=[],known_nonic31_binding_points=5,
 new_degree5_pair_tests=4095,new_quadratic_fibre_tests=45,seconds=time.monotonic()-tt,
 complete_marked_point_and_2division_census=True,all_diamond_orbits_reconstructed=True,all_primary_sums_and_memberships_verified=True,native_ideal_identities_separately_replayed=True)
(H/'independent_decic_replay.json').write_text(json.dumps(out,indent=2)+chr(10));print('D31A_C62_INDEPENDENT_REPLAY_COMPLETE',out,flush=True)
