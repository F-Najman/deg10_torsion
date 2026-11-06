from sage.all import *
from pathlib import Path
from collections import Counter,defaultdict
import ast,json,re,hashlib,time
C=Path(__file__).resolve().parents[1]/'continuation_20261031b/C66';rd=lambda n:ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'));f=rd('F2_proven_frame.m');nc=dict(rd('F2_NC_classes.m'));m=ast.parse((C/'inputs/class_arithmetic33.py').read_text());old={n.targets[0].id:ast.literal_eval(n.value)for n in m.body if isinstance(n,ast.Assign)};assert f[4]==old['rational_cusp_generators'];N=42373650;k=GF(2);R=PolynomialRing(k,'z');start=time.monotonic()
def xy(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;c1=2*b4/b3+a1;c2=b2-(b4/b3)**2-a1*b4/b3;b=-c2**3/b3**2;c=-(c1*c2-b3)/b3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return((s-r)/t,t/(s*s-s-r+1))
def vf(ef,x,y):
 def rat(q):return sum(y**i*k(a)for i,a in enumerate(q[0]))/sum(y**i*k(a)for i,a in enumerate(q[1]))
 return sum(rat(v)*x**i for i,v in enumerate(ef))
actual={};SS=Counter();checked=0
for d,rows in nc.items():
 seen=set()
 for pd,cl,ch,pol,vals in rows:
  L=GF(2**d,name='a',modulus=R(pol));a=L.gen();cv=lambda v:sum(L(c)*a**i for i,c in enumerate(v));x,y=map(cv,vals);r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1);E=EllipticCurve(L,[1-c,-b,-b,0,0]);P=E(0,0);assert P.order()==33 and xy(E,P)==(x,y);assert sum(vf([coef],L(1),y)*x**i for i,coef in enumerate(f[0]))==0;assert all(vf(g,x,y)==0 for g in pd[2]);assert pd[0]==d and pd[1]==1
  orb=[(x**(2**i),y**(2**i))for i in range(d)];assert len(set(orb))==d;key=min(tuple(tuple(int(q)for q in v.polynomial().list())for v in p)for p in orb);assert key not in seen;seen.add(key);assert cl[2]==0
  if ch==[0] and d==10:SS['supersingular'if E.j_invariant()==0 else'ordinary']+=1
  checked+=1
 actual[d]=len(seen)
coef=ast.literal_eval(re.search(r'ExpectedCoefficients:=(.*?);',(C/'inputs/frobenius2_input33.m').read_text())[1]);F=R(coef);assert F.degree()==42 and F(1)==15*N;tr=F.reverse().power_series_ring(11)(F.reverse()).log() if False else None
ps=[];closed=[]
for n in range(1,11):
 s=-n*int(F[42-n])-sum(int(F[42-j])*ps[n-j-1]for j in range(1,n));ps.append(s);pt=2**n+1-s;v=(pt-sum(d*closed[d-1]for d in divisors(n)if d<n))//n;closed.append(v)
assert all(closed[d-1]==actual[d]+sum(z[0]==d for z in f[6])for d in range(1,11))
chi=[z[0]for z in f[5]];pack=lambda a,b:(int(a)%15,int(b)%N);add=lambda a,b:pack(a[0]+b[0],a[1]+b[1]);ch=lambda z:(chi[0]*z[0]+chi[1]*z[1])%3
cp=[Counter()for _ in range(11)];cp[0][(0,0)]=1
for d,sd,cl,cc in f[6]:
 z=pack(*cl[:2]);nxt=[Counter(c)for c in cp]
 for t in range(11):
  zz=z
  for m in range(1,(10-t)//d+1):
   for v,mul in cp[t].items():nxt[t+m*d][add(v,zz)]+=mul
   zz=add(zz,z)
 cp=nxt
nccomp={d:Counter(pack(*r[1][:2])for r in rows)for d,rows in nc.items()};tot=Counter(cp[10])
for d in range(5,11):
 for z,m in nccomp[d].items():
  for c,n in cp[10-d].items():tot[add(z,c)]+=m*n
for i,p in enumerate(nc[5]):
 for q in nc[5][i:]:tot[add(pack(*p[1][:2]),pack(*q[1][:2]))]+=1
H={z:n for z,n in tot.items()if ch(z)==0};h=Counter();gon={}
for z,n in H.items():assert(n+1)&n==0;h[int((n+1).bit_length()-1)]+=1
for z,n in H.items():
 if n>1:gon[z]=n
rep=json.loads((C/'F2_degree10_classification.json').read_text());assert len(H)==rep['rational_image_classes']==321210 and sum(H.values())==321230;assert h=={1:321200,2:10};assert gon=={tuple(z['class'][:2]):z['finite_effective_divisors']for z in rep['gonal_classes']}
known7={tuple(z[:2])for z in old['known_septic_plus_degree2cusp']};known8={tuple(z[:2])for z in old['known_octic_plus_cusp']};known9={tuple(z[:2])for z in old['known_nonic_classes']};signed=json.loads((C/'inputs/short_cusp_classes33.json').read_text());unk={tuple(z['finite_class'][:2])for z in signed['representatives']}
assert all(tuple(r[1][:2])in known7 for r in nc[7]if r[2]==[0]);assert all(tuple(r[1][:2])in known8 for r in nc[8]if r[2]==[0]);assert {tuple(r[1][:2])for r in nc[9]if r[2]==[0]}==known9|unk and known9.isdisjoint(unk)
out={'status':'PASS','native_points_bound_to_actual_33torsion':checked,'noncuspidal_counts':actual,'Frobenius_closed_counts':closed,'independent_rational_image_classes':len(H),'independent_effective_divisors':sum(H.values()),'independent_h0_histogram':dict(h),'independent_gonal_classes':[list(z)+[0]for z in sorted(gon)],'NC10_rational_image_ordinary_supersingular':dict(SS),'global_witness_Picard_basis_unchanged':True,'all_known_mixed_branches_bound':True,'algorithm':'generating-function convolution over each closed cusp, independent of recursive divisor enumeration','elapsed_seconds':time.monotonic()-start};(C/'F2_independent_verified.json').write_text(json.dumps(out,indent=2,default=int)+'\n');print('D31B_C66_F2_INDEPENDENT_COMPLETE',checked,len(H),dict(SS),flush=True)
