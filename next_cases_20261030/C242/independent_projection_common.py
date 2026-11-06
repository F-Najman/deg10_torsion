from pathlib import Path
import ast,hashlib,json,time
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
ROOT=Path('work/deg10/next_cases_20261030/C242');t0=time.monotonic();used={}
def rd(name):
 p=ROOT/name;s=p.read_text();used[name]=hashlib.sha256(s.encode()).hexdigest();return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
R=PolynomialRing(GF(5),'z');fields={}
def fld(d):
 if d not in fields:
  src=rd(f'source{d}.m');fields[d]=GF(5**d,'z',modulus=R(src[1]))
 return fields[d]
def key(vals,d):
 ans=[]
 for _ in range(d):ans.append(tuple(int(a) for v in vals for a in v._vector_()));vals=[v**5 for v in vals]
 return min(ans)
def source(d,ky):
 F=fld(d);z=F.gen();b=R(ky[:d])(z);c=R(ky[d:])(z);E=EllipticCurve(F,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b)
 assert 42*Q==E(0) and 21*Q!=E(0) and 14*Q!=E(0) and 6*Q!=E(0) and P!=21*Q
 return E,P,Q
mods=(2,2,12,36)
def same(u,v):return u[-1]==v[-1]==0 and all((a-b)%m==0 for a,b,m in zip(u,v,mods))
def bc(E,P,Q):
 x=Q[0]-P[0];y=Q[1];return [x**3/y**2,(E.a2()+3*P[0])*x*x/y**2]
def iso_bc(E,P,Q):
 T=14*Q
 def image(U):
  A=U+T;B=U-T;return (U[0]+A[0]+B[0]-2*T[0],U[1]+A[1]+B[1])
 xp,yp=image(P);xq,yq=image(Q);assert yp==0;x=xq-xp
 return [x**3/yq**2,(E.a2()+3*xp)*x*x/yq**2]
reports=[]
for tag in TAGS:
 ds=[4,5,6] if tag=='small' else [10] if tag.startswith('D10s') else [int(tag[1:])]
 ims=rd(f'filter14_images_{tag}.m');cache={tuple(im[0]):im for im in ims};assert len(cache)==len(ims);seen=set();rowscount=0;maps=0
 for d in ds:
  a=rd(f'filter14_D{d}_{tag}.m');assert a[0]==rd(f'source{d}.m')[1]
  for row in a[1]+a[2]:
   dd,i,ky,cls=row;assert dd==d and len(ky)==2*d;E,P,Q=source(d,ky)
   for no,cl in enumerate(cls):
    vals=bc(E,P,3*Q) if no==0 else iso_bc(E,P,Q);kk=(d,)+key(vals,d);assert kk in cache;im=cache[kk];seen.add(kk);maps+=1
    e=min(int(j) for j in divisors(d) if all(v**(5**int(j))==v for v in vals));assert im[1]==e and im[2]==d//e
    assert same(cl,im[3]) and same(im[3],[im[2]*x for x in im[6]])
   if len(cls)==1:assert d>=7 and (cls[0][2]%2 or cls[0][3]%2)
   rowscount+=1
  if d<=6:
   actual=set()
   for row in a[1]:
    E,P,Q=source(d,row[2]);twos=[P,21*Q,P+21*Q];V=22*Q
    orbit={key(bc(E,twos[i],twos[j]+u*V),d) for i in range(3) for j in range(3) if i!=j for u in [1,2,4,5,8,10]}
    assert tuple(row[2])==min(orbit) and actual.isdisjoint(orbit);actual|=orbit
   assert [row[2] for row in a[2]]==[list(v) for v in sorted(actual)] and len(actual)==[126,648,2080][d-4]
 assert seen==set(cache)
 report={'tag':tag,'rows':rowscount,'maps':maps,'unique_images':len(ims),'status':'PASS'};reports.append(report);print('C242_SAGE_PROJECTION_TAG_COMPLETE',report,time.monotonic()-t0,flush=True)
result={'status':'PASS','reports':reports,'seconds':time.monotonic()-t0,'sha256':used};(ROOT/f'independent_projection_{RUN_TAG}.json').write_text(json.dumps(result,indent=2)+'\n');print('C242_SAGE_ALL_PROJECTIONS_COMPLETE',RUN_TAG,flush=True)
