from pathlib import Path
import ast,json,time,math
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve
ROOT=Path('work/deg10/next_cases_20261030/C242');t0=time.monotonic()
def rd(path):return ast.literal_eval(path.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=rd(ROOT/'geometry/cyclic21_full_arithmetic.m');fields,remaining=a[0],a[3];R=PolynomialRing(GF(5),'z');fs={d:GF(5**d,'z',modulus=R(pol)) for d,pol in fields}
def key(v,d):
 out=[]
 for _ in range(d):out.append(tuple(int(x) for w in v for x in w._vector_()));v=[w**5 for w in v]
 return min(out)
autos=[]
for ip in range(3):
 for iq in range(3):
  if ip==iq:continue
  perm=[ip,iq,3-ip-iq];p=list(range(3));po=0
  while True:
   p=[perm[i] for i in p];po+=1
   if p==list(range(3)):break
  for u in [1,2,4,5,8,10]:
   uo=next(n for n in range(1,7) if pow(u,n,21) in [1,20]);autos.append((ip,iq,u,math.lcm(po,uo)))
assert len(autos)==36
# Exact geometric cusps; this also pins the level-action convention.
labels={(a,b,c,d) for a in range(42) for b in range(42) if math.gcd(math.gcd(a,b),42)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orb(v):
 a,b,c,d=v;return frozenset(((e*(a+k*b))%42,e*b%42,(c+k*d)%2,d) for e in [-1,1] for k in range(42))
cs=sorted(set(map(orb,labels)),key=lambda x:sorted(x));ix={v:i for i,C in enumerate(cs) for v in C};assert len(cs)==72
def fro(v):a,b,c,d=v;return (5*a%42,b,5*c%2,d)
froperm=[ix[fro(next(iter(C)))] for C in cs];rational=[i for i in range(72) if froperm[i]==i];assert len(rational)==18
def hlabel(v,h):
 a,b,c,d=v;ip,iq,u,o=h;ts=[(c,d),(a%2,b%2),((c+a)%2,(d+b)%2)];s,t=ts[iq];p,q=ts[ip];return ((21*s+22*u*a)%42,(21*t+22*u*b)%42,p,q)
cuspaction={}
for h in autos:
 pp=[ix[hlabel(next(iter(C)),h)] for C in cs]
 assert sorted(pp)==list(range(72)) and all(pp[froperm[i]]==froperm[pp[i]] for i in range(72))
 assert all(all(ix[hlabel(v,h)]==pp[i] for v in C) for i,C in enumerate(cs))
 cuspaction[h]=[i for i in rational if pp[i]==i]
curves={};cache={}
def image(term,h):
 d,ky=term;kk=(d,tuple(ky))
 if kk not in curves:
  F=fs[d];z=F.gen();b=R(ky[:d])(z);c=R(ky[d:])(z);E=EllipticCurve(F,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);assert 42*Q==E(0) and 21*Q!=E(0) and 14*Q!=E(0) and 6*Q!=E(0);curves[kk]=(E,c,[P,21*Q,P+21*Q],22*Q)
 ck=(kk,h)
 if ck not in cache:
  E,c,ts,V=curves[kk];ip,iq,u,o=h;P=ts[ip];Q=ts[iq]+u*V;x=Q[0]-P[0];y=Q[1];cache[ck]=(d,key([x**3/y**2,(c+3*P[0])*x*x/y**2],d))
 return cache[ck]
reports=[];left=[];removed=[]
for ci,terms in remaining:
 deg=sum(d for d,ky in terms);cd=10-deg;orig=Counter((d,tuple(ky)) for d,ky in terms);stabs=[h for h in autos if Counter(image(t,h) for t in terms)==orig];eligible=[]
 if cd<=1:
  for h in stabs:
   if h[3]%3==0:eligible.append({'h':h,'reason':'3-primary stabilizer; gonality and no low-degree CM j0 points'})
   elif cd==1 and not cuspaction[h]:eligible.append({'h':h,'reason':'stabilizer has no rational fixed cusp'})
 report={'case_id':ci,'noncuspidal_degrees':[d for d,ky in terms],'cusp_degree':cd,'stabilizers':[{'h':h,'fixed_rational_cusps':cuspaction[h]} for h in stabs],'cheap_exclusion':eligible};reports.append(report)
 (removed if eligible else left).append([ci,terms])
def unified_text(fields,cases):
 fs='[*'+','.join('[*'+str(d)+','+repr(pol)+'*]' for d,pol in fields)+'*]'
 cc='[*'+','.join('[*'+str(ci)+', [*'+','.join('[*'+str(d)+','+repr(ky)+'*]' for d,ky in terms)+'*]*]' for ci,terms in cases)+'*]'
 return '[*'+fs+','+cc+'*]\n'
(ROOT/'after_cheap_stabilizers_unified.m').write_text(unified_text(fields,left));out={'status':'PASS','input_cases':len(remaining),'removed':len(removed),'remaining':len(left),'removed_case_ids':[r[0] for r in removed],'remaining_case_ids':[r[0] for r in left],'geometric_cusps':72,'rational_cusps_at5':18,'reports':reports,'seconds':time.monotonic()-t0};(ROOT/'cheap_stabilizers.json').write_text(json.dumps(out,indent=2)+'\n');print('C242_CHEAP_STABILIZER_CERTIFICATE',len(remaining),len(removed),len(left),[r[0] for r in removed],time.monotonic()-t0,flush=True)
