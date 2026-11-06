from sage.all import *
from pathlib import Path
import ast,json,re,time
T=time.monotonic();R=Path('work/deg10/continuation_20261031b');Fp=GF(3);PR=PolynomialRing(Fp,'z')
def rd(p):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',p.read_text().replace('[*','[').replace('*]',']'))))
a=rd(R/'M4N28/fourmap_Y28.m');u=json.loads((R/'Y228_shared/union_normalization.json').read_text());bind={i:c for label,i,c in u['bindings']if label=='M4N28'};fields={d:GF(3**d,'a',modulus=PR(pol))for d,pol in a[0]}
def key(E,P,Q,d):
 r=P[0];x=Q[0]-r;y=Q[1];b=x**3/y**2;c=(E.a2()+3*r)*x*x/y**2;return min(tuple(int(z.polynomial()[i])for z in[b**(3**j),c**(3**j)]for i in range(d))for j in range(d))
checks=[]
for idx,ids,cusp in a[2]:
 rec=a[1][idx-1];assert rec[0]==idx and ids==[idx];_,d,mul,ky,*_=rec;K=fields[d];z=K.gen();b=PR(ky[:d])(z);c=PR(ky[d:])(z);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);stack=[(P,Q)];seen={}
 while stack:
  P,Q=stack.pop();tag=key(E,P,Q,d)
  if tag in seen:continue
  seen[tag]=(P,Q);stack.extend([(P+14*Q,Q),(P,Q+P),(P,3*Q),(P,5*Q)])
 canon=min(seen);uc=u['cases'][bind[idx]-1];assert uc==[bind[idx],cusp,[[d,mul,list(canon)]]];checks.append({'source_case':idx,'union_case':bind[idx],'orbit_size':len(seen),'canonical_key':list(canon)})
out={'status':'PASS','cases':len(checks),'all_simultaneous_orbits_independently_generated':True,'method':'BFS on marked elliptic points using four generating transformations; no 24-element enumeration reused','checks':checks,'seconds':time.monotonic()-T};(R/'M4N28/parent_independent_union.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_M428_INDEPENDENT_UNION_COMPLETE',len(checks),time.monotonic()-T,flush=True)
