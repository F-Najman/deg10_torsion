from pathlib import Path
import ast,json
from sage.all import GF,ZZ,PolynomialRing,EllipticCurve,matrix
P=Path('work/deg10/continuation_20261029b/C226')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=read(P/'filter26.m');fields={};R=PolynomialRing(GF(3),'z')
for d,pol,c in a[3]:fields[d]=GF(3**d,'z',modulus=R(pol))
own={(r[0],r[1]):r for r in a[1]};sm={(r[0],r[1]):r for r in a[2]}
def key(v,d):
 ims=[]
 for _ in range(d):ims.append(tuple(int(c) for z in v for c in z._vector_()));v=[z**3 for z in v]
 return min(ims)
def bc(E,T,Q):x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
def act(row,i,j,k):
 d,idx,ky=row;K=fields[d];z=K.gen();b=sum(K(v)*z**u for u,v in enumerate(ky[:d]));c=sum(K(v)*z**u for u,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);ts=[T,13*Q,T+13*Q];return key(bc(E,ts[i],ts[j]+14*k*Q),d)
out=[]
for desc in [[1,10,54],[2,5,1,5,19]]:
 rows=[own[desc[1],desc[2]]]+([sm[desc[3],desc[4]]] if desc[0]==2 else []);old=sorted(tuple(r[2]) for r in rows);stab=[]
 for i in range(3):
  for j in range(3):
   if i==j:continue
   for k in range(1,7):
    if sorted(act(r,i,j,k) for r in rows)==old:stab.append([i,j,k])
 ims=[]
 for r in rows:
  for kk,e,w,cl in a[4]:
   if kk[0]==1 and kk[1:]==[r[0]]+r[2]:ims.append([e,w])
 out.append({'descriptor':desc,'stabilizers_ts_i_ts_j_14kQ':stab,'first_forgetful_images':ims})
(P/'two_classes_stabilizers.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
