from pathlib import Path
from math import gcd
from collections import Counter
import ast,json
from sage.all import GF,PolynomialRing,EllipticCurve,FunctionField,QQ
P=Path(__file__).resolve().parent

def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
packet=rd(P/'case26_section_model.m');records,fail=rd(P/'case26_section_RR.m');assert not fail and len(records)==24
R5=PolynomialRing(GF(5),'a');L=GF(5**6,'a',modulus=R5(packet[9][0]));a=L.gen();key=packet[9][1];b=R5(key[:6])(a);c=R5(key[6:])(a);E=EllipticCurve(L,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert 42*Q==E(0) and 21*Q!=E(0) and 14*Q!=E(0) and 6*Q!=E(0) and P2!=21*Q

def bc(E,P,Q):
 x=Q[0]-P[0];y=Q[1];return[x**3/y**2,(E.a2()+3*P[0])*x*x/y**2]
H=bc(E,P2+21*Q,P2+4*Q);assert H==[b**25,c**25]
HM=bc(E,21*Q,P2+37*Q);assert HM==[b**625,c**625]
PP,QQ0=P2,Q
for _ in range(3):PP,QQ0=PP+21*QQ0,PP+4*QQ0
assert PP==P2 and QQ0==Q
S=6*Q;x,y=S[0],S[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2p=3*x+E.a2()-m*m;a3p=2*y;kap=a3p/a2p;b7=-a2p**3/a3p**2;c7=1-2*m*a2p/a3p;v=b7/c7;R2=21*Q;Z2=P2+R2;u=4*(P2[0]-R2[0])*(P2[0]-Z2[0])*(R2[0]-Z2[0])/(kap**6*v**4*(v-1)**4)
assert u==R5(packet[9][3])(a) and list(u.minpoly())==packet[9][2] and packet[9][2]==[2,4,1]
s=u*u+8;assert v**3-s*v*v+(s-3)*v+1==0
# Exact universal cyclic cubic/Kummer identity over both Q and F25.
RQ=PolynomialRing(QQ,'s');ss=RQ.gen();DD=ss*ss-3*ss+9;AA=2*ss**3-9*ss*ss+27*ss-27;assert AA*AA-4*DD**3==-27*DD*DD
k=GF(25,'theta',modulus=R5(packet[0]));A=FunctionField(k,'s');s=A.gen();R=PolynomialRing(A,'v');vvar=R.gen();p=vvar**3-s*vvar**2+(s-3)*vvar+1;K=A.extension(p,'v');vv=K.gen();v1=(vv-1)/vv;v2=1/(1-vv);D=s*s-3*s+9;C=2*s**3-9*s*s+27*s-27
zpoly=PolynomialRing(k,'z');z=zpoly.gen();roots=(z*z+z+1).roots(multiplicities=False);assert len(roots)==2
for zeta in roots:
 w=vv+zeta*v1+zeta*zeta*v2;wp=vv+zeta*zeta*v1+zeta*v2;f=(C+3*(zeta-zeta*zeta)*D)/2;assert w*wp==D and w+wp==3*vv-s
 if w**3!=f:w=wp
 assert w**3==f and 3*vv-s==w+D/w
# Complete full-level Tate cusp action and all compatible C4 orbit norms.
n=42;labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if(c*b-d*a)%2==1};left=set(labels);cusps=[]
while left:
 a,b,c,d=min(left);orb={(sg*(a+j*b)%n,sg*b%n,sg*(c+j*d)%2,sg*d%2) for sg in [-1,1] for j in range(n)};left-=orb;cusps.append(sorted(orb))
assert len(cusps)==72;ix={x:i for i,C in enumerate(cusps) for x in C};hp=[];fp=[]
for C0 in cusps:
 a,b,c,d=C0[0];hp.append(ix[((4*a+21*c)%42,(4*b+21*d)%42,(c+a)%2,(d+b)%2)]);fp.append(ix[(5*a%42,b,c,d)])
assert all(hp[hp[hp[i]]]==i and hp[i]!=i and hp[hp[i]]!=i for i in range(72));assert all(fp[hp[i]]==hp[fp[i]] for i in range(72))
def orbits(perms,n):
 unseen=set(range(n));out=[]
 while unseen:
  a=min(unseen);todo=[a];seen={a}
  while todo:
   x=todo.pop()
   for f in perms:
    y=f[x]
    if y not in seen:seen.add(y);todo.append(y)
  unseen-=seen;out.append(sorted(seen))
 return out
source=orbits([fp],72);assert Counter(map(len,source))=={1:18,2:9,3:6,6:3};globalperms=[[ix[(g*C[0][0]%42,C[0][1],g*C[0][2]%2,C[0][3])] for C in cusps] for g in range(1,42) if gcd(g,42)==1];assert Counter(map(len,orbits(globalperms,72)))=={1:18,2:9,6:6}
HC=orbits([hp],72);qi={x:i for i,C in enumerate(HC) for x in C};assert len(HC)==24;fq=[qi[fp[C[0]]] for C in HC];quot=orbits([fq],24);assert Counter(map(len,quot))=={1:6,2:3,3:2,6:1};closed={x:i for i,C in enumerate(quot) for x in C}
rat=[C for C in source if len(C)==1];cub=[C for C in source if len(C)==3];types=Counter();omegas=set()
for C3 in cub:
 for C1 in rat:
  q3={closed[qi[x]] for x in C3};q1={closed[qi[x]] for x in C1};assert len(q3)==len(q1)==1;assert len(quot[next(iter(q3))])==3 and len(quot[next(iter(q1))])==1;types[next(iter(q3)),next(iter(q1))]+=1
  vec=[0]*72
  for x in C3+C1:
   for y in[x,hp[x],hp[hp[x]]]:vec[y]+=1
  assert sum(vec)==12 and all(vec[hp[i]]==vec[i] for i in range(72));omegas.add(tuple(vec))
assert len(types)==len(omegas)==12 and set(types.values())=={9}
def deg(D):return sum(ep[1]*m for ep,m in D)
assert [deg(D) for D in packet[5]]==[0,-6,-6] and all(deg(D)==1 for D in packet[6]) and all(deg(D)==3 for D in packet[7]);assert len(packet[8])==2 and all(deg(r[1])==2 for r in packet[8])
for qid,ci,ri,dims,proof in records:
 assert dims==[1,0,0] and len(proof)==3
 for j,pr in enumerate(proof):assert pr[0]==dims[j] and deg(pr[1])==[6,0,0][j] and len(pr[2])==len(pr[3])==8 and len(pr[4])==dims[j]
out=dict(complete=True,degree=10,case=26,chosen_h='(P+21Q,P+4Q)',h_order=3,h_Frobenius5_phase=2,accepted_model_generator_Frobenius5_phase=4,source_cusp_signature={'1':18,'2':9,'3':6,'6':3},all_72_geometric_cusps_free=True,quotient_cusp_signature={'1':6,'2':3,'3':2,'6':1},compatible_source_C4=108,omega_types=12,source_options_per_type=9,Kummer_resolvent_identity_checked=True,quadratic_quotient_containing_places=2,character_RR_patterns=24,character_RR_spaces=72,all_character_dimensions=[1,0,0],absolute_rational_base_module_rank=8)
(P/'case26_section_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C242_CASE26_SOURCE_KUMMER_CUSPS_SECTION_PACKET_REPLAY_PASS',json.dumps(out),flush=True)
