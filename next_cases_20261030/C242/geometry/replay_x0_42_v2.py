from pathlib import Path
import ast,json
from math import gcd
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent

def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;assert all(v[-1]==0 for v in R+hs);self.H=hermite_normal_form(Matrix([r[:-1]for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
def normpoint(E,S):
 x,y=S[0],S[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2=3*x+E.a2()-m*m;a3=2*y;return -a2**3/a3**2,1-2*m*a2/a3,a3/a2
def project(E,Q):
 b6,c6,k6=normpoint(E,7*Q);b7,c7,k7=normpoint(E,6*Q);assert b6==c6+c6*c6;v=b7/c7;assert b7==v*v*(v-1) and c7==v*(v-1);t=(v**3-3*v+1)/(v*(v-1))-8;ratio=(k7*(v*v-v+1)/k6)**2;return[c6,t,ratio]
def iso2(E,T,Q):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;u=(x+A+B/x,Q[1]*(1-B/x**2));EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);return EP,EP(u)
pic=rd(P/'x0_42_proven_picard_v2.m');rels=pic[5];hg=pic[12];cc=pic[9];G=Quotient(rels,[]);H=Quotient(rels,hg);assert G.index==18432 and H.index==192
# Independent squarefree-level Tate labels and cusp widths.
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
n=42;labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1};left=set(labels);cusps=[]
while left:
 a,b=min(left);orb={(ss*(a+k*b)%n,ss*b%n) for ss in [-1,1] for k in range(n)};left-=orb;cusps.append(sorted(orb))
assert len(cusps)==48;ix={x:i for i,C in enumerate(cusps) for x in C};perms=[[ix[(u*C[0][0]%n,C[0][1])] for C in cusps] for u in range(1,n) if gcd(u,n)==1];frob=[ix[(5*C[0][0]%n,C[0][1])] for C in cusps];global_orbits=orbits(perms,48);finite_orbits=orbits([frob],48);assert Counter(map(len,global_orbits))=={1:12,2:6,6:4} and Counter(map(len,finite_orbits))=={1:12,2:6,3:4,6:2}
def cuspkey(g):
 g6=gcd(g,6);e=1 if g%7==0 else 7;zkeys={1:[e,0,0],2:[0,e,0],3:[-e,-e,-e],6:[0,0,e]};return tuple(zkeys[g6]+[(6//g6)*(1 if g%7==0 else -1)])
def orbkey(o):
 gs={gcd(cusps[i][0][1],n) for i in o};assert len(gs)==1;return cuspkey(gs.pop())
aud=rd(P/'x0_42_direct_cusp_pushforward.m');assert aud[7] is True and aud[4]==[6]*8
assert Counter((len(o),orbkey(o)) for o in global_orbits)==Counter((r[1],tuple(aud[0][r[2]-1])) for r in aud[2]);assert Counter((len(o),orbkey(o)) for o in finite_orbits)==Counter((r[1],tuple(r[3])) for r in aud[1])
hg_direct=[[r[1]*a for a in cc[r[2]-1]] for r in aud[2]];assert Quotient(rels,hg_direct).H==H.H
comp=[{(0,)*4}]+[set()for _ in range(6)]
for degree,v in zip(pic[8],cc):
 if degree>6:continue
 for d in range(degree,7):comp[d].update(H.key([a+b for a,b in zip(list(z)+[0],v)])for z in list(comp[d-degree]))
assert [len(x)for x in comp]==[1,6,12,16,20,24,24]
input=rd(P.parent/'after_cheap_stabilizers_unified.m');byid={r[0]:r for r in input[1]};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(f))for d,f in input[0]};data=rd(P/'x0_42_full_arithmetic.m');cache={};isogs=0
for mapno,d,key,e,w,cl,vals in data[1]:
 K=fields[d];z=K.gen();b=R(key[:d])(z);c=R(key[d:])(z);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b)
 if mapno==1:EP,U=E,Q
 elif mapno==2:EP,U=E,P2+Q
 elif mapno==3:EP,U=E,P2+22*Q
 elif mapno in [4,5]:EP,U=iso2(E,P2 if mapno==4 else P2+21*Q,Q);isogs+=1
 else:
  assert mapno==6;EP,U1=iso2(E,21*Q,P2);EP2,U2=iso2(E,21*Q,Q);assert EP==EP2;U=U1+U2;isogs+=1
 assert 42*U==EP(0) and 21*U!=EP(0) and 14*U!=EP(0) and 6*U!=EP(0);xy=project(EP,U);assert xy==[R(v)(z) for v in vals];de=min(int(j) for j in divisors(d) if all(v**(5**int(j))==v for v in xy));assert de==e and w==d//e;cache[mapno,d,tuple(key)]=cl
expected=[];seen=set()
for id,cd,checks,yes in data[2]:
 assert id not in seen;seen.add(id);row=byid[id];assert cd==10-sum(p[0]for p in row[1])and[q[0]for q in checks]==list(range(1,len(checks)+1))
 for mapno,cl,flag in checks:
  summed=[0]*5
  for d,key in row[1]:summed=[a+b for a,b in zip(summed,cache[mapno,d,tuple(key)])]
  assert G.key(cl)==G.key(summed)and bool(flag)==(H.key([-a for a in cl])in comp[cd])
 assert bool(yes)==all(r[2]for r in checks)and(len(checks)==6 or not yes)
 if yes:expected.append(row)
assert seen==set(byid)and expected==data[3]
out=dict(complete=True,degree=10,genus=5,target='X0(42)',finite_picard_order=G.index,source_rational_image_order=G.index//H.index,source_rational_image_index=H.index,global_cusp_orbit_norms_independently_checked=True,finite_cusp_orbits_checked=len(finite_orbits),global_cusp_orbits_checked=len(global_orbits),input_supports=len(byid),remaining=len(expected),projection_images=len(cache),independent_velu2_images=isogs,cusp_complement_counts=[len(x) for x in comp]);(P/'x0_42_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C242_D10_X0_42_INDEPENDENT_REPLAY_PASS',json.dumps(out),flush=True)
