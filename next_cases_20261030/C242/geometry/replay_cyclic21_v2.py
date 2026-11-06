from pathlib import Path
import ast,json
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
def tate(E,Q):
 x,y=Q[0],Q[1];a2,a4=E.a2(),E.a4();aa3=2*y;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4;aa1=2*aa4/aa3;aa2-=(aa4/aa3)**2;b=-aa2**3/aa3**2;c=-(aa1*aa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def iso2(E,T,Q):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;u=(x+A+B/x,Q[1]*(1-B/x**2));EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);return EP,EP(u)
pic=rd(P/'cyclic21_proven_picard.m');rels=pic[5];hg=pic[12];cc=pic[9];G=Quotient(rels,[]);H=Quotient(rels,hg);assert G.index==2184 and H.index==6
reconstructed=[[sum(n*cc[j][i]for j,n in enumerate(row))for i in range(2)]for row in pic[14]];assert all(G.key(a)==G.key(b)for a,b in zip(reconstructed,hg))
comp=[{(0,)}]+[set()for _ in range(6)]
for degree,v in zip(pic[8],cc):
 if degree>6:continue
 for d in range(degree,7):comp[d].update(H.key([a+b for a,b in zip(list(z)+[0],v)])for z in list(comp[d-degree]))
assert [len(x)for x in comp]==[1,1,1,2,2,2,2]
input=rd(P/'genus5_merged_arithmetic.m');byid={r[0]:r for r in input[3]};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(f))for d,f in input[0]};data=rd(P/'cyclic21_full_arithmetic.m');cache={};isogs=0
for mapno,d,key,e,w,cl,vals in data[1]:
 K=fields[d];z=K.gen();b=R(key[:d])(z);c=R(key[d:])(z);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b)
 if mapno==1:EP,U=E,2*Q
 else:EP,U=iso2(E,P2 if mapno==2 else(21*Q if mapno==3 else P2+21*Q),2*Q);isogs+=1
 assert 21*U==EP(0)and 7*U!=EP(0)and 3*U!=EP(0);xy=tate(EP,U);assert xy==[R(v)(z)for v in vals];de=min(int(j)for j in divisors(d)if all(v**(5**int(j))==v for v in xy));assert de==e and w==d//e;cache[mapno,d,tuple(key)]=cl
expected=[];seen=set()
for id,cd,checks,yes in data[2]:
 assert id not in seen;seen.add(id);row=byid[id];assert cd==10-sum(p[0]for p in row[1])and[q[0]for q in checks]==list(range(1,len(checks)+1))
 for mapno,cl,flag in checks:
  summed=[0]*2
  for d,key in row[1]:summed=[a+b for a,b in zip(summed,cache[mapno,d,tuple(key)])]
  assert G.key(cl)==G.key(summed)and bool(flag)==(H.key([-a for a in cl])in comp[cd])
 assert bool(yes)==all(r[2]for r in checks)and(len(checks)==4 or not yes)
 if yes:expected.append(row)
assert seen==set(byid)and expected==data[3]
out=dict(complete=True,degree=10,genus=5,finite_picard_order=2184,rational_image_index=6,input_supports=len(byid),remaining=len(expected),projection_images=len(cache),independent_velu2_images=isogs,cusp_complement_counts=[len(x)for x in comp]);(P/'cyclic21_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C242_D10_CYCLIC21_INDEPENDENT_REPLAY_PASS',json.dumps(out),flush=True)
