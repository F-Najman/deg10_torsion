from pathlib import Path
import ast,json
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent

def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1;assert all(v[-1]==0 for v in R+hs)
  self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T);self.index=abs(int(self.H.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.H[j,j])
   for i in range(j+1):w[i]-=q*int(self.H[i,j])
  return tuple(w)
def tate(E,Q):
 x,y=Q[0],Q[1];a2,a4=E.a2(),E.a4();aa3=2*y;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4;aa1=2*aa4/aa3;aa2-=(aa4/aa3)**2
 b=-aa2**3/aa3**2;c=-(aa1*aa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def iso2(E,T,Q):
 K=E.base_field();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;u=(x+A+B/x,Q[1]*(1-B/x**2));EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);return EP,EP(u)
pic=rd(P/'cyclic36_proven_picard.m');rels=pic[5];hg=pic[12];cdeg=pic[8];cc=pic[9];G=Quotient(rels,[]);H=Quotient(rels,hg);assert H.index==214272
reconstructed=[[sum(n*cc[j][i] for j,n in enumerate(row)) for i in range(8)] for row in pic[14]];assert all(G.key(a)==G.key(b) for a,b in zip(reconstructed,hg)) and len(reconstructed)==len(hg)
comp=[{(0,)*7}]+[set() for _ in range(7)]
for degree,vec in zip(cdeg,cc):
 if degree>7:continue
 for d in range(degree,8):comp[d].update(H.key([a+b for a,b in zip(list(z)+[0],vec)]) for z in list(comp[d-degree]))
assert [len(x) for x in comp]==[1,7,27,77,181,273,279,279]
allinput=rd(P/'genus9_stage2_arithmetic.m');allcases=allinput[3];byid={r[0]:r for r in allcases};R=PolynomialRing(GF(5),'z');fields={d:GF(5**d,'z',modulus=R(f)) for d,f in allinput[0]}
stages=[];alltested={};final=[]
for label in ['prefix','tail']:
 p=P/('cyclic36_'+label+'_arithmetic.m')
 if not p.exists():continue
 fieldsdata,images,tests,remaining=rd(p);cached={};isogs=0
 for mapno,d,ky,e,w,cl,values in images:
  K=fields[d];z=K.gen();b=R(ky[:d])(z);c=R(ky[d:])(z);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b)
  if mapno==1:EP,U=E,Q
  elif mapno==2:EP,U=E,T+Q
  elif mapno in [3,4]:EP,U=iso2(E,T if mapno==3 else T+18*Q,Q);isogs+=1
  else:raise AssertionError(mapno)
  assert 36*U==EP(0) and 18*U!=EP(0) and 12*U!=EP(0);xy=tate(EP,U);assert xy==[R(v)(z) for v in values]
  de=min(int(j) for j in divisors(d) if all(t**(5**int(j))==t for t in xy));assert de==e and w==d//e
  cached[mapno,d,tuple(ky)]=cl
 expected=[]
 for caseid,cd,checks,yes in tests:
  assert caseid not in alltested;row=byid[caseid];assert cd==10-sum(q[0] for q in row[1]);assert [q[0] for q in checks]==list(range(1,len(checks)+1))
  for mapno,cl,flag in checks:
   summed=[0]*8
   for d,ky in row[1]:summed=[a+b for a,b in zip(summed,cached[mapno,d,tuple(ky)])]
   assert G.key(cl)==G.key(summed) and bool(flag)==(H.key([-a for a in cl]) in comp[cd])
  assert bool(yes)==all(q[2] for q in checks) and (len(checks)==4 or not yes);alltested[caseid]=bool(yes)
  if yes:expected.append(row)
 assert expected==remaining;final+=remaining
 stages.append(dict(stage=label,tested=len(tests),remaining=len(remaining),images=len(images),independent_velu2_images=isogs))
 print('C236_D10_CYCLIC36_SAGE_STAGE_PASS',stages[-1],flush=True)
complete=set(alltested)==set(byid)
out=dict(complete=complete,finite_rational_image_index=H.index,cusp_complement_counts=[len(x) for x in comp],stages=stages,genus9_survivors=len(allcases),tested=len(alltested),remaining=final)
(P/'cyclic36_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_D10_CYCLIC36_SAGE_REPLAY_PASS',json.dumps(out),flush=True)
