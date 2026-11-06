from pathlib import Path
import ast,json,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
P=Path('work/deg10/C232_resolution_20261030');t0=time.monotonic()
def read(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def key(v,d):
 out=[]
 for _ in range(d):out.append(tuple(int(c) for x in v for c in x._vector_()));v=[x**3 for x in v]
 return min(out)
def bc(E,T,Q):
 x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
def tate(E,Q):
 x,y=Q[0],Q[1];aa3=2*y;aa2=3*x+E.a2();aa4=3*x*x+2*x*E.a2()+E.a4();a1=2*aa4/aa3;a2=aa2-(aa4/aa3)**2;b=-a2**3/aa3**2;c=-(a1*a2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return [(s-r)/t,t/(s*s-s-r+1)]
def iso2(E,T,Q):
 K=E.base_ring();r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();x=Q[0]-r;y=Q[1];EG=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);U=EG(x+A+B/x,y*(1-B/x**2));return EG,U
fp,cd,nc=read('frobenius.m');tr=[];closed=[]
for d in range(1,11):
 a=-d*fp[74-d]-sum(fp[74-j]*tr[d-j-1] for j in range(1,d));tr.append(a);b=3**d+1-a-sum(e*closed[e-1] for e in range(1,d) if d%e==0);assert b%d==0;closed.append(b//d)
assert [closed[d-1]-Counter(cd)[d] for d in range(1,11)]==nc
T=read('target_proven_picard.m');rels=T[5];hg=T[8];NN=len(rels[0])-1;M=hermite_normal_form(Matrix([r[:-1] for r in rels+hg]).T);assert abs(int(M.det()))==68

def qkey(v):
 assert v[-1]==0;w=list(v[:-1])
 for j in range(NN-1,-1,-1):
  q=w[j]//int(M[j,j])
  for i in range(j+1):w[i]-=q*int(M[i,j])
 return tuple(w)
def add(u,v):return qkey([a+b for a,b in zip(u,v)]+[0])
allowed=[set() for _ in range(5)];allowed[0]={(0,)*NN}
for w,v in zip(T[9],T[10]):
 for d in range(w,5):allowed[d]|={add(x,qkey(v)) for x in allowed[d-w]}
assert list(map(len,allowed))==[1,1,3,3,4]
filters=[read('filter'+x+'.m') for x in ['probe6','789','10']];tests={};maps=0;isogs=0;left=[]
for images,rows,remaining in filters:
 for row in rows:assert (row[0],row[1]) not in tests;tests[row[0],row[1]]=(row,images)
 left+=remaining
R=PolynomialRing(GF(3),'x');coverage=[]
for d in range(6,11):
 s=read('source'+str(d)+'.m');assert s[0]==d and s[2];K=GF(3**d,'z',modulus=R(s[1]));z=K.gen();seen=set();total=Counter();n=0
 for i,(dd,ky) in enumerate(s[6],1):
  b=sum(K(v)*z**j for j,v in enumerate(ky[:d]));c=sum(K(v)*z**j for j,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);assert 2*P==E(0) and 32*Q==E(0) and 16*Q!=E(0) and P!=16*Q
  actual=min(int(e) for e in divisors(d) if b**(3**int(e))==b and c**(3**int(e))==c);assert actual==dd
  ims={key(bc(E,P+bb*(16*Q),cc*P+aa*Q),d) for bb in range(2) for cc in range(2) for aa in [1,3,5,7,9,11,13,15]};assert tuple(ky)==min(ims) and seen.isdisjoint(ims);seen.update(ims);total[dd]+=len(ims);n+=dd==d
  row,images=tests[d,i];assert row[:3]==[d,i,ky];flags=[]
  variants=[(E,Q),(E,P+Q),iso2(E,P,Q),iso2(E,P+16*Q,Q)]
  for var,ti,e,mult,cl,yes in row[3]:
   EG,U=variants[var-1];assert 32*U==EG(0) and 16*U!=EG(0);xy=tate(EG,U);im=images[ti-1];assert im[0]==[d]+list(key(xy,d)) and im[1]==e and im[2]==mult and im[3]==cl
   ee=min(int(v) for v in divisors(d) if all(x**(3**int(v))==x for x in xy));assert e==ee and mult*e==d;assert bool(yes)==(qkey([-v for v in cl]) in allowed[10-d]);flags.append(bool(yes));maps+=1;isogs+=var>=3
  assert bool(row[4])==all(flags)
 assert [total[int(e)] for e in divisors(d)]==s[5] and total[d]==nc[d-1];coverage.append({'degree':d,'closed_points':total[d],'full_automorphism_orbits':n});print('C232_INDEPENDENT_SIEVE_DEGREE',d,total[d],n,time.monotonic()-t0,flush=True)
assert len(tests)==305 and len(left)==1 and left[0][:2]==[6,8]
out={'status':'PASS','classification_step':'Necessary four-map test; one residual orbit remains','closed_non_cuspidal_counts':nc,'coverage':coverage,'tested_orbits':len(tests),'checked_maps':maps,'independent_two_isogenies':isogs,'target_rational_image_index':68,'allowed_cusp_complement_counts':list(map(len,allowed)),'residual_supports':left}
(P/'independent_sieve.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_INDEPENDENT_SIEVE_COMPLETE',maps,isogs,len(left),time.monotonic()-t0,flush=True)
