from pathlib import Path
import ast,json,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
SAVED=Path('work/deg10/C234_resolution_20261030');t0=time.monotonic()
def read(n):return ast.literal_eval((SAVED/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
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
 a=-d*fp[98-d]-sum(fp[98-j]*tr[d-j-1] for j in range(1,d));tr.append(a);b=3**d+1-a-sum(e*closed[e-1] for e in range(1,d) if d%e==0);assert b%d==0;closed.append(b//d)
assert [closed[d-1]-Counter(cd)[d] for d in range(1,11)]==nc
T=read('target_proven_picard.m');rels=T[5];hg=T[8];NN=len(rels[0])-1;M=hermite_normal_form(Matrix([r[:-1] for r in rels+hg]).T);assert abs(int(M.det()))==128

def qkey(v):
 assert v[-1]==0;w=list(v[:-1])
 for j in range(NN-1,-1,-1):
  q=w[j]//int(M[j,j])
  for i in range(j+1):w[i]-=q*int(M[i,j])
 return tuple(w)
def add(u,v):return qkey([a+b for a,b in zip(u,v)]+[0])
assert all(qkey(v)==(0,)*NN for v in T[10])
GMat=hermite_normal_form(Matrix([r[:-1] for r in rels]).T)
def gkey(v):
 assert v[-1]==0;w=list(v[:-1])
 for j in range(NN-1,-1,-1):
  q=w[j]//int(GMat[j,j])
  for i in range(j+1):w[i]-=q*int(GMat[i,j])
 return tuple(w)
filters=[read('filter'+str(i)+'.m') for i in range(3)];R=PolynomialRing(GF(3),'x');sources={};small={};fields={};coverage=[]
for d in range(4,11):
 src=read('source'+str(d)+'.m');K=GF(3**d,'z',modulus=R(src[1]));z=K.gen();fields[d]=K;sources[d]=[];seen=set();ownkeys=set();total=Counter()
 for dd,ky in src[6]:
  b=sum(K(v)*z**i for i,v in enumerate(ky[:d]));c=sum(K(v)*z**i for i,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);assert 34*Q==E(0) and 17*Q!=E(0) and 2*Q!=E(0) and P!=17*Q
  actual=min(int(e) for e in divisors(d) if b**(3**int(e))==b and c**(3**int(e))==c);assert actual==dd;ts=[P,17*Q,P+17*Q];V=18*Q;ims={key(bc(E,ts[i],ts[j]+a*V),d) for i in range(3) for j in range(3) if i!=j for a in range(1,9)};assert tuple(ky)==min(ims) and seen.isdisjoint(ims);seen.update(ims);total[dd]+=len(ims)
  if dd==d:sources[d].append([d,len(sources[d])+1,ky]);ownkeys.update(ims)
 assert [total[int(e)] for e in divisors(d)]==src[5] and total[d]==nc[d-1]
 if d<=6:small[d]=[[d,i,list(ky)] for i,ky in enumerate(sorted(ownkeys),1)]
 coverage.append({'degree':d,'closed_points':int(total[d]),'full_level_orbits':len(sources[d])});print('C234_SOURCE_INDEPENDENT',d,total[d],len(sources[d]),flush=True)
expected=[[1,d,r[1]] for d in range(4,11) for r in sources[d]]
for d,e in [(4,4),(4,5),(4,6),(5,5)]:
 for x in sources[d]:
  for y in small[e]:expected.append([2,d,x[1],e,y[1]])
assert len(expected)==521
maps=0;tested={};remaining=[]
for source,sm,cover,images,tests,left in filters:
 assert source==[r for d in range(4,11) for r in sources[d]];assert {d:rows for d,rows in sm}==small
 imdict={tuple(im[0]):im for im in images};assert len(imdict)==len(images)
 for k,desc,cd,checks,passed in tests:
  assert k not in tested and desc==expected[k-1];tested[k]=True;sp=[sources[desc[1]][desc[2]-1]]
  if desc[0]==2:sp.append(small[desc[3]][desc[4]-1])
  assert cd==10-sum(r[0] for r in sp);flags=[]
  for variant,cl,yes in checks:
   sumcl=[0]*len(cl)
   for d,idx,ky in sp:
    K=fields[d];z=K.gen();b=sum(K(v)*z**i for i,v in enumerate(ky[:d]));c=sum(K(v)*z**i for i,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b)
    if variant==1:EG,U=E,Q
    elif variant==2:EG,U=E,P+Q
    elif variant==3:EG,U=iso2(E,P,Q)
    else:EG,U=iso2(E,P+17*Q,Q)
    assert 34*U==EG(0) and 17*U!=EG(0) and 2*U!=EG(0);xy=tate(EG,U);im=imdict[tuple([d]+list(key(xy,d)))];e=im[1];w=im[2];assert d==e*w;ee=min(int(e) for e in divisors(d) if all(x**(3**int(e))==x for x in xy));assert ee==e;sumcl=[a+b for a,b in zip(sumcl,im[3])];maps+=1
   assert gkey(sumcl)==gkey(cl);accepted=qkey(cl)==(0,)*NN;assert bool(yes)==accepted;flags.append(accepted)
  assert bool(passed)==all(flags)
 remaining+=left
assert set(tested)==set(range(1,522)) and not remaining
out={'status':'PASS','source_genus':49,'source_non_cuspidal_counts':nc,'coverage':coverage,'support_patterns':len(expected),'checked_normed_summands':maps,'target_rational_image_index':128,'every_target_cusp_class_rational':True,'remaining_non_cuspidal_patterns':0,'pure_cusp_branch_requires_separate_short_unit_argument':True};(SAVED/'independent_sieve.json').write_text(json.dumps(out,indent=2)+'\n');print('C234_INDEPENDENT_SIEVE_COMPLETE',maps,time.monotonic()-t0,flush=True)
