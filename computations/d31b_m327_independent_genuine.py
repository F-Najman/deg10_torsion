from sage.all import *
from pathlib import Path
from collections import Counter
import ast,re,json,time,itertools
ROOT=Path('work/deg10/continuation_20261031b/M3N27');T=time.monotonic()
def rd(name):
 s=(ROOT/name).read_text();s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
k=GF(2);R=PolynomialRing(k,'u');n=27
src=rd('physical_source5.m');pri=rd('quotient_A7_RR.m');bind=rd('quotient_A7_v5_RR.m')
assert len(pri)==len(bind)==6 and all(a==b[:11] for a,b in zip(pri,bind))
K=GF(2**10,'z',modulus=R(src[0]));z=K.gen();zet=R(src[1])(z)
def bc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();v3=2*y+a3+a1*x;v2=3*x+a2;v4=3*x*x+2*x*a2+a4-a1*y;t1=2*v4/v3+a1;t2=v2-(v4/v3)**2-a1*v4/v3;return [-t2**3/v3**2,-(t1*t2-v3)/v3]
def xy(E,P):
 b,c=bc(E,P);r=b/c;s=c*c/(b-c);u=r*s-2*r+1;return [(s-r)/u,u/(s*s-s-r+1)]
def normal3(E,P,Q):
 a1,a2,a3,a4,a6=E.a_invariants();r,tt=P.xy();ss=(3*r*r+2*a2*r+a4-a1*tt)/(2*tt+a1*r+a3);aa1=a1+2*ss;aa3=a3+a1*r+2*tt;xx=Q[0]-r;yy=Q[1]-tt-ss*xx;sc=yy/xx**2;t=-sc**2*xx;q=(sc*aa1-t-2)/(t-1);assert sc**3*aa3==q*t*t+(1-q)*t;return [q,t]
def fk(v):
 d=v[0].parent().degree();return min(tuple(int((a**(4**j)).polynomial()[i]) for a in v for i in range(d)) for j in range(d//2))
def ik(v):
 A=v[0].parent();d=int(lcm([a.minpoly().degree() for a in v]));found=False
 choices=[[1,0,0],[0,1,0],[0,0,1]]+[list(a) for a in itertools.product([0,1],repeat=3) if sum(a)>1]
 for bits in choices:
  g=sum((A(a)*b for a,b in zip(bits,v)),A(0));pol=g.minpoly()
  if pol.degree()==d:found=True;break
 assert found
 def vec(a):return vector(k,[A(a).polynomial()[i] for i in range(A.degree())])
 M=matrix(k,[vec(g**i) for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vec(a)))) for a in v];return [d,bits,list(map(int,pol.list())),co]
# Independent exhaustive elliptic rational-point enumeration over every F4^d, d<=5.
allcurves=rd('inputs/all_curves.m');counts=[];orbits=set();allsource=set();allcyclic=set();curves_checked=0
for dd,pol,curves in allcurves:
 if dd%2:continue
 A=GF(2**dd,'a',modulus=R(pol));a=A.gen();keys=set();bkeys=set();rr=R([1,1,1]).roots(A,multiplicities=False);omega=rr[0] if dd!=10 else R(src[1])(a)
 for jj,aa,card in curves:
  if card%81:continue
  E=EllipticCurve(A,[R(v)(a) for v in aa]);assert E.cardinality()==card;curves_checked+=1
  pts=E.points();pp=[P for P in pts if P!=E(0) and 3*P==E(0)]
  if len(pp)!=8:continue
  qq=[Q for Q in pts if 27*Q==E(0) and 9*Q!=E(0)]
  for Q in qq:
   qk=fk(bc(E,Q));bkeys.add(qk);orb=min(fk(bc(E,j*Q)) for j in range(1,27) if gcd(j,27)==1);orbits.add(orb)
   for P in pp:
    if P.weil_pairing(9*Q,3)!=omega:continue
    keys.add(fk(normal3(E,P,Q)))
 counts.append(len(keys));print('M327_INDEPENDENT_SOURCE_DEGREE',dd,len(keys),len(bkeys),flush=True)
 if dd==10:allsource=keys;allcyclic=bkeys
assert counts==[0,0,0,0,162]
assert allsource==set(map(tuple,src[3])) and allcyclic==set(map(tuple,src[4]))
assert orbits=={tuple(v[3]) for v in src[5]} and len(orbits)==6
# Decode source constant orientation saved by the actual producer.
L=GF(2**60,'v');S=PolynomialRing(L,'X');X=S.gen();root=S(src[0]).roots(multiplicities=False)[0]
emb=lambda v:sum((L(v.polynomial()[i])*root**i for i in range(10)),L(0))
checks=[]
for case,(sr,pr) in enumerate(zip(src[5],bind),1):
 ky=sr[3];b=R(ky[:10])(z);c=R(ky[10:])(z);E0=EllipticCurve(K,[1-c,-b,-b,0,0]);Q0=E0(0,0);w0=R(pr[11])(z);assert w0*w0+w0+1==0
 E=EllipticCurve(L,list(map(emb,E0.a_invariants())));Q=E(list(map(emb,Q0.xy())));w=emb(w0)
 assert ik([w]+xy(E,Q))==ik([w0]+[R(co)(z) for co in pr[6]])
 assert ik([w]+xy(E,7*Q))==ik([R(pr[12])(z)]+[R(co)(z) for co in pr[7]])
 xs=E.division_polynomial(7).roots(multiplicities=False);assert len(xs)==24
 points=[P for x0 in xs for P in E.lift_x(x0,all=True)];assert len(points)==48
 kernels={}
 for P in points:
  assert 7*P==E(0) and P!=E(0);f=prod(X-(j*P)[0] for j in range(1,4));kernels[tuple(f.list())]=P
 assert len(kernels)==8;actual=Counter()
 for co,P in kernels.items():
  f=S(list(co));phi=E.isogeny(f);ET=phi.codomain();QT=phi(Q);assert phi.degree()==7 and 27*QT==ET(0) and 9*QT!=ET(0)
  nonzero=[j*P for j in range(1,7)]
  vx=Q[0]+sum(((Q+A)[0]-A[0] for A in nonzero),L(0));vy=Q[1]+sum(((Q+A)[1]-A[1] for A in nonzero),L(0));assert (vx,vy)==QT.xy()
  actual[json.dumps(ik([w]+xy(ET,QT)))]+=1
 expected=Counter();norms=[]
 for orbit,norm,pd,vals,pol,wco in pr[8]:
  A=GF(2**(len(pol)-1),'b',modulus=R(pol));a=A.gen();vv=[R(co)(a) for co in vals];wa=R(wco)(a);assert wa*wa+wa+1==0
  key=ik([wa]+vv);assert key[0]==2*pd[0] and 5*orbit==norm*pd[0];expected[json.dumps(key)]+=orbit;norms.append([orbit,norm,pd[0]])
 assert actual==expected
 checks.append({'case':case,'geometric_kernels':8,'kernel_orbits':norms,'source_GF4_constant':pr[11],'relative_GF4_keys_match':True,'manual_Velu_checks':8})
 print('M327_INDEPENDENT_GENUINE7_CASE_PASS',case,norms,flush=True)
report={'status':'PASS','source_NC_closed_counts_F4':counts,'actual_source_keys':len(allsource),'cyclic_projection_keys':len(allcyclic),'cyclic_diamond_orbits':6,'eligible_curve_models_checked':curves_checked,'geometric7kernels_checked':48,'manual_Velu_checks':48,'constant_field_orientation_explicit':True,'v4_v5_all_prime_and_divisor_records_identical':True,'cases':checks,'elapsed_seconds':time.monotonic()-T}
(ROOT/'independent_genuine7.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_M327_INDEPENDENT_SOURCE_GENUINE7_COMPLETE',time.monotonic()-T,flush=True)
