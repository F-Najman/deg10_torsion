from sage.all import *
from pathlib import Path
from collections import Counter
import ast,re,json,time,itertools
ROOT=Path('work/deg10/continuation_20261031b/M3N33');T=time.monotonic()
def rd(name):
 s=(ROOT/name).read_text();s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
k=GF(2);R=PolynomialRing(k,'u');n=33
src=rd('physical_source5.m');pri=rd('rational_trace_A5_RR.m');extra=rd('rational_trace_3isog_A5_RR.m')
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

allcurves=rd('inputs/all_curves.m');counts=[];orbits=set();allsource=set();allcyclic=set();curves_checked=0
for dd,pol,curves in allcurves:
 if dd%2:continue
 A=GF(2**dd,'a',modulus=R(pol));a=A.gen();keys=set();bkeys=set();rr=R([1,1,1]).roots(A,multiplicities=False);omega=rr[0] if dd!=10 else R(src[1])(a)
 for jj,aa,card in curves:
  if card%99:continue
  E=EllipticCurve(A,[R(v)(a) for v in aa]);assert E.cardinality()==card;curves_checked+=1;pts=E.points();pp=[P for P in pts if P!=E(0) and 3*P==E(0)]
  if len(pp)!=8:continue
  qq=[Q for Q in pts if 33*Q==E(0) and 11*Q!=E(0) and 3*Q!=E(0)]
  for Q in qq:
   bkeys.add(fk(bc(E,Q)));orbits.add(min(fk(bc(E,j*Q)) for j in range(1,33) if gcd(j,33)==1))
   for P in pp:
    if P.weil_pairing(11*Q,3)==omega:keys.add(fk(normal3(E,P,Q)))
 counts.append(len(keys));print('M333_INDEPENDENT_SOURCE_DEGREE',dd,len(keys),len(bkeys),flush=True)
 if dd==10:allsource=keys;allcyclic=bkeys
assert counts==[0,0,0,0,264] and allsource==set(map(tuple,src[3])) and allcyclic==set(map(tuple,src[4]))
assert len(allcyclic)==88 and orbits=={tuple(v[3]) for v in src[5]} and len(orbits)==12
L=GF(2**40,'v');S=PolynomialRing(L,'X');X=S.gen();root=S(src[0]).roots(multiplicities=False)[0]
emb=lambda v:sum((L(v.polynomial()[i])*root**i for i in range(10)),L(0))
selected=[];checks=[]
for case,sr in enumerate(src[5],1):
 pr=pri[case-1];assert pr[0]==case
 ky=sr[3];b=R(ky[:10])(z);c=R(ky[10:])(z);Eorig=EllipticCurve(K,[1-c,-b,-b,0,0]);Qorig=Eorig(0,0);assert 33*Qorig==Eorig(0) and 11*Qorig!=Eorig(0) and 3*Qorig!=Eorig(0)
 used3=False;kernelx=None
 if pr[9]!=0:
  good=[r for r in extra if r[0]==case and r[12]==0];assert good
  rr=good[0];kernelx=R(rr[1])(z);P3=Eorig.lift_x(kernelx);assert 3*P3==Eorig(0)
  phi3=Eorig.isogeny(PolynomialRing(K,'T').gen()-kernelx);E0=phi3.codomain();Q0=phi3(Qorig);assert phi3.degree()==3 and 33*Q0==E0(0) and 11*Q0!=E0(0) and 3*Q0!=E0(0)
  nonzero=[P3,2*P3];vx=Qorig[0]+sum(((Qorig+A)[0]-A[0] for A in nonzero),K(0));vy=Qorig[1]+sum(((Qorig+A)[1]-A[1] for A in nonzero),K(0));assert (vx,vy)==Q0.xy()
  EP=EllipticCurve(K,[R(v)(z) for v in rr[2]]);QP=EP([R(v)(z) for v in rr[3]])
  assert bc(EP,QP)==bc(E0,Q0)
  pr=[rr[0],rr[4],rr[5],rr[6],rr[7],rr[8],rr[9],rr[10],rr[11],rr[12],rr[13]];used3=True
 else:E0=Eorig;Q0=Qorig
 assert pr[9]==0;selected.append({'case':case,'used3isogeny':used3,'extra_row':(extra.index(rr)+1) if used3 else None})
 E=EllipticCurve(L,list(map(emb,E0.a_invariants())));Q=E(list(map(emb,Q0.xy())));actual=Counter()
 assert ik([L(0)]+xy(E,Q))==ik([K(0)]+[R(v)(z) for v in pr[6]])
 assert ik([L(0)]+xy(E,5*Q))==ik([K(0)]+[R(v)(z) for v in pr[7]])
 xx=E.division_polynomial(5).roots(multiplicities=False);assert len(xx)==12;points=[P for x0 in xx for P in E.lift_x(x0,all=True)];assert len(points)==24;kernels={}
 for P in points:
  assert 5*P==E(0) and P!=E(0);f=(X-P[0])*(X-(2*P)[0]);kernels[tuple(f.list())]=P
 assert len(kernels)==6
 for co,P in kernels.items():
  phi=E.isogeny(S(list(co)));ET=phi.codomain();QT=phi(Q);assert phi.degree()==5 and 33*QT==ET(0) and 11*QT!=ET(0) and 3*QT!=ET(0)
  nonzero=[j*P for j in range(1,5)];vx=Q[0]+sum(((Q+A)[0]-A[0] for A in nonzero),L(0));vy=Q[1]+sum(((Q+A)[1]-A[1] for A in nonzero),L(0));assert (vx,vy)==QT.xy();actual[json.dumps(ik([L(0)]+xy(ET,QT)))]+=1
 expected=Counter();norms=[]
 for orbit,norm,pd,vals,pol in pr[8]:
  A=GF(2**(len(pol)-1),'b',modulus=R(pol));a=A.gen();vv=[R(v)(a) for v in vals];key=ik([A(0)]+vv);assert key[0]==pd[0] and 10*orbit==norm*pd[0];expected[json.dumps(key)]+=orbit;norms.append([orbit,norm,pd[0]])
 assert expected==actual
 checks.append({'case':case,'rational3isogeny_used':used3,'source_point_F2_degree':pr[2][0],'source_norm':pr[3],'actual5kernels':6,'kernel_orbits_and_norms':norms,'manual_Velu_maps':6+int(used3)})
 print('M333_INDEPENDENT_GENUINE5_CASE_PASS',case,used3,norms,flush=True)
report={'status':'PASS','source_NC_closed_counts_F4':counts,'source_places':264,'cyclic33_projection_places':88,'cyclic33_diamond_orbits':12,'source_elliptic_models_checked':curves_checked,'selected_tests':selected,'actual5kernels':72,'actual3kernels_for_selected_lifts':3,'all_marked_Velu_maps':75,'all_QQtrace_norms_and_keys_bound':True,'cases':checks,'elapsed_seconds':time.monotonic()-T}
(ROOT/'independent_genuine35.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_M333_INDEPENDENT_SOURCE_GENUINE35_COMPLETE',time.monotonic()-T,flush=True)
