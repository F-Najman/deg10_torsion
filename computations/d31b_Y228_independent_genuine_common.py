from sage.all import *
from pathlib import Path
from collections import Counter
import ast,re,json,time
T=time.monotonic();ROOT=Path('work/deg10/continuation_20261031b/Y228_shared');Fp=GF(3);RP=PolynomialRing(Fp,'z')
def rd(n):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',(ROOT/n).read_text().replace('[*','[').replace('*]',']'))))
uc=json.loads((ROOT/'union_normalization.json').read_text());packet=rd('parent_native_all45_and_neighbors_v2.m');neighbors=rd('m428_neighbor_A5.m');fields={d:pol for d,pol in rd('m428_neighbor_A5.m')[:0]}
# Union field polynomials were retained in the independent normalization input.
a=ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',(ROOT.parent/'M4N28/fourmap_Y28.m').read_text().replace('[*','[').replace('*]',']'))));fields={int(d):GF(3**int(d),'a',modulus=RP(pol))for d,pol in a[0]}
def ik(v):
 A=v[0].parent();d=int(lcm([x.minpoly().degree()for x in v]));found=False
 for bits in[[1,0],[0,1],[1,1],[1,2]]:
  g=sum((A(a)*x for a,x in zip(bits,v)),A(0));pol=g.minpoly()
  if pol.degree()==d:found=True;break
 assert found
 def vec(x):return vector(Fp,[A(x).polynomial()[i]for i in range(A.degree())])
 M=matrix(Fp,[vec(g**i)for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vec(x))))for x in v];return[d,bits,list(map(int,pol.list())),co]
def normal(E,P,Q):
 assert E.a1()==0 and E.a3()==0 and P[1]==0
 x=Q[0]-P[0];y=Q[1];return[x**3/y**2,(E.a2()+3*P[0])*x*x/y**2]
def verify(E,P,Q,source,dia,hec):
 K=E.base_field();e=K.degree();assert 2*P==E(0)and P!=E(0)and 28*Q==E(0)and 14*Q!=E(0)and 4*Q!=E(0)and P!=14*Q
 assert ik(normal(E,P,Q))==source and ik(normal(E,P,5*Q))==dia
 m=int(4*lcm([int(h[0])for h in hec]));L=GF(3**(e*m),'v');S=PolynomialRing(L,'X');X=S.gen();roots=S(K.modulus().list()).roots(multiplicities=False);assert len(roots)==e;z=roots[0]
 emb=lambda v:sum((L(v.polynomial()[i])*z**i for i in range(e)),L(0))
 A=EllipticCurve(L,list(map(emb,E.a_invariants())));PP=A(list(map(emb,P.xy())));QQ=A(list(map(emb,Q.xy())));xx=A.division_polynomial(5).roots(multiplicities=False);assert len(xx)==12;points=[V for x in xx for V in A.lift_x(x,all=True)];assert len(points)==24;kernels={}
 for V in points:
  assert 5*V==A(0)and V!=A(0);f=(X-V[0])*(X-(2*V)[0]);kernels[tuple(f.list())]=V
 assert len(kernels)==6;actual=Counter()
 for co,V in kernels.items():
  phi=A.isogeny(S(list(co)));B=phi.codomain();P2=phi(PP);TQ=phi(QQ);assert phi.degree()==5 and 2*P2==B(0)and P2!=B(0)and P2!=14*TQ and 28*TQ==B(0)and 14*TQ!=B(0)and 4*TQ!=B(0)
  nonzero=[j*V for j in range(1,5)]
  for Z,W in[(PP,P2),(QQ,TQ)]:
   assert W.xy()==(Z[0]+sum(((Z+Vv)[0]-Vv[0]for Vv in nonzero),L(0)),Z[1]+sum(((Z+Vv)[1]-Vv[1]for Vv in nonzero),L(0)))
  actual[json.dumps(ik(normal(B,P2,TQ)))]+=1
 expected=Counter()
 for ell,n,key in hec:
  assert e*ell==n*key[0];expected[json.dumps(key)]+=ell
 assert actual==expected
 return {'source_field_degree':e,'genuine5kernels':6,'manual_marked_Velu_maps':12,'splitting_field_degree':int(e*m)}
checks=[];ncomp=0
for row in packet[0]:
 cid=int(row[0]);
 if (cid-1)%NSHARDS!=SHARD:continue
 for ob in row[2]:
  e,mul,ky,source,dia,hec=ob;K=fields[int(e)];z=K.gen();b=RP(ky[:e])(z);c=RP(ky[e:])(z);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);ans=verify(E,P,Q,source,dia,hec);checks.append({'case':cid,'multiplicity':mul,**ans});ncomp+=1
 print('Y228_INDEPENDENT_GENUINE_CASE',cid,ncomp,time.monotonic()-T,flush=True)
for r in packet[1]:
 cid,ki,source,dia,hec,h0=r
 if (int(cid)-1)%NSHARDS!=SHARD:continue
 nr=[v for v in neighbors if v[0]==cid and v[1]==ki][0];K=fields[10];z=K.gen();E=EllipticCurve(K,[RP(v)(z)for v in nr[3]]);P=E([RP(v)(z)for v in nr[4]]);Q=E([RP(v)(z)for v in nr[5]]);Ukey=uc['cases'][int(cid)-1][2][0][2];b=RP(Ukey[:10])(z);c=RP(Ukey[10:])(z);E0=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P0=E0(0,0);Q0=E0(b,b);U=E0([RP(v)(z)for v in nr[2]]);assert 2*U==P0;phi=E0.isogeny(P0);B=phi.codomain();PP=phi(U);QQ=phi(Q0)
 for Z,W in[(U,PP),(Q0,QQ)]:assert W.xy()==(Z[0]+(Z+P0)[0]-P0[0],Z[1]+(Z+P0)[1]-P0[1])
 assert normal(B,PP,QQ)==normal(E,P,Q);ans=verify(E,P,Q,source,dia,hec);checks.append({'neighbor_case':cid,**ans});print('M428_INDEPENDENT_NEIGHBOR_PASS',cid,time.monotonic()-T,flush=True)
out={'status':'PASS','shard':SHARD,'component_tests':ncomp,'checks':checks,'seconds':time.monotonic()-T};(ROOT/f'independent_genuine_s{SHARD}.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_Y228_INDEPENDENT_GENUINE_COMPLETE',SHARD,len(checks),time.monotonic()-T,flush=True)
