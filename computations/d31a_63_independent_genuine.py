from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C63';start=time.monotonic();k=GF(2);R2=PolynomialRing(k,'T')
def read(name):return ast.literal_eval((C/name).read_text().replace('[*','[').replace('*]',']'))
def E9at(t):
 c=t*t*(t-1);b=c*(t*t-t+1);return EllipticCurve(t.parent(),[1-c,-b,-b,0,0])
def normalized(E,R,S):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=R.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;r=b4/b3;n2=b2-r*r-a1*r;scale=b3/n2
 b=-n2**3/b3**2;c=-((2*r+a1)*n2-b3)/b3;t=c*c/(b-c)
 assert c==t*t*(t-1) and b==c*(t*t-t+1)
 xx=(S[0]-x)/scale**2;yy=(S[1]-y-r*(S[0]-x))/scale**3;En=E9at(t);Q=En(xx,yy);assert 7*Q==En(0) and Q!=En(0)
 return [t,xx]
def key2(vals):
 L=vals[0].parent();d=lcm([v.minpoly().degree()for v in vals]);found=False
 for bits in[[1,0],[0,1],[1,1]]:
  g=bits[0]*vals[0]+bits[1]*vals[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found;n=L.degree();vv=lambda x:vector(k,[L(x).polynomial()[i]for i in range(n)]);M=matrix(k,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(v))))for v in vals]]
def frobkey(vals,e):return min(tuple(int(v.polynomial()[i])for v in[vals[0]**(2**j),vals[1]**(2**j)]for i in range(e))for j in range(e))
ans=[];census=[]
for e in[7,10]:
 model=read('good%d_model.m'%e);assert model[2]==e;poly=model[4];K=GF(2**e,'a',modulus=R2(poly));a=K.gen();allpairs=set();orbits=set();tot=0
 for t in K:
  try:E=E9at(t)
  except (ArithmeticError,ValueError,ZeroDivisionError):continue
  if E.cardinality()%63:continue
  assert E.cardinality()==(126 if e==7 else 1008)
  roots=E.division_polynomial(7).roots(multiplicities=False)
  for xx in roots:
   points=E.lift_x(xx,all=True);assert len(points)==2 and 7*points[0]==E(0)
   key=frobkey([t,xx],e);allpairs.add(key);R=E(0,0);S=points[0];ks={frobkey(normalized(E,u*R,u*S),e)for u in range(1,32)if gcd(u,63)==1};assert len(ks)==9;orbits.add(min(ks));tot+=1
 assert len(allpairs)==27 and len(orbits)==3 and tot==(189 if e==7 else 270)
 supplied=set()
 for ti,xi in model[3]:
  t=R2(ti)(a);xx=R2(xi)(a);E=E9at(t);R=E(0,0);S=E.lift_x(xx);supplied.add(min(frobkey(normalized(E,u*R,u*S),e)for u in range(1,32)if gcd(u,63)==1))
 assert supplied==orbits;census.append({'degree':e,'geometric_pairs':tot,'closed_places':27,'diamond_orbits':3})
 L=GF(2**(24*e),'w');RL=PolynomialRing(L,'x');x=RL.gen();root=RL(poly).roots(multiplicities=False)[0]
 for ci,(ti,xi)in enumerate(model[3],1):
  t=R2(ti)(root);xx=R2(xi)(root);E=E9at(t);R=E(0,0);S=E.lift_x(xx);assert 9*R==E(0)and 3*R!=E(0)and 7*S==E(0)
  sk=key2(normalized(E,R,S));dk=key2(normalized(E,5*R,5*S));assert sk[0]==e and dk[0]==e
  xs=E.division_polynomial(5).roots(multiplicities=False);pts=[Q for xx5 in xs for Q in E.lift_x(xx5,all=True)];assert len(xs)==12 and len(pts)==24 and all(5*Q==E(0)and Q!=E(0)for Q in pts)
  kernels={tuple(((x-Q[0])*(x-(2*Q)[0])).list()):Q for Q in pts};assert len(kernels)==6
  counts=Counter();obs={}
  for co,Q in kernels.items():
   phi=E.isogeny(RL(list(co)));EE=phi.codomain();RR=phi(R);SS=phi(S);assert phi.degree()==5 and 9*RR==EE(0)and 3*RR!=EE(0)and 7*SS==EE(0)
   nz=[j*Q for j in range(1,5)]
   for A,AA in[(R,RR),(S,SS)]:
    velu=(A[0]+sum(((A+Z)[0]-Z[0]for Z in nz),L(0)),A[1]+sum(((A+Z)[1]-Z[1]for Z in nz),L(0)));assert velu==AA.xy()
   ky=key2(normalized(EE,RR,SS));tag=json.dumps(ky);counts[tag]+=1;obs[tag]=ky
  rows=[[int(counts[tag]),obs[tag]]for tag in sorted(counts)];ans.append({'degree':e,'case':ci,'source_key':sk,'diamond5_key':dk,'target_keys':rows,'actual_geometric5_kernels':6,'manual_Velu_marked_coordinate_checks':12});print('D31A_C63_GENUINE_CASE',e,ci,[v[0]for v in rows],time.monotonic()-start,flush=True)
result={'status':'PASS','independent_source_census':census,'source_cases':ans,'actual_geometric5_kernels':36,'manual_Velu_marked_coordinate_checks':72,'method':'Independent Sage exhaustive E9(t), marked 7-point normalization, all six actual 5-kernels and manual four-point Velu sums on both marked points','elapsed_seconds':time.monotonic()-start}
(C/'independent_genuine_keys.json').write_text(json.dumps(result,indent=2)+'\n');print('D31A_C63_INDEPENDENT_GENUINE_COMPLETE',36,72,time.monotonic()-start,flush=True)
