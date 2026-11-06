from sage.all import *
from pathlib import Path
from collections import Counter,defaultdict
import ast,json,time
ROOT=Path('work/deg10/continuation_20261031a/C70');START=time.monotonic();Fp=GF(3);RP=PolynomialRing(Fp,'T')
src=ast.literal_eval((ROOT/'source70_p3.m').read_text().replace('[*','[').replace('*]',']'))
primary=json.loads((ROOT/'A2_all.json').read_text())[2]
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y
 n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def xy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return[(s-r)/t,t/(s*s-s-r+1)]
def key2(v):
 K=v[0].parent();n=K.degree();d=lcm([a.minpoly().degree()for a in v]);found=False
 for bits in[[1,0],[0,1],[1,1],[1,2]]:
  g=bits[0]*v[0]+bits[1]*v[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found
 vv=lambda a:vector(Fp,[K(a).polynomial()[i]for i in range(n)])
 M=matrix(Fp,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(a))))for a in v]]
def frobkey(v,d):return min(tuple(int(a.polynomial()[i])for a in[v[0]**(3**j),v[1]**(3**j)]for i in range(d))for j in range(d))
fields={};output=[];closed=defaultdict(int);nker=0
for d,pol in src[0]:
 if d<6:continue
 K=GF(3**d,'a',modulus=RP(pol));L=GF(3**(2*d),'u');RL=PolynomialRing(L,'x');a=RL(pol).roots(multiplicities=False)[0]
 fields[d]=(K,L,a)
for (ci,d,maps),(dd,ky) in zip(primary,src[2]):
 assert ci==len(output)+1 and d==dd
 K,L,a=fields[d];z=K.gen();b=RP(ky[:d])(z);c=RP(ky[d:])(z);E0=EllipticCurve(K,[1-c,-b,-b,0,0]);P0=E0(0,0)
 assert 70*P0==E0(0)and all((70//q)*P0!=E0(0)for q in[2,5,7])
 ks={frobkey(pointbc(E0,u*P0),d)for u in range(1,36)if gcd(u,70)==1};assert min(ks)==tuple(ky);closed[d]+=len(ks)
 emb=lambda v:sum((L(v.polynomial()[i])*a**i for i in range(d)),L(0))
 E=EllipticCurve(L,list(map(emb,E0.a_invariants())));P=E(list(map(emb,P0.xy())));Q=2*P;stages=[]
 for stage in range(len(maps)):
  if stage==0:EE=E;QQ=Q
  else:
   R=PolynomialRing(L,'v');v=R.gen();T=35*P;ph0=E.isogeny(v-T[0]);EE=ph0.codomain();QQ=ph0(Q)
   assert (Q[0]+(Q+T)[0]-T[0],Q[1]+(Q+T)[1]-T[1])==QQ.xy()
  assert 35*QQ==EE(0)and 7*QQ!=EE(0)and 5*QQ!=EE(0)
  sourcekey=key2(xy(EE,QQ));diamondkey=key2(xy(EE,2*QQ));assert d==maps[stage][1]*sourcekey[0]and 2*d==maps[stage][3]*diamondkey[0]
  roots=EE.division_polynomial(2).roots(multiplicities=False);assert len(roots)==3
  R=PolynomialRing(L,'v');v=R.gen();actual=Counter();keys={}
  for sx in roots:
   ai=EE.a_invariants();T=EE(sx,-(ai[0]*sx+ai[2])/2);assert T!=EE(0)and 2*T==EE(0)
   ph=EE.isogeny(v-sx);ET=ph.codomain();QT=ph(QQ)
   assert(QQ[0]+(QQ+T)[0]-T[0],QQ[1]+(QQ+T)[1]-T[1])==QT.xy()and 35*QT==ET(0)
   key=key2(xy(ET,QT));s=json.dumps(key);actual[s]+=1;keys[s]=key;nker+=1
  expected=Counter((row[0],row[1],row[2][0])for row in maps[stage][4])
  obs=[]
  for s,e in actual.items():
   key=keys[s];assert(d*e)%key[0]==0;mult=d*e//key[0];obs.append([e,mult,key[0],key])
  assert Counter((a,b,c)for a,b,c,k in obs)==expected
  stages.append({'map':stage+1,'source_key':sourcekey,'diamond2_key':diamondkey,'targets':obs,'h0':maps[stage][5]})
 output.append({'case':ci,'degree':d,'source_orbit_size':len(ks),'maps':stages})
 if ci%100==0:print('D31A_C70_GENUINE_SOURCE_KERNEL_PROGRESS',ci,nker,time.monotonic()-START,flush=True)
assert len(output)==773 and[closed[d]for d in range(1,11)]==src[1]and nker==2328
out={'status':'PASS','prime':3,'level':70,'complete_source_closed_counts':src[1],'source_diamond_orbits':773,'all_source_canonical_keys_verified':True,'tested_projection_maps':776,'actual_cyclic2_kernels':nker,'manual_Velu_coordinates_verified':True,'norm_degree_orbit_multisets_match_native_producer':True,'cases':output,'elapsed_seconds':time.monotonic()-START}
(ROOT/'independent_genuine_isogenies.json').write_text(json.dumps(out,indent=2)+'\n')
print('D31A_C70_ALL_GENUINE_ISOGENY_CONTROLS_COMPLETE',773,776,nker,time.monotonic()-START,flush=True)
