from sage.all import *
from pathlib import Path
from collections import Counter,defaultdict
import json,time
T=time.monotonic();ROOT=Path('work/deg10/continuation_20261031a/C56');Fp=GF(3);PR=PolynomialRing(Fp,'z')
def rd(s):return json.loads((ROOT/s).read_text().replace('[*','[').replace('*]',']'))
a=rd('fourmap_filter.m');target=rd('target28_proven_picard.m');src=rd('source56_p3.m');fields={d:GF(3**d,'a',modulus=PR(p))for d,p in a[0]if d>=5}
def bc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def xy(E,P):
 b,c=bc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def ky(v,d):return min(tuple(int(a.polynomial()[i])for a in[v[0]**(3**j),v[1]**(3**j)]for i in range(d))for j in range(d))
def isok(E,P):return 28*P==E(0)and 14*P!=E(0)and 4*P!=E(0)
seen=defaultdict(set);observations=[]
for rec in a[1]:
 si,d,key,cls,obs=rec;K=fields[d];b=PR(key[:d])(K.gen());c=PR(key[d:])(K.gen());E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0);assert 56*P==E(0)and 28*P!=E(0)and 8*P!=E(0)
 if si<=len(src[2]):
  ks={ky(bc(E,u*P),d)for u in range(1,29)if gcd(u,56)==1};assert min(ks)==tuple(key)and not(seen[d]&ks);seen[d]|=ks
 R=PolynomialRing(K,'X');X=R.gen();T2=28*P;phi=E.isogeny(X-T2[0]);EE=phi.codomain();Q=phi(P);assert isok(EE,Q)and(Q[0],Q[1])==(P[0]+(P+T2)[0]-T2[0],P[1]+(P+T2)[1]-T2[1]);dual=phi.dual();assert dual(Q)==2*P
 ai=EE.a_invariants();twos=[EE(t,-(ai[0]*t+ai[2])/2)for t in EE.division_polynomial(2).roots(multiplicities=False)];assert len(twos)==3;ds=[S for S in twos if dual(S)==E(0)];assert len(ds)==1;P2=ds[0];assert P2!=14*Q
 rows=[]
 for j in range(4):
  if j<2:ET=EE;QT=Q if j==0 else P2+Q
  else:
   S=P2 if j==2 else P2+14*Q;psi=EE.isogeny(X-S[0]);ET=psi.codomain();QT=psi(Q);assert QT.xy()==(Q[0]+(Q+S)[0]-S[0],Q[1]+(Q+S)[1]-S[1])
  assert isok(ET,QT);v=xy(ET,QT);e=lcm([z.minpoly().degree()for z in v]);ck=[d]+list(ky(v,d));o=obs[j];image=a[4][o[1]-1];assert o[0]==j+1 and image[0]==ck and image[1]==cls[j]and image[2]==o[2]==e and o[3]*e==d;rows.append([j+1,ck,int(e),d//int(e)])
 observations.append([si,rows])
 if si%100==0:print('D31A_C56_INDEPENDENT_FOURMAP_PROGRESS',si,time.monotonic()-T,flush=True)
assert[ len(seen[d])for d in range(1,11)]==src[1]
assert {tuple(v[2])for v in a[1][len(src[2]):]}==seen[5]
L=matrix(ZZ,target[5]+target[12]+[[0,0,0,0,1]]).row_module();B=L.basis_matrix();assert abs(B.det())==4;BI=B.inverse()
def qkey(v):return tuple(int(4*z)%4 for z in vector(ZZ,v)*BI)
allowed=[set()for _ in range(6)];allowed[0].add(qkey([0]*5))
for d,v in zip(target[8],target[9]):
 for k in range(d,6):
  allowed[k]|={tuple((x+y)%4 for x,y in zip(qkey(v),z))for z in allowed[k-d]}
assert list(map(len,allowed))==[1,1,1,3,3,3]
left=[]
for row in a[2]:
 ti,ids,cd,flags=row;assert sum(a[1][i-1][1]for i in ids)+cd==10
 actual=[]
 for j in range(4):
  v=sum((vector(ZZ,a[1][i-1][3][j])for i in ids),vector(ZZ,[0]*5));actual.append(qkey(-v)in allowed[cd])
 assert actual==flags
 if all(actual):left.append(row)
assert left==a[3]and len(left)==135
out={'status':'PASS','degree':10,'group':'C56','source_diamond_orbits':len(src[2]),'all_source_closed_counts_independently_verified':src[1],'physical_source_map_records':len(observations),'genuine_marked_projection_maps':4*len(observations),'explicit_Velu_coordinate_checks':3*len(observations),'all_target_keys_and_norms_match':True,'full_integer_quotient_sieve_replayed':True,'tested_support_patterns':len(a[2]),'remaining_support_patterns':len(left),'remaining_pattern_types':[[list(k[0]),k[1],n]for k,n in Counter((tuple(a[1][i-1][1]for i in v[1]),v[2])for v in left).items()],'observations':observations,'elapsed_seconds':time.monotonic()-T}
(ROOT/'independent_fourmap.json').write_text(json.dumps(out,indent=2)+'\n');print('D31A_C56_INDEPENDENT_FOURMAP_COMPLETE',len(observations),len(a[2]),len(left),time.monotonic()-T,flush=True)
