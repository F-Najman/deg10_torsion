from sage.all import *
from pathlib import Path
from collections import Counter
import ast,re,json,time
T=time.monotonic();ROOT=Path('work/deg10/continuation_20261031b/M4N28');Fp=GF(3);PR=PolynomialRing(Fp,'z')
def rd(s):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',(ROOT/s).read_text().replace('[*','[').replace('*]',']'))))
a=rd('fourmap_Y28.m');target=rd('inputs/target28_proven_picard.m');fields={d:GF(3**d,'a',modulus=PR(p))for d,p in a[0]}
def bc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def xy(E,P):
 b,c=bc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def ky(v,d):return min(tuple(int(a.polynomial()[i])for a in[v[0]**(3**j),v[1]**(3**j)]for i in range(d))for j in range(d))
def isok(E,P):return 28*P==E(0)and 14*P!=E(0)and 4*P!=E(0)
relations=matrix(ZZ,target[5]).row_module();observations=[];counts=Counter()
for rec in a[1]:
 si,d,mul,key,cls,obs,flags=rec;K=fields[d];b=PR(key[:d])(K.gen());c=PR(key[d:])(K.gen());E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert isok(E,Q)and P2!=E(0)and 2*P2==E(0)and P2!=14*Q
 roots=E.division_polynomial(2).roots(multiplicities=False);assert len(roots)==3
 if d!=5:assert d%2==0 and all((roots[i]-roots[j]).is_square()for i in range(3)for j in range(i))
 assert mul==(2 if d==5 else 1);counts[d]+=1;rows=[];S=PolynomialRing(K,'X');X=S.gen()
 for j in range(4):
  if j<2:ET=E;QT=Q if j==0 else P2+Q
  else:
   P=P2 if j==2 else P2+14*Q;phi=E.isogeny(X-P[0]);ET=phi.codomain();QT=phi(Q);assert QT.xy()==(Q[0]+(Q+P)[0]-P[0],Q[1]+(Q+P)[1]-P[1])
  assert isok(ET,QT);v=xy(ET,QT);e=lcm([z.minpoly().degree()for z in v]);ck=[d]+list(ky(v,d));o=obs[j];image=a[4][o[1]-1];assert o[0]==j+1 and image[0]==ck and vector(ZZ,[mul*z-y for z,y in zip(image[1],cls[j])])in relations and image[2]==o[2]==e and o[3]*e==mul*d;rows.append([j+1,ck,int(e),mul*d//int(e)])
 observations.append([si,rows])
assert dict(counts)=={5:2,6:3,8:11,10:158}
L=matrix(ZZ,target[5]+target[12]+[[0,0,0,0,1]]).row_module();B=L.basis_matrix();assert abs(B.det())==4;BI=B.inverse()
def qkey(v):return tuple(int(4*z)%4 for z in vector(ZZ,v)*BI)
allowed=[set()for _ in range(6)];allowed[0].add(qkey([0]*5))
for d,v in zip(target[8],target[9]):
 for k in range(d,6):allowed[k]|={tuple((x+y)%4 for x,y in zip(qkey(v),z))for z in allowed[k-d]}
left=[]
for rec in a[1]:
 si,d,mul,key,cls,obs,flags=rec;actual=[qkey([-z for z in v])in allowed[10-d*mul]for v in cls];assert actual==flags
 if all(actual):left.append([si,[si],10-d*mul])
assert left==a[2]and len(left)==19
out={'status':'PASS','source_models':174,'eligible_source_orbits_by_degree':dict(counts),'genuine_marked_maps':696,'manual_Velu_maps':348,'integer_quotient_index':4,'complete_flags_replayed':True,'remaining_cases':19,'observations':observations,'seconds':time.monotonic()-T};(ROOT/'parent_independent_fourmap_v3.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_M428_INDEPENDENT_FOURMAP_COMPLETE',174,696,19,time.monotonic()-T,flush=True)
