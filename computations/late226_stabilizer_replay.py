from pathlib import Path
import ast,json
from sage.all import GF,QQ,ZZ,PolynomialRing,NumberField,EllipticCurve,matrix
P=Path('work/deg10/continuation_20261029b/C226')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=read(P/'filter26.m');row=next(r for r in a[1] if r[:2]==[9,52]);pol=next(r[1] for r in a[3] if r[0]==9);R=PolynomialRing(GF(3),'a');K=GF(3**9,'a',modulus=R(pol));z=K.gen();ky=row[2];b=sum(K(v)*z**i for i,v in enumerate(ky[:9]));c=sum(K(v)*z**i for i,v in enumerate(ky[9:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);ts=[T,13*Q,T+13*Q];V=14*Q;coeff=[(1,0),(0,13),(1,13)]
def bc(T,Q):x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
found=[]
for i in range(3):
 for j in range(3):
  if i==j:continue
  for k in range(1,7):
   pair=bc(ts[i],ts[j]+k*V);powers=[s for s in range(9) if pair==[b**(3**s),c**(3**s)]]
   if not powers:continue
   M=matrix(ZZ,[[coeff[i][0],coeff[j][0]],[coeff[i][1],coeff[j][1]+14*k]])
   def scalar(A,s):return all((A[r,t]-(s if r==t else 0))%(2 if r==0 else 26)==0 for r in range(2) for t in range(2))
   if (scalar(M**3,1) or scalar(M**3,-1)) and not(scalar(M,1) or scalar(M,-1)):found.append({'matrix':list(map(list,M.rows())),'frobenius_powers':powers,'cube_scalar':1 if scalar(M**3,1) else -1})
assert found
cm=read(P/'j0_twists.m');S=PolynomialRing(QQ,'z');F=NumberField(S(cm[0]),'z');U=PolynomialRing(F,'x');
def en(v):return sum(QQ(n)/QQ(d)*F.gen()**i for i,(n,d) in enumerate(v))
f=U([en(v) for v in cm[1]]);assert f==EllipticCurve(F,[0,0,0,0,-1]).division_polynomial(13);prod=U(1);degrees=[]
for d,m,co in cm[2]:
 h=U([en(v) for v in co]);assert h.degree()==d and h.is_irreducible();prod*=h**m;degrees += [d]*m
assert prod.monic()==f.monic() and sorted(degrees)==[6,6,72]
report={'status':'PASS','degree_9_descriptor':[1,9,52],'order_three_stabilizers':found,'j0_division_polynomial_factor_degrees':degrees,'no_factor_degree_divides_5':True}
(P/'stabilizer_j0_replay.json').write_text(json.dumps(report,indent=2,default=int)+'\n');print(json.dumps(report,indent=2,default=int),flush=True)
