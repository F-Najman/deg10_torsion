from pathlib import Path
import ast,json,hashlib,time
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
H=Path(__file__).resolve().parents[1];B=H/'next_cases_20261030/C69';start=time.monotonic()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
data=rd(B/'cover_results.m');model=rd(B/'cover_model.m');jf=rd(B/'base_j_functions.m');R=PolynomialRing(GF(2),'z');cache={};checked=0;norms=0;divpolys=0;seen=[]
def FF(pol):
 key=tuple(pol)
 if key not in cache:cache[key]=GF(2**(len(pol)-1),'a',modulus=R(pol))
 return cache[key]
def re(v,x):
 S=PolynomialRing(x.parent(),'q');return S(v[0])(x)/S(v[1])(x)
def ke(v,x,y):return sum((re(c,x)*y**i for i,c in enumerate(v)),x.parent().zero())
def ae(v,x,y,t):return sum((ke(c,x,y)*t**i for i,c in enumerate(v)),x.parent().zero())
def qcurve(E,xs):
 a1,a2,a3,a4,a6=E.a_invariants();b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6
 tv=sum(6*x*x+b2*x+b4 for x in xs);wv=sum(10*x**3+2*b2*x*x+3*b4*x+b6 for x in xs)
 return EllipticCurve(E.base_ring(),[a1,a2,a3,a4-5*tv,a6-b2*tv-7*wv])
def key3(v):
 F=v[0].parent();ans=[]
 for _ in range(F.degree()):
  ans.append(tuple(int(c)for a in v for c in a._vector_()));v=[a*a for a in v]
 return min(ans)
for d,idx,sourcekey,packs in data[5]:
 assert len(packs)==6;fields=[p[0][0]for p in packs];assert all(v==fields[0]for v in fields);F=FF(fields[0]);z=F.gen();S=PolynomialRing(F,'X');polys=[];ourkeys=[];nativekeys=[];es=[];first=None
 for (pol,coeffs),e,nkey in packs:
  vals=[S(c)(z)for c in coeffs];assert len(vals)==24;E=EllipticCurve(F,vals[:5]);P3=E(vals[5],vals[6]);P23=E(vals[7],vals[8]);assert P3!=E(0)and 3*P3==E(0)and P23!=E(0)and 23*P23==E(0)
  if first is None:first=vals[:9]
  else:assert first==vals[:9]
  ker=S(vals[9:12]);assert ker.degree()==2 and ker.is_monic();roots=ker.roots(multiplicities=False);assert len(roots)==2;pts=[p for x in roots for p in E.lift_x(x,all=True)];assert len(pts)==4 and all(5*p==E(0)and p!=E(0)for p in pts)
  target=qcurve(E,roots);assert list(target.a_invariants())==vals[12:17]
  for point,coordinates in [(P3,vals[17:19]),(P23,vals[19:21])]:
   xx=point[0]+sum((point+t)[0]-t[0]for t in pts);yy=point[1]+sum((point+t)[1]-t[1]for t in pts);assert [xx,yy]==coordinates;target(xx,yy)
  U=target(vals[17],vals[18]);V=target(vals[19],vals[20]);assert 3*U==target(0)and 23*V==target(0)
  projected=vals[21:24];x=projected[0];yyexpr=model[3][1][0];assert len(yyexpr)==2;y=(projected[1]-re(yyexpr[0],x))/re(yyexpr[1],x)
  p3expr=model[3][2];t=(projected[2]-ke(p3expr[0],x,y))/ke(p3expr[1],x,y)
  assert sum(re(c,x)*y**i for i,c in enumerate(model[0]))==0
  assert sum(ke(c,x,y)*t**i for i,c in enumerate(model[1]))==0
  assert all(ae(c,x,y,t)==v for c,v in zip(model[3],projected))
  assert ke(jf[0],x,y)==target.j_invariant();a1,a2,a3,a4,a6=target.a_invariants();xx,yy=U[0],U[1];aa3=2*yy+a3+a1*xx;aa2=3*xx+a2;aa4=3*xx*xx+2*xx*a2+a4-a1*yy;slope=aa4/aa3;assert aa2-slope*slope-a1*slope==0;tt=aa3/(a1+2*slope)**3;assert tt==t
  x23=[(j*V)[0]for j in range(1,12)];assert len(set(x23))==11;E23=qcurve(target,x23);assert ke(jf[1],x,y)==E23.j_invariant()
  ee=int(lcm([v.minimal_polynomial().degree()for v in projected]));assert ee==e;es.append(e);ourkeys.append(key3(projected));nativekeys.append(tuple(nkey));polys.append(ker);checked+=1
 assert len(set(polys))==6;prod=S.one()
 for p in polys:prod*=p
 psi=S(E.division_polynomial(5));assert prod.degree()==psi.degree()==12;assert prod.monic()==psi.monic(), (d,idx,prod,psi);divpolys+=1
 for i in range(6):
  for j in range(6):assert (ourkeys[i]==ourkeys[j])==(nativekeys[i]==nativekeys[j])
 total=0
 for k in set(ourkeys):
  ii=ourkeys.index(k);n=ourkeys.count(k);e=es[ii];assert d*n%e==0;total+=(d*n//e)*e;norms+=1
 assert total==6*d;seen.append([d,idx,len(set(ourkeys))])
assert checked==54 and divpolys==9
out={'status':'PASS','all_nine_base_source_classes_checked':True,'independent_Velu5_isogenies':checked,'complete_five_division_polynomials_recovered':divpolys,'independent_Velu23_quotient_j_checks':checked,'finite_moduli_coordinates_and_degrees_verified':checked,'positive_norm_summands':norms,'classes':seen,'input_sha256':{str(p.relative_to(H)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[B/'cover_results.m',B/'cover_model.m',B/'base_j_functions.m']},'seconds':time.monotonic()-start}
(B/'independent_moduli_isogenies.json').write_text(json.dumps(out,indent=2)+'\n');print('CONT69_INDEPENDENT_VELU_MODULI_NORMS_PASS',checked,norms,out['seconds'],flush=True)
