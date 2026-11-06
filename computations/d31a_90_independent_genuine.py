from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time,hashlib
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C90';start=time.monotonic();k=GF(2);R2=PolynomialRing(k,'T')
def read(n):return ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'))
source=read('source45.m');primary=sorted(read('A7_RR_part1.m')+read('A7_RR_part2.m'),key=lambda r:r[0])
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def unitxy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def key2(vals):
 L=vals[0].parent();d=lcm([v.minpoly().degree()for v in vals]);found=False
 for bits in[[1,0],[0,1],[1,1]]:
  g=bits[0]*vals[0]+bits[1]*vals[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found;n=L.degree();vv=lambda v:vector(k,[L(v).polynomial()[i]for i in range(n)]);M=matrix(k,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(v))))for v in vals]]
def frobkey(vals,d):return min(tuple(int(v.polynomial()[i])for v in[vals[0]**(2**j),vals[1]**(2**j)]for i in range(d))for j in range(d))

curves=read('inputs/all_curves.m');units=[u for u in range(1,23)if gcd(u,45)==1];d=10;pol=source[0][0][1];K=GF(2**d,'a',modulus=R2(pol));a=K.gen();keys=set();orbs=set();levels=0;eligible=0
for dd,poly,rows in curves:
 if dd!=10:continue
 for jj,aa,card in rows:
  if card%45:continue
  E=EllipticCurve(K,[R2(v)(a)for v in aa]);assert E.cardinality()==card;pts=[P for P in E.points()if P!=E(0)and 45*P==E(0)and 15*P!=E(0)and 9*P!=E(0)]
  if not pts:continue
  eligible+=1;levels+=len(pts)//2;seen=set()
  for P in pts:
   if P in seen:continue
   seen.update(u*P for u in range(1,45)if gcd(u,45)==1)
   ks={frobkey(pointbc(E,u*P),d)for u in units};keys.update(ks);orbs.add(min(ks))
assert len(keys)==108 and len(orbs)==10 and levels==1080 and eligible==70;assert sorted(orbs)==sorted(tuple(ky)for dd,ky in source[2])
print('D31A_C90_INDEPENDENT_SOURCE_PASS',108,10,flush=True)
w=json.loads((C/'inputs/known_decic45.json').read_text());fp=[k(QQ(v))for v in w['field']['poly_coeffs']];f=R2(fp);assert f.is_irreducible();root=PolynomialRing(K,'X')(fp).roots(multiplicities=False)[0];b=sum((K(QQ(v))*root**i for i,v in enumerate(w['curve']['tate_b_coeffs'])),K(0));c=sum((K(QQ(v))*root**i for i,v in enumerate(w['curve']['tate_c_coeffs'])),K(0));E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0);assert 45*P==E(0)and 15*P!=E(0)and 9*P!=E(0)
ks={frobkey(pointbc(E,u*P),d)for u in units};knownci=next(i for i,(dd,ky)in enumerate(source[2],1)if tuple(ky)==min(ks));assert knownci in[1,10];print('D31A_C90_KNOWN_DECIC_SOURCE_CASE',knownci,'orbit',len(ks),flush=True)
(C/'known_lift_binding.json').write_text(json.dumps({'status':'PASS','known_decic_source_case':knownci,'known_diamond_orbit_size':len(ks)},indent=2)+'\n')
red=[];bound=0
for p0 in prime_range(3,252):
 if 45%p0==0:continue
 Fp=GF(p0);Rp=PolynomialRing(Fp,'v');fp=Rp([Fp(QQ(v))for v in w['field']['poly_coeffs']])
 for r in fp.roots(multiplicities=False):
  if fp.derivative()(r)==0:continue
  b=Rp([Fp(QQ(v))for v in w['curve']['tate_b_coeffs']])(r);c=Rp([Fp(QQ(v))for v in w['curve']['tate_c_coeffs']])(r)
  try:E=EllipticCurve(Fp,[1-c,-b,-b,0,0])
  except ArithmeticError:continue
  n=E.cardinality();assert n%45==0;bound=gcd(bound,n);red.append([int(p0),int(r),int(n),int(bound)])
  if bound==45:break
 if bound==45:break
assert bound==45
out=[]
for ci,(d,ky)in enumerate(source[2],1):
 pol=next(pol for dd,pol in source[0]if dd==d);K=GF(2**d,'a',modulus=R2(pol));a=K.gen();b=R2(ky[:d])(a);c=R2(ky[d:])(a);E0=EllipticCurve(K,[1-c,-b,-b,0,0]);card=E0.cardinality();q=2**d;trace=q+1-card;M=matrix(GF(7),[[0,-q],[1,trace]]);Mm=identity_matrix(GF(7),2)
 for m in range(1,337):
  Mm*=M
  if Mm.is_one():break
 else:raise AssertionError('Frobenius order bound')
 L=GF(2**(d*m),'w');RL=PolynomialRing(L,'x');x=RL.gen();root=RL(pol).roots(multiplicities=False)[0];b=R2(ky[:d])(root);c=R2(ky[d:])(root);E=EllipticCurve(L,[1-c,-b,-b,0,0]);P=E(0,0);xs=E.division_polynomial(7).roots(multiplicities=False);pts=[T for xx in xs for T in E.lift_x(xx,all=True)];assert len(xs)==24 and len(pts)==48
 kernels={tuple(prod(x-(j*T)[0]for j in range(1,4)).list()):T for T in pts};assert len(kernels)==8;actual=Counter();obs={}
 for coeff,T in kernels.items():
  assert 7*T==E(0)and T!=E(0);phi=E.isogeny(RL(list(coeff)));EE=phi.codomain();Q=phi(P);assert 45*Q==EE(0)and 15*Q!=EE(0)and 9*Q!=EE(0)
  nz=[j*T for j in range(1,7)];velu=(P[0]+sum(((P+Z)[0]-Z[0]for Z in nz),L(0)),P[1]+sum(((P+Z)[1]-Z[1]for Z in nz),L(0)));assert velu==Q.xy();kyy=key2(unitxy(EE,Q));tag=json.dumps(kyy);actual[tag]+=1;obs[tag]=kyy
 expected=Counter();pobs=[]
 for orbit,n,pd,coords,polT in primary[ci-1][4]:
  A=GF(2**(len(polT)-1),'v',modulus=R2(polT));v=A.gen();vv=[sum((A(c)*v**i for i,c in enumerate(row)),A(0))for row in coords];key=key2(vv);assert d*orbit==n*pd[0]and key[0]==pd[0];expected[json.dumps(key)]+=orbit;pobs.append([int(orbit),int(n),key])
 assert actual==expected;out.append({'case':ci,'degree':d,'source_key':key2(unitxy(E,P)),'diamond7_key':key2(unitxy(E,7*P)),'target_keys':pobs,'native_h0':primary[ci-1][5]});print('D31A_C90_GENUINE_SEVEN_CASE',ci,d,m,[v[0]for v in pobs],time.monotonic()-start,flush=True)
result={'status':'PASS','independent_source_closed_places':108,'independent_source_diamond_orbits':10,'eligible_curve_count':70,'known_decic_source_case':knownci,'known_decic_full_torsion_order':45,'known_decic_odd_prime_bounds':red,'source_cases':out,'actual_geometric7_kernels':80,'manual_Velu_marked_coordinate_checks':80,'elapsed_seconds':time.monotonic()-start}
(C/'independent_genuine_and_known_lift.json').write_text(json.dumps(result,indent=2)+'\n');print('D31A_C90_INDEPENDENT_GENUINE_AND_LIFT_COMPLETE',80,time.monotonic()-start,flush=True)
