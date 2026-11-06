from pathlib import Path
import ast,json,time
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
S=Path('work/deg10/C232_resolution_20261030');t0=time.monotonic()
def rd(n):return ast.literal_eval((S/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def dec(K,v):return sum(K(x)*K.gen()**i for i,x in enumerate(v))
def bc(E,P,Q):
 x=Q[0]-P[0];y=Q[1];assert E.a1()==E.a3()==0;return [x**3/y**2,(E.a2()+3*P[0])*x*x/y**2]
R=PolynomialRing(GF(3),'x');source=rd('source6.m');K=GF(3**6,'z',modulus=R(source[1]));z=K.gen();ky=rd('filterprobe6.m')[2][0][2];b=dec(K,ky[:6]);c=dec(K,ky[6:]);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P=E(0,0);Q=E(b,b);assert 32*Q==E(0) and 16*Q!=E(0) and P!=16*Q
reports=[]
for ell in PRIMES:
 a=rd('hecke'+str(ell)+'_geometry.m');gprod=PolynomialRing(K,'X').one();count=0;rows=[]
 for row in a[6]:
  e,pol,kg,gcoeff,egcoeff,qq,pp,values,degree,mult,place,norm=row;L=GF(3**(len(pol)-1),'w',modulus=R(pol));w=L.gen();kroot=dec(L,kg);assert sum(L(x)*kroot**i for i,x in enumerate(source[1]))==0
  emb=K.hom([kroot],L);EL=EllipticCurve(L,[0,emb(c),0,emb(b*(1-b-c)),0]);PL=EL(0,0);QL=EL(emb(b),emb(b));LX=PolynomialRing(L,'X');g=LX([dec(L,v) for v in gcoeff]);assert g.degree()==(ell-1)//2
  phi=EL.isogeny(g,check=True);assert phi.degree()==ell;EG=phi.codomain();PP=phi(PL);QQ=phi(QL);actual=bc(EG,PP,QQ);target=[dec(L,v) for v in values];assert actual==target
  EGM=EllipticCurve(L,[dec(L,v) for v in egcoeff]);QM=EGM([dec(L,v) for v in qq]);PM=EGM([dec(L,v) for v in pp]);assert 32*QM==EGM(0) and 16*QM!=EGM(0) and 2*PM==EGM(0) and PM!=16*QM;assert bc(EGM,PM,QM)==target
  ng=LX.one()
  for i in range(e):ng*=LX([v**(3**(6*i)) for v in g.list()])
  NK=PolynomialRing(K,'X')([dec(K,v) for v in norm]);assert ng==LX([emb(v) for v in NK.list()]);gprod*=NK;count+=e;assert degree*mult==6*e and degree==lcm([x.minpoly().degree() for x in target]);rows.append({'kernel_orbit_degree':e,'target_residue_degree':degree,'norm_multiplicity':mult,'exact_normalized_moduli_pair':True})
 div=E.division_polynomial(ell);div/=div.leading_coefficient();assert count==ell+1 and gprod==div;reports.append({'prime':ell,'kernel_count':count,'norms_partition_full_division_polynomial':True,'Sage_isogenies_agree_with_exported_moduli':True,'rows':rows});print('C232_INDEPENDENT_HECKE_PRIME',ell,count,time.monotonic()-t0,flush=True)
(S/('independent_hecke_'+TAG+'.json')).write_text(json.dumps({'status':'PASS','reports':reports},indent=2)+'\n');print('C232_INDEPENDENT_HECKE_COMPLETE',PRIMES,time.monotonic()-t0,flush=True)
