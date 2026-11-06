from sage.all import *
from pathlib import Path
from collections import Counter
import json,time
P=Path('work/deg10/continuation_20261031b/C56');D=P/'kummer_rr';T=time.monotonic();Fp=GF(3);R=PolynomialRing(Fp,'t');cids=[271,327,491,578,65,818,837,978,1009]
def rd(p):return json.loads(p.read_text().replace('[*','[').replace('*]',']'))
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;s=b4/b3;n1=2*s+a1;n2=b2-s*s-a1*s;return[-n2**3/b3**2,1-n1*n2/b3]
def key(v):
 K=v[0].parent();n=K.degree();d=lcm([z.minpoly().degree()for z in v]);bits=next(z for z in[[1,0],[0,1],[1,1],[1,2]]if(z[0]*v[0]+z[1]*v[1]).minpoly().degree()==d);g=bits[0]*v[0]+bits[1]*v[1];M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(n)]for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vector(Fp,[K(z).polynomial()[j]for j in range(n)]))))for z in v];return[int(d),bits,list(map(int,g.minpoly().list())),co]
def physical(v):
 K=GF(3**(len(v[0])-1),'a',modulus=R(v[0]));a=K.gen();E=EllipticCurve(K,[R(c)(a)for c in v[1]]);Q=E([R(c)(a)for c in v[2]]);assert 56*Q==E(0)and 28*Q!=E(0)and 8*Q!=E(0)and key(bc(E,Q))==v[3];return K,E,Q
def manual(E,Q,phi):
 ker=phi.kernel_polynomial();RR=ker.parent();z=RR.gen();m=ker.degree();ss=[E.base_ring()(m)]
 for k in range(1,m):ss.append(-sum((ker[m-i]*ss[k-i]for i in range(1,k)),E.base_ring()(0))-k*ker[m-k])
 def tr(v):v=RR(v)%ker;return sum((v[i]*ss[i]for i in range(m)),E.base_ring()(0))
 x,y=Q.xy();a2=E.a2();g=z**3+a2*z*z+E.a4()*z+E.a6();inv=(x-z).inverse_mod(ker);px=2*(y*y+g)*inv**2-2*a2-2*x-4*z;py=-2*y+2*y*(2*x+a2+z)*inv-2*y*(y*y+3*g)*inv**3;return[x+tr(px),y+tr(py)]
rows={}
for ci in cids:
 q=11 if ci==65 else 5;op=rd(D/f'direct_case{ci}_A{q}_operator_snapshot.m')
 for ix in op[2]:
  row=op[0][ix-1];idx=(int(row[0]),int(row[6]))
  if idx in rows:assert rows[idx]==row
  else:rows[idx]=row
assert{v[0]for v in rows if v[1]==5}=={271,327,491,578,1,818,837,3,834,865}and{v[0]for v in rows if v[1]==11}=={65};out=[];nker=0
for(id,q),row in sorted(rows.items()):
 d=row[1];K,E,Q=physical(row[8]);assert K.degree()==d and row[3][0]==d and row[3][4]==key(bc(E,Q));KD,ED,PD=physical(row[9]);assert row[5][4]==key(bc(E,q*Q))==key(bc(ED,PD));orbit=lcm([h[0]for h in row[7]]);L=GF(3**(d*orbit),'u');RL=PolynomialRing(L,'x');u=RL(row[8][0]).roots(multiplicities=False)[0];emb=lambda z:sum((L(z.polynomial()[i])*u**i for i in range(d)),L(0));ai=list(map(emb,E.a_invariants()));xp,yp=map(emb,Q.xy());EL=EllipticCurve(L,[0,ai[1]+ai[0]**2/4,0,ai[3]+ai[0]*ai[2]/2,ai[4]+ai[2]**2/4]);QL=EL(xp,yp+(ai[0]*xp+ai[2])/2);assert key(bc(EL,QL))==row[3][4];phis=EL.isogenies_prime_degree(q);assert len(phis)==q+1;kers={tuple(ph.kernel_polynomial().list())for ph in phis};assert len(kers)==q+1;unseen=set(kers);actual=Counter();obs=[]
 while unseen:
  coeff=min(unseen,key=str);orb=[];z=coeff
  while z not in orb:assert z in kers;orb.append(z);z=tuple(c**(3**d)for c in z)
  assert z==coeff and set(orb)<=unseen;unseen-=set(orb);ph=EL.isogeny(RL(list(coeff)));assert ph.scaling_factor()==1;QT=ph(QL);assert list(QT.xy())==manual(EL,QL,ph);assert 56*QT==ph.codomain()(0)and 28*QT!=ph.codomain()(0)and 8*QT!=ph.codomain()(0);kk=key(bc(ph.codomain(),QT));assert d*len(orb)%kk[0]==0;actual[json.dumps(kk)]+=len(orb);obs.append([len(orb),d*len(orb)//kk[0],kk[0],kk])
 expected=Counter()
 for e,n,desc,phys in row[7]:
  KT,ET,PT=physical(phys);kk=key(bc(ET,PT));assert kk==desc[4]and kk[0]==desc[0]and d*e==n*kk[0]and KT.degree()==d*e;expected[json.dumps(kk)]+=e
 assert actual==expected and sum(h[0]for h in row[7])==q+1;nker+=q+1;out.append({'source_ID':id,'hecke':q,'source_degree':d,'source_key':row[3][4],'diamond_key':row[5][4],'kernel_orbits':obs,'geometric_kernels':q+1});print('D31B_C56_SIEVE_GENUINE_POINT_PASS',id,q,len(phis),time.monotonic()-T,flush=True)
assert nker==72 and len(out)==11;result={'status':'PASS','source_level':56,'prime':3,'source_contexts':11,'cyclic5_kernels':60,'cyclic11_kernels':12,'total_geometric_kernels':72,'all_source_and_diamond_keys_match':True,'all_positive_physical_Tatebc_keys_and_norms_match':True,'complete_Sage_kernel_census_and_independent_pairwise_Velu_traces':True,'observations':out,'elapsed_seconds':time.monotonic()-T};(P/'sieve_independent_genuine_isogenies.json').write_text(json.dumps(result,indent=2)+'\n');print('D31B_C56_SIEVE_GENUINE_ALL_COMPLETE',11,72,time.monotonic()-T,flush=True)
