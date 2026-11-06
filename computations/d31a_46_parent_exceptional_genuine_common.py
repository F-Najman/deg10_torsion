from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time
C=Path('work/deg10/continuation_20261031a/C46');start=time.monotonic();Fp=GF(3);RP=PolynomialRing(Fp,'t')
def read(name):return ast.literal_eval((C/name).read_text().replace('[*','[').replace('*]',']'))
raw=read(f'parent_exceptional_hecke_s{SHARD}.m');native=read(f'parent_exceptional_native_s{SHARD}.m');assert len(raw)==len(native)
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def unitxy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1
 return[(s-r)/t,t/(s*s-s-r+1)]
def key2(vals):
 L=vals[0].parent();d=lcm([v.minpoly().degree()for v in vals]);found=False
 for bits in[[1,0],[0,1],[1,1],[1,2]]:
  g=bits[0]*vals[0]+bits[1]*vals[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found
 n=L.degree();vv=lambda v:vector(Fp,[L(v).polynomial()[i]for i in range(n)]);M=matrix(Fp,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(v))))for v in vals]]
result=[];nker=0;ncomponents=0
for crow,nrow in zip(raw,native):
 typ,idx,attempts=crow;ntyp,nidx,nattempts=nrow;assert(typ,idx)==(ntyp,nidx)and len(attempts)==len(nattempts);ar=[]
 for (q,comps,h0),(nq,ncomps,nh0)in zip(attempts,nattempts):
  assert(q,h0)==(nq,nh0)and len(comps)==len(ncomps);cr=[]
  for row,nr in zip(comps,ncomps):
   d,mul,pol,aa,qq,pp,dd,hec=row;nd,nmul,kp,kd,nt=nr;assert(d,mul)==(nd,nmul)
   K=GF(3**d,'a',modulus=RP(pol));a=K.gen();E0=EllipticCurve(K,[RP(v)(a)for v in aa]);Q0=E0([RP(v)(a)for v in qq]);assert 46*Q0==E0(0)and 23*Q0!=E0(0)and 2*Q0!=E0(0)
   assert key2(unitxy(E0,Q0))==kp and key2(unitxy(E0,q*Q0))==kd
   card=E0.cardinality();trace=3**d+1-card;M=matrix(GF(q),[[0,-3**d],[1,trace]]);B=identity_matrix(GF(q),2)
   for m in range(1,q*q+q+1):
    B*=M
    if B.is_one():break
   else:raise AssertionError('Frobenius order not found')
   L=GF(3**(d*m),'w');RL=PolynomialRing(L,'x');x=RL.gen();z=RL(pol).roots(multiplicities=False)[0];E=EllipticCurve(L,[RP(v)(z)for v in aa]);Q=E([RP(v)(z)for v in qq]);xs=E.division_polynomial(q).roots(multiplicities=False);pts=[T for sx in xs for T in E.lift_x(sx,all=True)];assert len(xs)==(q*q-1)//2 and len(pts)==q*q-1
   kernels={tuple(prod((x-(j*T)[0]for j in range(1,(q-1)//2+1))).list()):T for T in pts};assert len(kernels)==q+1
   actual=Counter()
   for coeff,T in kernels.items():
    assert q*T==E(0)and T!=E(0);ph=E.isogeny(RL(list(coeff)));ET=ph.codomain();QT=ph(Q);assert 46*QT==ET(0)and 23*QT!=ET(0)and 2*QT!=ET(0)
    nz=[j*T for j in range(1,q)];velu=(Q[0]+sum(((Q+Z)[0]-Z[0]for Z in nz),L(0)),Q[1]+sum(((Q+Z)[1]-Z[1]for Z in nz),L(0)));assert velu==QT.xy();actual[json.dumps(key2(unitxy(ET,QT)))]+=1;nker+=1
   expected=Counter()
   for (e,n,pd),(ne,nn,key)in zip(hec,nt):
    assert(e,n)==(ne,nn)and d*e==n*pd[0]and key[0]==pd[0];expected[json.dumps(key)]+=e
   assert actual==expected;ncomponents+=1;cr.append({'degree':int(d),'multiplicity':int(mul),'q':int(q),'torsion_splitting_degree':int(m),'geometric_kernels':int(q+1),'source_key':kp,'diamond_key':kd,'all_target_keys_and_norms_match':True})
  ar.append({'q':int(q),'native_h0':int(h0),'components':cr})
 assert attempts[-1][2]==0;result.append({'type':typ,'index':int(idx),'attempts':ar});print('C46_PARENT_GENUINE_CASE',SHARD,typ,idx,nker,time.monotonic()-start,flush=True)
out={'status':'PASS','case_count':len(result),'source_components':ncomponents,'geometric_kernels':nker,'manual_Velu_coordinates_verified':True,'all_native_prime_keys_and_norms_match':True,'coincident_targets_aggregated':True,'cases':result,'elapsed_seconds':time.monotonic()-start}
(C/f'parent_exceptional_genuine_s{SHARD}.json').write_text(json.dumps(out,indent=2)+'\n');print('C46_PARENT_EXCEPTIONAL_GENUINE_COMPLETE',SHARD,len(result),nker,time.monotonic()-start,flush=True)
