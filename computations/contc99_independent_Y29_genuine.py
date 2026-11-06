"""Physical point9/cyclic11 keys for the genus-29 residual obstruction."""
from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,time
ROOT=Path('work/deg10/continuation_20261030c/C99');START=time.monotonic();F2=GF(2);R2=PolynomialRing(F2,'T')
def read(n):return ast.literal_eval((ROOT/n).read_text().replace('[*','[').replace('*]',']'))
source=read('source99.m');pol=next(v[1]for v in source[0]if v[0]==10)
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y
 n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return [-n2**3/b3**2,-(n1*n2-b3)/b3]
def key2(v):
 K=v[0].parent();n=K.degree();d=lcm([a.minpoly().degree()for a in v]);found=False
 for bits in[[1,0],[0,1],[1,1]]:
  g=bits[0]*v[0]+bits[1]*v[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found
 vv=lambda a:vector(F2,[K(a).polynomial()[i]for i in range(n)])
 M=matrix(F2,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(a))))for a in v]]
def Ykey(E,P):
 assert 99*P==E(0)and 33*P!=E(0)and 9*P!=E(0)
 Q=11*P;S=9*P;assert Q.order()==9 and S.order()==11
 b,c=pointbc(E,Q);t=c*c/(b-c)
 assert c==t*t*(t-1)and b==c*(t*t-t+1)
 RP=PolynomialRing(E.base_ring(),'v');v=RP.gen();ker=prod(v-(j*S)[0]for j in range(1,6));assert ker.degree()==5 and ker.is_squarefree()
 phi=E.isogeny(ker);assert phi.degree()==11 and phi(S)==phi.codomain()(0)
 return key2([t,phi.codomain().j_invariant()])
L=GF(2**20,'u');RP=PolynomialRing(L,'x');x=RP.gen();a=RP(pol).roots(multiplicities=False)[0]
d,ky=source[2][2];assert d==10;b=R2(ky[:d])(a);c=R2(ky[d:])(a);E=EllipticCurve(L,[1-c,-b,-b,0,0]);P=E(0,0)
sources={str(u):Ykey(E,u*P)for u in[1,2,5,10]};assert sources['10']==sources['1']and len({json.dumps(v)for v in sources.values()})==3
xs=E.division_polynomial(5).roots(multiplicities=False);pts=[T for xx in xs for T in E.lift_x(xx,all=True)]
assert len(xs)==12 and len(pts)==24 and all(T.order()==5 for T in pts)
kernels={tuple(((x-T[0])*(x-(2*T)[0])).list()):T for T in pts};assert len(kernels)==6
counts=Counter();obs={}
for coeff,T in kernels.items():
 ker=RP(list(coeff));phi=E.isogeny(ker);ET=phi.codomain();PT=phi(P);assert PT.order()==99
 nonzero=[j*T for j in range(1,5)]
 assert (P[0]+sum(((P+S)[0]-S[0]for S in nonzero),L(0)),P[1]+sum(((P+S)[1]-S[1]for S in nonzero),L(0)))==PT.xy()
 pair=[Ykey(ET,PT),Ykey(ET,2*PT)];k=json.dumps(pair);counts[k]+=1;obs[k]=pair
rows=[[counts[k],obs[k][0],obs[k][1]]for k in sorted(counts)]
assert sorted(v[0]for v in rows)==[1,1,2,2]and sorted(v[1][0]for v in rows)==[10,10,20,20]
result={'status':'PASS','case':3,'model':'X1(9) fibre X0(11)','coordinate_keys':'[t9,j(E/C11)], t9=c9^2/(b9-c9)','source_phase_keys':sources,'paired_Hecke_keys_weighted_by_actual_kernel_orbits':rows,'actual_geometric_cyclic5_kernels':6,'explicit_5_Velu_coordinate_checks':6,'actual_cyclic11_isogenies_constructed':16,'diamond10_equals_identity_on_moduli':True,'elapsed_seconds':time.monotonic()-START}
(ROOT/'independent_Y29_moduli_keys.json').write_text(json.dumps(result,indent=2)+'\n')
print('C99_INDEPENDENT_Y29_GENUINE_MARKED_KEYS_COMPLETE',3,6,[v[0]for v in rows],time.monotonic()-START,flush=True)
