from pathlib import Path
import sys,json,time,os
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from sage.all import GammaH,ModularSymbols,QQ,ZZ,matrix,prime_range
from mdsage.cuspidal_classgroup import cuspidal_integral_structure
P=Path('work/deg10/continuation_20261029b');t=time.monotonic()
M=ModularSymbols(GammaH(54,[19]),2,0)
S=cuspidal_integral_structure(M)
assert S.dimension()==20
C=M.star_involution().matrix().restrict(S).change_ring(ZZ)
I=matrix.identity(ZZ,20)
assert C*C==I
plus=(C-I).transpose().row_module();minus=(C+I).transpose().row_module();packet=[]
print('INTEGRAL_HOMOLOGY_READY',M.dimension(),S.dimension(),time.monotonic()-t,flush=True)
for q in prime_range(5,40):
 if 54%q==0:continue
 T=M.hecke_matrix(q).restrict(S).change_ring(ZZ)
 D=M.diamond_bracket_matrix(q).restrict(S).change_ring(ZZ)
 assert C*T==T*C and C*D==D*C
 e=1 if q%3==1 else -1
 plus+=(T-I-q*D).transpose().row_module()
 minus+=(T-e*(I+q*D)).transpose().row_module()
 if plus.rank()<20 or minus.rank()<20:continue
 np=abs(plus.matrix().det());nm=abs(minus.matrix().det())
 row={'q':int(q),'plus_order':int(np),'minus_order':int(nm),'total_v3':int(np.valuation(3)+nm.valuation(3)),'T':[list(map(int,r)) for r in T.rows()],'D':[list(map(int,r)) for r in D.rows()]}
 packet.append(row)
 (P/'homology318_sage.json').write_text(json.dumps({'integral_congruence_subgroup':'GammaH(54,[19])','complex_conjugation':[list(map(int,r)) for r in C.rows()],'rows':packet},indent=2)+'\n')
 print('REAL_EISENSTEIN_TORSION_BOUND',q,np,nm,row['total_v3'],time.monotonic()-t,flush=True)
 if row['total_v3']==10:break
print('HOMOLOGY318_SAGE_COMPLETE',time.monotonic()-t,flush=True)
