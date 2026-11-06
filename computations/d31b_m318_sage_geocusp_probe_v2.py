from sage.all import *
from pathlib import Path
import time,json
T=time.monotonic();M=ModularSymbols(GammaH(54,[19]),2,sign=0);print('M318_SAGE_MODSYMS',M.dimension(),time.monotonic()-T,flush=True);C=M.cuspidal_submodule();E=M.eisenstein_submodule();print('M318_SAGE_CUSP_EIS',C.dimension(),E.dimension(),flush=True)
for N,ob in [('M',M),('C',C),('E',E)]:
 print(N,[s for s in dir(ob)if any(x in s for x in ['integral','projection','cusps','boundary','manin'])],flush=True)
print('M318_CUSPIDAL_LATTICE',C.integral_structure().rank(),flush=True);print('M318_RELATIVE_LATTICE',M.integral_structure().rank(),flush=True)
B=matrix(QQ,list(C.basis_matrix().rows())+list(E.basis_matrix().rows()));assert B.is_square()and B.is_invertible();pr=B.inverse()[:,:C.dimension()];L=C.integral_structure().basis_matrix();ML=M.integral_structure().basis_matrix();print('LATTICE_SHAPES',L.dimensions(),ML.dimensions(),pr.dimensions(),flush=True)
# Both integral structures are embedded in the ambient modular symbol space.
CL=L*pr;assert CL.is_square();P=ML*pr*CL.inverse();den=P.denominator();A=matrix(ZZ,den*P);H=matrix(ZZ,list(A.rows())+list(den*identity_matrix(ZZ,C.dimension()).rows()));S=H.row_module(ZZ).basis_matrix();index=den**C.dimension()/abs(S.det());invs=[int(den//v)for v in S.smith_form()[0].diagonal()if den//v>1];rep={'status':'PASS','GammaH':[54,[19]],'relative_dimension':M.dimension(),'cuspidal_dimension':C.dimension(),'geometric_cuspidal_order':int(index),'invariants':sorted(invs),'denominator':int(den),'elapsed_seconds':time.monotonic()-T};Path('work/deg10/continuation_20261031b/M3N18/sage_geometric_cuspidal_group.json').write_text(json.dumps(rep,indent=2)+'\n');print('D31B_M318_SAGE_GEOMETRIC_CUSPIDAL_COMPLETE',rep,flush=True)
