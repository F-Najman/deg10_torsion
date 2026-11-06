from pathlib import Path
import sys,json,time
from sage.all import GammaH,ModularSymbols,ZZ,identity_matrix,prod
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.cuspidal_classgroup import cuspidal_integral_structure
P=Path('work/deg10/C232_resolution_20261030');tt=time.monotonic();Gamma=GammaH(64,[15,-1]);M=ModularSymbols(Gamma,2,sign=0);S=cuspidal_integral_structure(M);assert S.rank()==26
star=(M.star_involution().matrix().restrict(S)-1).change_ring(ZZ);kill=star.transpose().row_module();records=[];mats=[]
for q in [3,5,7,11,13,17,19,23,29,31]:
 A=(M.hecke_matrix(q)-M.diamond_bracket_matrix(q)-q).restrict(S).change_ring(ZZ);kill+=A.transpose().row_module();mat=kill.basis_matrix();order=abs(int(mat.det())) if mat.is_square() else None;records.append({'prime':q,'rank':int(kill.rank()),'upper_bound_order':order});mats.append({'prime':q,'matrix':[[int(x) for x in row] for row in A.rows()]});print('C232_QUOTIENT13_HECKE_BOUND',q,order,time.monotonic()-tt,flush=True)
 if order==512000:break
bound=abs(int(kill.basis_matrix().det()));inv=[int(x) for x in kill.basis_matrix().elementary_divisors() if x!=1];assert bound%512000==0
out={'status':'PASS','level':64,'Delta':[1,15,17,31,33,47,49,63],'genus':13,'known_global_rational_cuspidal_order':512000,'rational_torsion_upper_bound':bound,'index_upper_bound':bound//512000,'finite_bound_invariants':inv,'full_integral_homology_rank':26,'uses_integral_boundary_kernel':True,'Xmu_operator_convention':'T_q - diamond(q) - q; Fricke gives the Q-isomorphic X_Delta model','integral_lattice_basis':[[[int(x.numerator()),int(x.denominator())] for x in row] for row in S.basis_matrix().rows()],'star_minus_identity':[[int(x) for x in row] for row in star.rows()],'hecke_matrices':mats,'steps':records,'final_relation_lattice':[[int(x) for x in row] for row in kill.basis_matrix().rows()]}
(P/'quotient13_rational_torsion_bound.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_QUOTIENT13_TORSION_BOUND_COMPLETE',bound,bound//512000,time.monotonic()-tt,flush=True)
