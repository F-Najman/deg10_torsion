from pathlib import Path
import sys,json,time
from functools import reduce
from sage.all import GammaH,ModularSymbols,ZZ,QQ,Matrix,Integers,Cusp,infinity
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.cuspidal_classgroup import cuspidal_integral_structure,intersection
from mdsage.maartens_sage_functions import integral_period_mapping,modular_unit_lattice,rational_modular_unit_lattice,galois_action,galois_orbits
P=Path('work/deg10/C232_resolution_20261030');tt=time.monotonic();Gamma=GammaH(64,[15,-1]);M=ModularSymbols(Gamma,2,sign=0);S=cuspidal_integral_structure(M);n=S.rank();assert n==26
star=(M.star_involution().matrix().restrict(S)-1).change_ring(ZZ);kill=star.transpose().row_module();ops=[]
for p in [3,5]:
 A=(M.hecke_matrix(p)-M.diamond_bracket_matrix(p)-p).restrict(S).change_ring(ZZ);kill+=A.transpose().row_module();ops.append(A)
K=kill.basis_matrix().transpose();pm=integral_period_mapping(M);periods=Matrix([M.coordinate_vector(M([c,infinity]))*pm for c in Gamma.cusps() if not c.is_infinity()]);MH=K.inverse().row_module(ZZ);CG=periods.stack(Matrix.identity(QQ,n)).row_module(ZZ);I=MH.intersection(CG);index=I.index_in(MH);assert index==2
L,D=modular_unit_lattice(Gamma);Lrat,Drat=rational_modular_unit_lattice(Gamma);C=D.quotient(L);ambient=ZZ**Gamma.ncusps();kernels=[];perms=[];cusps=Gamma.cusps()
for u in Integers(64).unit_gens():
 perm=[cusps.index(Gamma.reduce_cusp(galois_action(c,u,64))) for c in cusps];mat=Matrix([ambient.gen(i) for i in perm]);f=C.hom([C(i.lift()*(mat-1)) for i in C.gens()]);kernels.append(f.kernel());perms.append({'unit':int(u),'permutation':perm})
fixed=reduce(intersection,kernels);assert fixed.V()==Drat+L;assert fixed.cardinality()==512000;assert Drat.quotient(Lrat).cardinality()==512000

def ser(A):return [[[int(x.numerator()),int(x.denominator())] for x in row] for row in A.rows()]
def ints(A):return [[int(x) for x in row] for row in A.rows()]
out={'status':'PASS','method':'Integral homology with Hecke operators at 3 and 5 and complex conjugation; exact period and Galois lattices','homology_rank':26,'relation_lattice_order':abs(int(K.det())),'rational_cuspidal_order':512000,'rational_torsion_index_divides':int(index),'full_homology_basis':ser(S.basis_matrix()),'star_minus_identity':ints(star),'hecke_T_minus_diamond_minus_prime':[ints(A) for A in ops],'relation_matrix':ints(K),'periods_of_cusp_differences':ser(periods),'Hecke_bound_lattice':ser(MH.basis_matrix()),'geometric_cuspidal_lattice':ser(CG.basis_matrix()),'intersection_lattice':ser(I.basis_matrix()),'modular_unit_lattice':ints(L.basis_matrix()),'degree_zero_cusp_divisor_lattice':ints(D.basis_matrix()),'rational_cusp_divisor_lattice':ints(Drat.basis_matrix()),'Galois_fixed_cuspidal_lattice':ints(fixed.V().basis_matrix()),'Galois_permutations':perms,'global_cusp_orbits':[[cusps.index(c) for c in orb] for orb in galois_orbits(Gamma)]};(P/'integral_torsion_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_QUOTIENT13_INTEGRAL_CERTIFICATE_COMPLETE',index,fixed.cardinality(),time.monotonic()-tt,flush=True)
