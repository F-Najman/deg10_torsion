from pathlib import Path
import json
from functools import reduce
from sage.all import Matrix,QQ,ZZ
P=Path('work/deg10/C232_resolution_20261030');a=json.loads((P/'integral_torsion_certificate.json').read_text());n=a['homology_rank'];assert n==26

def mat(v):return Matrix(QQ,[[QQ(x)/QQ(y) for x,y in row] for row in v])
star=Matrix(ZZ,a['star_minus_identity']);ops=[Matrix(ZZ,v) for v in a['hecke_T_minus_diamond_minus_prime']];kill=star.transpose().row_module()
for A in ops:kill+=A.transpose().row_module()
K=kill.basis_matrix().transpose();assert K==Matrix(ZZ,a['relation_matrix']);MH=K.inverse().row_module(ZZ);assert MH==mat(a['Hecke_bound_lattice']).row_module(ZZ);periods=mat(a['periods_of_cusp_differences']);CG=periods.stack(Matrix.identity(QQ,n)).row_module(ZZ);assert CG==mat(a['geometric_cuspidal_lattice']).row_module(ZZ);I=MH.intersection(CG);assert I==mat(a['intersection_lattice']).row_module(ZZ);assert I.index_in(MH)==2
L=Matrix(ZZ,a['modular_unit_lattice']).row_module();D=Matrix(ZZ,a['degree_zero_cusp_divisor_lattice']).row_module();Dr=Matrix(ZZ,a['rational_cusp_divisor_lattice']).row_module();nc=D.degree();V=ZZ**nc;C=D.quotient(L);kernels=[]
for r in a['Galois_permutations']:
 permutation=r['permutation'];assert sorted(permutation)==list(range(nc));m=Matrix([V.gen(i) for i in permutation]);f=C.hom([C(v.lift()*(m-1)) for v in C.gens()]);kernels.append(f.kernel())
def intersection(A,B):assert A.W()==B.W();return A.V().intersection(B.V()).quotient(A.W())
fixed=reduce(intersection,kernels);assert fixed.V()==Dr+L;assert fixed.cardinality()==512000;assert fixed.V()==Matrix(ZZ,a['Galois_fixed_cuspidal_lattice']).row_module();assert Dr.intersection(L).index_in(Dr)==512000
out={'status':'PASS','integral_homology_rank':n,'Hecke_bound_order':abs(int(K.det())),'Hecke_bound_mod_geometric_cuspidals_order':2,'Galois_fixed_geometric_cuspidal_order':512000,'Galois_fixed_classes_equal_rational_cuspidal_divisors':True,'rational_Jacobian_index_over_rational_cuspidals_divides':2};(P/'independent_integral_lattices.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_INDEPENDENT_INTEGRAL_LATTICES_COMPLETE',out,flush=True)
