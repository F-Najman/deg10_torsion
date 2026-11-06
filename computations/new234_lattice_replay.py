from pathlib import Path
import sys,ast,json,time
from sage.all import GammaH,ZZ,QQ,Matrix,pari
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.maartens_sage_functions import integral_period_mapping,galois_orbits
P=Path('work/deg10/C234_resolution_20261030');tt=time.monotonic();a=json.loads((P/'rational_cusp_unit_lattice.json').read_text());G=GammaH(68,[35,-1]);M=G.modular_symbols();pm=integral_period_mapping(M);cusps=G.cusps();rat=a['rational_cusp_indices'];periods=Matrix(QQ,[M.coordinate_vector(M([cusps[i],cusps[rat[0]]]))*pm for i in rat]);DR=Matrix(ZZ,[[int(i==j)-int(i==0) for i in range(24)] for j in range(1,24)]);V=DR*periods;den=V.denominator();U=(V.stack(Matrix.identity(QQ,98))*den).change_ring(ZZ);ker=U.left_kernel().basis_matrix();coef=ker[:,:23];exact=(coef*DR).row_module();B=Matrix(ZZ,a['restricted_modular_unit_lattice']);assert exact==B.row_module();gram=B*B.transpose();result=pari(gram).qfminim(60);assert int(result[0])==0
squares=[]
for n1 in range(6,11):
 for n2 in range((24-n1+2)//3,n1+1):squares.append((11-n1)**2+n1-1+(11-n2)**2+n2-1)
assert max(squares)==60
out={'status':'PASS','rank':23,'rational_cusps':24,'integral_homology_rank':98,'period_kernel_rebuilt_exactly':True,'vectors_of_norm_at_most_60':0,'fibre_bound':60,'period_matrix':[[[int(x.numerator()),int(x.denominator())] for x in row] for row in periods.rows()],'relative_degree_zero_basis':[[int(x) for x in row] for row in DR.rows()]};(P/'independent_cusp_lattice.json').write_text(json.dumps(out,indent=2)+'\n');print('C234_INDEPENDENT_CUSP_LATTICE_COMPLETE',int(result[0]),time.monotonic()-tt,flush=True)
