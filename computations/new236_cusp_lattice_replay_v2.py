from pathlib import Path
import json,sys,ast,time
from sage.all import ZZ,QQ,Matrix,GammaH,pari,vector
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.maartens_sage_functions import integral_period_mapping
P=Path('work/deg10/C236_resolution_20261030/cuspidal_boundary');a=json.loads((P/'lattice.json').read_text());tt=time.monotonic();G=GammaH(72,[37,-1]);M=G.modular_symbols();pm=integral_period_mapping(M);cs=G.cusps();periods=Matrix([M.coordinate_vector(M([c,cs[-1]]))*pm for c in cs]);assert cs[-1].is_infinity()
A=Matrix(ZZ,[[int(j in o) for j in range(64)] for o in a['Frobenius_orbits']]);B=Matrix(ZZ,a['Frobenius_fixed_unit_lattice']);assert B.row_module().is_submodule(A.row_module());assert all(sum(row)==0 for row in B.rows())
# Recompute the complete period kernel on Frobenius-fixed cusp divisors.
C=A*periods;degs=Matrix(QQ,A.nrows(),1,[len(o) for o in a["Frobenius_orbits"]]);T=C.augment(degs).stack(Matrix.identity(QQ,C.ncols()).augment(Matrix(QQ,C.ncols(),1)));K=(T*T.denominator()).change_ring(ZZ).kernel().basis_matrix();coef=K[:,:A.nrows()];rebuilt=(coef*A).row_module();assert rebuilt==B.row_module()
BB=B.LLL();Q=BB*BB.transpose();res=pari(Q).qfminim(105,None,2);V=Matrix(ZZ,res[2]);assert int(res[0])==2*V.ncols();count=0;hits=[];hist={}
for col in V.columns():
 v=col*BB;norm=int(v*v)
 if norm>104:continue
 count+=1;pos=sum(max(0,int(x)) for x in v);hist[pos]=hist.get(pos,0)+1
 if pos==10:hits.append([int(x) for x in v])
assert count==24552 and not hits
# Bound for the two rational fibres with largest distinct rational support.
bounds=[]
for n1 in range(4,11):
 n2=(24-n1+4)//5;assert n2<=n1;bounds.append(((11-n1)**2+n1-1)+((11-n2)**2+n2-1))
assert max(bounds)==104
out={'status':'PASS','period_kernel_rebuilt':True,'rank':35,'geometric_cusps':64,'prime':5,'norm_bound':104,'PARI_short_vector_pairs':count,'positive_degree_histogram':hist,'positive_degree_10_vectors':0,'two_fibre_norm_bounds':bounds};(P/'independent_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C236_CUSP_INDEPENDENT_REPLAY_COMPLETE',count,hist,time.monotonic()-tt,flush=True)
