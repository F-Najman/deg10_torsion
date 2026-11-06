from pathlib import Path
import sys,json,time
from sage.all import GammaH,ZZ,Matrix,Cusp,Integers
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.maartens_sage_functions import modular_unit_lattice,galois_action
P=Path('work/deg10/C236_resolution_20261030/cuspidal_boundary');tt=time.monotonic();G=GammaH(72,[37,-1]);assert G.genus()==41;cusps=G.cusps();assert len(cusps)==64
perm=[cusps.index(G.reduce_cusp(galois_action(c,Integers(72)(5),72))) for c in cusps];seen=set();orbs=[]
for i in range(64):
 if i in seen:continue
 o=[];j=i
 while j not in o:o.append(j);seen.add(j);j=perm[j]
 orbs.append(o)
assert sorted(map(len,orbs))==[1]*24+[2]*8+[6]*4
A=Matrix(ZZ,[[int(j in o) for j in range(64)] for o in orbs]);L,D=modular_unit_lattice(G);LC=L.intersection(A.row_module());assert LC.rank()==35;B=LC.basis_matrix();assert B*Matrix(ZZ,64,1,[1]*64)==0
# Coordinates in the disjoint Frobenius-orbit basis are integral.
coeff=A.solve_left(B);assert coeff in Matrix(ZZ,coeff.nrows(),coeff.ncols()).parent();coeff=coeff.change_ring(ZZ)
out={'status':'PASS','level':72,'Delta':[1,35,37,71],'genus':41,'prime':5,'geometric_cusps':64,'closed_cusp_degrees':list(map(len,orbs)),'Frobenius_orbits':orbs,'geometric_modular_unit_lattice':[[int(x) for x in r] for r in L.basis_matrix().rows()],'Frobenius_fixed_unit_lattice':[[int(x) for x in r] for r in B.rows()],'closed_unit_coordinates':[[int(x) for x in r] for r in coeff.rows()],'norm_bound':104}
(P/'lattice.json').write_text(json.dumps(out,indent=2)+'\n');(P/'lattice.m').write_text('CuspUnitBasis:='+repr(out['Frobenius_fixed_unit_lattice'])+';\nCuspOrbits:='+repr(orbs)+';\n');print('C236_CUSP_LATTICE_COMPLETE',LC.rank(),B.ncols(),time.monotonic()-tt,flush=True)
