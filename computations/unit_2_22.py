import sys,os,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'shared/mdsage'))
from sage.all import *
from mdsage.maartens_sage_functions import Gamma11,has_modular_unit_of_degree,rational_modular_unit_lattice,galois_orbits
m,n,d=2,22,10
G=Gamma11(m,n//m)
print('DECIC_UNIT_START',m,n,'genus',G.genus(),'relative_degree',d,flush=True)
ok,v=has_modular_unit_of_degree(G,d,l2_step=2)
if ok:
 L,D=rational_modular_unit_lattice(G)
 assert v in L and sum(v)==0 and sum(abs(x) for x in v)==2*d
 out={'m':m,'n':n,'relative_degree':d,'base_degree':1 if m<=2 else 2,'divisor':[int(x) for x in v],'cusps':[str(c) for c in G.cusps()],'lattice_basis':[[int(x) for x in row] for row in L.basis_matrix()],'degree_verified':True}
 Path(__file__).resolve().parents[1].joinpath('occurrences/unit_2_22.json').write_text(json.dumps(out,indent=2)+'\n')
 print('DECIC_UNIT_SUCCESS',m,n,'degree',d,'divisor',list(v),flush=True)
else: print('DECIC_UNIT_NOT_FOUND',m,n,'This does not exclude degree-d points',flush=True)
print('DECIC_UNIT_SEARCH_COMPLETE',flush=True)
