import sys,json,hashlib
from pathlib import Path
H=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(H/'shared/mdsage'))
from sage.all import *
from mdsage.maartens_sage_functions import Gamma11,rational_modular_unit_lattice,galois_orbits
out=[]
for p in sorted((H/'occurrences').glob('unit_*.json')):
 d=json.loads(p.read_text());m,n=d['m'],d['n'];G=Gamma11(m,n//m)
 assert d['cusps']==[str(c) for c in G.cusps()]
 L,D=rational_modular_unit_lattice(G);v=vector(ZZ,d['divisor'])
 assert v in L and sum(v)==0 and sum(max(0,x) for x in v)==d['relative_degree']
 assert d['base_degree']*(d['relative_degree'])==10
 assert d['lattice_basis']==[[int(x) for x in row] for row in L.basis_matrix()]
 print('DECIC_UNIT_REPLAY_PASS',m,n,'degree',d['relative_degree'],flush=True);out.append([m,n])
print('DECIC_UNIT_REPLAY_COMPLETE count='+str(len(out)),flush=True)
