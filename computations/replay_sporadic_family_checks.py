import sys,json
from pathlib import Path
H=Path(__file__).resolve().parents[1];sys.path.insert(0,str(H/'shared/mdsage'))
from sage.all import *
from mdsage.maartens_sage_functions import Gamma11,rational_modular_unit_lattice
for m,n in [(1,28),(1,30),(2,18)]:
 d=json.loads((H/f'notes/sporadic_check_unit_{m}_{n}.json').read_text());G=Gamma11(m,n//m);L,D=rational_modular_unit_lattice(G);v=vector(ZZ,d['divisor']);assert v in L and sum(v)==0 and sum(max(0,x) for x in v)==10
 print('DEGREE10_INFINITE_FAMILY_REPLAY_PASS',m,n,flush=True)
print('SPORADIC_FAMILY_CHECKS_COMPLETE',flush=True)
