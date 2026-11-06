from pathlib import Path
import sys,json,time
from sage.all import GammaH
sys.path.insert(0,str(Path('work/deg10/shared/mdsage').resolve()))
from mdsage.cuspidal_classgroup import upper_bound_index_cusps_in_JG_torsion,cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup
from mdsage.maartens_sage_functions import rational_cuspidal_classgroup
P=Path('work/deg10/C232_resolution_20261030');G=GammaH(64,[15,-1]);tt=time.monotonic();C=rational_cuspidal_classgroup(G);print('C232_QUOTIENT13_ANALYTIC_CUSPS',C.invariants(),C.cardinality(),time.monotonic()-tt,flush=True);assert C.cardinality()==512000
idx=upper_bound_index_cusps_in_JG_torsion(G,512000,bound=20);print('C232_QUOTIENT13_CUSP_INDEX_BOUND',idx,time.monotonic()-tt,flush=True);extra=cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(G);print('C232_QUOTIENT13_CUSP_DESCENT_QUOTIENT',extra,time.monotonic()-tt,flush=True)
out={'status':'PASS','rational_cuspidal_invariants':[int(x) for x in C.invariants()],'rational_cuspidal_order':int(C.cardinality()),'index_bound_over_geometric_cuspidals':int(idx),'Galois_invariant_cuspidal_classes_mod_rational_cuspidal_divisors':[int(x) for x in extra]};(P/'quotient13_cuspidal_saturation.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_QUOTIENT13_CUSPIDAL_SATURATION_COMPLETE',time.monotonic()-tt,flush=True)
