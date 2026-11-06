from pathlib import Path
import json,time
from sage.all import GammaH,ModularSymbols,ZZ
P=Path('work/deg10/C236_resolution_20261030/geometry');t=time.monotonic();Gamma=GammaH(72,[35,-1]);M=ModularSymbols(Gamma,2,sign=0).cuspidal_subspace();assert M.dimension()==82
classes=[];seen=set();H={1,35,37,71}
for a in range(1,72):
 if Gamma.base_ring if False else False:pass
 if any(a%p==0 for p in [2,3]) or a in seen:continue
 classes.append(a);seen|={a*h%72 for h in H}
rows=[]
for d in [8,9,72]:
 W=M.atkin_lehner_operator(d).matrix()
 for a in classes:
  A=M.diamond_bracket_matrix(a)*W;involution=A*A==1;tr=int(A.trace());row={'Atkin_Lehner_divisor':d,'diamond_twist':a,'homology_trace':tr,'is_involution':involution,'fixed_points_if_involution':2-tr if involution else None,'genus_quotient_if_involution':int((82+tr)/4) if involution else None};rows.append(row);print(row,flush=True)
(P/'fricke_traces.json').write_text(json.dumps({'status':'PASS','genus':41,'rows':rows},indent=2)+'\n');print('C236_FRICKE_TRACE_PROBE_COMPLETE',time.monotonic()-t,flush=True)
