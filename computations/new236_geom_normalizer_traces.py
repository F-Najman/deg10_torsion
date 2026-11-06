from pathlib import Path
import json,time
from sage.all import GammaH,ModularSymbols,ZZ,xgcd
P=Path('work/deg10/C236_resolution_20261030/geometry');t=time.monotonic();Gamma=GammaH(72,[35,-1]);MA=ModularSymbols(Gamma,2,sign=0);S=MA.cuspidal_subspace();SB=S.free_module();assert S.dimension()==82;Ws={d:S.atkin_lehner_operator(d).matrix() for d in [1,8,9,72]};rows=[];autos=[]
for b in [0,1]:
 for c in [0,1]:
  for d in [1,5,7,11,13,17]:
   lower=18*c if c else 36;gg,e,z=xgcd(ZZ(d),ZZ(lower));assert gg==1;upper=-z;tweak=(b-upper)%2;upper+=tweak*d;e+=tweak*lower;assert d*e-upper*lower==1 and upper%2==b;W=[2*d,upper,4*lower,2*e];A=MA._action_on_modular_symbols(W).restrict(SB);assert A**12==1
   if b==c==0:assert A==S.diamond_bracket_matrix(d) or A==S.diamond_bracket_matrix(d).inverse()
   autos.append(A);params=[b,c,d]
   for div,F in Ws.items():
    B=A*F;inv=B*B==1;tr=int(B.trace());row={'level_parameters':params,'Atkin_Lehner_divisor':div,'homology_trace':tr,'is_involution':inv,'fixed_points_if_involution':2-tr if inv else None};rows.append(row)
    if inv:print(row,flush=True)
assert len({tuple(A.list()) for A in autos})==24
(P/'normalizer_traces.json').write_text(json.dumps({'status':'PASS','genus':41,'automorphism_count':24,'rows':rows},indent=2)+'\n');print('C236_FULL_NORMALIZER_TRACE_PROBE_COMPLETE',time.monotonic()-t,flush=True)
