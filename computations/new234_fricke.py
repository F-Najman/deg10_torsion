from sage.all import GammaH,ModularSymbols
import json,time
from pathlib import Path
T=time.monotonic();G=GammaH(68,[35,-1]);M=ModularSymbols(G,2,sign=0).cuspidal_subspace();print('GENUS',G.genus(),M.dimension(),flush=True)
for q in [4,17,68]:
 try:
  W=M.atkin_lehner_operator(q).matrix();tr=W.trace();print('FRICKE',q,W*W==1,tr,'fixed',100+tr,'genus',(98+tr)/4,time.monotonic()-T,flush=True)
 except Exception as e:print(type(e).__name__,str(e),flush=True)
