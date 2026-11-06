from sage.all import *
from pathlib import Path
import ast,re,json
H=Path('work/deg10/continuation_20261031b/M3N18')
def rd(n):
 s=(H/n).read_text();s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
a=rd('cusp_automorphisms.m');tau=rd('tau_cusp_permutation_p7.m')[1];G=PermutationGroup(a[2]+[tau]);re={'order':int(G.order()),'involutions':sum(g.order()==2 for g in G),'old27':int(PermutationGroup(a[2]).order()),'cusp_permutation_generators':a[2]+[tau]};print(json.dumps(re));(H/'full_normalizer_order.json').write_text(json.dumps(re,indent=2)+'\n');print('D31B_M318_AUT_ORDER_COMPLETE')
