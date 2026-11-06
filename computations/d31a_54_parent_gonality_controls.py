"""Independent coset-permutation genus and Castelnuovo--Severi controls."""
from math import gcd
from fractions import Fraction
from pathlib import Path
import json
N=54
def genus(H):
 def rep(x):
  c,d=x
  return min(((a*c)%N,(a*d)%N)for a in H)
 reps=sorted({rep((c,d))for c in range(N)for d in range(N)if gcd(gcd(c,d),N)==1})
 S={x:rep((x[1],-x[0]))for x in reps}
 ST={x:rep((x[1],x[1]-x[0]))for x in reps}
 T={x:rep((x[0],x[0]+x[1]))for x in reps}
 assert set(S.values())==set(ST.values())==set(T.values())==set(reps)
 assert all(S[S[x]]==x and ST[ST[ST[x]]]==x for x in reps)
 unused=set(reps);widths=[]
 while unused:
  x=min(unused);z=x;w=0
  while z in unused:
   unused.remove(z);w+=1;z=T[z]
  assert z==x
  widths.append(w)
 e2=sum(S[x]==x for x in reps);e3=sum(ST[x]==x for x in reps)
 g=1+Fraction(len(reps),12)-Fraction(e2,4)-Fraction(e3,3)-Fraction(len(widths),2)
 assert g.denominator==1
 return {'index':len(reps),'cusps':len(widths),'cusp_widths':sorted(widths),'elliptic2':e2,'elliptic3':e3,'genus':int(g)}
H={pow(17,e,54)for e in range(6)};assert H=={1,17,19,35,37,53}
x=genus({1,53});y=genus(H);assert(x['genus'],y['genus'],x['index'],y['index'])==(52,10,972,324)
lo=Fraction(325*972,32768);assert 9<lo<10
cs=3*y['genus']+(3-1)*(10-1);assert cs==48<x['genus']and 10%3!=0
out={'status':'PASS','source':x,'quotient':y,'diamond_subgroup':sorted(H),'map_degree':3,'retained_analytic_bound':[lo.numerator,lo.denominator],'CS_upper_bound_for_degree10_map':cs,'gonality_exceeds10':True,'method':'Exact modular coset permutations, retained analytic bound, and Castelnuovo-Severi; no cusp jets required'}
Path('work/deg10/continuation_20261031a/C54/parent_gonality_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print('C54_PARENT_GONALITY_COMPLETE',x['genus'],y['genus'],cs)
