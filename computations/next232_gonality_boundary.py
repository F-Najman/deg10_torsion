from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path('work/deg10/C232_resolution_20261030');n=32
labels={(a,r,s,t) for a in range(n) for r in range(n) if gcd(gcd(a,r),n)==1 for s in range(2) for t in range(2) if (s*r-t*a)%2==1}
def orb(x):
 a,r,s,t=x;return frozenset(((e*(a+k*r))%n,e*r%n,e*(s+k*t)%2,e*t%2) for e in [-1,1] for k in range(n))
cs=sorted(set(map(orb,labels)),key=lambda x:sorted(x));ix={x:i for i,c in enumerate(cs) for x in c};assert len(cs)==56
fro=[]
for c in cs:
 im={ix[(3*a%n,r,s,t)] for a,r,s,t in c};assert len(im)==1;fro.append(im.pop())
def cycles(perm):
 unseen=set(range(len(perm)));out=[]
 while unseen:
  x=min(unseen);c=[]
  while x not in c:c.append(x);unseen.remove(x);x=perm[x]
  out.append(c)
 return out
assert Counter(map(len,cycles(fro)))=={1:16,2:8,4:2,8:2}

from sage.all import GammaH
X=GammaH(64,[33,-1]);Y=GammaH(64,[15,-1]);assert X.genus()==37 and Y.genus()==13 and X.index()==2*Y.index()
assert sorted(Y._list_of_elements_in_H())==[1,15,17,31,33,47,49,63]
sigma=[]
for C in cs:
 ans={ix[(11*a%n,r,s,t)] for a,r,s,t in C};assert len(ans)==1;sigma.append(ans.pop())
f7=[]
for i in range(56):
 j=i
 for _ in range(7):j=fro[j]
 f7.append(j)
assert f7==sigma and pow(3,7,32)==11
bounds=[2*13+(d-1) for d in range(1,12)];assert max(bounds)==36<37
out={'status':'PASS','source_congruence_subgroup':{'level':64,'Delta':[1,31,33,63],'index':int(X.index()),'genus':int(X.genus())},'quotient_congruence_subgroup':{'level':64,'Delta':[1,15,17,31,33,47,49,63],'index':int(Y.index()),'genus':int(Y.genus())},'quotient_degree':2,'literature_quotient_Q_gonality_lower_bound':6,'source_Q_gonality_lower_bound':12,'Castelnuovo_Severi_bounds_for_degrees_1_through_11':bounds,'cusps':56,'cusp_Frobenius_signature':dict(Counter(map(len,cycles(fro)))),'Hecke_boundary_prime':11,'boundary_sigma11_is_Frobenius3_power':7,'boundary_sigma11_fixes_every_F3_closed_cusp':True}
(P/'gonality_boundary_verification.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_GONALITY_AND_BOUNDARY_COMPLETE',out,flush=True)
