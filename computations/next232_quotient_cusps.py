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
rows=[]
for b in range(2):
 for c in range(2):
  for d in range(1,16,2):
   if (d*d+16*b*c)%32!=1 or (b,c,d)==(0,0,1):continue
   gg=[]
   for C in cs:
    ans={ix[((16*c*s+d*a)%n,(16*c*t+d*r)%n,(s+b*a)%2,(t+b*r)%2)] for a,r,s,t in C};assert len(ans)==1;gg.append(ans.pop())
   assert all(gg[gg[i]]==i and gg[fro[i]]==fro[gg[i]] for i in range(56));fixed=sum(i==gg[i] for i in range(56));genus=(76-fixed)//4;assert 4*genus==76-fixed
   oc=cycles(gg);oi={i:j for j,C in enumerate(oc) for i in C};qfro=[oi[fro[C[0]]] for C in oc];rat=sum(i==qfro[i] for i in range(len(oc)));gon=(rat+3)//4
   rows.append({'parameters':[b,c,d],'fixed_geometric_cusps':fixed,'quotient_genus':genus,'quotient_cusp_signature':dict(Counter(map(len,cycles(qfro)))),'rational_cusps_over_F3':rat,'quotient_gonality_lower_bound_by_point_count':gon,'CS_excludes_independent_degree10_map':37>2*genus+9,'cusp_permutation':gg})
(P/'involution_quotients.json').write_text(json.dumps({'status':'PASS','geometric_cusps':56,'Frobenius_signature':dict(Counter(map(len,cycles(fro)))),'rows':rows},indent=2)+'\n')
for r in rows:print({k:v for k,v in r.items() if k!='cusp_permutation'},flush=True)
print('C232_INVOLUTION_CUSP_QUOTIENTS_COMPLETE',flush=True)
