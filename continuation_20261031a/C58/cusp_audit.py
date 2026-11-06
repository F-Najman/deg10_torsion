"""Exact primitive-column cusp orbits and degree-ten norm-image coverage."""
from pathlib import Path
from math import gcd,prod
from collections import Counter
import json,hashlib
HERE=Path(__file__).resolve().parent
n=58
primitive={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x
 return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,primitive)),key=lambda x:sorted(x))
assert len(classes)>0
ix={x:i for i,c in enumerate(classes) for x in c};reps=[min(c) for c in classes]
units=[u for u in range(n) if gcd(u,n)==1]
gal={u:[ix[(u*a%n,b)] for a,b in reps] for u in units}
def diamond(u):return[ix[(pow(u,-1,n)*a%n,u*b%n)] for a,b in reps]
def cycles(p):
 left=set(range(len(p)));out=[]
 while left:
  first=min(left);i=first;o=[]
  while i not in o:o.append(i);left.remove(i);i=p[i]
  assert i==first;out.append(frozenset(o))
 return out
def global_orbits(action,count):
 left=set(range(count));out=[]
 while left:
  i=min(left);o=frozenset(action[u][i] for u in units);assert o<=left;left-=o;out.append(o)
 return out



sf=cycles(gal[3]);sg=global_orbits(gal,len(classes))
rat=[next(iter(o))for o in sf if len(o)==1]
width=lambda i:n//gcd(reps[i][1],n)
dd=diamond(3)
assert all(dd[dd[dd[i]]] in range(len(classes)) for i in range(len(classes)))
record=dict(status='PASS',level=n,geometric_cusps=len(classes),global_signature=dict(Counter(map(len,sg))),F3_signature=dict(Counter(map(len,sf))),
 global_orbits=[sorted(o)for o in sg],F3_orbits=[sorted(o)for o in sf],F3_permutation=gal[3],diamond3_permutation=dd,
 rational_cusp_widths=dict(Counter(width(i)for i in rat)),rational_cusps=rat,
 rational_diamond_cycles=[sorted(o)for o in cycles(dd) if o<=set(rat)])
(HERE/'cusp_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print('CUSP_AUDIT',n,record['global_signature'],record['F3_signature'],record['rational_cusp_widths'])
