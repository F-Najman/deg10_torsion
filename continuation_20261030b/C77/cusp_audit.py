"""Exact primitive-column cusp orbits and degree-ten norm-image coverage."""
from pathlib import Path
from math import gcd,prod
from collections import Counter
import json,hashlib
HERE=Path(__file__).resolve().parent
n=77
primitive={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
def orbit(x):
 a,b=x
 return frozenset(((s*(a+k*b))%n,s*b%n) for s in (-1,1) for k in range(n))
classes=sorted(set(map(orbit,primitive)),key=lambda x:sorted(x))
assert len(classes)==120
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
sf=cycles(gal[2]);sg=global_orbits(gal,len(classes));assert set(sf)!=set(sg)
assert Counter(map(len,sf))=={1:30,3:10,10:3,30:1};assert Counter(map(len,sg))=={1:30,6:5,10:3,30:1}
h=cycles(diamond(36));assert len(h)==24 and all(len(o)==5 for o in h)
qi={i:j for j,o in enumerate(h) for i in o}
qgal={u:[qi[gal[u][next(iter(o))]] for o in h] for u in units}
qf=cycles(qgal[2]);qg=global_orbits(qgal,len(h));assert set(qf)!=set(qg)
assert Counter(map(len,qf))=={1:6,2:3,3:2,6:1};assert Counter(map(len,qg))=={1:6,2:3,6:2}
closed={j:i for i,o in enumerate(qf) for j in o};ds=[len(o) for o in qf]
rows=[];generators=set()
for o in sf:
 targets={closed[qi[i]] for i in o};assert len(targets)==1
 at=next(iter(targets));e=ds[at];d=len(o);assert d%e==0
 v=tuple(d//e if i==at else 0 for i in range(len(qf)))
 if d<=10:generators.add(v)
 width=n//gcd(reps[next(iter(o))][1],n)
 rows.append(dict(source_degree=d,width=width,quotient_closed_index=at,quotient_degree=e,norm_weight=d//e))
bydegree=[set()for d in range(11)];bydegree[0].add((0,)*len(ds))
for d in range(11):
 for v in list(bydegree[d]):
  for w in generators:
   total=d+sum(a*b for a,b in zip(w,ds))
   if total<=10:bydegree[total].add(tuple(a+b for a,b in zip(v,w)))
assert len(bydegree[10])==4992
normal=[];extra=[]
for v in sorted(bydegree[10]):
 if any(v[i] for i,e in enumerate(ds)if e==2):
  assert sum(v)==5 and sum(v[i]>0 for i in range(len(v)))==1;extra.append(v)
 else:normal.append(v)
assert len(normal)==4989 and len(extra)==3
record=dict(status='PASS',level=77,source_geometric_cusps=120,
 source_global_signature=dict(Counter(map(len,sg))),source_F2_signature=dict(Counter(map(len,sf))),
 every_source_cusp_global_orbit_intact_at_2=False,
 H_generator=36,H_order=5,H_free_on_all_geometric_cusps=True,
 quotient_global_signature=dict(Counter(map(len,qg))),quotient_F2_signature=dict(Counter(map(len,qf))),
 every_quotient_cusp_global_orbit_intact_at_2=False,
 source_to_quotient_closed_norms=rows,quotient_closed_degrees=ds,
 degree10_pushforward_patterns=4992,normal_degree10_patterns=4989,
 additional_images_5_times_quadratic_cusp=extra,
 complete_degree10_pushforward_vectors=sorted(bydegree[10]),
 source_global_orbits=[sorted(o)for o in sg],source_Frob2_permutation=gal[2],
 quotient_global_orbits=[sorted(o)for o in qg],quotient_Frob2_permutation=qgal[2])
(HERE/'cusp_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print('C77_CUSP_AUDIT_PASS source120 Y24 patterns4992=4989+3 global_cubic_fragments_retained=true')
