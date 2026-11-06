from math import gcd
from collections import Counter
import json
from pathlib import Path
H=Path(__file__).resolve().parent
out=[]
for m,n,q in [(3,21,4),(3,24,7),(4,24,5),(4,16,5),(3,18,7),(6,12,7)]:
 unseen={(a,c) for a in range(n) for c in range(n) if gcd(gcd(a,c),n)==1};orbits=[];which={}
 while unseen:
  a,c=min(unseen);o={((sgn*(a+m*k*c))%n,(sgn*c)%n) for k in range(n//m) for sgn in [1,-1]};i=len(orbits);orbits.append(o)
  for p in o:which[p]=i
  unseen-=o
 def gal(i,r):
  a,c=min(orbits[i]);return which[(r*a%n,c)]
 group=[r for r in range(n) if gcd(r,n)==1 and (r-1)%m==0]
 # Check that the Galois action is independent of the cusp-label representative.
 for r in group:
  for i,o in enumerate(orbits):assert {which[(r*a%n,c)] for a,c in o}=={gal(i,r)}
 remain=set(range(len(orbits)));local=[]
 while remain:
  i=min(remain);o={i};j=gal(i,q)
  while j!=i:o.add(j);j=gal(j,q)
  local.append(o);remain-=o
 rat=[next(iter(o)) for o in local if len(o)==1]
 assert all(all(gal(i,r)==i for r in group) for i in rat)
 seen=set();global_orbits=[]
 for i in range(len(orbits)):
  if i not in seen:o={gal(i,r) for r in group};global_orbits.append(o);seen|=o
 out.append({'m':m,'n':n,'q':q,'geometric_cusps':len(orbits),'local_signature':dict(Counter(map(len,local))),'global_signature':dict(Counter(map(len,global_orbits))),'all_residue_rational_cusps_globally_rational':True,'Galois_multipliers':group,'global_orbits':[sorted(o) for o in global_orbits],'local_orbits':[sorted(o) for o in local]})
(H/'cusp_tate_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('DEGREE10_TATE_CUSP_REPLAY_PASS cases=6')
