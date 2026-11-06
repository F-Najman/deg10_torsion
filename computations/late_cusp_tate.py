from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path('work/deg10/continuation_20261029b');out=[]
for m,n in [(1,38),(1,40),(1,50),(2,38),(2,40),(2,50),(2,26)]:
 unseen={(a,c) for a in range(n) for c in range(n) if gcd(gcd(a,c),n)==1};orbits=[];which={}
 while unseen:
  a,c=min(unseen);orb={((s*(a+m*k*c))%n,(s*c)%n) for k in range(n//m) for s in [1,-1]};idx=len(orbits);orbits.append(orb)
  for z in orb:which[z]=idx
  unseen-=orb
 def gal(i,r):a,c=min(orbits[i]);return which[r*a%n,c]
 group=[r for r in range(n) if gcd(r,n)==1 and (r-1)%m==0]
 for r in group:
  for i,o in enumerate(orbits):assert {which[r*a%n,c] for a,c in o}=={gal(i,r)}
 remain=set(range(len(orbits)));local=[]
 while remain:
  i=min(remain);orb={i};j=gal(i,3)
  while j!=i:orb.add(j);j=gal(j,3)
  local.append(orb);remain-=orb
 global_orbits=[];seen=set()
 for i in range(len(orbits)):
  if i not in seen:o={gal(i,r) for r in group};global_orbits.append(o);seen|=o
 assert all(all(gal(next(iter(o)),r)==next(iter(o)) for r in group) for o in local if len(o)==1)
 if (m,n)==(2,38):assert Counter(map(len,local))=={1:27,9:3}
 if (m,n)==(2,40):assert Counter(map(len,local))=={1:16,2:12,4:10}
 if (m,n)==(2,50):assert Counter(map(len,local))=={1:30,4:6,10:3}
 out.append({'m':m,'n':n,'prime':3,'geometric_cusps':len(orbits),'local_signature':dict(Counter(map(len,local))),'global_signature':dict(Counter(map(len,global_orbits))),'every_residue_rational_cusp_globally_rational':True,'global_orbits':[sorted(o) for o in global_orbits],'local_orbits':[sorted(o) for o in local]})
(P/'cusp_tate.json').write_text(json.dumps(out,indent=2)+'\n')
print('ALL_TATE_CUSP_ACTIONS_VERIFIED',[(a['m'],a['n'],a['local_signature'],a['global_signature']) for a in out],flush=True)
