from pathlib import Path
from math import gcd
from collections import Counter
import json
N=48;p=5;q=13
pairs=[(a,b)for a in range(N)for b in range(N)if gcd(gcd(a,b),N)==1]
def can(v):
 a,b=v;return min(((s*(a+k*b))%N,(s*b)%N)for s in[1,-1]for k in range(N))
labels=sorted({can(v)for v in pairs});where={v:i for i,v in enumerate(labels)}
def sigma(u,i):a,b=labels[i];return where[can((u*a,b))]
def delta(u,i):a,b=labels[i];return where[can((u*a,pow(u,-1,N)*b))]
frob=[sigma(p,i)for i in range(len(labels))];seen=set();orbits=[]
for i in range(len(labels)):
 if i in seen:continue
 orb=[];j=i
 while j not in orb:orb.append(j);seen.add(j);j=frob[j]
 assert j==i;orbits.append(orb)
oi={j:i for i,o in enumerate(orbits)for j in o};cols=[];allcols=[]
for i,o in enumerate(orbits):
 v=Counter()
 for j in o:v[sigma(q,j)]+=1;v[delta(q,sigma(pow(q,-1,N),j))]+=q;v[j]-=1;v[delta(q,j)]-=q
 nz={j:n for j,n in v.items()if n};cv={oi[j]:nz[j]for j in nz}
 for j in nz:assert all(nz.get(k,0)==nz[j]for k in orbits[oi[j]])
 if cv:
  assert len(cv)==2 and sorted(cv.values())==[-12,12]
  widths={N//gcd(labels[j][1],N)for j in o};assert len(widths)==1
  cols.append([i,len(o),widths.pop(),sorted(cv.items())])
 allcols.append(sorted(cv.items()))
assert len(labels)==56 and Counter(len(o)for o in orbits)==Counter({1:16,2:12,4:4})
assert len(cols)==4 and sorted((r[1],r[2])for r in cols)==[(2,2),(2,2),(4,1),(4,1)]
out={'status':'PASS','level':N,'prime':p,'hecke':q,'geometric_cusps':len(labels),'closed_cusp_degrees':dict(Counter(map(len,orbits))),'primitive_labels':labels,'Frob5_orbits':orbits,'nonzero_A13_columns':cols,'all_closed_columns':allcols,'conclusion':'All A13 cusp images lie in the span of12 times the differences of the width2 quadratic pair and width1 quartic pair.'};P=Path(__file__).parent;(P/'independent_cusp_A13_labels.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C48_INDEPENDENT_CUSP_A13_LABELS_COMPLETE',len(labels),len(orbits),len(cols))
