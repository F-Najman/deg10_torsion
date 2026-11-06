from pathlib import Path
from math import gcd
from collections import Counter
import json
N=56;p=3;q=11
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
oi={j:i for i,o in enumerate(orbits)for j in o};cols=[]
for i,o in enumerate(orbits):
 if len(o)>2:continue
 v=Counter()
 for j in o:v[sigma(q,j)]+=1;v[j]-=1;v[delta(q,sigma(pow(q,-1,N),j))]-=q;v[delta(q,j)]-=0
 # T_q=sigma_q+q delta_q sigma_q^-1; subtract1+qdelta_q.
 v=Counter()
 for j in o:v[sigma(q,j)]+=1;v[delta(q,sigma(pow(q,-1,N),j))]+=q;v[j]-=1;v[delta(q,j)]-=q
 nz={j:n for j,n in v.items()if n};assert not nz,(i,len(o),nz)
 cols.append([i,len(o)])
out={'status':'PASS','level':N,'prime':p,'hecke':q,'geometric_cusps':len(labels),'closed_cusp_degrees':dict(Counter(map(len,orbits))),'tested_cusps_degree_at_most2':len(cols),'all_A11_columns_zero_literally':True,'primitive_labels':labels,'Frob3_orbits':orbits,'tested_closed_columns':cols};P=Path(__file__).parent;(P/'cusp_A11_degree2_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C56_CUSP_A11_COMPLETE',len(labels),len(cols),out['closed_cusp_degrees'])
