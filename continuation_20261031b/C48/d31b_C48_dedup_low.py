from sage.all import *
from pathlib import Path
from collections import defaultdict,Counter
import json,time
P=Path('work/deg10/continuation_20261031b/C48');T=time.monotonic();Fp=GF(5);PR=PolynomialRing(Fp,'t')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
a=rd('low_patterns_remaining.m');fs={d:GF(5**d,'a',modulus=PR(p))for d,p in a[0]if d>=3};units=[u for u in range(1,25)if gcd(u,48)==1]
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3;return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def key(v,d):return min(tuple(int(z.polynomial()[i])for z in[v[0]**(5**j),v[1]**(5**j)]for i in range(d))for j in range(d))
active={i for r in a[2]for i in r[2]};actions={}
for i in sorted(active):
 rec=a[1][i-1];d,k=rec[1:3];K=fs[d];b=PR(k[:d])(K.gen());c=PR(k[d:])(K.gen());E=EllipticCurve(K,[1-c,-b,-b,0,0]);Q=E(0,0);assert 48*Q==E(0)and 24*Q!=E(0)and 16*Q!=E(0)
 actions[i]=[tuple([d]+list(key(bc(E,u*Q),d)))for u in units]
clusters={};out=[];maprows=[]
for row in a[2]:
 can=min(tuple(sorted(actions[i][j]for i in row[2]))for j in range(8))
 if can not in clusters:clusters[can]=len(out)+1;out.append([row[0],row[1],row[2],row[3]])
 maprows.append([row[0],clusters[can]])
summary={'status':'PASS','original_cases':len(a[2]),'source_diamond_orbits':len(out),'active_source_points':len(active),'orbit_types':[[list(k),v]for k,v in Counter(tuple(r[1])for r in out).items()],'case_to_orbit':maprows,'all_source_exact_orders_verified':True,'elapsed_seconds':time.monotonic()-T};(P/'low_source_orbits.json').write_text(json.dumps(summary,indent=2)+'\n')
def mag(v):
 if isinstance(v,(int,Integer)):return str(v)
 f=all(isinstance(z,(int,Integer))for z in v);return('['if f else'[*')+','.join(mag(z)for z in v)+(']'if f else'*]')
(P/'low_hecke_orbits.m').write_text(mag([a[0],a[1],out])+'\n');print('D31B_C48_LOW_ORBITS_COMPLETE',len(a[2]),len(out),summary['orbit_types'],time.monotonic()-T,flush=True)
