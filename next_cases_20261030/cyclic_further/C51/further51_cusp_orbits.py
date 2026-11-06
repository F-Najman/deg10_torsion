#!/usr/bin/env python3
from pathlib import Path
from math import gcd
from collections import Counter
import json
H=Path(__file__).resolve().parent;N=51

def cusp(a,c):
 c%=N;g=gcd(c,N);return min((a%g,c),((-a)%g,(-c)%N))
def sigma(P,u):return cusp(u*P[0],P[1])
def diamond(P,u):return cusp(pow(u,-1,N)*P[0],u*P[1])
def orbits(values,action):
 left=set(values);answer=[]
 while left:
  first=min(left);x=first;cycle=[]
  while x not in cycle:cycle.append(x);x=action(x)
  assert x==first;part=frozenset(cycle);assert part<=left;left-=part;answer.append(part)
 return answer
labels={cusp(a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1}
frobenius=orbits(labels,lambda p:sigma(p,2));assert len(labels)==64 and Counter(map(len,frobenius))=={1:16,2:8,8:4}
units=[u for u in range(N)if gcd(u,N)==1]
global_rational={p for p in labels if all(sigma(p,u)==p for u in units)}
local_rational={next(iter(o))for o in frobenius if len(o)==1}
assert global_rational==local_rational and len(global_rational)==16
result={d:sorted(map(len,orbits({o for o in frobenius if len(o)==d},lambda o:frozenset(diamond(p,5)for p in o))))for d in [1,2,8]}
assert result[1]==[16]and result[2]==[8]
out={'N':51,'geometric_cusps':64,'global_rational_cusps':16,'Frobenius2_signature':{'1':16,'2':8,'8':4},'diamond_generator':5,'diamond_orbit_lengths_on_closed_cusps':result,'rational_cusps_transitive':True,'quadratic_cusps_transitive':True,'status':'PASS'}
(H/'cusp_orbits_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('DECIC51_CUSP_ORBITS_PASS',json.dumps(out,sort_keys=True))
