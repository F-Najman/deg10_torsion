from sage.all import *
from pathlib import Path
from collections import Counter,defaultdict
import json,time
ROOT=Path('work/deg10/continuation_20261031b/C56');T=time.monotonic();Fp=GF(3);R=PolynomialRing(Fp,'t')
def rd(n):return json.loads((ROOT/'inputs'/n).read_text().replace('[*','[').replace('*]',']'))
a=rd('fourmap_filter.m');fs={d:GF(3**d,'a',modulus=R(p))for d,p in a[0]if d>=5}
def bc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3;return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def key(v,d):return min(tuple(int(z.polynomial()[i])for z in[v[0]**(3**j),v[1]**(3**j)]for i in range(d))for j in range(d))
h={};lookup={}
for rec in a[1]:
 si,d,k,cls,obs=rec;K=fs[d];b=R(k[:d])(K.gen());c=R(k[d:])(K.gen());E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0);h[si]=(d,key(bc(E,29*P),d));lookup[si]=(d,tuple(k))
rows=[];inv=[];non=[]
for row in a[3]:
 ci,ids,cd,_=row;v=Counter(lookup[i]for i in ids);hv=Counter(h[i]for i in ids);yes=v==hv;rows.append([ci,ids,cd,yes]);(inv if yes else non).append(row)
out={'status':'PASS','involution':29,'tested_patterns':len(rows),'invariant_patterns':len(inv),'noninvariant_patterns':len(non),'invariant_types':[[list(k[0]),k[1],v]for k,v in Counter((tuple(a[1][i-1][1]for i in r[1]),r[2])for r in inv).items()],'noninvariant_types':[[list(k[0]),k[1],v]for k,v in Counter((tuple(a[1][i-1][1]for i in r[1]),r[2])for r in non).items()],'rows':rows,'noninvariant_cases':non,'elapsed_seconds':time.monotonic()-T}
(ROOT/'involution_action.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C56_INVOLUTION_COMPLETE',out['invariant_types'],out['noninvariant_types'],time.monotonic()-T,flush=True)
