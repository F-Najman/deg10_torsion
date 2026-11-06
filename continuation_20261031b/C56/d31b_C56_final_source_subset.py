from sage.all import *
from pathlib import Path
import json,time
P=Path('work/deg10/continuation_20261031b/C56');T=time.monotonic();a=json.loads((P/'inputs/fourmap_filter.m').read_text().replace('[*','[').replace('*]',']'));ids={65,271,327,491,578,492,737,818,837,871,925,978,1009,1023,1078};cs=[r for r in a[3]if r[0]in ids];assert len(cs)==len(ids);active=sorted({i for r in cs for i in r[1]});records=[a[1][i-1][:3]for i in active];Fp=GF(3);PR=PolynomialRing(Fp,'t');fs={d:GF(3**d,'a',modulus=PR(p))for d,p in a[0]};us=[u for u in range(1,29)if gcd(u,56)==1]
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;s=b4/b3;n1=2*s+a1;n2=b2-s*s-a1*s;return[-n2**3/b3**2,1-n1*n2/b3]
def key(v,d):return min(tuple(int(z.polynomial()[i])for z in[v[0]**(3**j),v[1]**(3**j)]for i in range(d))for j in range(d))
acts={}
for i,d,ky in records:
 K=fs[d];b=PR(ky[:d])(K.gen());c=PR(ky[d:])(K.gen());E=EllipticCurve(K,[1-c,-b,-b,0,0]);Q=E(0,0);assert Q.order()==56;acts[i]=[tuple([d]+list(key(bc(E,u*Q),d)))for u in us]
clusters={};out=[];mapping=[]
for r in cs:
 can=min(tuple(sorted(acts[i][j]for i in r[1]))for j in range(len(us)))
 if can not in clusters:clusters[can]=r[0];out.append([r[0],r[1],r[2]])
 mapping.append([r[0],clusters[can]])
def mag(v):
 if isinstance(v,(int,Integer)):return str(v)
 f=all(isinstance(z,(int,Integer))for z in v);return('['if f else'[*')+','.join(mag(z)for z in v)+(']'if f else'*]')
(P/'final_original_source_subset.m').write_text(mag([a[0],records,out,mapping])+'\n');(P/'final_original_source_subset.json').write_text(json.dumps({'status':'PASS','original_cases':len(cs),'diamond_orbits':len(out),'active_points':len(records),'case_to_orbit':mapping,'elapsed_seconds':time.monotonic()-T},indent=2)+'\n');print('D31B_C56_FINAL_ORIGINAL_SOURCE_COMPLETE',len(cs),len(out),len(records),mapping,flush=True)
