#!/usr/bin/env python3
"""Independent integer replay of C242 census bookkeeping and genus4 filters."""
import ast,hashlib,itertools,json,math,pathlib,re
from collections import Counter
ROOT=pathlib.Path(__file__).resolve().parent;USED={}
def rd(name):
 p=ROOT/name;txt=p.read_text();USED[name]=hashlib.sha256(txt.encode()).hexdigest();return ast.literal_eval(txt.replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def det(a):
 if len(a)==1:return a[0][0]
 return sum((-1)**j*a[0][j]*det([r[:j]+r[j+1:] for r in a[1:]]) for j in range(len(a)))
fp,cusps,nc=rd('frobenius.m');assert len(fp)==123 and fp[-1]==1
tr=[];closed=[]
for d in range(1,11):
 a=-d*fp[122-d]-sum(fp[122-j]*tr[d-j-1] for j in range(1,d));tr.append(a);b=5**d+1-a-sum(e*closed[e-1] for e in range(1,d) if d%e==0);assert b%d==0;closed.append(b//d)
assert [closed[d-1]-Counter(cusps)[d] for d in range(1,11)]==nc==[0,0,0,126,648,2080,11376,48816,211776,983196]
for row in json.loads((ROOT/'reuse_manifest.json').read_text()):
 assert hashlib.sha256(pathlib.Path(row['source']).read_bytes()).hexdigest()==row['sha256']==hashlib.sha256(pathlib.Path(row['copy']).read_bytes()).hexdigest()
merged=rd('source10.m');records=[];counts=[0]*4;seen=set()
for shard in range(4):
 a=rd(f'source10_s{shard}.m');assert a[0]==10 and a[1]==merged[1] and a[2] and a[3]==976887;counts=[x+y for x,y in zip(counts,a[7])]
 for row in a[8]:assert tuple(row[1]) not in seen;seen.add(tuple(row[1]));records.append(row)
assert records==merged[8] and counts==merged[7]==[0,0,648,983196]
target=rd('target14_proven_picard.m');assert target[5]==[[2,0,0,0,0],[0,2,0,0,0],[0,0,12,0,0],[0,0,0,36,0]]
def qkey(v):assert len(v)==5 and v[-1]==0;return (v[2]%2,v[3]%2)
assert all(qkey(v)==(0,0) for v in target[5]+target[8]+target[10]);lat=[v[:-1] for v in target[5]+target[8]];index=math.gcd(*(abs(det(list(rows))) for rows in itertools.combinations(lat,4)));assert index==4
own={};small={};fields={}
for d in range(4,10):
 tag='small' if d<=6 else f'D{d}';a=rd(f'filter14_D{d}_{tag}.m');fields[d]=a[0];own[d]=a[1]
 if d<=6:small[d]=a[2]
dict10={}
for shard in range(8):
 a=rd(f'filter14_D10_D10s{shard}.m')
 if 10 in fields:assert fields[10]==a[0]
 fields[10]=a[0]
 for row in a[1]:assert row[1] not in dict10 and (row[1]-1)%8==shard;dict10[row[1]]=row
n10=sum(row[0]==10 for row in merged[8]);assert sorted(dict10)==list(range(1,n10+1));own[10]=[dict10[i] for i in sorted(dict10)]
assert [len(own[d]) for d in range(4,10)]==[5,18,64,316,1372,5888]
assert [len(small[d]) for d in range(4,7)]==[126,648,2080]
for d in range(4,11):
 a=rd(f'source{d}.m');assert a[0]==d and a[1]==fields[d] and a[2];divisors=[e for e in range(1,d+1) if d%e==0];assert a[7]==[nc[e-1] for e in divisors]
 assert [row[2] for row in own[d]]==[ky for degree,ky in a[8] if degree==d]
 assert [row[1] for row in own[d]]==list(range(1,len(own[d])+1))
sid=0;survivors=[];stats=[]
def pattern(pat,rows):
 global sid
 trials=accepts=0
 for terms in rows:
  sid+=1;trials+=1
  if any(len(r[3])!=2 for r in terms):continue
  if any(sum(qkey(r[3][j])[k] for r in terms)%2 for j in range(2) for k in range(2)):continue
  accepts+=1;survivors.append([sid,[[r[0],r[2]] for r in terms]])
 stats.append([pat,trials,accepts])
for d in range(4,11):pattern([d],((r,) for r in own[d]))
for a,b in [(4,4),(4,5),(4,6),(5,5)]:pattern([a,b],itertools.product(own[a],small[b]))
assert sid==7663+n10+25934
out=rd('first_survivors_unified.m');details=rd('first_survivors_details.m');assert out==[[[d,fields[d]] for d in range(4,11)],survivors] and details[0]==stats and [r[0] for r in details[1]]==[r[0] for r in survivors]
result={'status':'PASS','source_non_cuspidal_counts':nc,'source10_record_count':len(records),'own_counts':{d:len(own[d]) for d in own},'target_quotient_index':index,'support_cases':sid,'survivors':len(survivors),'pattern_statistics':stats,'sha256':USED}
(ROOT/'independent_arithmetic_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
