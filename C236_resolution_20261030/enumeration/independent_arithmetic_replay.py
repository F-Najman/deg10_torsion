#!/usr/bin/env python3
"""Independent exact integer replay of the first C236 quotient sieve.

Magma supplies proven Picard presentations and exact place-class witnesses.
This program checks the resulting lattice quotient and repeats support arithmetic
without calling Magma.  It is not an independent proof of the Picard presentation.
"""
import ast, hashlib, itertools, json, math, pathlib, re
from collections import Counter
ROOT=pathlib.Path(__file__).resolve().parent
USED={}
def read(name):
    p=ROOT/name; txt=p.read_text(); USED[name]=hashlib.sha256(txt.encode()).hexdigest()
    txt=re.sub(r'\btrue\b','True',txt);txt=re.sub(r'\bfalse\b','False',txt)
    return ast.literal_eval(txt.replace('[*','[').replace('*]',']'))
def det3(a,b,c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
fp,cusps,nc=read('frobenius.m');assert len(fp)==83 and fp[-1]==1
tr=[];closed=[]
for d in range(1,11):
    a=-d*fp[82-d]-sum(fp[82-j]*tr[d-j-1] for j in range(1,d));tr.append(a)
    b=5**d+1-a-sum(e*closed[e-1] for e in range(1,d) if d%e==0)
    assert b%d==0;closed.append(b//d)
assert [closed[d-1]-Counter(cusps)[d] for d in range(1,11)]==nc
assert nc==[0,0,48,120,720,2904,10416,50334,218080,979512]
# Recheck the exact accepted-source copies and the three-way degree10 merge.
reuse=json.loads((ROOT/'reuse_manifest.json').read_text())
for item in reuse:
    source=pathlib.Path(item['source']);copy=pathlib.Path(item['copy'])
    assert hashlib.sha256(source.read_bytes()).hexdigest()==item['sha256']
    assert hashlib.sha256(copy.read_bytes()).hexdigest()==item['sha256']
merged=read('source10.m');records=[];counts=[0]*4;seen=set()
for shard in range(3):
    a=read(f'source10_s{shard}.m');assert a[0]==10 and a[2] and a[1]==merged[1] and a[3]==976887
    counts=[x+y for x,y in zip(counts,a[7])]
    for row in a[8]:
        assert tuple(row[1]) not in seen;seen.add(tuple(row[1]));records.append(row)
assert records==merged[8] and counts==merged[7]==[0,0,720,979512]
ad=read('target18_aut_actions.m'); target=read('target18_proven_picard.m')
assert ad[0]==[2,6] and ad[1]==[[2,0],[0,6]]
assert target[5]==[[6,0,0,0],[0,84,0,0],[0,0,252,0]]
def raw(v):
    assert len(v)==4 and v[3]==0
    return tuple(sum(v[i]*ad[2][i][j] for i in range(3))%[2,6][j] for j in range(2))
assert all(raw(r)==(0,0) for r in target[5]+target[8]+target[10])
assert (raw([0,1,-5,0]),raw([0,0,1,0]))==((1,0),(0,1))
lat=[r[:3] for r in target[5]+target[8]]
index=math.gcd(*(abs(det3(*rows)) for rows in itertools.combinations(lat,3)))
assert index==12  # Hence the displayed surjective Q map has exactly kernel H.
# Check serialized action matrices and their exact generating-place tests.
def action(v,mm):
    return tuple(sum(v[i]*mm[i][j] for i in range(2))%[2,6][j] for j in range(2))
for test in ad[4]:
    assert all(action(test[1],ad[3][i])==tuple(test[i+2]) for i in range(3))
assert len(ad[4])==48
own={};small={};fields={}
for d in range(3,8):
    a=read(f'small_qclasses{d}.m');fields[d]=a[0];own[d]=a[1];small[d]=a[2]
    base=read(f'filter18_D{d}_ownsmall.m')
    assert a[0]==base[0]
    assert len(base[1])==len(own[d])
    for r,s in zip(own[d],base[1]):
        assert r[:3]==s[:3] and r[3]==[list(raw(v)) for v in s[3]]
    assert len({tuple(r[2]) for r in small[d]})==len(small[d])
    if d<=6:
        independent=read(f'filter18_D{d}.m')
        old={tuple(r[2]):[list(raw(v)) for v in r[3]] for r in independent[2]}
        assert all(old[tuple(r[2])]==r[3] for r in small[d])
for d in [8,9]:
    a=read(f'filter18_D{d}_D{d}.m');fields[d]=a[0]
    own[d]=[[r[0],r[1],r[2],[list(raw(v)) for v in r[3]]] for r in a[1]]
dict10={}
for shard in range(8):
    a=read(f'filter18_D10_D10s{shard}.m')
    if 10 in fields:assert fields[10]==a[0]
    fields[10]=a[0]
    for r in a[1]:
        assert r[1] not in dict10 and (r[1]-1)%8==shard
        dict10[r[1]]=[r[0],r[1],r[2],[list(raw(v)) for v in r[3]]]
assert sorted(dict10)==list(range(1,40875))
own[10]=[dict10[i] for i in sorted(dict10)]
assert [len(own[d]) for d in range(3,11)]==[2,5,30,126,434,2115,9088,40874]
for d in range(3,11):
    a=read(f'source{d}.m');assert a[0]==d and a[2] and a[1]==fields[d]
    divisors=[e for e in range(1,d+1) if d%e==0]
    assert a[7]==[nc[e-1] for e in divisors]
    assert len({tuple(r[1]) for r in a[8]})==len(a[8])
    keys=[r[1] for r in a[8] if r[0]==d]
    assert [r[2] for r in own[d]]==keys
    assert all(len(r[2])==2*d for r in own[d])
    assert [r[1] for r in own[d]]==list(range(1,len(own[d])+1))
sid=0;survivors=[];stats=[]
def pattern(pat,rows):
    global sid
    trials=accepts=0
    for terms in rows:
        sid+=1;trials+=1
        if any(len(r[3])!=2 for r in terms):continue
        if any(sum(r[3][j][k] for r in terms)%[2,6][k] for j in range(2) for k in range(2)):continue
        accepts+=1;survivors.append([sid,[[r[0],r[2]] for r in terms]])
    stats.append([pat,trials,accepts])
for d in range(3,11):pattern([d],((r,) for r in own[d]))
for a,b in [[3,3],[3,4],[3,5],[3,6],[3,7],[4,4],[4,5],[4,6],[5,5]]:
    pattern([a,b],itertools.product(own[a],small[b]))
pattern([3,3,3],((r,s,t) for r in own[3] for s,t in itertools.combinations_with_replacement(small[3],2)))
pattern([3,3,4],itertools.product(own[3],small[3],small[4]))
assert sid==135282 and len(survivors)==1251
record=read('first_survivors_unified.m');detail=read('first_survivors_details.m')
assert record[0]==[[d,fields[d]] for d in range(3,11)]
assert survivors==record[1] and stats==detail[0]
assert [r[0] for r in detail[1]]==[r[0] for r in survivors]
result={'status':'PASS','scope':'independent integer quotient and support arithmetic; Picard geometry supplied by Magma exact identities',
 'frobenius_non_cuspidal_counts':nc,'source_copies_verified':len(reuse),'degree10_records':len(records),'quotient_invariants':[2,6],'lattice_index':index,'automorphism_place_tests':len(ad[4]),
 'source_own_counts':{d:len(own[d]) for d in range(3,11)},'small_counts':{d:len(small[d]) for d in range(3,8)},
 'support_cases':sid,'survivors':len(survivors),'pattern_statistics':stats,'sha256':USED}
(ROOT/'independent_arithmetic_replay.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
