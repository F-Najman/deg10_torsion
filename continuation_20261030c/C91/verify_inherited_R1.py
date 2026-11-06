"""Replay the degree-independent ordinary-support canonical rank certificate."""
from pathlib import Path
import gzip,json,hashlib,time
H=Path(__file__).resolve().parent
p=H/'inputs/composite91/canonical_R1_matrices.json.gz';data=json.loads(gzip.decompress(p.read_bytes()));assert len(data)==288
records=[];tt=time.time()
for a in data:
 assert len(a['rows'])==265 and a['columns']==590
 basis={}
 for row in a['rows']:
  v=int(row,16);assert v.bit_length()<=590
  while v:
   k=v.bit_length()-1
   if k in basis:v^=basis[k]
   else:basis[k]=v;break
 assert len(basis)==265
 records.append({k:a[k]for k in ['source','target_diamond_translate','convention']}|{'rank':len(basis)})
assert {(x['source'],x['target_diamond_translate'],x['convention'])for x in records}=={(s,t,c)for s in range(1,5)for t in range(36)for c in [1,2]}
x=dict(status='PASS',matrices=288,shape=[265,590],all_ranks=265,method='Independent integer-bit Gaussian elimination',input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),records=records,seconds=time.time()-tt)
(H/'inherited_R1_rank_replay.json').write_text(json.dumps(x,indent=2)+'\n')
print('C91_INHERITED_R1_RANK_REPLAY_PASS 288 matrices rank265')
