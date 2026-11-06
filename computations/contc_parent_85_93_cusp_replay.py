"""Independent relation kernels and full-rank-at-leaf mixed-cusp checks."""
from pathlib import Path
from math import comb
import json,time,subprocess,hashlib
H=Path(__file__).resolve().parents[1];C=H/'continuation_20261030c';t0=time.time();families=[];pure={}
def kernel(rows):
 piv={};rels=[]
 for i,row in enumerate(rows):
  v=row;w=1<<i
  while v:
   p=v.bit_length()-1
   if p in piv:a,b=piv[p];v^=a;w^=b
   else:piv[p]=(v,w);break
  if not v:rels.append(w)
 for r in rels:
  v=0
  for i,row in enumerate(rows):
   if r>>i&1:v^=row
  assert v==0
 assert len(piv)+len(rels)==len(rows)
 return rels
for N in[85,93]:
 P=C/f'C{N}';d=json.loads((P/'cusp_blocks.json').read_text());rat=[v[0]for v in d['rational']['local']];n=len(rat);rows=rat[0][:10]+sum((v[:5]for v in rat[1:]),[]);rel=kernel(rows);assert len(rel)==(13 if N==85 else 5)
 best=999;best_raw=999;v=0
 for j in range(1,1<<len(rel)):
  v^=rel[(j&-j).bit_length()-1];z=v;costs=[(z&1023).bit_length()];z>>=10
  for i in range(n-1):costs.append((z&31).bit_length());z>>=5
  assert not z;best=min(best,max(1,max(costs))+sum(costs[1:]));best_raw=min(best_raw,sum(costs))
 assert best>10 and best_raw>10
 pure[N]={'relation_dimension':len(rel),'all_nonzero_relations_checked':(1<<len(rel))-1,'minimum_normalized_prefix_cost':best,'minimum_raw_prefix_cost':best_raw,'all_pure_rational_rank10':True}
 def add(label,fix,rem):families.append((N,label,fix,rem,rat))
 if N==85:
  Q=d['nonrational']['q4'];add('onequartic_plus_rat6',Q[0][:4],6);add('doublequartic_plus_rat2',Q[0],2)
  for i in range(1,8):add(f'twoquartics_{i}_plus_rat2',Q[0][:4]+Q[i][:4],2)
  for name in['q8width5','q8width1']:add(name+'_plus_rat2',d['nonrational'][name][0],2)
 else:
  Q=d['nonrational']['q5'];add('onequintic_plus_rat5',Q[0][:5],5);add('doublequintic',Q[0],0)
  for i in range(1,6):add(f'twoquintics_{i}',Q[0][:5]+Q[i][:5],0)
lines=[str(len(families))]
def encoded(v):return ' '.join(str(v>>(64*i)&((1<<64)-1))for i in range(4))
for N,label,fix,rem,rat in families:
 assert len(fix)+rem==10 and all(v.bit_length()<=256 for v in fix+sum(rat,[]))
 lines.append(f'{len(fix)} {rem} {len(rat)}');lines.extend(encoded(v)for v in fix)
 for r in rat:lines.extend(encoded(v)for v in r[:rem])
p=C/'parent_mixed_rank_families.txt';p.write_text('\n'.join(lines)+'\n');cpp=H/'computations/contc_parent_mixed_ranks.cpp';exe=C/'parent_mixed_rank_native';subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
r=subprocess.run([str(exe),str(p)],capture_output=True,text=True,check=True);print(r.stdout,flush=True);rs=[list(map(int,v.split()))for v in r.stdout.splitlines()if v[:1].isdigit()];assert len(rs)==len(families)
for N in[85,93]:
 records=[]
 for i,count in rs:
  n,label,fix,rem,rat=families[i];assert count==comb(rem+len(rat)-1,len(rat)-1)
  if n==N:records.append({'family':label,'count':count,'all_full_rank':True})
 d={'status':'PASS','level':N,'method':'Separate integer relation-kernel elimination and Gray enumeration; independent C++ full rank at every mixed leaf, without incremental basis reuse','rational':pure[N],'mixed_families':records,'total_mixed_patterns':sum(v['count']for v in records),'source_blocks_sha256':hashlib.sha256((C/f'C{N}/cusp_blocks.json').read_bytes()).hexdigest(),'cpp_sha256':hashlib.sha256(cpp.read_bytes()).hexdigest(),'seconds':time.time()-t0}
 (C/f'C{N}/parent_cusp_rank_replay.json').write_text(json.dumps(d,indent=2)+'\n')
print('PARENT_C85_C93_CUSP_REPLAY_COMPLETE',len(families),sum(count for i,count in rs),round(time.time()-t0,2),flush=True)
