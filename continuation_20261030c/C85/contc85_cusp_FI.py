LEVEL=85;CASEPATH='work/deg10/continuation_20261030c/C85'
from pathlib import Path
import json,time
from sage.all import GF,Matrix
P=Path(CASEPATH);N=LEVEL;t0=time.monotonic();a=json.loads((P/'cusp_blocks.json').read_text());rat=[v[0]for v in a['rational']['local']];amb=a['Wrank'];g=len(rat)
def insert(B,v):
 while v:
  k=v.bit_length()-1
  if k in B:v^=B[k]
  else:B[k]=v;return True
 return False
def rank(rows):
 B={}
 for v in rows:insert(B,v)
 return len(B)
R=rat[0][:10]+sum((v[:5]for v in rat[1:]),[]);M=Matrix(GF(2),[[row>>j&1 for j in range(amb)]for row in R]);ker=[sum(int(c)<<i for i,c in enumerate(row))for row in M.left_kernel().basis()];best=1000;badrel=None
for bits in range(1,2**len(ker)):
 v=0
 for i,r in enumerate(ker):
  if bits>>i&1:v^=r
 pre=[];mask=(1<<10)-1;pre.append((v&mask).bit_length());v>>=10
 for k in range(g-1):pre.append((v&31).bit_length());v>>=5
 cost=max(1,max(pre))+sum(pre[1:])
 if cost<best:best=cost;badrel=pre
print('CONTC_RATIONAL_PREFIX_COST',N,len(R),M.rank(),len(ker),best,badrel,flush=True)
records=[];failures=[]
def testfamily(label,fixed,total):
 base={}
 for v in fixed:insert(base,v)
 assert len(base)==len(fixed),(label,'fixed rank',len(base),len(fixed))
 count=[0];bad=[]
 def visit(j,left,B,choices):
  if j==g:
   if left==0:
    count[0]+=1
    if len(B)!=10:
     if len(bad)<10:bad.append([choices,len(B)])
   return
  cur=B.copy()
  for m in range(left+1):
   visit(j+1,left-m,cur,choices+[m])
   if m<left:insert(cur,rat[j][m])
 visit(0,total,base,[])
 row=dict(label=label,configurations=count[0],failures=bad);records.append(row);failures.extend(bad);print('CONTC_MIXED_CUSP_FAMILY',N,label,count[0],bad[:1],time.monotonic()-t0,flush=True)
if N==85:
 Q=a['nonrational']['q4'];testfamily('onequartic_plus_rat6',Q[0][:4],6)
 testfamily('doublequartic_plus_rat2',Q[0][:8],2)
 for k in range(1,8):testfamily('twoquartics_'+str(k)+'_plus_rat2',Q[0][:4]+Q[k][:4],2)
 for name in ['q8width5','q8width1']:testfamily(name+'_plus_rat2',a['nonrational'][name][0][:8],2)
else:
 Q=a['nonrational']['q5'];testfamily('onequintic_plus_rat5',Q[0][:5],5)
 testfamily('doublequintic',Q[0][:10],0)
 for k in range(1,6):testfamily('twoquintics_'+str(k),Q[0][:5]+Q[k][:5],0)
out=dict(N=N,complete=True,operator_dimension=amb,rational_prefix_rows=len(R),rational_prefix_rank=M.rank(),rational_relation_dimension=len(ker),minimum_rational_prefix_cost=best,minimum_cost_witness=badrel,rational_all_degree10_rank10=best>10,mixed_families=records,mixed_failures=failures,all_tested_ranks10=best>10 and not failures,seconds=time.monotonic()-t0)
(P/'cusp_FI_primary.json').write_text(json.dumps(out,indent=2)+'\n');print('CONTC_CUSP_PRIMARY_COMPLETE',N,out['all_tested_ranks10'],time.monotonic()-t0,flush=True)
