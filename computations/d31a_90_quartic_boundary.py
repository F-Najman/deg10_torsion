from pathlib import Path
import re,ast,json,time
H=Path(__file__).resolve().parents[1];O=H/'continuation_20261031a/C90';s=(O/'inputs/C90/cartier_kernels_D9.m').read_text();t0=time.time()
v=lambda k:ast.literal_eval(re.findall(r"\b"+k+r":=(.*?);",s,re.S)[0])
ds=v('ds');assert len(ds)==23
blocks=[]
for i,d in enumerate(ds,1):
 m=8//d;a=re.findall(r"localjets\[<"+str(i)+','+str(m)+r">\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);",s,re.S);assert len(a)==1
 nr,nc=int(a[0][0]),int(a[0][1]);es=ast.literal_eval(a[0][2]);assert nr==39 and nc==m*d
 blocks.append([sum(es[nc*r+j]<<r for r in range(nr))for j in range(nc)])
count=0;bad=[];total_expected=0
# Every path builds a restriction basis by independent bit elimination.
def rec(i,left,pivs,vs):
 global count
 if i==len(ds):
  if left:return
  if not any(m and d==4 for m,d in zip(vs,ds)):return
  count+=1
  if any(m and d==4 for m,d in zip(vs,ds)) and len(pivs)!=10:bad.append([vs,len(pivs)])
  return
 q=pivs.copy();degree=ds[i]
 for m in range(min(left//degree,len(blocks[i])//degree)+1):
  rec(i+1,left-m*degree,q,vs+[m])
  if (m+1)*degree>left or (m+1)*degree>len(blocks[i]):break
  for row in blocks[i][m*degree:(m+1)*degree]:
   x=row
   while x:
    j=x.bit_length()-1
    if j in q:x^=q[j]
    else:q[j]=x;break
rec(0,10,{},[])
missing=[{'cusp':i+1,'degree':d,'multiplicity':10//d}for i,d in enumerate(ds)if 10//d>8//d]
result={'status':'PASS'if not bad else'RESIDUAL','tested_existing_degree10_patterns':count,'rank_defects':bad,'filter':'All configurations containing at least one quartic cusp; other supports handled by global lift and gonality18',
 'missing_pure_support_patterns':missing,'elapsed_seconds':time.time()-t0}
(O/'quartic_boundary_verified.json').write_text(json.dumps(result,indent=2)+'\n');print('D31A_C90_QUARTIC_BOUNDARY_COMPLETE',count,'defects',len(bad),'missing',len(missing),'seconds',time.time()-t0,flush=True)
