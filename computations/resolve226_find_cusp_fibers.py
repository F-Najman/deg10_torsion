from pathlib import Path
import json,ast,itertools
from sage.all import ZZ,matrix,vector
P=Path('work/deg10/C226_resolution_20261030');d=json.loads((P/'optimized_cusp_divisors.json').read_text());B=matrix(ZZ,d['principal_relation_basis']);out=[]
for idx,desc,co,candidates in d['classes']:
 A=vector(ZZ,[max(0,x) for x in co]);found={tuple(A)};rels=[]
 for size in [1,2,3]:
  for inds in itertools.combinations(range(17),size):
   for signs in itertools.product([-1,1],repeat=size):
    v=sum((s*B.row(i) for i,s in zip(inds,signs)),vector(ZZ,[0]*18));C=A-v
    if min(C)>=0 and v:
     found.add(tuple(C));rels.append([list(map(int,v)),list(map(int,C))])
  if len(found)>1:break
 assert len(found)>1,idx
 print('C226_CUSPIDAL_FIBERS',idx,'number',len(found),'positive_degree',sum(A),flush=True)
 out.append([idx,desc,co,[list(v) for v in sorted(found)],rels])
def enc(x):
 if isinstance(x,list):
  a,b=('[',']') if all(isinstance(v,int) for v in x) else ('[*','*]');return a+','.join(enc(v) for v in x)+b
 return str(x)
(P/'cuspidal_pencil_fibers.m').write_text(enc([d['labels'],out])+'\n');(P/'cuspidal_pencil_fibers.json').write_text(json.dumps({'status':'PASS','labels':d['labels'],'classes':out},indent=2)+'\n')
