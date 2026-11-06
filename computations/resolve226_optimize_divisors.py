from pathlib import Path
import ast,json,random
from sage.all import ZZ,QQ,matrix,vector
H=Path('work/deg10');P=H/'C226_resolution_20261030';C=H/'continuation_20261029b/C226'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
labels,base,hgens,rels,ker,rows=read(C/'cusp_preimages.m');ib=labels.index(base);actual=[]
for k in ker:
 v=vector(ZZ,k);v[ib]-=sum(v);actual.append(v)
B=matrix(ZZ,actual).row_module().basis_matrix().LLL();assert B.nrows()==17 and all(sum(b)==0 for b in B.rows());BS,_=B.gram_schmidt();rng=random.Random(2261030)
def babai(w):
 v=vector(QQ,w)
 for b,s in reversed(list(zip(B.rows(),BS.rows()))):v-=round(v.dot_product(s)/s.dot_product(s))*b
 return vector(ZZ,v)
def score(v):return (sum(max(0,x) for x in v),max(abs(x) for x in v),sum(x*x for x in v))
out=[]
for idx,desc,cl,old in rows:
 co=vector(ZZ,old);co[ib]+=10-sum(co);best=babai(co);found={tuple(best)}
 for j in range(10000):
  target=vector(QQ,[rng.randint(-4,4) for _ in range(18)]);target-=vector(QQ,[sum(target)/18]*18);v=babai(co-target)+vector(ZZ,[0]*18) if False else None
  # Round an affine translate to a nearby lattice point; retain the original coset.
  d=co-target;near=vector(ZZ,[0]*18)
  for b,s in reversed(list(zip(B.rows(),BS.rows()))):
   q=round(d.dot_product(s)/s.dot_product(s));d-=q*b;near+=q*b
  v=co-near
  for _ in range(3):
   changed=False
   for b in B.rows():
    choices=[v,v+b,v-b];w=min(choices,key=score)
    if score(w)<score(v):v=w;changed=True
   if not changed:break
  if score(v)<score(best):best=v;print('C226_CUSP_REP_IMPROVEMENT',idx,j,score(best),list(best),flush=True)
  if score(v)[0]<=15:found.add(tuple(v))
 assert sum(best)==10 and co-best in B.row_module();assert best*matrix(ZZ,hgens)-vector(ZZ,cl) in matrix(ZZ,rels).row_module()
 candidates=sorted(found,key=score)[:30];out.append([idx,desc,list(best),[list(v) for v in candidates]])
 print('C226_OPTIMAL_FOUND',idx,score(best),list(best),'candidates',len(candidates),flush=True)
def enc(x):
 if isinstance(x,list):
  a,b=('[',']') if all(isinstance(v,int) for v in x) else ('[*','*]');return a+','.join(enc(v) for v in x)+b
 return str(x)
(P/'optimized_cusp_divisors.m').write_text(enc([labels,out])+'\n');(P/'optimized_cusp_divisors.json').write_text(json.dumps({'status':'PASS','labels':labels,'principal_relation_basis':[list(map(int,b)) for b in B.rows()],'classes':out},indent=2,default=int)+'\n')
print('C226_OPTIMIZED_DIVISORS_COMPLETE',flush=True)
