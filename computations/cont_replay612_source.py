from pathlib import Path
from collections import defaultdict,Counter
import ast,json,itertools
from sage.all import ZZ,QQ,matrix,vector,block_matrix,identity_matrix
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=read(D/'diagonal612_D5_source_checkpoint.m')
sigs,ds,invs,relations,cgs,basis,images,source,scov=a
assert Counter(ds)=={1:16,2:8}
assert scov==[[3,96,144],[4,468,16],[5,3264,1]] and len(source)==24576
R7,R13=[matrix(ZZ,[r[:-1] for r in R]) for R in relations]
C7,C13=[matrix(ZZ,[r[:-1] for r in C]) for C in cgs]
assert all(v[-1]==0 for data in relations+cgs+[basis,images] for v in data)
B7=R7.stack(C7).row_module().basis_matrix();B13=R13.stack(C13).row_module().basis_matrix()
order7,order13=abs(R7.det()),abs(R13.det());assert order7==95551488 and order13==61917364224
assert order7//abs(B7.det())==order13//abs(B13.det())==5971968
# Exact graph of the single global cusp subgroup; both projections have its full order.
R=block_matrix(ZZ,[[R7,0],[0,R13]])
Graph=R.stack(C7.augment(C13)).row_module().basis_matrix();GraphInv=Graph.inverse()
assert abs(R.det())//abs(Graph.det())==5971968
# Reconstruct the actual T13 map and its kernel independently of invariant-factor guesses.
B=matrix(ZZ,[v[:-1] for v in basis]);A=matrix(ZZ,[v[:-1] for v in images]);r=R7.ncols()
Rel=block_matrix(ZZ,[[R7,0],[0,R7]])
AGraph=Rel.stack(B.augment(A)).row_module().basis_matrix();AGraphInv=AGraph.inverse()
assert abs(Rel.det())//abs(AGraph.det())==order7
assert abs(R7.stack(B).row_module().basis_matrix().det())==1
assert order7//abs(R7.stack(A).row_module().basis_matrix().det())==4
assert all((vector(ZZ,list(c)+[0]*r)*AGraphInv).denominator()==1 for c in C7.rows())
Sm,U,V=B7.smith_form();Vi=V.inverse();quot=[abs(Sm[i,i]) for i in range(r)]
assert __import__('math').prod(quot)==16
cosets=[vector(ZZ,t)*Vi for t in itertools.product(*[range(int(q)) for q in quot])]
allowed=[v for v in cosets if (vector(ZZ,list(v)+[0]*r)*AGraphInv).denominator()==1]
assert len(allowed)==4 and all((2*v*B7.inverse()).denominator()==1 for v in allowed)
assert any((v*B7.inverse()).denominator()!=1 for v in allowed)
assert all(Sm[i,i] for i in range(r))
# Canonical representatives use the actual relation lattices, including their map bases.
canonical_bases={id(R):R.row_module().basis_matrix() for R in (R7,R13)}
def key(v,R):
 assert v[-1]==0;w=vector(ZZ,v[:-1])
 for row in canonical_bases[id(R)].rows():
  j=next(j for j,c in enumerate(row) if c);assert row[j]>0;w-=(w[j]//row[j])*row
 return tuple(w)
def weighted(ws,d):
 if not ws:
  if d==0:yield ()
  return
 for n in range(d//ws[0]+1):
  for tail in weighted(ws[1:],d-n*ws[0]):yield (n,)+tail
expected={d:set(weighted(ds,5-d)) for d in (3,4,5)};seen=set();pred=defaultdict(list);pointclasses={};passed=0
for ind,(d,i,v,cl,ok,pr) in enumerate(source,1):
 ident=(d,i,tuple(v));assert ident not in seen and tuple(v) in expected[d];seen.add(ident)
 assert 1<=i<={3:96,4:468,5:3264}[d]
 truth=(vector(ZZ,cl[:-1]+[0]*r)*AGraphInv).denominator()==1;assert ok==truth
 base=[cl[j]-sum(v[k]*cgs[0][k][j] for k in range(24)) for j in range(len(cl))];ky=key(base,R7)
 if (d,i) in pointclasses:assert pointclasses[d,i]==ky
 else:pointclasses[d,i]=ky
 if ok:
  assert (vector(ZZ,[2*c for c in cl[:-1]]+pr[:-1])*GraphInv).denominator()==1
  pred[key(pr,R13)].append(ind);passed+=1
 else:assert not pr
assert len(seen)==24576 and len(pointclasses)==96+468+3264
# Newton identities verify both complete place censuses and native group orders.
closed={}
for q,f,order in read(D/'frob612_arithmetic.txt'):
 assert len(f)==19 and f[-1]==1
 tr=[0];cl=[0]
 for k in range(1,6):
  tr.append(-k*f[18-k]-sum(f[18-j]*tr[k-j] for j in range(1,k)))
  value=q**k+1-tr[k]-sum(d*cl[d] for d in range(1,k) if k%d==0);assert value%k==0;cl.append(value//k)
 closed[q]=cl[1:];assert sum(f)==order
 if q in (7,13):assert order=={7:order7,13:order13}[q]
assert closed[7]==[16,8,96,468,3264]
print('DECIC612_SOURCE_INDEPENDENT_REPLAY_COMPLETE',len(source),passed,'kernel_over_cusps_C2_squared',flush=True)
