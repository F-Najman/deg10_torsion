from pathlib import Path
from collections import defaultdict,Counter
import ast,json,itertools
from sage.all import ZZ,QQ,matrix,vector,block_matrix,identity_matrix
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=read(D/'diagonal_6_12_D5_arithmetic.txt')
sigs,ds,invs,relations,cgs,basis,images,source,scov,rawds,rawcusp,targets,matches,counts,coverage,allcl=a
assert Counter(ds)=={1:16,2:8} and Counter(rawds)=={1:32}
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
rawplaces={};nc={};cov={};repindices={}
for d,keys,labels,reps in coverage:
 assert len(keys)==len(labels)==len(set(map(tuple,keys)))
 assert set(labels)==set(range(1,len(reps)+1)) and all(labels[j-1]==i for i,j in enumerate(reps,1))
 n=closed[13][d-1]-(32 if d==1 else 0);assert len(keys)==n
 rawplaces[d]=n;cov[d]=len(reps);repindices[d]=reps
assert [x[0] for x in counts]==[2,3,4,5]
comps={d:set(weighted(rawds,5-d)) for d in range(2,6)}
assert all([d,rawplaces[d],cov[d],len(comps[d])] in counts for d in range(2,6))
assert len(allcl)==2 and len(allcl[0])==rawplaces[2] and len(allcl[1])==rawplaces[3]
repc={};seen=set();recomputed=[]
for k,(desc,cl) in enumerate(targets,1):
 typ=desc[0];ident=tuple(tuple(x) if isinstance(x,list) else x for x in desc);assert ident not in seen;seen.add(ident)
 if typ in range(2,6):
  d,i,v=desc;assert 1<=i<=cov[d] and tuple(v) in comps[d]
  base=[cl[j]-sum(v[t]*rawcusp[t][j] for t in range(32)) for j in range(len(cl))];ky=key(base,R13)
  if (d,i) in repc:assert repc[d,i]==ky
  else:repc[d,i]=ky
  if d in (2,3):assert ky==key(allcl[d-2][repindices[d][i-1]-1],R13)
 elif typ==22:
  _,i,j,v=desc;assert 1<=i<=cov[2] and 1<=j<=rawplaces[2] and tuple(v) in set(weighted(rawds,1))
  known=list(repc[2,i])+[0];known=[known[t]+allcl[0][j-1][t]+sum(v[z]*rawcusp[z][t] for z in range(32)) for t in range(len(cl))];assert key(known,R13)==key(cl,R13)
 else:
  assert typ==23;_,i,j=desc;assert 1<=i<=cov[2] and 1<=j<=rawplaces[3]
  known=list(repc[2,i])+[0];known=[known[t]+allcl[1][j-1][t] for t in range(len(cl))];assert key(known,R13)==key(cl,R13)
 hits=pred.get(key([2*c for c in cl],R13),[])
 if hits:recomputed.append([k,hits])
assert matches==recomputed
assert len(targets)==sum(cov[d]*len(comps[d]) for d in range(2,6))+cov[2]*(rawplaces[2]*32+rawplaces[3])
report={'group':[6,12],'relative_degree':5,'native_source_divisors':len(source),'source_Hecke_survivors':passed,'Hecke_kernel_over_cusp_invariants':[2,2],'global_cusp_order':5971968,'target_representatives':len(targets),'diagonal_matches':len(matches),'target_counts':counts,'exact_Hecke_kernel_common_cusp_graph_and_complete_mixed_pattern_replay':'PASS'}
rrfile=D/'M6N12_D5_RR19_arithmetic.txt'
if rrfile.exists():
 _,_,lifts=read(D/'M6N12_D5_cusp_lifts.m');rr,left=read(rrfile);matched={m[0] for m in matches};covered=set();dimensions=Counter();systems=0
 for (ids,v,unused),record in zip(lifts,rr):
  k,ids2,v2,dim,rows,positive=record;assert ids==ids2 and v==v2
  cusp=[sum(v[i]*cgs[1][i][j] for i in range(24)) for j in range(len(cgs[1][0]))]
  for i in ids:assert i in matched and key(cusp,R13)==key([2*c for c in targets[i-1][1]],R13);covered.add(i)
  dimensions[dim]+=1
  for u,degrees,mults,even in rows:
   assert len(degrees)==len(mults) and sum(d*m for d,m in zip(degrees,mults))==10
   assert even==all(m%2==0 for m in mults)
  if not positive:
   canonical={tuple([0]*j+[1]+list(t)) for j in range(dim) for t in itertools.product(range(19),repeat=dim-j-1)}
   assert {tuple(u) for u,degrees,mults,even in rows}==canonical and not any(row[-1] for row in rows)
  else:assert rows[-1][-1] and k in left
  systems+=len(rows)
 assert len(rr)==len(lifts) and covered==matched
 report.update({'RR_classes':len(lifts),'RR_dimensions':dict(dimensions),'RR_members_tested':systems,'remaining_parity_classes':left,'independent_lift_and_complete_parity_replay':'PASS'})
(D/'M6N12_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC612_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
