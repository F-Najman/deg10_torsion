from pathlib import Path
from collections import defaultdict,Counter
import ast,json
from sage.all import ZZ,QQ,matrix,vector,block_matrix,identity_matrix
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
a=read(D/'diagonal_4_16_D5_arithmetic.txt');sigs,ds,inv5,r5,inv9,r9,cg,source,targets,matches,counts,coverage,rawcusp,fps,sourcecoverage=a
assert Counter(ds)=={1:16,2:4,4:4} and counts==[[3,320,5,308],[4,1808,33,24],[5,11904,186,1]]
assert sourcecoverage==[[3,64,140],[4,160,16],[5,512,1]]
assert len(source)==12032 and len(targets)==2518 and len(matches)==158
for R,cs in [(r5,cg[0]),(r9,cg[1])]:
 assert all(r[-1]==0 for r in R) and all(v[-1]==0 for v in cs)
R5=matrix(ZZ,[r[:-1] for r in r5 if any(r)]);R9=matrix(ZZ,[r[:-1] for r in r9 if any(r)])
C5=matrix(ZZ,[r[:-1] for r in cg[0]]);C9=matrix(ZZ,[r[:-1] for r in cg[1]])
B5=R5.stack(C5).row_module().basis_matrix();B9=R9.stack(C9).row_module().basis_matrix()
ordc5=abs(R5.det())//abs(B5.det());ordc9=abs(R9.det())//abs(B9.det());assert ordc5==ordc9==327680000
assert abs(B5.det())==32
assert all((vector(ZZ,2*e)*B5.inverse()).denominator()==1 for e in identity_matrix(ZZ,R5.ncols()).rows())
Relations=block_matrix(ZZ,[[R5,0],[0,R9]])
Graph=Relations.stack(C5.augment(C9)).row_module().basis_matrix();GraphInv=Graph.inverse()
assert abs(Relations.det())//abs(Graph.det())==327680000
# Thus both graph projections identify exactly the actual global cusp group.
R9basis=R9.row_module().basis_matrix()
def key(v):
 assert v[-1]==0;w=vector(ZZ,v[:-1])
 for row in R9basis.rows():
  j=next(j for j,c in enumerate(row) if c);b=row[j];assert b>0
  w-=(w[j]//b)*row
 return tuple(w)
def weighted(ws,d):
 if not ws:
  if d==0:yield ()
  return
 for n in range(d//ws[0]+1):
  for tail in weighted(ws[1:],d-n*ws[0]):yield (n,)+tail
expected={d:set(weighted(ds,5-d)) for d in (3,4,5)}
source_seen=set();predicted=defaultdict(list);pointclasses={}
for ind,(d,i,v,cl,pred) in enumerate(source,1):
 assert tuple(v) in expected[d] and 1<=i<={3:64,4:160,5:512}[d]
 ident=(d,i,tuple(v));assert ident not in source_seen;source_seen.add(ident)
 assert cl[-1]==pred[-1]==0
 vv=vector(ZZ,[2*c for c in cl[:-1]]+pred[:-1]);assert (vv*GraphInv).denominator()==1
 predicted[key(pred)].append(ind)
 base=tuple(cl[j]-sum(v[k]*cg[0][k][j] for k in range(24)) for j in range(len(cl)))
 # Check a fixed place has the same class for every cusp complement, modulo the actual relation lattice.
 h=R5.row_module().basis_matrix();w=vector(ZZ,base[:-1])
 for row in h.rows():
  j=next(j for j,c in enumerate(row) if c);w-=(w[j]//row[j])*row
 if (d,i) in pointclasses:assert pointclasses[d,i]==tuple(w)
 else:pointclasses[d,i]=tuple(w)
assert len(source_seen)==sum({3:64,4:160,5:512}[d]*len(expected[d]) for d in expected)
refcounts={d:(n,r,c) for d,n,r,c in counts};cov={}
for d,keys,labels,reps in coverage:
 n,r,c=refcounts[d];assert len(keys)==len(set(map(tuple,keys)))==len(labels)==n and len(reps)==r
 assert set(labels)==set(range(1,r+1)) and all(labels[j-1]==i for i,j in enumerate(reps,1));cov[d]=(r,n)
ratindices={j for d,i,v,cl in targets if d==4 for j,t in enumerate(v) if t};assert len(ratindices)==24
rawds=[1 if j in ratindices else 2 for j in range(32)];comps={d:set(weighted(rawds,5-d)) for d in (3,4,5)}
seen=set();recomputed=[]
for k,(d,i,v,cl) in enumerate(targets,1):
 assert (d,i,tuple(v)) not in seen and tuple(v) in comps[d] and 1<=i<=cov[d][0];seen.add((d,i,tuple(v)))
 ky=key([2*c for c in cl]);ps=predicted.get(ky,[])
 if ps:recomputed.append([k,ps])
assert recomputed==matches and len(seen)==2518
# The signed lifts are genuine common-cusp preimages of every remaining doubled class.
_,_,lifts=read(D/'M4N16_D5_cusp_lifts.m');rr,left=read(D/'M4N16_D5_RR13_arithmetic.txt');assert not left and len(lifts)==len(rr)==158
matched={m[0] for m in matches};covered=set();dimensions=Counter();systems=0
for (ids,v,unused),r in zip(lifts,rr):
 k,ids2,v2,dim,rows,positive=r;assert ids==ids2 and v==v2 and not positive
 cusp=[sum(v[i]*cg[1][i][j] for i in range(24)) for j in range(len(cg[1][0]))]
 for i in ids:
  assert i in matched and key(cusp)==key([2*c for c in targets[i-1][3]]);covered.add(i)
 assert len(rows)==(13**dim-1)//12
 canonical={tuple([0]*j+[1]+list(t)) for j in range(dim) for t in __import__('itertools').product(range(13),repeat=dim-j-1)}
 assert {tuple(u) for u,degrees,mults,even in rows}==canonical
 for u,degrees,mults,even in rows:
  assert len(degrees)==len(mults) and sum(d*m for d,m in zip(degrees,mults))==10
  assert not even and any(m%2 for m in mults)
 dimensions[dim]+=1;systems+=len(rows)
assert covered==matched and dimensions=={1:157,2:1} and systems==171
report={'group':[4,16],'absolute_degree':10,'relative_degree':5,'source_divisors':12032,'target_divisor_representatives':2518,'diagonal_matches':158,'global_cusp_order':327680000,'source_J_over_cusp_exponent':2,'independent_actual_relation_graph_and_complete_match_replay':'PASS','RR_dimensions':dict(dimensions),'rational_members_tested_at_13':171,'even_effective_divisors':0,'independent_lift_and_parity_replay':'PASS'}
(D/'M4N16_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC416_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
