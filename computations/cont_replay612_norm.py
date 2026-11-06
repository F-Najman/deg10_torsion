from pathlib import Path
from collections import defaultdict,Counter
import ast,json,itertools,math
from sage.all import ZZ,matrix,vector,block_matrix
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def completed(stem):
 s=json.loads((H/f'jobs/{stem}.status.json').read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0,stem
 return s
sig,ds,invs,relations,cgs,basis,images,source,coverage=read(D/'diagonal612_D5_source_checkpoint_v5.m')
assert coverage==[[3,96,144],[4,468,16],[5,3264,1]] and len(source)==24576 and Counter(ds)=={1:16,2:8}
R=matrix(ZZ,[r[:-1] for r in relations[0]]);C=matrix(ZZ,[c[:-1] for c in cgs[0]]);r=R.ncols();assert abs(R.det())==95551488
RB=R.row_module().basis_matrix();CB=R.stack(C).row_module().basis_matrix();assert abs(R.det())//abs(CB.det())==5971968
B=matrix(ZZ,[v[:-1] for v in basis]);A=matrix(ZZ,[v[:-1] for v in images]);RR=block_matrix(ZZ,[[R,0],[0,R]])
AG=RR.stack(B.augment(A)).row_module().basis_matrix();AI=AG.inverse();assert abs(RR.det())//abs(AG.det())==95551488
assert abs(R.stack(B).row_module().basis_matrix().det())==1 and abs(R.det())//abs(R.stack(A).row_module().basis_matrix().det())==4
assert all((vector(ZZ,list(c)+[0]*r)*AI).denominator()==1 for c in C.rows())
Sm,U,V=CB.smith_form();Vi=V.inverse();qq=[abs(Sm[i,i]) for i in range(r)];assert math.prod(qq)==16
cosets=[vector(ZZ,t)*Vi for t in itertools.product(*[range(int(q)) for q in qq])];allowed=[v for v in cosets if (vector(ZZ,list(v)+[0]*r)*AI).denominator()==1]
assert len(allowed)==4 and all((2*v*CB.inverse()).denominator()==1 for v in allowed)
def key(v):
 assert v[-1]==0;w=vector(ZZ,v[:-1])
 for row in RB.rows():
  j=next(j for j,c in enumerate(row) if c);w-=(w[j]//row[j])*row
 return tuple(w)
def weighted(ws,d):
 if not ws:
  if d==0:yield ()
  return
 for n in range(d//ws[0]+1):
  for tail in weighted(ws[1:],d-n*ws[0]):yield (n,)+tail
expected={d:set(weighted(ds,5-d)) for d in (3,4,5)};seen=set();point={};passed=set()
for ind,(d,i,v,cl,ok,pr) in enumerate(source,1):
 ident=(d,i,tuple(v));assert ident not in seen and tuple(v) in expected[d] and 1<=i<={3:96,4:468,5:3264}[d];seen.add(ident)
 truth=(vector(ZZ,cl[:-1]+[0]*r)*AI).denominator()==1;assert ok==truth
 base=[cl[j]-sum(v[t]*cgs[0][t][j] for t in range(24)) for j in range(len(cl))]
 if (d,i) in point:assert point[d,i]==key(base)
 else:point[d,i]=key(base)
 if ok:passed.add(ind)
assert len(point)==3828 and len(passed)==3552
# These eight successful native jobs recompute every actual point class; no terminated run is promoted.
recert=set()
for sh in range(8):
 stem=f'cont_source612_recert_v2_shard{sh}';completed(stem)
 for d,i,cl in read(D/(stem+'_arithmetic.txt')):
  assert (d,i) not in recert and point[d,i]==key(cl);recert.add((d,i))
assert recert==set(point)
sigs,deg,lifts=read(D/'M6N12_D5_all_source_lifts.m');assert sigs==sig and deg==ds and len(lifts)==3360
covered=set();unique=set()
for ids,v,unused in lifts:
 vv=[sum(v[t]*cgs[0][t][j] for t in range(24)) for j in range(len(cgs[0][0]))];ky=key(vv);assert ky not in unique;unique.add(ky)
 for i in ids:
  assert i in passed and i not in covered and key([2*x for x in source[i-1][3]])==ky;covered.add(i)
assert covered==passed
# Independent monic polynomial square test over each odd prime field.
def square(a,p):
 assert a and a[-1]==1 and all(0<=x<p for x in a)
 if (len(a)-1)%2:return False
 m=(len(a)-1)//2;b=[0]*(m+1);b[m]=1;half=pow(2,-1,p)
 for j in range(m-1,-1,-1):b[j]=(a[m+j]-sum(b[i]*b[m+j-i] for i in range(j+1,m))) * half%p
 cc=[0]*(2*m+1)
 for i,x in enumerate(b):
  for j,y in enumerate(b):cc[i+j]=(cc[i+j]+x*y)%p
 return cc==a

def stage(p,shards,required):
 done=set();left=set();dims=Counter();members=0;norm_rejected=0;full=0
 for sh in range(shards):
  stem=f'cont_RR612_norm{p}_shard{sh}';completed(stem);rr,remaining=read(D/(stem+'_arithmetic.txt'));computed=[]
  for k,v,dim,rows,positive in rr:
   assert k in required and k not in done and lifts[k-1][1]==v;done.add(k);dims[dim]+=1
   if isinstance(positive,str):assert positive=='dimension cap' and dim>4 and not rows;computed.append(k);continue
   seen=set()
   for u,pol,passed_norm,degs,mults,even,cusp in rows:
    assert len(u)==dim and all(0<=x<p for x in u) and tuple(u) not in seen;seen.add(tuple(u));assert next(x for x in u if x)!=0 and next(x for x in u if x)==1
    assert len(pol)<=11 and square(pol,p)==passed_norm
    if not passed_norm:norm_rejected+=1;assert not even
    else:assert degs;full+=1
    if degs:
     assert len(degs)==len(mults) and all(d>0 for d in degs) and all(m>0 for m in mults) and sum(d*m for d,m in zip(degs,mults))==10
     assert even==all(m%2==0 for m in mults)
    if passed_norm and even and not cusp:assert positive
   if positive:assert rows[-1][2] and rows[-1][5] and not rows[-1][6];computed.append(k)
   else:
    canon={tuple([0]*j+[1]+list(t)) for j in range(dim) for t in itertools.product(range(p),repeat=dim-j-1)}
    assert seen==canon and len(rows)==(p**dim-1)//(p-1) and not any(row[2] and row[5] and not row[6] for row in rows)
   members+=len(rows)
  assert computed==remaining;left.update(remaining)
 assert done==set(required)
 return left,{'prime':p,'input_classes':len(required),'survivors':len(left),'RR_dimensions':dict(dims),'projective_members_checked':members,'norm_polynomials_rejecting_even_divisors':norm_rejected,'norm_square_members_with_full_support_checked':full,'independent_polynomial_squares_and_complete_RR_coverage':'PASS'}
left,first=stage(13,12,set(range(1,3361)));assert len(left)==192
stages=[first]
for p,sh in [(19,8),(31,8),(37,8),(43,8)]:
 if not (H/f'jobs/cont_RR612_norm{p}_shard0.status.json').exists():break
 left,rp=stage(p,sh,left);stages.append(rp)
 if not left:break
report={'group':[6,12],'absolute_degree':10,'relative_degree':5,'genus':9,'native_source_point_classes_recomputed':3828,'complete_source_divisors':24576,'source_patterns_surviving_Hecke':3552,'doubled_global_cusp_classes':3360,'Hecke_kernel_over_cusp_invariants':[2,2],'global_cusp_order':5971968,'native_source_recertification_and_independent_relation_replay':'PASS','stages':stages,'remaining_classes':sorted(left),'classification_case':'EXCLUDED' if not left else 'UNRESOLVED'}
(D/'M6N12_norm_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC612_NORM_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
