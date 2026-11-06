from pathlib import Path
import ast,json,itertools,time
P=Path('work/deg10/C226_resolution_20261030');tt=time.monotonic()
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
labels,hgens,glob,rels=read(P/'all_cusp_classes.m');ivs=[r[i] for i,r in enumerate(rels)];assert len(ivs)==6
hgens=[tuple(int(x)%d for x,d in zip(h[:-1],ivs)) for h in hgens];glob=[(v,tuple(int(x)%d for x,d in zip(h[:-1],ivs))) for v,h in glob]
opt=json.loads((P/'optimized_cusp_divisors.json').read_text());assert opt['labels']==labels

def add(a,b):return tuple((x+y)%d for x,y,d in zip(a,b,ivs))
def keyclass(co):
 ans=(0,)*6
 for x,g in zip(co,hgens):ans=add(ans,tuple(x*y%d for y,d in zip(g,ivs)))
 return ans

def entries(gs,D):
 mults=[[tuple(k*x%d for x,d in zip(g,ivs)) for k in range(D+1)] for g in gs];co=[0]*len(gs)
 def rec(i,left,cl):
  if i==len(gs):yield D-left,cl,tuple(co);return
  for k in range(left+1):co[i]=k;yield from rec(i+1,left-k,add(cl,mults[i][k]))
 yield from rec(0,D,(0,)*6)
out=[]
for idx,desc,co,alternatives in opt['classes']:
 A=[max(0,x) for x in co];D=sum(A);target=keyclass(A);table={};cnt=0
 for deg,cl,v in entries(hgens[:9],D):
  key=(deg,)+cl
  if key in table:table[key].append(v)
  else:table[key]=[v]
  cnt+=1
 print('MITM_FIRST_HALF',idx,cnt,len(table),time.monotonic()-tt,flush=True)
 offsets=[]
 for ns in itertools.product(range(D//6+1),repeat=3):
  if sum(ns)*6>D:continue
  cl=target
  for n,(value,g) in zip(ns,glob):cl=tuple((a-n*b)%d for a,b,d in zip(cl,g,ivs))
  offsets.append((ns,D-6*sum(ns),cl))
 found=[]
 for deg,cl,v in entries(hgens[9:],D):
  for ns,rd,t in offsets:
   if deg>rd:continue
   key=(rd-deg,)+tuple((a-b)%d for a,b,d in zip(t,cl,ivs))
   for w in table.get(key,[]):
    rat=list(w+v)
    if rat==A and not any(ns):continue
    assert sum(rat)+6*sum(ns)==D
    ans=[0]*6
    for n,g in zip(rat,hgens):ans=[(a+n*b)%d for a,b,d in zip(ans,g,ivs)]
    for n,(value,g) in zip(ns,glob):ans=[(a+n*b)%d for a,b,d in zip(ans,g,ivs)]
    assert tuple(ans)==target;found.append([rat,list(ns)])
 print('MITM_CUSPIDAL_FIBERS',idx,len(found),time.monotonic()-tt,flush=True);out.append([idx,desc,co,found])
 del table
(P/'cuspidal_pencil_fibers_mitm.json').write_text(json.dumps({'status':'PASS','labels':labels,'global_degree6_cusp_labels':[g[0] for g in glob],'classes':out,'elapsed_seconds':time.monotonic()-tt},indent=2)+'\n')
print('C226_ALL_CUSPIDAL_PENCIL_FIBERS_COMPLETE',flush=True)
