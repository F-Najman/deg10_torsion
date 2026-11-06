from pathlib import Path
from collections import Counter,defaultdict
import ast,json,time,math
C=Path(__file__).resolve().parents[1]/'continuation_20261031b/C66';start=time.monotonic();read=lambda n:ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'));f=read('F2_proven_frame.m');nc=dict(read('F2_NC_classes.m'));N=42373650
assert f[3]==[[15,0,0],[0,N,0]]or f[3]==[[15,0,0],[0,N,0],[0,0,0]]
chi=[a[0]for a in f[5]];assert len(chi)==3
pack=lambda a,b:(a%15)+15*(b%N)
def un(z):return z%15,z//15
def add(a,b):aa,ab=un(a);ba,bb=un(b);return pack(aa+ba,ab+bb)
def ch(z):a,b=un(z);return(chi[0]*a+chi[1]*b)%3
def cl(v):assert v[2]==0;return pack(v[0],v[1])
cs=[(r[0],cl(r[2]))for r in f[6]];qs=[(r[0],cl(r[2]))for r in f[7]];assert Counter(d for d,z in cs)=={1:10,2:5,5:2,10:1}
assert all(ch(z)==0 for d,z in qs)
def cusps(ds,total):
 ans=[]
 def rec(i,left,z,vals):
  if i==len(ds):
   if not left:ans.append((z,tuple(vals)))
   return
  d,c=ds[i];zz=z
  for m in range(left//d+1):rec(i+1,left-m*d,zz,vals+[m]);zz=add(zz,c)
 rec(0,total,0,[]);return ans
comp={d:cusps(cs,d)for d in range(6)};pure=cusps(cs,10);qpure=cusps(qs,10)
counts=Counter();examples=defaultdict(list);kindcounts=defaultdict(Counter);qcounts=Counter();qex={}
def insert(z,kind,record):
 if ch(z):return
 counts[z]+=1;kindcounts[z][kind]+=1
 if len(examples[z])<9:examples[z].append(record)
for z,v in pure:insert(z,'cusps',{'type':'cusps','multiplicities':v})
for z,v in qpure:qcounts[z]+=1;qex.setdefault(z,v)
for d in range(5,11):
 for i,row in enumerate(nc[d]):
  z=cl(row[1]);assert ch(z)==row[2][0]
  for zz,v in comp[10-d]:insert(add(z,zz),str(d),{'type':'NC+cusp','degree':d,'point_index':i+1,'cusp_multiplicities':v})
for i in range(len(nc[5])):
 for j in range(i,len(nc[5])):insert(add(cl(nc[5][i][1]),cl(nc[5][j][1])),'5+5',{'type':'NC5+NC5','point_indices':[i+1,j+1]})
for z,n in counts.items():assert(n+1)&n==0,(z,n)
gonal=[z for z,n in counts.items()if n>=3]
summary=[]
for z in sorted(gonal):summary.append({'class':list(un(z))+[0],'finite_effective_divisors':counts[z],'h0_F2':int(math.log2(counts[z]+1)),'global_effective_cuspidal_representatives':qcounts[z],'types':dict(kindcounts[z]),'examples':examples[z]})
noncusp=[]
for z in counts:
 if any(k!='cusps'for k in kindcounts[z]):noncusp.append({'class':list(un(z))+[0],'finite_effective_divisors':counts[z],'types':dict(kindcounts[z]),'examples':examples[z]})
report={'status':'PASS','finite_rational_image_effective_degree10_divisors':sum(counts.values()),'rational_image_classes':len(counts),'pure_cusp_divisors_enumerated':len(pure),'global_effective_cusp_divisors_enumerated':len(qpure),'h0_histogram':dict(Counter(int(math.log2(n+1))for n in counts.values())),'gonal_classes':summary,'noncuspidal_classes':noncusp,'elapsed_seconds':time.monotonic()-start}
(C/'F2_degree10_classification.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_C66_F2_DEGREE10_CLASSES_COMPLETE classes',len(counts),'gonal',len(gonal),'h0hist',report['h0_histogram'],'seconds',time.monotonic()-start,flush=True)
for r in summary:print('C66_GONAL_CLASS',r['class'],r['finite_effective_divisors'],r['global_effective_cuspidal_representatives'],r['types'],flush=True)
print('NONCUSP_TYPE_COUNTS',Counter(k for r in noncusp for k in r['types']),flush=True)
