from pathlib import Path
from collections import Counter
from math import gcd
import ast,json,re,subprocess,sys,time
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C63';start=time.time()
def read(n):return ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'))
raw=read('native_replay.m');ind=json.loads((C/'independent_genuine_keys.json').read_text());assert ind['status']=='PASS';total=0
for e,rows in raw:
 for ci,sk,dk,targets,h0s in rows:
  r=next(r for r in ind['source_cases']if r['degree']==e and r['case']==ci);assert sk==r['source_key']and dk==r['diamond5_key']
  expected=Counter({json.dumps(ky):mult for mult,ky in r['target_keys']});actual=Counter()
  for orbit,n,ky in targets:
   assert e*orbit==n*ky[0];actual[json.dumps(ky)]+=orbit
  assert actual==expected and sum(actual.values())==6 and all(x==0 for x in h0s);total+=len(h0s)
assert total==21
# Independent primitive Tate-column cusp arithmetic, including Gal5 orientation.
n=63;labels={(a,b)for a in range(n)for b in range(n)if gcd(gcd(a,b),n)==1}
def orb(a,b):return frozenset((sgn*(a+i*b)%n,sgn*b%n)for sgn in[-1,1]for i in range(n))
classes=sorted(set(orb(a,b)for a,b in labels),key=lambda x:sorted(x));ix={x:i for i,o in enumerate(classes)for x in o};h=[ix[(8*min(o)[0]%n,8*min(o)[1]%n)]for o in classes];assert all(h[h[i]]==i and h[i]!=i for i in range(len(h)))
qs=sorted(set(tuple(sorted((i,h[i])))for i in range(len(h))));qi={i:j for j,o in enumerate(qs)for i in o}
def gp(u):return [qi[ix[(u*min(classes[o[0]])[0]%n,min(classes[o[0]])[1])]]for o in qs]
g2=gp(2);g5=gp(5);unit=[j for j,o in enumerate(qs)if min(classes[o[0]])[1]==0];assert len(unit)==9
left=set(unit);f=[]
while left:
 x=min(left);o=[];y=x
 while y notin o:o.append(y);left.remove(y);y=g2[y]
 assert y==x and len(o)==3;f.append(frozenset(o))
assert len(f)==3;perm=[f.index(frozenset(g5[x]for x in o))for o in f];assert sorted(perm)==[0,1,2]and all(perm[i]!=i and perm[perm[perm[i]]]==i for i in range(3))
# Fresh high kernels each live in their own basis; no jets are mixed.
def rows(M):
 nr,nc,es=M;assert len(es)==nr*nc;return [sum(es[nc*i+j]<<j for j in range(nc))for i in range(nr)]
def rank(rs):
 p={}
 for x in rs:
  while x:
   j=x.bit_length()-1
   if j in p:x^=p[j]
   else:p[j]=x;break
 return len(p)
seen=set()
for part in range(1,5):
 model=read('model_%d.m'%part);cm=rows(model[2]);U=rows(model[3]);assert len(cm)==49 and rank([cm[i]^(1<<i)for i in range(49)])==45 and rank(U)==45
 assert rank(U+[cm[i]^(1<<i)for i in range(49)])==45
 for i,d,m,W,rr,finite,gen in read('high_kernels_%d.m'%part):
  w=rows(W);assert i notin seen and d*m==10 and rr==10 and rank(w)==len(w);seen.add(i)
  assert rank(U+w)-rank(w)==10
assert len(seen)==12
s=(C/'inputs/C63/cartier_kernels_D9.m').read_text();ds=ast.literal_eval(re.findall(r'\bds:=(.*?);',s,re.S)[0]);blocks=[]
for i,d in enumerate(ds,1):
 m=9//d;nr,nc,ss=re.findall(r'localjets\[<'+str(i)+','+str(m)+r'>\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S)[0];nr=int(nr);nc=int(nc);es=ast.literal_eval(ss);blocks.append([sum(es[nc*r+j]<<r for r in range(nr))for j in range(nc)])
count=0
def enum(i,left,rs):
 global count
 if i==len(ds):
  if not left:assert len(rs)==10 and rank(rs)==10;count+=1
  return
 for m in range(min(left//ds[i],len(blocks[i])//ds[i])+1):enum(i+1,left-m*ds[i],rs+blocks[i][:m*ds[i]])
enum(0,10,[]);assert count==239226
r=subprocess.run([sys.executable,str(C/'inputs/C63/replay_rr_kernels.py'),str(C/'inputs/C63/cartier_kernels_D9.m')],text=True,capture_output=True);assert r.returncode==0,r.stdout+r.stderr
(C/'old_kernel_replay.log').write_text(r.stdout+r.stderr)
res={'status':'PASS','native_RR_tests_bound_to_genuine_isogenies':21,'actual_geometric5_kernels':36,'manual_Velu_checks':72,'independent_full_matrix_boundary_patterns':239226,'independent_new_high_kernel_rank_checks':12,'total_pure_cusp_patterns':239238,'all_cusp_ranks':10,'nonglobal_width1_cubic_cusps':3,'Gal5_on_three_closed_cusps':'3-cycle; all six ordered pairs tested','old_RR_kernel_equalities_independently_replayed':121,'elapsed_seconds':time.time()-start}
(C/'independent_final_controls.json').write_text(json.dumps(res,indent=2)+'\n');print('D31A_C63_INDEPENDENT_BIND_BOUNDARY_COMPLETE',21,239238,time.time()-start,flush=True)
