from pathlib import Path
from collections import Counter
import ast,json
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
N,d,g,ds,keys,blocks,bad,tested,total=read(D/'double_cusp_3_21_D5_dt_arithmetic.txt')
assert (N,d,g)==(21,5,25) and Counter(ds)=={1:24,3:8} and tested==total==100680 and not bad
assert len(keys)==128
# Pack each GF(4) column as two binary bit planes; elimination operates on exact field elements.
mask=(1<<g)-1
cols={}
for (i,m),B in zip(keys,blocks):
 assert len(B)==g and all(len(row)==2*m*ds[i-1] for row in B)
 cc=[]
 for j in range(len(B[0])):
  lo=sum((B[k][j]&1)<<k for k in range(g));hi=sum(((B[k][j]>>1)&1)<<k for k in range(g));cc.append(lo|(hi<<g))
 cols[i,m]=cc
 # Every longer block begins with the independently checked shorter restrictions.
 if m>1:assert cc[:2*(m-1)*ds[i-1]]==cols[i,m-1]
def mul(w,a):
 lo=w&mask;hi=w>>g
 if a==1:return w
 if a==2:return hi|((lo^hi)<<g)
 if a==3:return (lo^hi)|(lo<<g)
 return 0
inv=[0,1,3,2]
def rank(cc):
 piv={}
 for w in cc:
  while w:
   bits=(w&mask)|(w>>g);k=bits.bit_length()-1;a=((w>>k)&1)|(((w>>(g+k))&1)<<1)
   if k in piv:w^=mul(piv[k],a)
   else:piv[k]=mul(w,inv[a]);break
 return len(piv)
def weighted(weights,n):
 if not weights:
  if n==0:yield []
  return
 for m in range(n//weights[0]+1):
  for tail in weighted(weights[1:],n-m*weights[0]):yield [m]+tail
count=0
for v in weighted(ds,5):
 cc=[c for i,m in enumerate(v,1) if m for c in cols[i,m]]
 assert len(cc)==10 and rank(cc)==10;count+=1
assert count==100680
c=read(D/'cont_local_3_21_q4_census.m');m,n,q,cusp,ns,f,full=c
assert (m,n,q)==(3,21,4) and ns==[0,0,0,144,216]
tr=[0]
for k in range(1,11):tr.append(-k*f[2*g-k]-sum(f[2*g-j]*tr[k-j] for j in range(1,k)))
cl=[0]
for k in range(1,6):
 z=4**k+1-tr[2*k]-sum(j*cl[j] for j in range(1,k) if k%j==0);assert z%k==0;cl.append(z//k)
assert cl[1:]==full==[ns[k-1]+cusp.count(k) for k in range(1,6)]
pointkeys,labels,reps,tests,remaining=read(D/'sieve_3_21_D5_arithmetic.txt')
assert len(pointkeys)==len(set(map(tuple,pointkeys)))==216 and Counter(labels)=={1:72,2:72,3:72}
assert tests==[[1,13,0],[2,13,0],[3,13,0]] and not remaining
old=read(D/'degree_independent_321/sieve_3_21_D4_arithmetic.txt');assert len(old[0])==144 and Counter(old[1])=={1:72,2:72} and old[3]==[[1,13,0],[2,13,0]] and not old[4]
tate=json.loads((D/'cusp_tate_replay.json').read_text());tate=[r for r in tate if (r['m'],r['n'])==(3,21)][0];assert tate['local_signature']=={'1':24,'3':8} and tate['all_residue_rational_cusps_globally_rational']
report={'group':[3,21],'absolute_degree':10,'relative_degree':5,'noncuspidal_counts':ns,'cusp_divisors_tested':count,'doubled_cusp_rank':10,'local_matrices':len(keys),'quartic_Hecke_orbits':2,'quintic_Hecke_orbits':3,'all_Hecke_images_nonprincipal':True,'independent_replay':'PASS'}
(D/'M3N21_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC_M3N21_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
