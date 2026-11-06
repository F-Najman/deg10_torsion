from pathlib import Path
from math import gcd,isqrt
from collections import Counter
from fractions import Fraction
import re,json,hashlib,time
P=Path(__file__).resolve().parent

def phi(n):return sum(gcd(a,n)==1 for a in range(1,n+1))
def projective_order(a,n):
 u=1
 for k in range(1,phi(n)+1):
  u=u*a%n
  if u in (1,n-1):return k
 raise AssertionError('not a unit')
def orbits(perms,n):
 unseen=set(range(n));out=[]
 while unseen:
  x=min(unseen);todo=[x];seen={x}
  while todo:
   y=todo.pop()
   for perm in perms:
    z=perm[y]
    if z not in seen:seen.add(z);todo.append(z)
  unseen-=seen;out.append(sorted(seen))
 return out
def rank_mod(A,p):
 A=[list(map(lambda x:x%p,row)) for row in A];m=len(A);n=len(A[0]);r=0
 for j in range(n):
  pivot=next((i for i in range(r,m) if A[i][j]),None)
  if pivot is None:continue
  A[r],A[pivot]=A[pivot],A[r];u=pow(A[r][j],-1,p);A[r]=[x*u%p for x in A[r]]
  for i in range(m):
   if i!=r and A[i][j]:
    u=A[i][j];A[i]=[(x-u*y)%p for x,y in zip(A[i],A[r])]
  r+=1
  if r==m:break
 return r

def tate(n):
 labels={(a,b) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1}
 unseen=set(labels);cs=[]
 while unseen:
  a,b=min(unseen);o={(s*(a+k*b)%n,s*b%n) for k in range(n) for s in (-1,1)}
  assert o<=labels;unseen-=o;cs.append(sorted(o))
 ix={x:i for i,C in enumerate(cs) for x in C};assert len(ix)==len(labels)
 perms=[]
 for u in range(1,n):
  if gcd(u,n)!=1:continue
  perm=[]
  for C in cs:
   ims={ix[(u*a%n,b)] for a,b in C};assert len(ims)==1;perm.append(ims.pop())
  assert sorted(perm)==list(range(len(cs)));perms.append(perm)
 frob=[ix[(2*C[0][0]%n,C[0][1])] for C in cs]
 global_orbits=orbits(perms,len(cs));finite_orbits=orbits([frob],len(cs))
 assert {tuple(x) for x in global_orbits}=={tuple(x) for x in finite_orbits}
 return {'geometric_cusps':len(cs),'tate_label_count':len(labels),'global_signature':dict(Counter(map(len,global_orbits))),'finite_signature':dict(Counter(map(len,finite_orbits))),'every_global_cusp_orbit_unchanged_at2':True,'global_orbits':global_orbits,'frobenius2_permutation':frob,'cusp_representatives':[C[0] for C in cs]}

def waterhouse(n):
 rows=[]
 for f in range(1,11):
  q=2**f;b=isqrt(4*q);ss={0}
  if f%2:ss|={2**((f+1)//2),-2**((f+1)//2)}
  else:ss|={2**(f//2),-2**(f//2),2**(f//2+1),-2**(f//2+1)}
  cand=[v for v in range(q+1-b,q+2+b) if v%n==0]
  accepted=[v for v in cand if (q+1-v)%2 or q+1-v in ss]
  assert not accepted,(n,f,accepted)
  rows.append({'f':f,'q':q,'hasse_interval':[q+1-b,q+1+b],'candidate_orders':cand,'candidate_traces':[q+1-v for v in cand],'supersingular_traces':sorted(t for t in ss if abs(t)<=b),'eligible_orders':accepted})
 return rows

def winding(n,a):
 text=(P/'inputs'/('winding%d.log'%n)).read_text()
 pat=r'WINDING N=(\d+) order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)'
 rows=[tuple(map(int,m)) for m in re.findall(pat,text)];assert all(r[0]==n for r in rows)
 order=projective_order(2,n);assert order==36
 exp=next(e for e in range(36) if pow(2,e,n) in (a,n-a))
 assert sorted(r[1] for r in rows)==[1,2,3,4,6,9,12,18,36]
 assert all(r[5]==phi(r[1]) for r in rows)
 genus=sum(r[3]*r[5] for r in rows);active=[r for r in rows if exp%r[1]];assert {r[1] for r in active}=={9,18,36};assert all(r[3]==r[4] for r in active)
 kp=next((P/'inputs').glob('krylov_N%d_chi*_order36.json'%n));k=json.loads(kp.read_text());d=k['d'];p=k['p'];A=[k['A'][d*i:d*(i+1)] for i in range(d)];v=k['v'];K=[]
 for _ in range(d):K.append(v);v=[sum(v[i]*A[i][j] for i in range(d))%p for j in range(d)]
 assert rank_mod(K,p)==d and d==next(r[3] for r in rows if r[1]==36)
 return {'full_genus_dimension':genus,'all_character_orders':[r[1] for r in rows],'active_character_orders':[r[1] for r in active],'active_relative_dimensions':[r[3] for r in active],'finite_image_abelian_dimension':sum(r[3]*r[5] for r in active),'diamond_generator_order':36,'diamond_exponent':exp,'diamond_order':projective_order(a,n),'krylov_prime':p,'krylov_dimension':d,'krylov_rank':d,'complete_exact_winding_records_reused':True}

out={}
for n,a,expected_sig in [(111,10,{1:36,2:18,36:2}),(95,11,{1:36,4:9,18:2,36:1})]:
 tc=tate(n);assert tc['global_signature']==expected_sig;local=waterhouse(n);w=winding(n,a)
 mu=n*n
 for p in [p for p in range(2,n+1) if n%p==0 and all(p%d for d in range(2,isqrt(p)+1))]:mu=mu*(p*p-1)//(p*p)
 assert mu%2==0;mu//=2;bound=Fraction(325*mu,32768);assert bound>40 and 10%w['diamond_order']
 expected_genus=1+mu//12-tc['geometric_cusps']//2;assert w['full_genus_dimension']==expected_genus
 out[str(n)]={'complete':True,'degree':10,'diamond':a,'PSL2_index':mu,'gonality_lower_bound':[bound.numerator,bound.denominator],'comparison_degree':40,'no_nonCM_exception_needed_because_all_reduction_is_cuspidal':True,'tate':tc,'waterhouse':local,'winding':w}
summary={'complete':True,'degrees_reproved':[10],'levels':[111,95],'no_new_magma_jobs_started':True,'results':out}
(P/'independent_replay.json').write_text(json.dumps(summary,indent=2)+'\n')
message='C111_C95_DYADIC_TATE_WATERHOUSE_WINDING_KRYLOV_REPLAY_PASS '+json.dumps({n:{'cusps':v['tate']['global_signature'],'diamond_order':v['winding']['diamond_order'],'genus':v['winding']['full_genus_dimension'],'finite_image_dimension':v['winding']['finite_image_abelian_dimension'],'bound':v['gonality_lower_bound']} for n,v in out.items()})
(P/'replay.log').write_text(message+'\n');print(message)
