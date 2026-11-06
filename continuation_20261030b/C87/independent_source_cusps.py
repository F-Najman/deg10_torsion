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
 pass
 return {'geometric_cusps':len(cs),'tate_label_count':len(labels),'global_signature':dict(Counter(map(len,global_orbits))),'finite_signature':dict(Counter(map(len,finite_orbits))),'every_global_cusp_orbit_unchanged_at2':{tuple(x) for x in global_orbits}=={tuple(x) for x in finite_orbits},'global_orbits':global_orbits,'frobenius2_permutation':frob,'cusp_representatives':[C[0] for C in cs]}

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

import ast
from sage.all import GF,PolynomialRing,EllipticCurve
H=P.parents[1];raw=H/'next_cases_20261030/dyadic_local/all_curves.m';check=json.loads((H/'next_cases_20261030/dyadic_local/independent_point_counts.json').read_text());assert check['status']=='PASS' and hashlib.sha256(raw.read_bytes()).hexdigest()==check['input_sha256'];data=ast.literal_eval(raw.read_text().replace('[*','[').replace('*]',']'));rows=[];closed=[];counts={}
for d,pol,curves in data:
 R=PolynomialRing(GF(2),'z');F=GF(2**d,'a',modulus=R(pol));z=F.gen();fv=lambda v:R(v)(z);points=0;types=Counter();nc=0
 for j,aa,n in curves:
  if n%87:continue
  E=EllipticCurve(F,[fv(v)for v in aa]);assert E.cardinality()==n and any(j) and n%29==0 and n%(29*29)!=0;nc+=1
  pp=E.division_polynomial(3);pts=[p for xx in pp.roots(multiplicities=False)for p in E.lift_x(xx,all=True)];assert all(p!=E(0) and 3*p==E(0) for p in pts);n3=len(pts)+1;assert n3 in[3,9];types[(n,n3)]+=1;points+=14*(n3-1)
 v=points-sum(e*closed[e-1]for e in range(1,d)if d%e==0);assert v%d==0;closed.append(v//d);rows.append([d,nc,points,v//d,[[list(k),v]for k,v in types.items()]])
assert closed==[0,0,0,0,0,0,0,0,112,224]
tc=tate(87);assert tc['global_signature']=={1:28,2:14,28:2} and tc['finite_signature']=={1:28,2:14,28:2};assert tc['every_global_cusp_orbit_unchanged_at2'];bound=Fraction(7*3360,800);assert bound==Fraction(147,5) and bound>20
res={'status':'PASS','level':87,'genus':225,'source_closed_counts':closed,'independent_division3_counts':rows,'global_and_finite_cusp_signature':tc['global_signature'],'all_source_cusps_globally_intact':True,'gonality_lower_bound':[bound.numerator,bound.denominator],'comparison_degree':20,'source_data_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'tate_data':tc}
(P/'independent_source_cusps.json').write_text(json.dumps(res,indent=2)+'\n');print('B87_INDEPENDENT_SOURCE_CUSPS_PASS',closed,flush=True)
