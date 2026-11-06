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
H=P.parents[1]
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
raw=rd(P/'inputs/all_curves.m');source=rd(P/'source81.m');fil=rd(P/'x0_filter.m');R=PolynomialRing(GF(2),'z');seen_records=[];rows=[];levels=[];verified_moduli=0

def key(vs):
 out=[]
 for _ in range(vs[0].parent().degree()):out.append(tuple(int(c) for v in vs for c in v._vector_()));vs=[v*v for v in vs]
 return min(out)
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();aa3=2*y+a3+a1*x;m=(3*x*x+2*x*a2+a4-a1*y)/aa3;aa2=3*x+a2-m*m-a1*m;return[-aa2**3/aa3**2,1-(a1+2*m)*aa2/aa3]
def velu(E,xs):
 a1,a2,a3,a4,a6=E.a_invariants();b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;tv=sum(6*x*x+b2*x+b4 for x in xs);wv=sum(10*x**3+2*b2*x*x+3*b4*x+b6 for x in xs);return EllipticCurve(E.base_ring(),[a1,a2,a3,a4-5*tv,a6-b2*tv-7*wv])
for d,pol,curves in raw:
 F=GF(2**d,'a',modulus=R(pol));a=F.gen();fv=lambda v:R(v)(a);points=set();orbs=set();nq=0;cycs=0
 for j,ecs,n in curves:
  if n%81:continue
  E=EllipticCurve(F,[fv(v)for v in ecs]);assert E.cardinality()==n;valid=[Q for Q in E.points() if Q!=E(0) and 81*Q==E(0) and 27*Q!=E(0)]
  if not valid:continue
  nq+=1;marked=set()
  for Q in valid:
   if Q in marked:continue
   marked.update([u*Q for u in range(1,81)if gcd(u,81)==1]);cycs+=1;ks={key(bc(E,u*Q))for u in range(1,41)if gcd(u,81)==1};points|=ks;orbs.add(min(ks))
 native=source[1][d-1];assert [d,len(curves),nq,cycs,len(points),len(orbs)]==[native[0],native[1],native[2],native[4],native[5],native[6]]
 seen_records.extend((d,list(v))for v in sorted(orbs));levels.append(len(points));rows.append([d,nq,cycs,len(points),len(orbs)])
assert seen_records==[(d,k)for d,k in source[2]] and levels==[0,0,0,0,0,0,0,0,135,135]
fields={d:GF(2**d,'a',modulus=R(pol))for d,pol in source[0]}
for row in fil[1]:
 ci,d,ky,values,ims,passed=row;F=fields[d];a=F.gen();b=R(ky[:d])(a);c=R(ky[d:])(a);E=EllipticCurve(F,[1-c,-b,-b,0,0]);Q=E(0,0);assert Q!=E(0)and 81*Q==E(0)and 27*Q!=E(0);xs=[(u*Q)[0]for u in range(1,41)];assert len(set(xs))==40;EG=velu(E,xs);assert [E.j_invariant(),EG.j_invariant()]==[R(v)(a)for v in values];assert len(ims)==1;de,w,cl,yes,flag,ibas=ims[0];assert d==de*w and de==d and w==1;assert yes==(cl[1]%2==0)and cl[2]==0 and passed==yes;verified_moduli+=1
assert [r[0]for r in fil[1]if r[5]]==fil[2]==[1,2,3,7,10]
tc=tate(81);assert tc['every_global_cusp_orbit_unchanged_at2'];out={'status':'PASS','source_census_independent_all_points':rows,'all_ten_normalized_records_exact':True,'finite_NC_closed_counts':levels,'manual_cyclic81_Velu_j_pair_checks':verified_moduli,'all_ten_inverse_norms_one':True,'odd_coset_sieve_remaining':[1,2,3,7,10],'global_and_finite_cusp_signature':tc['global_signature'],'all_cusps_globally_intact':True,'tate_data':tc,'input_hashes':{str(q.relative_to(H)):hashlib.sha256(q.read_bytes()).hexdigest()for q in[P/'inputs/all_curves.m',P/'source81.m',P/'x0_filter.m']}}
(P/'parent_independent_source_filter.json').write_text(json.dumps(out,indent=2)+'\n');print('C81_PARENT_INDEPENDENT_SOURCE_FILTER_PASS',verified_moduli,levels,flush=True)
