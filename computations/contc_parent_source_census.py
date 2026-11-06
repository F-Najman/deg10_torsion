"""Independent exact modular point/cusp census for the new decic frontier."""
from pathlib import Path
from math import gcd
from collections import Counter
from fractions import Fraction
import ast,json,hashlib,time
from sage.all import GF,PolynomialRing,EllipticCurve,divisors,moebius,euler_phi
H=Path(__file__).resolve().parents[1]; OUT=H/'continuation_20261030c'; t0=time.time()
raw=H/'next_cases_20261030/dyadic_local/all_curves.m'
base=json.loads((raw.parent/'independent_point_counts.json').read_text()); assert base['status']=='PASS' and base['all_finite_j_values_and_twists_present']
assert hashlib.sha256(raw.read_bytes()).hexdigest()==base['input_sha256']
data=ast.literal_eval(raw.read_text().replace('[*','[').replace('*]',']'));levels=[63,65,85,91,93,99,105]
counts={n:[]for n in levels};closed={n:[]for n in levels};rows={n:[]for n in levels};R=PolynomialRing(GF(2),'z')
for d,pol,curves in data:
 F=GF(2**d,'a',modulus=R(pol));z=F.gen();fv=lambda v:R(v)(z)
 totals={n:Fraction(0)for n in levels};types={n:Counter()for n in levels};eligible={n:[]for n in levels}
 for jj,aa,card in curves:
  ls=[n for n in levels if card%n==0]
  if not ls:continue
  E=EllipticCurve(F,[fv(v)for v in aa]);assert E.cardinality()==card
  inv=list(map(int,E.abelian_group().invariants()));assert len(inv)<=2
  inv=[1]*(2-len(inv))+inv;a,b=inv;assert a*b==card and b%a==0
  aut=len(E.automorphisms()) if not any(jj) else 2
  for n in ls:
   exact=sum(int(moebius(n//k))*gcd(a,k)*gcd(b,k)for k in divisors(n))
   if not exact:continue
   totals[n]+=Fraction(exact,aut);types[n][(card,a,b,aut,exact)]+=1
   j=fv(jj);jdeg=next(k for k in divisors(d)if j**(2**k)==j)
   eligible[n].append({'j':jj,'a_invariants':aa,'cardinality':card,'group':[a,b],'automorphisms':aut,'points_exact_order':exact,'j_field_degree':int(jdeg)})
 for n in levels:
  assert totals[n].denominator==1,(n,d,totals[n]);p=int(totals[n]);counts[n].append(p)
  v=p-sum(e*closed[n][e-1]for e in range(1,d)if d%e==0);assert v%d==0;(closed[n]).append(v//d)
  rows[n].append({'degree':d,'curve_type_counts':[[list(k),v]for k,v in sorted(types[n].items())],'eligible_curves':eligible[n],'rational_nc_points':p,'nc_closed_points':v//d,'field_polynomial':pol})
 print('CONT_C_PARENT_CENSUS_DEGREE',d,{n:closed[n][-1]for n in levels},flush=True)

def orbits(perms,n):
 unseen=set(range(n));ans=[]
 while unseen:
  x=min(unseen);todo=[x];seen={x}
  while todo:
   a=todo.pop()
   for p in perms:
    y=p[a]
    if y not in seen:seen.add(y);todo.append(y)
  unseen-=seen;ans.append(sorted(seen))
 return ans

def cusp_data(n):
 labels={(a,b)for a in range(n)for b in range(n)if gcd(gcd(a,b),n)==1};unseen=set(labels);cs=[]
 while unseen:
  a,b=min(unseen);o={(sgn*(a+k*b)%n,sgn*b%n)for k in range(n)for sgn in [-1,1]};assert o<=labels;unseen-=o;cs.append(sorted(o))
 ix={x:i for i,c in enumerate(cs)for x in c};assert len(ix)==len(labels)
 def perm(u):
  ans=[]
  for c in cs:
   ims={ix[(u*a%n,b)]for a,b in c};assert len(ims)==1;ans.append(ims.pop())
  assert sorted(ans)==list(range(len(cs)));return ans
 glob=orbits([perm(u)for u in range(1,n)if gcd(u,n)==1],len(cs));fin=orbits([perm(2)],len(cs))
 return {'global_signature':dict(Counter(map(len,glob))),'finite_signature':dict(Counter(map(len,fin))),'all_orbits_intact':{tuple(o)for o in glob}=={tuple(o)for o in fin},'geometric_cusps':len(cs),'global_orbits':glob,'finite_orbits':fin,'cusp_representatives':[c[0]for c in cs]}
res={'status':'PASS','method':'Independent Sage exact finite elliptic groups and point-order formula, weighted by rational automorphisms; primitive Tate-column cusp orbit arithmetic','input_sha256':base['input_sha256'],'degrees':list(range(1,11)),'cases':{},'seconds':time.time()-t0}
for n in levels:
 idx=Fraction(n*n,2)
 for p in __import__('sage.all',fromlist=['prime_divisors']).prime_divisors(n):idx*=Fraction(int(p*p-1),int(p*p))
 assert idx.denominator==1
 gbound=Fraction(325,32768)*idx
 res['cases'][str(n)]={'level':n,'nc_closed_counts':closed[n],'nc_rational_counts':counts[n],'fields':rows[n],'cusps':cusp_data(n),'projective_index':int(idx),'gonality_lower_bound':[gbound.numerator,gbound.denominator]}
(OUT/'parent_source_census.json').write_text(json.dumps(res,indent=2)+'\n');print('CONT_C_PARENT_SOURCE_CENSUS_COMPLETE', {n:closed[n]for n in levels},flush=True)
