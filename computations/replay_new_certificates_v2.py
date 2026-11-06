from pathlib import Path
import json,re,hashlib,math
from fractions import Fraction
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,matrix,vector,euler_phi,is_prime
H=Path(__file__).resolve().parents[1]
# Replay the rational matrices independently of Magma.
total=0
for p in sorted((H/'exclusions').glob('image121_chi*.json')):
 d=json.loads(p.read_text());R=PolynomialRing(QQ,'t');rat=lambda q:QQ(q[0])/q[1]
 f=R([rat(q) for q in d['field_polynomial']]);assert f.is_irreducible();K=NumberField(f,'a');a=K.gen()
 el=lambda row:sum((rat(q)*a**i for i,q in enumerate(row)),K(0))
 A=matrix(K,d['d'],d['d'],[el(row) for row in d['A']]);z=A.det().norm()
 assert z==d['norm'] and ZZ(z)%2==1 and d['norm_power']==K.degree()//euler_phi(d['chi_order'])
 total+=int(euler_phi(d['chi_order']))*d['d'];print('DECIC_121_MATRIX_NORM_REPLAY_PASS',d['chi_order'],'norm_is_odd',flush=True)
assert total==520
# Every serialized modular Krylov matrix is replayed in its finite field.
ms=list((H/'exclusions').glob('krylov_*.json'))
for p in sorted(ms):
 d=json.loads(p.read_text());assert is_prime(d['p']);F=GF(d['p']);n=d['d'];A=matrix(F,n,n,d['A']);v=vector(F,d['v']);vs=[]
 for _ in range(n):vs.append(v);v=v*A
 rank=matrix(F,vs).rank()
 if rank<n:
  logs='\n'.join(q.read_text() for q in (H/'jobs').glob('winding*.log'))
  pat=r'WINDING_KRYLOV_MATRIX[^\n]*'+re.escape(p.name)+r'\nWINDING_KRYLOV N=\d+ order=\d+ span_rank='+str(rank)+r' dimension='+str(n)+r' verified=false'
  assert re.search(pat,logs),p.name
  print('DECIC_WINDING_KRYLOV_DIAGNOSTIC_REPLAY',p.name,'rank',rank,'dimension',n,'not_used_as_full_rank_certificate',flush=True)
 else: print('DECIC_WINDING_KRYLOV_REPLAY_PASS',p.name,flush=True)
# The dimension sum includes each character and every oldform multiplicity.
D=json.loads((H/'exclusions/global3_additional.json').read_text())
for r in D['rows']:
 n=r['N'];idx=Fraction(n*n,2)
 for p in range(2,n+1):
  if n%p==0 and is_prime(p):idx*=Fraction(p*p-1,p*p)
 phi=lambda a:int(euler_phi(a));g=1+idx/12-Fraction(sum(phi(d)*phi(n//d) for d in range(1,n+1) if n%d==0),4)
 assert idx==r['index'] and g==r['genus'];assert 10*r['kq']<Fraction(r['gonality_lower_bound'])==Fraction(325,32768)*idx
 o=1
 while pow(r['a'],o,n) not in [1,n-1]:o+=1
 assert o==r['diamond_order'] and 10%o and r['minimum_CM_degree']>10
 if 'character_data' in r:assert sum(phi(z[0])*z[2] for z in r['character_data'])==g
 print('DECIC_ADDITIONAL_ARITHMETIC_PASS',n,flush=True)
print('DECIC_NEW_CERTIFICATES_REPLAY_COMPLETE matrices='+str(len(ms))+' image121_dimension='+str(total),flush=True)
