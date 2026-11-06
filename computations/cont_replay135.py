from pathlib import Path
import json,re,math
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,matrix,vector,euler_phi,is_prime,divisors
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029';total=0;rows=[]
for p in sorted(D.glob('image135_a16_q31_*.json')):
 d=json.loads(p.read_text());rat=lambda r:QQ(r[0])/r[1];R=PolynomialRing(QQ,'t');f=R([rat(c) for c in d['field_polynomial']]);assert f.is_irreducible();K=NumberField(f,'a');a=K.gen()
 el=lambda cs:sum((rat(c)*a**i for i,c in enumerate(cs)),K(0))
 A=matrix(K,d['d'],d['d'],[el(cs) for cs in d['A']]);z=A.det().norm();assert z==d['norm'] and ZZ(z)%2==1
 assert d['norm_power']==K.degree()//euler_phi(d['chi_order'])
 total+=int(euler_phi(d['chi_order']))*d['d'];rows.append([d['chi_order'],d['d'],str(z)])
 print('DECIC135_MATRIX_REPLAY_PASS',d['chi_order'],flush=True)
assert total==480 and len(rows)==6
p=D/'krylov_N135_chi4_order36.json';d=json.loads(p.read_text());assert is_prime(d['p']);F=GF(d['p']);n=d['d'];A=matrix(F,n,n,d['A']);v=vector(F,d['v']);vs=[]
for _ in range(n):vs.append(v);v=v*A
assert matrix(F,vs).rank()==n==16
s=(H/'jobs/cont_winding135.log').read_text();cs=re.findall(r'WINDING N=135 order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',s)
assert sum(int(euler_phi(int(o)))*int(d) for o,c,d,w,k in cs)==529
assert [(int(o),int(c)) for o,c,d,w,k in cs if d!=w]==[(1,1)]
N=135;idx=QQ(N*N)/2
for p in ZZ(N).prime_divisors():idx*=1-QQ(1)/(p*p)
assert idx==7776 and 60<325*idx/32768==QQ(78975)/1024
assert next(o for o in range(1,100) if pow(16,o,135) in (1,134))==9
assert math.gcd(9,10)==1
report={'N':135,'diamond':16,'diamond_order':9,'image_dimension':480,'matrix_rows':rows,'genus':529,'winding_dimension_sum':529,'only_nonfull_winding_component':'trivial character','k2':6,'degree_times_k2':60,'gonality_lower_bound':'78975/1024','minimum_CM_degree':72,'independent_matrix_and_arithmetic_replay':'PASS'}
(D/'C135_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('DECIC135_INDEPENDENT_REPLAY_COMPLETE',flush=True)
