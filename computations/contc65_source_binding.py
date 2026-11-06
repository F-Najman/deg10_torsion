from sage.all import *
from pathlib import Path
import ast,json,sys
sys.path.insert(0,'work/deg10/continuation_20261030c/C65/inputs')
from unit_arithmetic import lit,decode_record
ROOT=Path('work/deg10/continuation_20261030c/C65')
src=ast.literal_eval((ROOT/'source65.m').read_text().replace('[*','[').replace('*]',']'))
records=lit((ROOT/'ordinary10_all.m').read_text(),'SourceRecords')
R=PolynomialRing(GF(2),'T');pol=next(v[1] for v in src[0] if v[0]==10);K=GF(2**10,'a',modulus=R(pol));a=K.gen()
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y
 n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3
 return [-n2**3/b3**2,-(n1*n2-b3)/b3]
def code(x):return [int(x.polynomial()[j]) for j in range(10)]
def key(v):return min(tuple(code(v[0]**(2**j))+code(v[1]**(2**j))) for j in range(10))
allkeys=set();sizes=[]
for i,rec in enumerate(records):
 d,ky=src[2][8+i];assert d==10
 b=R(ky[:10])(a);c=R(ky[10:])(a);E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0)
 data=decode_record(rec);assert data['field_poly']==pol
 cv=lambda x:R([(x>>j)&1 for j in range(10)])(a)
 EE=EllipticCurve(K,list(map(cv,data['ainv'])));PP=EE(list(map(cv,data['point'])))
 assert pointbc(EE,PP)==[b,c] and 65*PP==EE(0) and 5*PP!=EE(0) and 13*PP!=EE(0)
 ks={key(pointbc(E,u*P)) for u in range(1,33) if gcd(u,65)==1}
 assert min(ks)==tuple(ky) and not allkeys.intersection(ks)
 allkeys.update(ks);sizes.append(len(ks))
assert sizes==[24,24,24,12,24] and len(allkeys)==108
(ROOT/'independent_source_binding.json').write_text(json.dumps({'status':'PASS','normalized_degree10_sources':5,'closed_places':108,'diamond_orbit_sizes':sizes,'all_original_Tate_bc_keys_and_canonical_level_points_match':True},indent=2)+'\n')
print('C65_DECIC_PHYSICAL_SOURCE_BINDING_PASS',sizes,108,flush=True)
