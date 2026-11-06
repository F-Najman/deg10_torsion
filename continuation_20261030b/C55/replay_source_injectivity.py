"""Independent integer replay of the rational cusp2-primary saturation."""
from pathlib import Path
import ast,json,hashlib
H=Path(__file__).resolve().parent
p=H/'inputs/N55_quotient_cusp_arithmetic.txt'
a=ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
assert a[0]==[5,5,10,1200,0]and a[3]==[5,10,1200]
gens=[(v[2]%2,v[3]%16)for v in a[2]];S={(0,0)}
while True:
 T=S|{((x+u)%2,(y+v)%16)for x,y in S for u,v in gens}
 if T==S:break
 S=T
assert len(S)==32
iv=[5,550,8972396739917886000]
def v2(n):
 r=0
 while n%2==0:r+=1;n//=2
 return r
assert [2**v2(n)for n in iv if n%2==0]==[2,16]
c=json.loads((H/'cusp_audit.json').read_text());assert c['every_source_cusp_global_orbit_intact_at_2']and c['every_quotient_cusp_global_orbit_intact_at_2']
r=json.loads((H/'source_torsion_injectivity.json').read_text());assert r['source_Q_torsion_invariants']==iv and r['quotient_cusp_projected_2primary_generators']==[list(v)for v in gens]
assert r['explicitly_generated_2primary_elements']==[list(v)for v in sorted(S)]
for row in r['evidence']:
 p=Path(row['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
print('C55_SOURCE_TORSION_SATURATION_INTEGER_REPLAY_PASS 32elements')
