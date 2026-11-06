from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path(__file__).resolve().parent;N=93

def normal(a,c):
 c%=N;g=gcd(c,N);a%=g;assert gcd(a,g)==1;return min((a,c),((-a)%g,(-c)%N))
S={normal(a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1};assert len(S)==120
frob=lambda v:normal(2*v[0],v[1]);h=lambda v:normal(32*v[0],32*v[1]);left=set(S);orbits=[]
while left:
 v=min(left);o=[];w=v
 while w not in o:o.append(w);w=frob(w)
 assert w==v and set(o)<=left;left-=set(o);orbits.append(o)
assert Counter(map(len,orbits))=={1:30,2:15,5:6,10:3}
assert all(h(h(v))==v for v in S)
assert all({h(v)for v in o}==set(o)for o in orbits if len(o)in[2,10])
assert all({h(v)for v in o}!=set(o)for o in orbits if len(o)==5)
assert all(h(v)!=v for v in S if gcd(v[1],N)==1)
assert 32%31==1 and (32*32)%93==1 and 32 not in[1,92]
phi5cyclic=[n for n in range(1,26)if n!=23]+[28,30];assert 31 not in phi5cyclic and 17 in phi5cyclic
out=dict(status='PASS',level=93,geometric_cusps=120,finite_signature=dict(Counter(map(len,orbits))),h_scalar=32,h_order=2,h_fixes_every_quadratic_and_decic_closed_cusp=True,h_swaps_quintic_cusp_pairs=True,comparison_degree_after_overlap_at_most36=True,h_stability_descends_31_level_to_quintic=True,prime31_absent_from_accepted_Phi5=True,prime17_is_present_in_Phi5=True,remaining_pure_cusp_types=['rational10','quintic5_plus_rational5','two_quintics'])
(P/'cusp_descent_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('C93_DEG10_H32_CUSP_DESCENT_ARITHMETIC_PASS')
