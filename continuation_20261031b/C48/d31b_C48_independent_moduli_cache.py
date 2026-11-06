from sage.all import *
from pathlib import Path
import json,time
P=Path('work/deg10/continuation_20261031b/C48');T=time.monotonic();Fp=GF(5);R=PolynomialRing(Fp,'t')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
lo=rd('low_geometry.m');hi=rd('high_character_merged.json');raw=rd('new_high_sources.m');labels=rd('geometric_action_labels.m');fs={d:GF(5**d,'a',modulus=R(pol))for d,pol in lo[0]}
def decode(d,ky):K=fs[d];return R(ky[:d])(K.gen()),R(ky[d:])(K.gen())
def key(v,d):
 b,c=v;out=[]
 for j in range(d):out.append(tuple(int(z.polynomial()[i])for z in[b,c]for i in range(d)));b=b**5;c=c**5
 return min(out)
def desc(E,S,Q):
 assert S[1]==0 and E.a1()==0 and E.a3()==0
 z=Q[0]-S[0];return[z**3/Q[1]**2,(3*S[0]+E.a2())*z**2/Q[1]**2]
def project(d,ky):
 b,c=decode(d,ky);K=fs[d];E=EllipticCurve(K,[1-c,-b,-b,0,0]);P0=E(0,0);S0=24*P0;assert S0!=E(0)and 2*S0==E(0)and 16*P0!=E(0);Z=P0+S0;r=S0[0];a1=1-c;A2=-b+a1*a1/4;A4=-a1*b/2;A6=b*b/4;v=3*r*r+2*A2*r+A4;ET=EllipticCurve(K,[0,A2,0,A4-5*v,A6-(4*A2+7*r)*v]);xx=P0[0]+Z[0]-S0[0];yy=P0[1]+Z[1]-S0[1]+(a1*xx-b)/2;Q=ET(xx,yy);S=ET(-A2-2*r,0);assert 24*Q==ET(0)and 12*Q!=ET(0)and 8*Q!=ET(0)and 2*S==ET(0)and S!=12*Q;return ET,S,Q
for ii,r in enumerate(lo[1],1):
 rid,d,ky,rep,ti,hh=r;E,S,Q=project(d,ky);eps,eta,u=labels[hh-1];vv=desc(E,S+eps*(12*Q),u*Q+eta*S);dt,kt,e,n=lo[2][ti-1];assert dt==d and e*n==d and key(vv,d)==tuple(kt);assert lcm([z.minpoly().degree()for z in vv])==e
 if ii%2000==0:print('D31B_C48_INDEPENDENT_LOW_CACHE_PROGRESS',ii,time.monotonic()-T,flush=True)
print('D31B_C48_INDEPENDENT_LOW_CACHE_PASS',len(lo[1]),time.monotonic()-T,flush=True)
# Independently expand each cached high target's complete rational automorphism orbit.
orbits={}
for ti,(d,ky,e,n)in enumerate(hi[2],1):
 b,c=decode(d,ky);K=fs[d];E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);S=E(0,0);Q=E(b,b);Q12=12*Q;QQ={u:u*Q for u in[1,5,7,11]}
 for eps,eta,u in labels:
  kyh=(d,key(desc(E,S+eps*Q12,QQ[u]+eta*S),d))
  if kyh in orbits:assert orbits[kyh]==ti
  else:orbits[kyh]=ti
 if ti%1000==0:print('D31B_C48_INDEPENDENT_HIGH_ORBITS_PROGRESS',ti,len(orbits),time.monotonic()-T,flush=True)
source={r[0]:r for r in raw[1]};kept=0
for r in hi[1]:
 if not r[5]:continue
 id,d,e,n,chi,keep,ti=r;_,dd,ky=source[id];assert dd==d;E,S,Q=project(d,ky);v=desc(E,S,Q);assert lcm([z.minpoly().degree()for z in v])==e and e*n==d;assert orbits[(d,key(v,d))]==ti;kept+=1
 if kept%2000==0:print('D31B_C48_INDEPENDENT_HIGH_CACHE_PROGRESS',kept,time.monotonic()-T,flush=True)
assert kept==12604
out={'status':'PASS','low_actual_source_points':len(lo[1]),'low_target_queries':len(lo[2]),'all_mixed_normalizing_automorphism_labels_bound_to_actual_projected_points':True,'high_retained_source_points':kept,'high_target_queries':len(hi[2]),'high_expanded_target_orbit_keys':len(orbits),'all_high_cluster_bindings_match_actual_projected_points':True,'method':'Actual source24P and manual two-point Velu formulas in normalized short coordinates; actual dual-kernel and marked24-level; exact rational automorphism action and full Frobenius orbit, independently computed in Sage','elapsed_seconds':time.monotonic()-T};(P/'independent_moduli_cache.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C48_INDEPENDENT_ALL_MODULI_CACHE_COMPLETE',len(lo[1]),kept,len(orbits),time.monotonic()-T,flush=True)
