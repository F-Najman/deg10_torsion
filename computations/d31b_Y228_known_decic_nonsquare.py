from sage.all import *
from pathlib import Path
import json,time
C=Path(__file__).resolve().parents[1]/'continuation_20261031b/Y228_shared';found={};start=time.monotonic()
for p in prime_range(5,402):
 if p==7:continue
 F=GF(p);R=PolynomialRing(F,'t');f=R([2,2,-1,-2,-1,1])
 for av in f.roots(multiplicities=False):
  if f.derivative()(av)==0:continue
  b=R([1032,-239,-222,-759,419])(av);c=R([42,-11,-8,-31,17])(av);E=EllipticCurve(F,[1-c,-b,-b,0,0]);
  if E.discriminant()==0:continue
  P=E(0,0)
  if P.order()!=28:continue
  T=14*P;a1,a2,a3,a4,a6=E.a_invariants();S=PolynomialRing(F,'X');X=S.gen();cub=4*X**3+(a1*a1+4*a2)*X**2+2*(a1*a3+2*a4)*X+(a3*a3+4*a6);h=cub//(X-T[0]);assert h.degree()==2
  for rv in h.roots(multiplicities=False):
   if h.derivative()(rv)==0:continue
   Q2=E(rv,-(a1*rv+a3)/2);assert Q2.order()==2 and Q2!=T
   for bb in [0,1]:
    PP=Q2+bb*T
    for cc in [0,1]:
     key=(bb,cc)
     if key in found:continue
     QQ=P+cc*PP;delta=QQ[0]-PP[0];assert delta!=0
     if not delta.is_square():found[key]={'p':int(p),'quintic_root':int(av),'quadratic_root':int(rv),'P2_shift':bb,'Q_shift':cc,'x_difference':int(delta),'Legendre_symbol':int(delta**((p-1)//2))}
 if len(found)==4:break
assert len(found)==4
out={'status':'PASS','criterion':'A point lies in the dual rational2-isogeny image only if x(Q)-x(P2) is a square; normalized Derickx B differs by a square','marking_cosets':4,'independent_prime_nonsquare_witnesses':[found[v]for v in sorted(found)],'elapsed_seconds':time.monotonic()-start};(C/'known_decic_nonsquare_witnesses.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_Y228_KNOWN_DECIC_NONSQUARE_COMPLETE',len(found),out['independent_prime_nonsquare_witnesses'],flush=True)
