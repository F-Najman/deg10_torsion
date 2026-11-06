from pathlib import Path
import json
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,EllipticCurve
H=Path(__file__).resolve().parents[1];d=json.loads((H/'exclusions/full11_generator.json').read_text());rat=lambda q:QQ(q[0])/q[1]
R=PolynomialRing(QQ,'t');f=R([rat(q) for q in d['field_polynomial']]);assert f==R.cyclotomic_polynomial(11);K=NumberField(f,'z');z=K.gen()
el=lambda row:sum((rat(q)*z**i for i,q in enumerate(row)),K(0))
A=EllipticCurve(K,[0,-1,1,0,0]);Q=A([el(row) for row in d['generator']]);assert 25*Q==A(0) and 5*Q!=A(0)
Ar=EllipticCurve(GF(23),[0,-1,1,0,0]);assert Ar.cardinality()==25 and f.change_ring(GF(23))(GF(23)(2))==0
cache=json.loads((H/'inputs/global_newforms_cache.json').read_text());used=[]
for lev in [1,11,121]:
 for a in cache[str(lev)]:
  assert a['analytic_rank_proved']
  if not a['is_cm']:assert a['analytic_rank']==0;used.append(a['label'])
assert A.j_invariant() not in ZZ
# Non-CM Dirichlet twists at characters of conductor 11 have level 121.
# Together with rank-zero Kato, the nonzero L-values prove rank A(K)=0.
print('DECIC_FULL11_RANK_INPUT_REPLAY_PASS levels=1,11,121 nonCM_forms='+str(len(used)),flush=True)
cusps=0;js=set()
for i in range(25):
 q=i*Q
 if q==A(0) or q[0]==0:cusps+=1;continue
 x=(q[1]+1)/q[0];y=1/q[0];assert x*x*y-x*y*y+y-1==0
 if x==0 or x*y==1:cusps+=1;continue
 r=(x*x*y-x*y+y-1)/(x*(x*y-1));s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1);aa=[1-c,-b,-b,K(0),K(0)]
 a1,a2,a3,a4,a6=aa;b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 if -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6==0:cusps+=1;continue
 E=EllipticCurve(K,aa);P=E([0,0]);assert 11*P==E(0) and P!=E(0)
 red=lambda a:sum((GF(23)(coef)*GF(23)(2)**j for j,coef in enumerate(K(a).list())),GF(23)(0))
 Er=EllipticCurve(GF(23),[red(a) for a in aa]);assert Er.cardinality()<121
 js.add(QQ(E.j_invariant()));print('DECIC_FULL11_CASE_REPLAY_PASS',i,'good_prime=23','order',Er.cardinality(),flush=True)
assert cusps==10 and js=={QQ(-121),QQ(-32768),QQ(-24729001)}
(H/'exclusions/full11_rank_input.json').write_text(json.dumps({'curve':[0,-1,1,0,0],'base_degree':10,'field':'Q(zeta_11)','nonCM_proved_rank_zero_forms':used,'rank_zero_argument':'All twists by characters of Gal(Q(zeta_11)/Q) are non-CM forms of level dividing 121, and all such forms have proved nonzero central L-values. Apply rank-zero Kato to these Dirichlet twists.','noncuspidal_j_invariants':sorted(map(int,js)),'torsion_order':25},indent=2)+'\n')
print('DECIC_CYCLOTOMIC11_REPLAY_COMPLETE points=25 cusps=10 noncusps=15 full11_excluded=true',flush=True)
