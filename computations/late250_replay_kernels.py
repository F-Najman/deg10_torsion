from pathlib import Path
import ast,json,time
from sage.all import GF,PolynomialRing,EllipticCurve,prod
P=Path('work/deg10/continuation_20261029b');tt=time.monotonic()
helper=Path('work/deg10/computations/late_replay_filters_v2.py').read_text().split('# The two low-genus filters')[0];exec(helper)
reports=[]
for desc,d,fp,ky,rows in read(P/'C250/kernel_certificates.m'):
 K=field(d,fp);E,T,Q=src(d,ky,K);RK=PolynomialRing(K,'x');norms=[];targets=[];orb=0
 for e,fl,kg,gg,aa,pp,ng in rows:
  m=len(fl)-1;L=GF(3**m,'w',modulus=R(fl));w=L.gen();dec=lambda a:sum(L(v)*w**j for j,v in enumerate(a));alpha=dec(kg);emb=K.hom([alpha],L);assert len(fp)-1==d and R(fp)(alpha)==0
  EL=EllipticCurve(L,[emb(a) for a in E.ainvs()]);U=EL(emb(Q[0]),emb(Q[1]));RL=PolynomialRing(L,'X');g=RL([dec(v) for v in gg]);assert g.degree()==3 and g.is_monic()
  phi=EL.isogeny(g,check=True);assert phi.degree()==7
  EG=EllipticCurve(L,[dec(v) for v in aa]);PG=EG(dec(pp[0]),dec(pp[1]));assert 50*PG==EG(0) and 25*PG!=EG(0) and 10*PG!=EG(0)
  assert EG.j_invariant()==phi.codomain().j_invariant() and tate(EG,PG)==tate(phi.codomain(),phi(U))
  norm=RL.one();h=g
  for i in range(e):norm*=h;h=RL([a**(3**d) for a in h.list()])
  assert h==g
  nk=RK([sum(K(v)*K.gen()**j for j,v in enumerate(z)) for z in ng]);assert RL([emb(a) for a in nk.list()])==norm
  norms.append(nk);orb+=e;td=deg(tate(EG,PG),m);assert d*e%td==0;targets.append([td,d*e//td])
 assert orb==8 and prod(norms).monic()==RK(E.division_polynomial(7)).monic()
 reports.append({'source':desc,'all_eight_kernels_covered_once':True,'every_normalized_isogeny_and_level_image_verified_in_Sage':True,'normed_target_degrees_and_multiplicities':targets,'positive_Hecke_divisor_degree':sum(a*b for a,b in targets)})
 assert reports[-1]['positive_Hecke_divisor_degree']==64
(P/'C250/independent_kernel_replay.json').write_text(json.dumps({'status':'PASS','records':reports,'elapsed_seconds':time.monotonic()-tt},indent=2)+'\n')
print('ALL_C250_KERNELS_AND_POINT_IMAGES_INDEPENDENTLY_VERIFIED',time.monotonic()-tt,flush=True)
