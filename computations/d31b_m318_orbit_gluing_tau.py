from sage.all import *
from pathlib import Path
import ast,json,re,itertools,time,hashlib
START=time.monotonic();N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing')
exec(compile((N/'framed_prefix.snapshot.py').read_text(),str(N/'framed_prefix.snapshot.py'),'exec'))
N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');START=time.monotonic()
AUT=rd('cusp_automorphisms.m');GENPAR=[[1,1,0,1],[1,0,1,1],[2,0,0,5]];GENPER=AUT[2];BASEPER=AUT[3];assert len(BASEPER)==27
TAUPER=rd('tau_cusp_permutation_p7.m')[1];assert TAUPER==rd('tau_cusp_permutation_p13.m')[1]
GENPAR=GENPAR+[[0,0,0,0]];GENPER=GENPER+[TAUPER]
closure={tuple(range(1,25))};queue=list(closure)
while queue:
 pp=queue.pop(0)
 for gg in GENPER:
  nn=tuple(gg[i-1]for i in pp)
  if nn not in closure:closure.add(nn);queue.append(nn)
assert len(closure)==54
ALLPER=BASEPER+[list(v)for v in sorted(closure-set(map(tuple,BASEPER)))];assert len(ALLPER)==54

def normpair(E,P,Q):
 a1,a2,a3,a4,a6=E.a_invariants();r,tt=P.xy();sl=(3*r*r+2*a2*r+a4-a1*tt)/(2*tt+a1*r+a3);a1p=a1+2*sl;a2p=a2-sl*a1+3*r-sl*sl;a3p=a3+r*a1+2*tt;a4p=a4-sl*a3+2*r*a2-(tt+r*sl)*a1+3*r*r-2*sl*tt;a6p=a6+r*a4+r*r*a2+r**3-tt*a3-r*tt*a1-tt*tt
 assert a2p==a4p==a6p==0
 xx=Q[0]-r;yy=Q[1]-tt-sl*xx;scale=yy/xx**2;t=-scale**2*xx;aa1=scale*a1p;aa3=scale**3*a3p;q=(aa1-t-2)/(t-1);assert aa3==q*t*t+(1-q)*t
 return[q,t]

def fieldkey(qt,p,d):return min(tuple(int((z**(p**j)).polynomial()[i])for z in qt for i in range(d))for j in range(d))

def actions(p,G,RR,CC):
 nn=RR.ncols();zero=vector(ZZ,[0]*nn);rlist=list(RR.rows());cc=[vector(ZZ,z)for z in CC];Rinv=RR.inverse()
 def red(v):
  v=vector(ZZ,v);c=v*Rinv;return v-vector(ZZ,[floor(x)for x in c])*RR
 def ky(v):return tuple(red(v))
 old={d:rd(f'inputs/M3N18_9multiple_{"fresh_p7"if p==7 else"aut27_p13"}_degree{d}_data.m')for d in[2,3]}
 fields={};lookup={};points={};allpoints=[]
 for d,data in old.items():
  K=GF(p**d,'a',modulus=PolynomialRing(GF(p),'x')(data[2]));fields[d]=K;aa=K.gen();Rpol=PolynomialRing(GF(p),'x');look={}
  for qt,cl in zip(data[4],data[3]):
   vals=[Rpol(v)(aa)for v in qt];k=fieldkey(vals,p,d);assert k not in look;look[k]=vector(ZZ,cl[:-1]);allpoints.append((d,vals,vector(ZZ,cl[:-1])))
  lookup[d]=look
 H,_=hnf(rlist+cc);selected=[]
 for d,qt,cl in allpoints:
  Hnew,_=hnf(list(H.rows())+[cl])
  if abs(Hnew.det())<abs(H.det()):selected.append((d,qt,cl));H=Hnew
  if H==identity_matrix(ZZ,nn):break
 assert H==identity_matrix(ZZ,nn)
 inp=rlist+cc+[r[2]for r in selected];HH,UU=hnf(inp);assert HH==identity_matrix(ZZ,nn)
 print('M318_AUT_FRAME_GENERATORS',p,nn,[(d,fieldkey(qt,p,d))for d,qt,cl in selected],flush=True)
 mats=[]
 for gen,perm in zip(GENPAR,GENPER):
  ac,bc,ct,dt=gen;cbase=cc[perm[base]-1];imgs=[zero for _ in rlist]+[cc[perm[i]-1]-ds[i]*cbase for i in range(len(cc))]
  for d,qt,cl in selected:
   q,t=qt;K=fields[d];E=EllipticCurve(K,[q*t-q+t+2,0,q*t*t-q*t+t,0,0]);P=E(0,0);Q=E(-t,t*t);assert 3*P==E(0)and P!=E(0)and 18*Q==E(0)and 9*Q!=E(0)and 6*Q!=E(0);
   if ac==0:
    phi=E.isogeny(9*Q);ET=phi.codomain();dual=phi.dual();ai=ET.a_invariants();kx=dual.kernel_polynomial().roots(multiplicities=False);assert len(kx)==1;TT=ET(kx[0],-(ai[0]*kx[0]+ai[2])/2);assert dual(TT)==E(0);PP=phi(P);QQ=2*phi(Q)+TT;assert P.weil_pairing(6*Q,3)==PP.weil_pairing(6*QQ,3);v=normpair(ET,PP,QQ)
   else:
    PP=ac*P+6*bc*Q;QQ=ct*P+dt*Q;assert P.weil_pairing(6*Q,3)==PP.weil_pairing(6*QQ,3);v=normpair(E,PP,QQ)
   k=fieldkey(v,p,d);assert k in lookup[d];imgs.append(lookup[d][k]-d*cbase)
  A=matrix(ZZ,UU*matrix(ZZ,imgs));assert all(ky(v*A)==ky(zero)for v in rlist);assert all(ky(v*A)==ky(w)for v,w in zip(inp,imgs));assert all(ky(v*(A**(2 if ac==0 else 3)-1))==ky(zero)for v in identity_matrix(ZZ,nn));mats.append(A)
 # Synchronize the full action by its actual faithful cusp permutation.
 perm0=tuple(range(1,25));D={perm0:identity_matrix(ZZ,nn)};queue=[perm0]
 while queue:
  per0=queue.pop(0);M=D[per0]
  for per1,A in zip(GENPER,mats):
   pp=tuple(per1[i-1]for i in per0);MM=matrix(ZZ,[red(v)for v in M*A])
   if pp in D:assert all(ky(v)==ky(w)for v,w in zip(D[pp],MM))
   else:D[pp]=MM;queue.append(pp)
 assert len(D)==54 and set(D)==set(map(tuple,ALLPER));ALL=[D[tuple(pp)]for pp in ALLPER]
 for A,pp in zip(ALL,ALLPER):
  cb=cc[pp[base]-1]
  assert all(ky(cc[i]*A)==ky(cc[pp[i]-1]-ds[i]*cb)for i in range(24))
 return ALL,{'prime':p,'dimension':nn,'extra_actual_point_generators':[[d,list(fieldkey(qt,p,d))]for d,qt,cl in selected],'full_actions':[[list(map(int,row))for row in A.rows()]for A in ALL],'all_relations_and_cusp_actions_exact':True,'all_actual_marked_generator_images_exact':True}

A7,rep7=actions(7,g7,R7,c7);A13,rep13=actions(13,g13,R13,c13)
(N/'full_aut_actions_tau.json').write_text(json.dumps({'status':'PASS','data':[rep7,rep13]},indent=2)+'\n')
print('M318_FULL_AUT_ACTIONS_TAU_COMPLETE',time.monotonic()-START,flush=True)

def relations(vs,RR):
 M=matrix(ZZ,vs);m=M.nrows();S=matrix(ZZ,list(M.rows())+list(RR.rows()));kk=S.left_kernel().basis_matrix();assert kk.nrows()==m;K=kk.matrix_from_columns(range(m)).hermite_form();assert K.nrows()==K.ncols()==m
 assert all(z in ZZ for z in (K*M*RR.inverse()).list())
 HH,_=hnf(list(RR.rows())+list(M.rows()));assert abs(K.det())==abs(RR.det()/HH.det())
 return K

sur=json.loads((P/'paired_final_survivors.json').read_text())['survivors'];out=[];remaining=[]
for rec in sur:
 t=vector(ZZ,rec['target_class']);targetlist=[vector(ZZ,z)for z in c13]+[t*A for A in A13];KT=relations(targetlist,R13);good=[]
 for sid,h in rec['source_matches']:
  desc,sc=source[sid-1];sv=vector(ZZ,sc)*A7[h-1]+5*vector(ZZ,c7[ALLPER[h-1][base]-1]);sourcelist=[vector(ZZ,z)for z in c7]+[sv*A for A in A7];KS=relations(sourcelist,R7);same=KS==KT
  out.append({'target_id':rec['target_id'],'target_descriptor':rec['target_descriptor'],'source_pattern':sid,'source_descriptor':desc,'normalizing_Aut':h,'source_subgroup_order':int(abs(KS.det())),'target_subgroup_order':int(abs(KT.det())),'relation_lattices_equal':bool(same),'source_relations_fail_at_target':sum(any(z not in ZZ for z in row)for row in KS*matrix(ZZ,targetlist)*R13.inverse()),'target_relations_fail_at_source':sum(any(z not in ZZ for z in row)for row in KT*matrix(ZZ,sourcelist)*R7.inverse())})
  if same:good.append([sid,h])
 if good:rr=dict(rec);rr['source_matches']=good;remaining.append(rr)
 print('M318_ORBIT_RELATION_CASE',rec['target_id'],len(rec['source_matches']),len(good),time.monotonic()-START,flush=True)
report={'status':'PASS','scope':'Necessary orbit relation test; not by itself an exclusion if survivors remain','Aut_order':54,'global_cusp_generators':24,'tested_pairs':len(out),'remaining_pairs':sum(len(v['source_matches'])for v in remaining),'remaining_target_patterns':len(remaining),'observations':out,'survivors':remaining,'mathematical_condition':'Equality of full integral relation lattices of24globalcuspclasses and54Jacobian-Aut translates of the same normalized degree5divisorclass; torsion specializes injectively at both split odd primes','elapsed_seconds':time.monotonic()-START}
(N/'orbit_relations_tau_report.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_M318_ORBIT_RELATIONS_TAU_COMPLETE',len(out),report['remaining_pairs'],report['remaining_target_patterns'],time.monotonic()-START,flush=True)
