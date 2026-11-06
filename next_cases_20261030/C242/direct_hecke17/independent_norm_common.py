from pathlib import Path
import ast,json,time,hashlib
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors,matrix,vector,sage_eval
ROOT=Path('work/deg10/next_cases_20261030/C242/direct_hecke17');t0=time.monotonic();used={}
def rd(p):
 ss=p.read_text();used[str(p)]=hashlib.sha256(ss.encode()).hexdigest();return ast.literal_eval(ss.replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
summary=json.loads((ROOT/'final_test_summary.json').read_text());selected=[r for r in summary['selected_nonprincipal_tests'] if r['tag']==TAG];ell=11 if TAG=='A11pure' else 17
proofs=rd(ROOT/TAG/'explicit_hecke_data.m');need=summary['needed_kernel_record_indices_zero_based'][TAG]
R5=PolynomialRing(GF(5),'z');fcache={}
def fld(pol):
 kk=tuple(pol)
 if kk not in fcache:fcache[kk]=GF(5**(len(pol)-1),'z',modulus=R5(pol))
 return fcache[kk]
def ev(K,cs):return R5(cs)(K.gen())
def curve(dat):
 pol,ai,pv,qv=dat;K=fld(pol);E=EllipticCurve(K,[ev(K,a) for a in ai]);pp=[ev(K,a) for a in pv];qq=[ev(K,a) for a in qv];assert pp[2]==qq[2]==1;return K,E,E(pp),E(qq)
def order42(E,Q):return 42*Q==E(0) and 21*Q!=E(0) and 14*Q!=E(0) and 6*Q!=E(0)
def degree(vs,d):return min(int(e) for e in divisors(d) if all(v**(5**int(e))==v for v in vs))
def xy(E,Q):
 x,y=Q[0],Q[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2=3*x+E.a2()-m*m;a3=2*y;b=-a2**3/a3**2;c=1-2*m*a2/a3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return [(s-r)/t,t/(s*s-s-r+1)]
Rxy=PolynomialRing(GF(5),2,names=('x','y'));xx,yy=Rxy.gens();model=Path('work/deg10/shared/mdmagma/models_X1_n/FFFc42.txt');eq=sage_eval(model.read_text().replace('\n',' '),locals={'x':xx,'y':yy})
def pk(vals):
 L=vals[0].parent();d=degree(vals,L.degree());a=None
 for aa in range(5):
  g=vals[0]+aa*vals[1]
  if degree([g],L.degree())==d:a=aa;break
 if a is None:g=vals[1];assert degree([g],L.degree())==d;a=5
 B=matrix(GF(5),[(g**j)._vector_() for j in range(d)]);cx=B.solve_left(vector(GF(5),vals[0]._vector_()));cy=B.solve_left(vector(GF(5),vals[1]._vector_()));return(tuple(map(int,g.minpoly())),tuple(map(int,cx)),tuple(map(int,cy)),a)
points={}
def geom(vals):
 k=pk(vals);d=len(k[0])-1;point={'polynomial':list(k[0]),'x_coefficients':list(k[1]),'y_coefficients':list(k[2]),'primitive_linear_combination':k[3],'degree':d,'unramified_y':bool(eq.derivative(yy)(*vals)!=0),'x_primitive':degree([vals[0]],vals[0].parent().degree())==d};sk=json.dumps(k,separators=(',',':'));points[sk]=point;return sk,d
count=0;weighted=0;cache={};source_reports=[]
def evaluate(target,images,d):
 global count,weighted
 K,E,Paux,Q=curve(target);assert K.degree()==d and order42(E,Q);orig=xy(E,Q);delta=xy(E,ell*Q);assert eq(*orig)==eq(*delta)==0;ko,de=geom(orig);kd,de2=geom(delta);assert de==de2 and d%de==0;D=Counter();D[ko]-=d//de;D[kd]-=ell*d//de;RK=PolynomialRing(K,'w');normprod=RK.one();positive_degree=0
 for e,srcL,gdat,outL,embed,expect,dd,w in images:
  L,EL,PL,QL=curve(srcL);LL,EG,PG,QG=curve(outL);assert L is LL and L.degree()==d*e;zz=ev(L,embed);em=K.hom([zz],L);assert [em(a) for a in E.a_invariants()]==list(EL.a_invariants()) and EL([em(a) for a in Q])==QL and EL([em(a) for a in Paux])==PL
  RL=PolynomialRing(L,'w');g=RL([ev(L,a) for a in gdat]);assert g.degree()==(ell-1)//2 and g.is_monic();iso=EL.isogeny(g,check=True);assert iso.degree()==ell and iso.codomain()==EG and iso(PL)==PG and iso(QL)==QG
  vals=xy(EG,QG);assert vals==[ev(L,a) for a in expect] and eq(*vals)==0;k,de=geom(vals);assert de==dd and d*e%de==0 and w==d*e//de;D[k]+=w;positive_degree+=w*de
  ng=RL.one();gg=g
  for _ in range(e):ng*=gg;gg=RL([a**(5**d) for a in gg])
  assert gg==g;B=matrix(GF(5),[(zz**j)._vector_() for j in range(d)]);co=[]
  for aa in ng:
   vv=B.solve_left(vector(GF(5),aa._vector_()));bb=ev(K,list(vv));assert em(bb)==aa;co.append(bb)
  normprod*=RK(co);count+=1;weighted+=e
 assert sum(r[0] for r in images)==ell+1 and normprod==RK(E.division_polynomial(ell)).monic() and positive_degree==(ell+1)*d
 D={k:int(v) for k,v in D.items() if v};assert sum(points[k]['degree']*v for k,v in D.items())==0;return D
for index in need:
 kk,src,target,de,w,images=proofs[index];mapno,d=kk[:2];ky=kk[2:];K,E,P,Q=curve(src);assert P==E(0,0) and P!=21*Q and order42(E,Q) and list(map(int,Q[0]._vector_()))==ky[:d] and list(map(int,E.a2()._vector_()))==ky[d:]
 L,EG,Qa,Qb=curve(target);assert K is L and Qa==Qb
 if mapno<=3:assert EG==E and Qa==([Q,P+Q,P+22*Q][mapno-1])
 else:
  T=[P,21*Q,P+21*Q][mapno-4];U=P+22*Q if mapno==5 else Q;iso2=E.isogeny(T,check=True);assert iso2.codomain()==EG and iso2(U)==Qa
 assert degree(xy(EG,Qa),d)==de and w==d//de;cache[tuple(kk)]=evaluate(target,images,d);source_reports.append({'record_index_zero_based':index,'map':mapno,'source_degree':d,'kernel_orbits':len(images),'geometric_kernels':sum(r[0] for r in images)})
 print('C242_FINAL_INDEPENDENT_NORM_POINT_PASS',TAG,index,mapno,d,len(images),round(time.monotonic()-t0,2),flush=True)
casegeoms={};expected=set()
for row in selected:
 ci=row['case_id'];mapno=row['selected_map'];D=Counter()
 for d,ky in row['source_terms']:
  kk=tuple([mapno,d]+ky);expected.add(kk)
  for k,v in cache[kk].items():D[k]+=v
 casegeoms[str(ci)]={k:int(v) for k,v in D.items() if v}
 assert sum(points[k]['degree']*v for k,v in casegeoms[str(ci)].items())==0
assert expected==set(cache)
out={'status':'PASS','tag':TAG,'ell':ell,'selected_cases':[r['case_id'] for r in selected],'source_summands':len(need),'kernel_orbits':count,'geometric_kernel_count':weighted,'all_Velu_maps_and_level_images_checked':True,'all_division_polynomial_factors_exhausted':True,'all_residue_norm_weights_checked':True,'source_reports':source_reports,'seconds':time.monotonic()-t0,'sha256':used};(ROOT/f'independent_norm_{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/f'divisor_geometry_{TAG}.json').write_text(json.dumps({'status':'PASS','points':points,'cases':casegeoms},indent=2)+'\n');print('C242_FINAL_NINE_INDEPENDENT_NORM_PASS',TAG,len(casegeoms),count,weighted,round(time.monotonic()-t0,2),flush=True)
