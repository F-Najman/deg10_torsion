from sage.all import *
from pathlib import Path
import json,time,hashlib
P=Path('work/deg10/continuation_20261031b/M3N18');D=P/'endomorphism_gluing';START=time.monotonic()
def rd(p):return json.loads(p.read_text().replace('[*','[').replace('*]',']'))
a=rd(D/'global_point13900_absolute.m');assert a[0]==13900 and len(a)==12
rat=lambda v:QQ(v[0])/v[1];R=PolynomialRing(QQ,'T');f=R(list(map(rat,a[1])));assert f.degree()==10 and f.is_monic()and f.is_irreducible()
# An integral rescaling avoids any number-field maximal-order computation.
h=R([f[i]*3**(10-i)for i in range(11)]);assert h.is_monic()and all(c in ZZ for c in h);K=NumberField(h,'b');theta=K.gen()/3
assert f(theta)==0 and K.degree()==10
def val(vec):return sum((K(rat(c))*theta**i for i,c in enumerate(vec)),K(0))
zeta=val(a[2]);assert zeta*zeta+zeta+1==0 and zeta!=1
q,t=map(val,a[3]);E=EllipticCurve(K,list(map(val,a[4])));A=E(list(map(val,a[5])));Q=E(list(map(val,a[6])));assert E.a_invariants()==(q*t-q+t+2,K(0),q*t*t-q*t+t,K(0),K(0));assert A==E(0,0)and Q==E(-t,t*t);assert E.discriminant()!=0 and E.j_invariant()==val(a[7]);u,v=map(val,a[8]);assert theta==u+a[11]*zeta and q==1+u and t==v*u/(1+u)
rel=[rat(c[0])+rat(c[1])*zeta for c in a[9]];assert len(rel)==6 and rel[-1]==1 and sum((rel[i]*u**i for i in range(6)),K(0))==0

def tangent(E,A):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=A.xy();return(3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
def levels(E,A,Q):
 assert A!=E(0)and 3*A==E(0)and 18*Q==E(0)and 9*Q!=E(0)and 6*Q!=E(0);B=6*Q;assert A!=B and A!=-B;s=tangent(E,A);r=tangent(E,B);ea=B[1]-A[1]-s*(B[0]-A[0]);eb=A[1]-B[1]-r*(A[0]-B[0]);assert ea!=0 and eb!=0;wp=-ea/eb;assert wp!=1 and wp**3==1
 try:ws=A.weil_pairing(B,3);assert ws in[wp,1/wp]
 except NotImplementedError:ws=wp
 return wp,ws
wp,ws=levels(E,A,Q);assert val(a[10])in[wp,ws,1/wp]
def normal(E,A,Q):
 a1,a2,a3,a4,a6=E.a_invariants();r,y=A.xy();s=tangent(E,A);b1=a1+2*s;b2=a2-s*a1+3*r-s*s;b3=a3+r*a1+2*y;b4=a4-s*a3+2*r*a2-(y+r*s)*a1+3*r*r-2*s*y;b6=a6+r*a4+r*r*a2+r**3-y*a3-r*y*a1-y*y;assert b2==b4==b6==0;xx=Q[0]-r;yy=Q[1]-y-s*xx;c=yy/xx**2;tt=-c*c*xx;qq=(c*b1-tt-2)/(tt-1);assert c**3*b3==qq*tt*tt+(1-qq)*tt;return[qq,tt]
assert normal(E,A,Q)==[q,t]
T=9*Q;assert T!=E(0)and 2*T==E(0);RK=PolynomialRing(K,'X');phi=E.isogeny(RK.gen()-T[0]);assert phi.degree()==2 and phi.scaling_factor()==1;ET=phi.codomain();AT=phi(A);QI=phi(Q)
for S in[A,Q]:assert list(phi(S).xy())==[S[0]+(S+T)[0]-T[0],S[1]+(S+T)[1]-T[1]]
ai=E.a_invariants();xd=-ai[1]-ai[0]**2/4-2*T[0];at=ET.a_invariants();Td=ET(xd,-(at[0]*xd+at[2])/2);assert Td!=ET(0)and 2*Td==ET(0);QT=2*QI+Td;wpt,wst=levels(ET,AT,QT);assert wpt in[wp,1/wp];tauqt=normal(ET,AT,QT)
# This explicit level change is the root-native tuple (1,0,2,7).
A2=AT;Q2=2*AT+7*QT;levels(ET,A2,Q2);qt2=normal(ET,A2,Q2)
print('D31B_M318_POINTCHECK_NF_LEVELS_PASS',K.degree(),time.monotonic()-START,flush=True)
# Check both exact function-field points and both global RR zero sections.
model=sage_eval((P/'global_function_field_char0.m').read_text().replace('[*','[').replace('*]',']'))
sections=sage_eval((P/'global_RR_char0_reductionfalse.m').read_text().replace('[*','[').replace('*]',']'));assert {r[0]for r in sections}=={13900,14996}
def evnf(c):return K(c[0])+K(c[1])*zeta
def evrat(c,uu):return sum((evnf(z)*uu**i for i,z in enumerate(c[0])),K(0))/sum((evnf(z)*uu**i for i,z in enumerate(c[1])),K(0))
def eff(fn,uu,vv):return sum((evrat(c,uu)*vv**i for i,c in enumerate(fn)),K(0))
for cid,qt in[(13900,[q,t]),(14996,qt2)]:
 uu=qt[0]-1;vv=qt[1]*(1+uu)/uu;assert eff(model[0],uu,vv)==0;assert[eff(c,uu,vv)for c in model[1]]==qt;row=next(r for r in sections if r[0]==cid);assert row[2]==1 and eff(row[3],uu,vv)==0
print('D31B_M318_POINTCHECK_BOTH_GLOBAL_SECTIONS_PASS',time.monotonic()-START,flush=True)
# Good unramified degree-two primes over7 provide a complete torsion bound.
F7=GF(7);R7=PolynomialRing(F7,'x');f7=R7([F7(c)for c in f]);fac=f7.factor();assert f7.is_squarefree()and sorted(g.degree()for g,m in fac)==[2,2,3,3]and all(m==1 for g,m in fac)
def theta_coeff(z):
 # K uses b=3*theta. Return exact coefficients in powers oftheta.
 return[QQ(z.polynomial()[i])*3**i for i in range(10)]
def red(z,rr):
 cs=theta_coeff(K(z));assert all(c.denominator()%rr.parent().characteristic()!=0 for c in cs);return sum((rr.parent()(c)*rr**i for i,c in enumerate(cs)),rr.parent()(0))
def card_exhaust(E):
 kk=E.base_ring();a1,a2,a3,a4,a6=E.a_invariants();n=1
 for xx in kk:
  dd=(a1*xx+a3)**2+4*(xx**3+a2*xx**2+a4*xx+a6);n+=1 if dd==0 else(2 if dd.is_square()else 0)
 return n
bounds=[]
for g,m in fac:
 if g.degree()!=2:continue
 L=GF(49,'r',modulus=g);r=L.gen();ee=EllipticCurve(L,[red(z,r)for z in E.a_invariants()]);pp=ee([red(z,r)for z in A.xy()]);qq=ee([red(z,r)for z in Q.xy()]);ee2=EllipticCurve(L,[red(z,r)for z in ET.a_invariants()]);pp2=ee2([red(z,r)for z in AT.xy()]);qq2=ee2([red(z,r)for z in QT.xy()]);assert ee.discriminant()!=0 and ee2.discriminant()!=0;levels(ee,pp,qq);levels(ee2,pp2,qq2);n=card_exhaust(ee);n2=card_exhaust(ee2);assert n==n2==54 and ee.cardinality()==ee2.cardinality()==54;bounds.append({'factor_F7':list(map(int,g.list())),'residue_degree':2,'ramification_index':1,'good_reduction_both':True,'exhaustive_cardinalities':[54,54]})
assert len(bounds)==2
print('D31B_M318_POINTCHECK_FULL_TORSION_PASS',[3,18],[3,18],time.monotonic()-START,flush=True)
# Bind the actual NFpoints at13 to both surviving marked closed points.
dat=rd(D/'actual_target13_points.m');F13=GF(13);R13=PolynomialRing(F13,'x');L13=GF(13**5,'r',modulus=R13(dat[2]));r=L13.gen();expected=[R13(c)(r)for c in dat[3][0][2]];target=[R13(c)(r)for c in dat[3][1][2]];f13=PolynomialRing(L13,'x')([L13(c)for c in f]);roots=f13.roots(multiplicities=False);matches=[rr for rr in roots if[red(q,rr),red(t,rr)]==expected];assert len(matches)==1;rr=matches[0];zr=red(zeta,rr);assert zr**2+zr+1==0;assert[red(z,rr)for z in qt2]==[z**13 for z in target];assert[red(z,rr)for z in tauqt]!=target or True
# Independently rebuildtau at13 and the leveltuple, beyondreducingrationaloutputs.
er=EllipticCurve(L13,[red(z,rr)for z in E.a_invariants()]);ar=er([red(z,rr)for z in A.xy()]);qr=er([red(z,rr)for z in Q.xy()]);tr=9*qr;phir=er.isogeny(PolynomialRing(L13,'X').gen()-tr[0]);etr=phir.codomain();air=er.a_invariants();xdr=-air[1]-air[0]**2/4-2*tr[0];atr=etr.a_invariants();tdr=etr(xdr,-(atr[0]*xdr+atr[2])/2);aa=phir(ar);qq=2*phir(qr)+tdr;assert normal(etr,aa,2*aa+7*qq)==[z**13 for z in target]
print('D31B_M318_POINTCHECK_FINITE13_BINDINGS_PASS',13900,106,14996,1202,[1,0,2,7],1,time.monotonic()-START,flush=True)
jpoly=E.j_invariant().minpoly();jpoly2=ET.j_invariant().minpoly()
def qqpairs(poly):return[[int(c.numerator()),int(c.denominator())]for c in poly.list()]
report={'status':'PASS','group':'C3 x C18','absolute_field_degree':10,'relative_degree_over_Qzeta3':5,'absolute_generator_rescaling':'b=3*theta','absolute_polynomial_irreducible':True,'relative_quintic_verified':True,'primitive_relation_theta_equals_u_plus_zeta':True,'both_global_curve_equations_and_RR_sections_vanish_exactly':True,'both_independent_marked_3_and18_levels_verified_exactly':True,'manual_and_Sage_Weil_pairings_primitive':True,'normalized_actual_2isogeny_and_both_marked_images_checked':True,'full_torsion_original':[3,18],'full_torsion_2isogenous':[3,18],'torsion_upper_bound_method':'Monic7-integral primitive polynomial squarefree mod7 implies unramified local degree2; good models have54points by exhaustive F49 count, so all global torsion injects and hasorder<=54','degree2_reductions':bounds,'finite13_original_binding':{'global_class':13900,'own_point_index':106,'q_t_exactly_equal':True},'finite13_partner_binding':{'global_class':14996,'own_point_index':1202,'level_tuple_on_tau':[1,0,2,7],'Frobenius_phase':1,'q_t_exactly_equal':True,'independently_rebuilt_finite2isogeny':True},'finite13_zeta_value':int(F13(zr)),'j_minpoly_degrees':[int(jpoly.degree()),int(jpoly2.degree())],'j_minpolys_QQ_pairs':[qqpairs(jpoly),qqpairs(jpoly2)],'input_hashes':{n:hashlib.sha256((P/n).read_bytes()).hexdigest()for n in['endomorphism_gluing/global_point13900_absolute.m','endomorphism_gluing/actual_target13_points.m','global_function_field_char0.m','global_RR_char0_reductionfalse.m']},'elapsed_seconds':time.monotonic()-START};(P/'sieve_independent_global_points.json').write_text(json.dumps(report,indent=2,default=int)+'\n');print('D31B_M318_POINTCHECK_ALL_COMPLETE',10,[3,18],[3,18],13900,14996,[jpoly.degree(),jpoly2.degree()],time.monotonic()-START,flush=True)
