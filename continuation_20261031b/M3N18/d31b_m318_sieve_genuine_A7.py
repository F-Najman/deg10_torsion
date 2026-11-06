from sage.all import *
from pathlib import Path
from collections import Counter
import json,time
P=Path('work/deg10/continuation_20261031b/M3N18');START=time.monotonic();Fp=GF(13);RP=PolynomialRing(Fp,'t')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
def freeze(v):return tuple(map(freeze,v))if isinstance(v,list)else v
def normal(E,A,Q):
 a1,a2,a3,a4,a6=E.a_invariants();r,y=A.xy();s=(3*r*r+2*a2*r+a4-a1*y)/(2*y+a1*r+a3);b1=a1+2*s;b2=a2-s*a1+3*r-s*s;b3=a3+r*a1+2*y;b4=a4-s*a3+2*r*a2-(y+r*s)*a1+3*r*r-2*s*y;b6=a6+r*a4+r*r*a2+r**3-y*a3-r*y*a1-y*y;assert b2==b4==b6==0;xx=Q[0]-r;yy=Q[1]-y-s*xx;u=yy/xx**2;t=-u*u*xx;q=(u*b1-t-2)/(t-1);assert u**3*b3==q*t*t+(1-q)*t;return[q,t]
def key(v):
 K=v[0].parent();n=K.degree();d=lcm([z.minpoly().degree()for z in v]);bits=next(z for z in[[1,0],[0,1]]+[[1,i]for i in range(1,13)]if(z[0]*v[0]+z[1]*v[1]).minpoly().degree()==d);g=bits[0]*v[0]+bits[1]*v[1];M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(n)]for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vector(Fp,[K(z).polynomial()[j]for j in range(n)]))))for z in v];return[int(d),bits,list(map(int,g.minpoly().list())),co]
def decode(fx):
 K=GF(13**(len(fx[0])-1),'a',modulus=RP(fx[0]));aa=K.gen();return K,[RP(v)(aa)for v in fx[1]]
def levels(E,A,Q):
 assert A!=E(0)and 3*A==E(0)and 18*Q==E(0)and 9*Q!=E(0)and 6*Q!=E(0);w=A.weil_pairing(6*Q,3);assert w!=1 and w**3==1 and w**13==w;return w

def manual(E,Q,ph):
 ker=ph.kernel_polynomial();R=ker.parent();x=R.gen();m=ker.degree();ss=[E.base_ring()(m)]
 for j in range(1,m):ss.append(-sum((ker[m-i]*ss[j-i]for i in range(1,j)),E.base_ring()(0))-j*ker[m-j])
 def tr(v):v=R(v)%ker;return sum((v[i]*ss[i]for i in range(m)),E.base_ring()(0))
 xx,yy=Q.xy();a2=E.a2();g=x**3+a2*x*x+E.a4()*x+E.a6();iv=(xx-x).inverse_mod(ker);vx=2*(yy*yy+g)*iv**2-2*a2-2*xx-4*x;vy=-2*yy+2*yy*(2*xx+a2+x)*iv-2*yy*(yy*yy+3*g)*iv**3;return[xx+tr(vx),yy+tr(vy)]
packets=[rd(f'final_Hecke7_v2_s{i}.m')for i in range(3)];cases=sorted([r for p in packets for r in p[2]if r[4]==0]);assert len(cases)==43;native={}
for s in[0,1]:
 for r in rd(f'sieve_native_A7_s{s}.m')[0]:assert r[0]not in native;native[r[0]]=r
assert set(native)=={r[0]for r in cases};contexts={};kappas=set();out=[];references=0;pointuses=0;payloadorbits=0
for case in cases:
 ci,descriptor,q,parts,hh=case;assert q==7 and hh==0 and native[ci][1]==descriptor and native[ci][3]==0 and len(parts)==len(native[ci][2]);caseobs=[]
 for part,nat in zip(parts,native[ci][2]):
  K,v=decode(part[1]);qq,tt=v[:2];E=EllipticCurve(K,v[2:7]);A=E(v[7:9]);Q=E(v[9:11]);ww=levels(E,A,Q);kappas.add(int(Fp(ww)));assert normal(E,A,Q)==[qq,tt];ks=key([qq,tt]);d=ks[0];assert d==K.degree()==part[0][0]==nat[0]and ks==nat[1];kd=key(normal(E,7*A,7*Q));_,vd=decode(part[3]);assert kd==nat[2]==key(vd);csk=freeze(ks)
  if csk not in contexts:
   oe=lcm([h[2]for h in part[4]]);L=GF(13**(d*oe),'u');RL=PolynomialRing(L,'x');theta=RL(part[1][0]).roots(multiplicities=False)[0];vals=[RP(z)(theta)for z in part[1][1]];ainv=vals[2:7];EL=EllipticCurve(L,ainv);AL=EL(vals[7:9]);QL=EL(vals[9:11]);ES=EllipticCurve(L,[0,ainv[1]+ainv[0]**2/4,0,ainv[3]+ainv[0]*ainv[2]/2,ainv[4]+ainv[2]**2/4]);AS=ES(AL[0],AL[1]+(ainv[0]*AL[0]+ainv[2])/2);QS=ES(QL[0],QL[1]+(ainv[0]*QL[0]+ainv[2])/2);phis=ES.isogenies_prime_degree(7);assert len(phis)==8;kers={tuple(ph.kernel_polynomial().list())for ph in phis};assert len(kers)==8;contexts[csk]={'L':L,'RL':RL,'sourcevals':vals,'EL':EL,'AL':AL,'QL':QL,'ES':ES,'AS':AS,'QS':QS,'kers':kers,'d':d,'degree':d*oe,'uses':0}
  ctx=contexts[csk];ctx['uses']+=1;L=ctx['L'];RL=ctx['RL'];vals=ctx['sourcevals'];covered=set();actual=Counter();expected=Counter();observ=[];assert len(part[4])==len(nat[3]);references+=2
  for witness,nativepos in zip(part[4],nat[3]):
   pd,n,e,fx=witness;assert nativepos[:2]==[e,n]and len(fx[0])-1==d*e and len(fx[1])==24
   # Align the physical witness field with the chosen source, using a gcd.
   f=RL(fx[0]);f=f.gcd(RL(fx[1][9])+vals[1]);f=f.gcd(RL(fx[1][2])-vals[2]);assert f.degree()==e;theta=f.roots(multiplicities=False)[0];vv=[RP(z)(theta)for z in fx[1]];assert vv[2:11]==vals[2:11];assert all(z**(13**d)==z for z in vv[2:11]);ET=EllipticCurve(L,vv[11:16]);AT=ET(vv[16:18]);QT=ET(vv[18:20]);wt=levels(ET,AT,QT);assert wt==L(Fp(ww));targetqt=normal(ET,AT,QT);assert targetqt==vv[:2];ker=RL(vv[20:24]);assert ker.degree()==3 and ker.is_monic();kk=tuple(ker.list());assert kk in ctx['kers'];orb=[];kk1=kk
   while kk1 not in orb:assert kk1 in ctx['kers'];orb.append(kk1);kk1=tuple(z**(13**d)for z in kk1)
   assert kk1==kk and len(orb)==e and not(covered&set(orb));covered|=set(orb);ph=ctx['ES'].isogeny(ker);assert ph.scaling_factor()==1;AP=ph(ctx['AS']);QP=ph(ctx['QS']);assert list(AP.xy())==manual(ctx['ES'],ctx['AS'],ph)and list(QP.xy())==manual(ctx['ES'],ctx['QS'],ph);assert levels(ph.codomain(),AP,QP)==L(Fp(ww));assert normal(ph.codomain(),AP,QP)==targetqt;kt=key(targetqt);assert kt==nativepos[2]and kt[0]==pd[0]and d*e==n*kt[0];actual[json.dumps(kt)]+=e;expected[json.dumps(nativepos[2])]+=e;observ.append([e,n,kt[0],kt]);payloadorbits+=1;references+=1
  assert covered==ctx['kers']and actual==expected and sum(z[0]for z in observ)==8;pointuses+=1;caseobs.append({'degree':d,'source_key':ks,'diamond7_key':kd,'kernel_orbits':observ,'all8_actual_kernels_accounted_for':True})
 out.append({'case':ci,'descriptor':descriptor,'source_parts':caseobs,'native_RR_h0':0});print('D31B_M318_SIEVE_GENUINE_CASE_PASS',ci,len(out),pointuses,len(contexts),time.monotonic()-START,flush=True)
assert len(out)==43 and pointuses==55 and references==292 and len(kappas)==1
report={'status':'PASS','group':'C3 x C18','prime':13,'hecke':7,'excluded_supports':43,'physical_source_context_uses':55,'distinct_closed_marked_sources':len(contexts),'geometric7_kernels_per_source':8,'geometric7_kernels_over_distinct_sources':8*len(contexts),'explicit_kernel_orbit_payloads_checked':payloadorbits,'native_prime_references_bound':references,'chosen_Weil_pairing_in_F13':list(kappas)[0],'full_cyclic_kernel_census_matches_actual_Frobenius_orbits':True,'both_independent3_and18_marked_points_checked':True,'stored_kernel_polynomials_and_source_field_embeddings_checked':True,'independent_pairwise_Velu_traces_on_both_marked_points':True,'all_native_source_diamond_and_target_keys_and_norms_match':True,'observations':out,'elapsed_seconds':time.monotonic()-START};(P/'sieve_independent_genuine_A7.json').write_text(json.dumps(report,indent=2,default=int)+'\n');print('D31B_M318_SIEVE_GENUINE_ALL_COMPLETE',43,55,len(contexts),8*len(contexts),references,time.monotonic()-START,flush=True)
