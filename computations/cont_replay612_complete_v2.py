from pathlib import Path
import runpy,ast,json,itertools
from sage.all import QQ,ZZ,matrix,vector,PolynomialRing,NumberField
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
ctx=runpy.run_path(str(H/'computations/cont_replay612_norm.py'));read=ctx['read'];sig,ds,lifts=ctx['sig'],ctx['ds'],ctx['lifts'];C,key=ctx['C'],ctx['key'];S=set(ctx['left']);assert len(S)==192
U=matrix(ZZ,sig).transpose().row_module();assert U.rank()==23
base=ds.index(1);ub=U.basis_matrix();proj=ub.matrix_from_columns([j for j in range(24) if j!=base]);assert abs(proj.det())==5971968
old=read(D/'degree_independent_612_L8/M6N12_global_L8_cusp_classes_root.m');reduced=read(D/'degree_independent_612_L8/M6N12_reduced_L8_cusp_coefficients_magma_root.m');oldrr=read(D/'degree_independent_612_L8/M6N12_L8_char0_six_pencils_arithmetic_root.txt')
assert sig==old[0] and ds==old[1] and len(oldrr)==6 and all(row[2]==2 and row[3]==8 for row in oldrr)
def cl(v):return key([sum(v[j]*ctx['cgs'][0][j][i] for j in range(24)) for i in range(len(ctx['cgs'][0][0]))])
lookup={cl(lifts[k-1][1]):k for k in S};matches={}
for idx,c8 in reduced[:6]:
 assert min(c8)>=0 and sum(c8[i]*ds[i] for i in range(24))==8 and idx in {row[1] for row in oldrr}
 for j,d in enumerate(ds):
  if d!=1:continue
  eff=c8[:];eff[j]+=2;ky=cl(eff)
  if ky not in lookup:continue
  k=lookup[ky];assert k not in matches;vv=vector(ZZ,lifts[k-1][1]);vv[base]+=10-sum(vv[i]*ds[i] for i in range(24));assert vv-vector(ZZ,eff) in U;matches[k]=[idx,j+1]
assert set(matches)==set(range(769,865)) and len(matches)==96
# The remaining classes form a single actual 96-map orbit; every image class was freshly recomputed.
pars,orbits=read(D/'M6N12_quintic_class_orbits_v2.m');assert len(pars)==96 and len(orbits)==1
norm=lambda p:min(tuple(x%q for x,q in zip(p,[6,6,6,12])),tuple(-x%q for x,q in zip(p,[6,6,6,12])))
def product(x,y):
 a,b,c,d=x;e,f,g,h=y;return norm([a*e+2*b*g,a*f+b*h,c*e+d*g,2*c*f+d*h])
ps=set(map(norm,pars));assert len(ps)==96 and all(a%2==d%2==1 and (a*d-2*b*c)%6==1 for a,b,c,d in pars)
assert all(product(x,y) in ps for x,y in itertools.product(pars,repeat=2))
k,i,images,classes,keys,actuals=orbits[0];assert k==2791 and len(images)==len(classes)==len(keys)==len(actuals)==96 and len(set(classes))==len(set(map(tuple,keys)))==96
assert set(classes)==S-set(matches)
for j,c,actual in zip(images,classes,actuals):assert key(actual)==ctx['point'][5,j] and key([2*x for x in actual])==cl(lifts[c-1][1])
# The effective representative differs from the signed divisor by global modular-unit divisors.
red=json.loads((D/'M6N12_effective_quintic_lifts_v2.json').read_text())[0];assert red['index']==2791
signed=vector(ZZ,lifts[2790][1]);signed[base]+=10-sum(signed[j]*ds[j] for j in range(24));eff=vector(ZZ,red['effective']);assert vector(ZZ,red['signed'])==signed and eff-signed in U and min(eff)>=0 and eff*vector(ZZ,ds)==10
assert read(D/'M6N12_char0_coordinate_identity.txt')==[6,True,True]
cert=read(D/'M6N12_char0_norm_certificate.m');idx,degree,coeff,defdata,fdata,hdata,Hdata,Gdata,leaddata,fibers,found=cert;assert idx==2791 and degree==10 and coeff==list(eff) and found
Rz=PolynomialRing(QQ,'z');F=NumberField(Rz([1,-1,1]),'z');z=F.gen();Ru=PolynomialRing(F,'u');u=Ru.gen();B=Ru.fraction_field();Y=PolynomialRing(B,'v');v=Y.gen();Ra=PolynomialRing(F,'a');a=Ra.gen()
def df(row):return F(QQ(row[0][0])/QQ(row[0][1]))+F(QQ(row[1][0])/QQ(row[1][1]))*z
def dp(row,R):return R([df(c) for c in row])
def dr(row):return B(dp(row[0],Ru))/B(dp(row[1],Ru))
pol=Y([dr(c) for c in defdata]);fun=Y([dr(c) for c in fdata]);hL=dr(hdata);assert pol.degree()==6;pol=pol.monic()
rows=[]
for j in range(6):
 rem=(fun*v**j)%pol;row=rem.list()+[B(0)]*(6-len(rem.list()));rows.append(row)
M=matrix(B,rows);cp=M.charpoly('t');assert cp.degree()==6
H=[Ra(0) for j in range(11)]
for j in range(7):
 r=hL*cp[j];assert r.denominator().degree()==0;hp=Ru([c/r.denominator()[0] for c in r.numerator().list()]);assert hp.degree()<=10
 for k,c in enumerate(hp.list()):H[k]+=c*a**j
assert H==[dp(c,Ra) for c in Hdata]
Q=Ra.fraction_field();cc=[Q(c)/H[10] for c in H];bb=[Q(0)]*6;bb[5]=Q(1)
for j in range(4,-1,-1):bb[j]=(cc[5+j]-sum(bb[t]*bb[5+j-t] for t in range(j+1,5)))/2
G=Ra(0)
for j in range(5):G=G.gcd((cc[j]-sum(bb[t]*bb[j-t] for t in range(j+1))).numerator())
G=G.monic();assert G==dp(Gdata,Ra) and H[10]==dp(leaddata,Ra) and G==a-(90+QQ(15)/2*z)
candidate=(G*H[10]).monic();seen=set();good=[]
for rootdata,polydata,passed,degs,mults in fibers:
 root=df(rootdata);assert root not in seen and candidate(root)==0;seen.add(root)
 while candidate(root)==0:candidate,rem=candidate.quo_rem(a-root);assert rem==0
 spec=Ru([c(root) for c in H]);spec=spec.monic();assert spec==dp(polydata,Ru) and spec.is_square()==passed
 if passed:assert degs==[5] and mults==[2];good.append((root,spec))
assert candidate==1 and len(good)==1
root,spec=good[0];g5=spec.sqrt().monic();assert g5.degree()==5 and g5.is_irreducible();K=F.extension(g5,'w');w=K.gen();wv=(4*z+1)*w**4+(-z+4)*w**3+(-2*z+10)*w**2+(-13*z+2)*w+2*z-5

def ev(c):
 c=B(c);return sum((K(a)*w**j for j,a in enumerate(c.numerator().list())),K(0))/sum((K(a)*w**j for j,a in enumerate(c.denominator().list())),K(0))
def evaly(p):return sum((ev(c)*wv**j for j,c in enumerate(p.list())),K(0))
assert evaly(pol)==0 and evaly(fun)==K(root)
pu=Y([c.derivative() for c in pol.list()]);fu=Y([c.derivative() for c in fun.list()]);py=evaly(pol.derivative());assert py!=0
deriv=evaly(fu)-evaly(fun.derivative())*evaly(pu)/py;assert deriv==0
# A degree-ten map with a degree-five unramified-u point of zero multiplicity at least two has fiber exactly 2P.
report={'curve':'X1(6,12) over Q(zeta6)','absolute_degree':10,'relative_degree':5,'source_divisors':24576,'doubled_cusp_classes':3360,'classes_after_complete_F13_sieve':192,'classes_excluded_by_a_forced_cusp_basepoint':96,'remaining_class_orbits':1,'remaining_class_orbit_size':96,'independently_recomputed_characteristic_polynomial_and_square_constraints':'PASS','unique_finite_noncuspidal_even_fiber_parameter':'90 + 15*z/2','independent_quintic_point_and_implicit_derivative':'PASS','closed_quintic_points_up_to_full_level_automorphisms':1,'full_torsion_group':[6,12],'classification_case':'ONE SPORADIC OCCURRENCE; COMPLETE FOR CURVES CONTAINING C6 x C12'}
(D/'M6N12_complete_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');(D/'M6N12_basepoint_matches.json').write_text(json.dumps(matches,indent=2)+'\n');print('DECIC612_COMPLETE_INDEPENDENT_REPLAY_PASS',report,flush=True)
