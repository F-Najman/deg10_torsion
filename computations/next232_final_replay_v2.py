from pathlib import Path
from math import gcd
from collections import Counter
import json
P=Path('work/deg10/C232_resolution_20261030');n=32
labels={(a,r,s,t) for a in range(n) for r in range(n) if gcd(gcd(a,r),n)==1 for s in range(2) for t in range(2) if (s*r-t*a)%2==1}
def orb(x):
 a,r,s,t=x;return frozenset(((e*(a+k*r))%n,e*r%n,e*(s+k*t)%2,e*t%2) for e in [-1,1] for k in range(n))
cs=sorted(set(map(orb,labels)),key=lambda x:sorted(x));ix={x:i for i,c in enumerate(cs) for x in c};assert len(cs)==56
fro=[]
for c in cs:
 im={ix[(3*a%n,r,s,t)] for a,r,s,t in c};assert len(im)==1;fro.append(im.pop())
def cycles(perm):
 unseen=set(range(len(perm)));out=[]
 while unseen:
  x=min(unseen);c=[]
  while x not in c:c.append(x);unseen.remove(x);x=perm[x]
  out.append(c)
 return out
assert Counter(map(len,cycles(fro)))=={1:16,2:8,4:2,8:2}

import ast
from sage.all import GF,PolynomialRing,EllipticCurve
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
DATA=Path('work/deg10/C232_resolution_20261030')
def rd(n):return ast.literal_eval((DATA/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def act(a0):
 out=[]
 for C in cs:
  a,r,s,t=min(C);out.append(ix[(a0*a%32,a0*r%32,s,t)])
 return out
h=act(15);loc=cycles(fro);li={i:j for j,C in enumerate(loc) for i in C};gc=[]
for C in cs:
 a,r,s,t=min(C);gc.append(ix[(5*a%32,r,s,t)])
assert all(h[h[i]]==i for i in range(56))
def val(i):
 a,r,s,t=min(cs[i]);u=min(r,32-r);w=len(cs[i])//2;v=[-2*u,-2*u,32-4*u] if t==0 else [u,16-u,0];assert all(w*x%32==0 for x in v);return tuple(w*x//32 for x in v)
v=[val(i) for i in range(56)];cg=rd('hecke5_cusp_geometry.m');assert Counter(tuple(x[:1]+x[2:]) for x in cg[2])==Counter(v[C[0]] for C in loc if len(C)<=4)
assert all(gc[C[0]] in C for C in loc if len(C)==1)
assert all(h[i]==i for C in loc if len(C) in [2,4] for i in C)
oc=cycles(h);oi={i:j for j,C in enumerate(oc) for i in C};qfro=[oi[fro[C[0]]] for C in oc];qloc=cycles(qfro);pred=[]
for C in oc:
 i=C[0];a,r,s,t=min(cs[i]);u2=min(2*r%32,32-(2*r%32));width=len(cs[i])//2;e=2 if h[i]==i else 1;raw=[-2*u2,-2*u2,32-4*u2] if t==0 else [u2,16-u2,0];arr=[x*width//(32*e) for x in raw];assert all(x*width%(32*e)==0 for x in raw);rad=(u2 if t else -u2) if r%2 else 0;assert rad*width%(64*e)==0;arr.append(rad*width//(64*e));assert (v[i][0]+v[h[i]][0])%e==0;arr.append((v[i][0]+v[h[i]][0])//e);pred.append(tuple(arr))
qdata=rd('quotient13_filter.m');assert Counter(map(tuple,qdata[0]))==Counter(pred[C[0]] for C in qloc);assert len(qdata[1])==17
# Recompute the two disjoint coset sets by integer Hermite normal forms.
a=rd('quotient13_proven_picard.m');rels=a[5];hg=qdata[2];nn=len(rels[0])-1;HH=hermite_normal_form(Matrix([r[:-1] for r in rels+hg]).T);assert abs(int(HH.det()))==32

def qkey(vec):
 assert vec[-1]==0;w=list(vec[:-1])
 for j in range(nn-1,-1,-1):
  q=w[j]//int(HH[j,j])
  for i in range(j+1):w[i]-=q*int(HH[i,j])
 return tuple(w)
rows,common=rd('quotient13_orbit_filter.m');assert len(rows)==32 and common==[];R=PolynomialRing(GF(3),'x');src=rd('source6.m');K=GF(3**6,'aa',modulus=R(src[1]));aa=K.gen();ky=rd('filterprobe6.m')[2][0][2];b=sum(K(x)*aa**i for i,x in enumerate(ky[:6]));c=sum(K(x)*aa**i for i,x in enumerate(ky[6:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);S32=E(b,b);assert 32*S32==E(0) and 16*S32!=E(0)
sets=[];seen=None;hphase=None
for row in rows:
 pars,vals,cl,fbindex,checks,poss,cum=row;bb,dd,u=pars;PS=P2+bb*(16*S32);SS=dd*P2+u*S32;V=2*SS;T=16*SS;x=V[0]-PS[0];y=V[1];alpha=x/y;b16=x**3/y**2;c16=(c+3*PS[0])*x*x/y**2;bt=3*T[0]**2+2*c*T[0]+E.a4();ww=alpha*((SS[0]-T[0])**2-bt)/(2*SS[1]);actual=[b16,c16,ww];assert all(v**27==v for v in actual);assert vals==[[int(x) for x in v._vector_()] for v in actual];EB=EllipticCurve(K,[0,c16,0,b16*(1-b16-c16),0]);QV=EB(b16,b16);assert ww*ww==b16-(8*QV)[0]
 out=set()
 for ci,native,good in checks:
  z=qkey([x+y for x,y in zip(cl,a[9][ci-1])]);yes=qkey([2*x for x in z]+[0])==(0,)*nn;assert bool(good)==yes
  if yes:out.add(z)
 assert len(out)==len(poss);sets.append(out);seen=out if seen is None else seen&out;assert len(seen)==len(cum)
assert len(sets[0])==len(sets[1])==2 and not (sets[0]&sets[1]) and seen==set();assert rows[0][0]==[0,0,1] and rows[1][0]==[0,0,3]
# The final source support cannot have a rational cusp in its complement.
tar=rd('target_proven_picard.m');TH=hermite_normal_form(Matrix([r[:-1] for r in tar[5]+tar[8]]).T);tn=TH.rows
def tkey(vec):
 assert vec[-1]==0;w=list(vec[:-1])
 for j in range(tn-1,-1,-1):
  q=w[j]//int(TH[j,j])
  for i in range(j+1):w[i]-=q*int(TH[i,j])
 return tuple(w)
def add(x,y):return tkey([a+b for a,b in zip(x,y)]+[0])
allowed=[set() for _ in range(5)];allowed[0]={(0,)*tn}
for degree,cl in zip(tar[9],tar[10]):
 for d in range(degree,5):allowed[d]|={add(x,tkey(cl)) for x in allowed[d-degree]}
last=rd('filterprobe6.m')[1][-1];needs4=[]
for var,ti,e,mult,cl,yes in last[3]:
 z=tkey([-x for x in cl]);assert z in allowed[4]
 if z not in allowed[3]:needs4.append(var)
assert needs4
# The involution acts by Frobenius cubed on the noncuspidal source point.
U=15*S32;x=U[0];y=U[1];pair=[x**3/y**2,c*x*x/y**2];assert pair==[b**27,c**27]
out={'status':'PASS','full_source_cusps':56,'source_cusp_Frobenius_signature':dict(Counter(map(len,loc))),'all_16_residue_rational_cusps_globally_rational':True,'all_degree_2_and_4_source_cusps_fixed_geometrically_by_h':True,'target_maps_requiring_degree4_cusp_class':needs4,'source_involution_phase':3,'quotient_cusps':40,'quotient_local_signature':dict(Counter(map(len,qloc))),'quotient_global_cusp_blocks':17,'quotient_cuspidal_image_index':32,'all_32_quotient_moduli_images_checked':True,'first_two_coset_set_sizes':[2,2],'first_two_coset_sets_disjoint':True,'first_two_parameters':[rows[0][0],rows[1][0]],'common_cosets':0,'excluded_groups':[[2,32],[1,64]]};(DATA/'independent_final_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C232_INDEPENDENT_FINAL_REPLAY_COMPLETE',out,flush=True)
