from sage.all import *
from pathlib import Path
from collections import defaultdict,Counter
import ast,json,time
D=Path(__file__).resolve().parents[1];C=D/'continuation_20261031b/C66';rd=lambda n:ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'));M=rd('gonal_model_F5.m');F=rd('gonal_fibre_filter.m');frame=rd('F2_proven_frame.m');sig,units,vecs=rd('pencil_cusp_vectors_F5.m');assert units==M[2] and len(vecs)==20
wanted={tuple(z['class'])for z in json.loads((C/'F2_degree10_classification.json').read_text())['gonal_classes']};lookup={tuple(z[3]):z for z in frame[7]};classes=[]
for v in vecs:
 z=[sum(v[i]*lookup[tuple(sig[i])][2][j]for i in range(17))for j in range(3)];z[0]%=15;z[1]%=42373650;assert z[2]==0;classes.append(tuple(z))
assert set(classes)==wanted and len(set(classes[:10]))==len(set(classes[10:]))==5 and set(classes[:10]).isdisjoint(classes[10:]);print('D31B_C66_GONAL_CLASS_COVERAGE_COMPLETE 10 2',flush=True)
k=GF(5);R=PolynomialRing(k,'u');start=time.monotonic();fib=[];kernels=0;points=set();details=[]
def xy(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;c1=2*b4/b3+a1;c2=b2-(b4/b3)**2-a1*b4/b3;b=-c2**3/b3**2;c=-(c1*c2-b3)/b3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return((s-r)/t,t/(s*s-s-r+1))
def degreejoint(vals):
 d=vals[0].parent().degree()
 return next(i for i in divisors(d)if all(v**(5**i)==v for v in vals))
def valrat(er,y):return sum(k(c)*y**i for i,c in enumerate(er[0]))/sum(k(c)*y**i for i,c in enumerate(er[1]))
def valef(ef,x,y):return sum(valrat(c,y)*x**i for i,c in enumerate(ef))
def field(cv):return GF(5**(len(cv)-1),name='z',modulus=R(cv))
for bi,aa,cm,bounds,records in F:
 assert cm in(0,10);assert max(bounds)<=4;counts=[0]*120;byplace=defaultdict(list);source_mass=0
 if cm==10:assert not records
 for rec in records:
  pd,m,d,e,pol,ep,rv,eq,vals,bs=rec;assert e==1;L=field(pol);z=L.gen();cv=lambda v:sum(L(c)*z**i for i,c in enumerate(v));E=EllipticCurve(L,[cv(v)for v in ep[0]]);P=E([cv(v)for v in ep[1]]);r=cv(rv);a1,a2,a3,a4,a6=E.a_invariants();S=E(r,-(a1*r+a3)/2);assert S.order()==2 and P.order()==33;source66=17*P+S;assert source66.order()==66 and degreejoint(xy(E,source66))==d*e
  x,y=xy(E,P);assert degreejoint((x,y))==d;assert valef(M[1][bi-1],x,y)==L(aa);assert valef(M[4][0],x,y)==E.a1();assert all(valef(v,x,y)==a for v,a in zip(M[4],E.a_invariants()));assert sum(valrat(c,y)*x**i for i,c in enumerate(M[0]))==0
  assert pd[0]==d and pd[1]==1;assert len(pd[2])==2 and all(valef(g,x,y)==0 for g in pd[2]);key=json.dumps(pd);byplace[key].append(r);points.add((bi,aa,key));
  phi=E.isogeny(S,algorithm='velu');ES=phi.codomain();Q=phi(P);EN=EllipticCurve(L,[cv(v)for v in eq[0]]);QN=EN([cv(v)for v in eq[1]]);assert QN.order()==33;assert xy(ES,Q)==xy(EN,QN)
  check=[]
  for gf in M[1]:
   for a in units:
    xx,yy=xy(ES,a*Q);check.append(valef(gf,xx,yy))
  assert check==[cv(v)for v in vals];newbs=[]
  for i,v in enumerate(check):
   if v**5==v:b=int(k(v));counts[i*6+b]+=m*d*e;newbs.append([i+1,b])
  assert bs==newbs;kernels+=1;source_mass+=m*d*e;details.append({'base':bi,'value':aa,'source_degree':d,'kernel_orbit':e,'source66_degree':d*e,'target_moduli_degree':degreejoint(xy(ES,Q)),'rational_targets':bs})
 for key,roots in byplace.items():
  assert len(roots)==3 and len(set(roots))==3
 assert counts==bounds;assert cm==10 or source_mass==30
 fib.append({'base':bi,'value':aa,'cuspidal_mass_on_X133':cm,'marked_source_mass':source_mass,'maximum_paired_fibre_degree':max(bounds),'kernel_records':len(records)})
out={'status':'PASS','global_gonal_classes':10,'diamond_orbits':2,'transported_pencil_classes':[list(v)for v in classes],'fibres':fib,'genuine_marked_2kernels':kernels,'distinct_source33_places':len(points),'source66_degrees_verified':True,'all_TwoGenerators_bound_to_actual_moduli':True,'all_paired_fibre_degrees_at_most':4,'details':details,'elapsed_seconds':time.monotonic()-start};(C/'gonal_genuine_independent.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C66_GONAL_GENUINE_INDEPENDENT_COMPLETE',kernels,len(points),'maxdegree4',flush=True)
