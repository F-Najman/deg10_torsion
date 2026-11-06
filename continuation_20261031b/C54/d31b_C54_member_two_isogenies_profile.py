from sage.all import *
from pathlib import Path
from itertools import product
from collections import Counter
import json,time
P=Path('work/deg10/continuation_20261031b/C54');T=time.monotonic();Fp=GF(5);R=PolynomialRing(Fp,'t');a=json.loads((P/'class_members_profile.m').read_text().replace('[*','[').replace('*]',']'))
def mag(v):
 if isinstance(v,(int,Integer)):return str(v)
 f=all(isinstance(z,(int,Integer))for z in v);return('['if f else'[*')+','.join(mag(z)for z in v)+(']'if f else'*]')
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;s=b4/b3;n1=2*s+a1;n2=b2-s*s-a1*s;return[-n2**3/b3**2,1-n1*n2/b3]
def key(v):
 K=v[0].parent();n=K.degree();d=lcm([a.minpoly().degree()for a in v]);bits=next(z for z in[[1,0],[0,1],[1,1],[1,2],[1,3],[1,4]]if(z[0]*v[0]+z[1]*v[1]).minpoly().degree()==d);g=bits[0]*v[0]+bits[1]*v[1];M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(n)]for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vector(Fp,[K(z).polynomial()[j]for j in range(n)]))))for z in v];return[int(d)]+bits+list(map(int,g.minpoly().list()))+sum(co,[])
images=[];imidx={};branches={};pointdata=[]
for pi,d,ky in a[0]:
 pol=ky[3:4+d];K=GF(5**d,'a',modulus=R(pol));g=K.gen();b=R(ky[4+d:4+2*d])(g);c=R(ky[4+2*d:4+3*d])(g);E=EllipticCurve(K,[1-c,-b,-b,0,0]);Q=E(0,0);assert Q.order()==27;factors=E.division_polynomial(2).factor();assert sum(f.degree()*m for f,m in factors)==3 and all(m==1 for f,m in factors);pbs=[]
 for fi,(f,m)in enumerate(factors,1):
  e=f.degree();L=GF(5**(d*e),'u');RL=PolynomialRing(L,'x');u=RL(pol).roots(multiplicities=False)[0];emb=lambda z:sum((L(z.polynomial()[i])*u**i for i in range(d)),L(0));EL=EllipticCurve(L,list(map(emb,E.a_invariants())));QL=EL(list(map(emb,Q.xy())));sx=RL(list(map(emb,f.list()))).roots(multiplicities=False)[0];ai=EL.a_invariants();S=EL(sx,-(ai[0]*sx+ai[2])/2);assert S.order()==2;phi=EL.isogeny(RL.gen()-sx);ET=phi.codomain();QT=phi(QL);assert QT.order()==27 and list(QT.xy())==[QL[0]+(QL+S)[0]-S[0],QL[1]+(QL+S)[1]-S[1]];tk=key(bc(ET,QT));kt=tuple(tk)
  if kt not in imidx:imidx[kt]=len(images)+1;images.append([len(images)+1,tk])
  assert d*e%tk[0]==0;pbs.append([fi,int(e),imidx[kt],int(d*e//tk[0])])
 branches[pi]=pbs;pointdata.append([pi,d,pbs])
# Every factor is a Frobenius orbit of the geometric two-kernel fibre.
def choices(p,m):
 bs=branches[p];out=[]
 for ns in product(*(range(m//b[1]+1)for b in bs)):
  if sum(n*b[1]for n,b in zip(ns,bs))==m:out.append([[b[2],n*b[3]]for n,b in zip(ns,bs)if n])
 return out
lifts=[]
for mi,ci,cl,nc,cd in a[1]:
 opts=[choices(p,m)for p,m in nc]
 for sels in product(*opts):
  cc=Counter()
  for ss in sels:
   for i,n in ss:cc[i]+=n
  assert sum(images[i-1][1][0]*n for i,n in cc.items())==10-cd;lifts.append([len(lifts)+1,mi,ci,cl,cd,sorted([i,n]for i,n in cc.items())])
(P/'member_two_isogenies_profile.m').write_text(mag([images,pointdata,lifts])+'\n');(P/'member_two_isogenies_profile.json').write_text(json.dumps({'status':'PASS','target_points':len(a[0]),'actual_kernel_orbits':sum(len(x[2])for x in pointdata),'distinct_target_images':len(images),'RR_members':len(a[1]),'possible_full_source_NC_lifts':len(lifts),'cusp_degrees':dict(Counter(r[4]for r in lifts)),'elapsed_seconds':time.monotonic()-T},indent=2)+'\n');print('D31B_C54_MEMBER_TWO_ISOGENIES_PROFILE_COMPLETE',len(a[0]),len(images),len(a[1]),len(lifts),time.monotonic()-T,flush=True)
