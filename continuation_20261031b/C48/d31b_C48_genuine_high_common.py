from sage.all import *
from pathlib import Path
from collections import Counter
import json,time
P=Path('work/deg10/continuation_20261031b/C48');T=time.monotonic();Fp=GF(5);PR=PolynomialRing(Fp,'t')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
src=rd('high_hecke_sources.m');fieldpol=dict(src[0]);sources={r[0]:r for r in src[1]}
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;s=b4/b3;n1=2*s+a1;n2=b2-s*s-a1*s;return[-n2**3/b3**2,1-n1*n2/b3]
def xy(E,Q):
 b,c=bc(E,Q);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def key(v):
 K=v[0].parent();n=K.degree();d=lcm([a.minpoly().degree()for a in v]);bits=next(z for z in[[1,0],[0,1],[1,1],[1,2],[1,3],[1,4]]if(z[0]*v[0]+z[1]*v[1]).minpoly().degree()==d);g=bits[0]*v[0]+bits[1]*v[1];M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(n)]for i in range(d)]);co=[list(map(int,M.transpose().solve_right(vector(Fp,[K(a).polynomial()[j]for j in range(n)]))))for a in v];return[int(d),bits,list(map(int,g.minpoly().list())),co]
def manual(E,Q,phi):
 ker=phi.kernel_polynomial();R=ker.parent();z=R.gen();m=ker.degree();s=[E.base_ring()(m)]
 for k in range(1,m):s.append(-sum((ker[m-i]*s[k-i]for i in range(1,k)),E.base_ring()(0))-k*ker[m-k])
 def tr(v):v=R(v)%ker;return sum((v[i]*s[i]for i in range(m)),E.base_ring()(0))
 xp,yp=Q.xy();a2=E.a2();g=z**3+a2*z*z+E.a4()*z+E.a6();inv=(xp-z).inverse_mod(ker);px=2*(yp*yp+g)*inv**2-2*a2-2*xp-4*z;py=-2*yp+2*yp*(2*xp+a2+z)*inv-2*yp*(yp*yp+3*g)*inv**3;return[xp+tr(px),yp+tr(py)]

from sys import argv
q=int(argv[1]);sh=int(argv[2]);nsh=int(argv[3]);native=rd(f'parent_native_high_q{q}.m');out=[];fieldcache={}
for ii,row in enumerate(native):
 if ii%nsh!=sh:continue
 id,d,ky=sources[row[0]];assert row[1]==d and row[2]==q;e=lcm([v[0]for v in row[7]]);cachekey=(d,e)
 if cachekey not in fieldcache:
  L=GF(5**(d*e),'a');RL=PolynomialRing(L,'z');u=RL(fieldpol[d]).roots(multiplicities=False)[0];fieldcache[cachekey]=(L,RL,u)
 L,RL,u=fieldcache[cachekey];b=PR(ky[:d])(u);c=PR(ky[d:])(u);aa1=1-c;aa2=-b;aa3=-b;E=EllipticCurve(L,[0,aa2+aa1**2/4,0,aa1*aa3/2,aa3**2/4]);Q=E(0,aa3/2);assert Q.order()==48;sk=key(xy(E,Q));dk=key(xy(E,q*Q));assert sk==row[3]and dk==row[5]and d==row[4]*sk[0]and q*d==row[6]*dk[0]
 phis=E.isogenies_prime_degree(q);assert len(phis)==q+1;kers={tuple(ph.kernel_polynomial().list())for ph in phis};assert len(kers)==q+1;unseen=set(kers);actual=Counter();obs=[]
 while unseen:
  kk=min(unseen,key=str);orbit=[];z=kk
  while z not in orbit:orbit.append(z);assert z in kers;z=tuple(a**(5**d)for a in z)
  assert z==kk and set(orbit)<=unseen;unseen-=set(orbit);ph=E.isogeny(RL(list(kk)));assert ph.scaling_factor()==1;QQ=ph(Q);assert list(QQ.xy())==manual(E,Q,ph);assert QQ.order()==48;tk=key(xy(ph.codomain(),QQ));actual[json.dumps(tk)]+=len(orbit);assert d*len(orbit)%tk[0]==0;obs.append([len(orbit),d*len(orbit)//tk[0],tk[0],tk])
 expected=Counter()
 for oe,n,tk in row[7]:assert d*oe==n*tk[0];expected[json.dumps(tk)]+=oe
 assert actual==expected;assert row[8]==0;out.append({'source':id,'degree':d,'hecke':q,'source_key':sk,'diamond_key':dk,'target_orbits':obs,'all_geometric_kernels':q+1,'native_RR':row[8]});print('D31B_C48_GENUINE_HIGH_POINT_PASS',q,id,len(out),time.monotonic()-T,flush=True)
result={'status':'PASS','level':48,'prime':5,'hecke':q,'shard':sh,'shards':nsh,'cases':out,'method':'Complete Sage prime-kernel census, exact Frobenius permutation, normalized Velu isogenies and independent pairwise trace formulas, all source and norm-weighted marked-target keys bound to fresh native primes','elapsed_seconds':time.monotonic()-T};(P/f'independent_genuine_high_q{q}_s{sh}.json').write_text(json.dumps(result,indent=2)+'\n');print('D31B_C48_GENUINE_HIGH_COMPLETE',q,sh,len(out),time.monotonic()-T,flush=True)
