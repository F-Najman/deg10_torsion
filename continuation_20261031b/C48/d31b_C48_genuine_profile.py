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
outs=[]
for q,sh in[(7,0),(13,0)]:
 primary=rd(f'hecke{q}_s{sh}.m');r=next((r for r in primary[0]if max(v[0]for v in r[7])==14),primary[0][0]);id,d,ky=sources[r[0]];e=lcm([v[0]for v in r[7]]);L=GF(5**(d*e),'a');RL=PolynomialRing(L,'z');u=RL(fieldpol[d]).roots(multiplicities=False)[0];b=PR(ky[:d])(u);c=PR(ky[d:])(u);aa1=1-c;aa2=-b;aa3=-b;E=EllipticCurve(L,[0,aa2+aa1**2/4,0,aa1*aa3/2,aa3**2/4]);Q=E(0,aa3/2);assert Q.order()==48;print('D31B_C48_GENUINE_PROFILE_FIELD',q,id,d,e,L.degree(),time.monotonic()-T,flush=True)
 phis=E.isogenies_prime_degree(q);assert len(phis)==q+1;print('D31B_C48_GENUINE_PROFILE_KERNELS',q,len(phis),time.monotonic()-T,flush=True);targets=Counter()
 for ph in phis:
  QQ=ph(Q);assert list(QQ.xy())==manual(E,Q,ph);targets[json.dumps(key(xy(ph.codomain(),QQ)))]+=1
 outs.append([q,id,d,key(xy(E,Q)),key(xy(E,q*Q)),[[n,json.loads(k)]for k,n in targets.items()]]);print('D31B_C48_GENUINE_PROFILE_PASS',q,id,len(targets),time.monotonic()-T,flush=True)
(P/'genuine_profile.json').write_text(json.dumps({'status':'PASS','observations':outs,'elapsed_seconds':time.monotonic()-T},indent=2)+'\n');print('D31B_C48_GENUINE_PROFILE_COMPLETE',2,time.monotonic()-T,flush=True)
