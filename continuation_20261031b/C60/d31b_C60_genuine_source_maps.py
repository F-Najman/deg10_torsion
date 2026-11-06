from sage.all import *
from pathlib import Path
import json,time
P=Path('work/deg10/continuation_20261031b/C60');T=time.monotonic();Fp=GF(7);PR=PolynomialRing(Fp,'z')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
pts,cands,fibs=rd('mixed_fibres.json');native={}
for sh in range(3):
 a=rd(f'three_maps_s{sh}_v2.m')
 for rec in a[0]:native[rec[0]]=(rec,a[1])
def bc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;n1=2*b4/b3+a1;n2=b2-(b4/b3)**2-a1*b4/b3;return[-n2**3/b3**2,-(n1*n2-b3)/b3]
def xy(E,Q):
 b,c=bc(E,Q);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def key(v):
 K=v[0].parent();n=K.degree();d=int(lcm([z.minpoly().degree()for z in v]));ok=False
 for bits in[[1,0],[0,1]]+[[1,j]for j in range(1,7)]:
  g=bits[0]*v[0]+bits[1]*v[1];mp=g.minpoly()
  if mp.degree()==d:ok=True;break
 assert ok
 V=lambda z:vector(Fp,[K(z).polynomial()[i]for i in range(n)])
 M=matrix(Fp,[V(g**i)for i in range(d)]);assert M.rank()==d;co=[list(map(int,M.transpose().solve_right(V(z))))for z in v]
 return[d]+bits+list(map(int,mp.list()))+sum(co,[])
fs={};out=[]
for ci,r in enumerate(pts,1):
 d,bits,pol,co=r;ks=tuple(pol)
 if ks not in fs:fs[ks]=GF(7**d,'a',modulus=PR(pol))
 K=fs[ks];b=PR(co[0])(K.gen());c=PR(co[1])(K.gen());E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert 2*P2==E(0)and P2!=E(0)and 30*Q==E(0)and all((30//p)*Q!=E(0)for p in[2,3,5])and P2!=15*Q
 assert key([b,c])==[d]+bits+pol+sum(co,[])
 rec,ims=native[ci];assert rec[1]==d;rows=[]
 for j,QQ in enumerate([Q,P2+16*Q,P2+Q]):
  assert 30*QQ==E(0)and all((30//p)*QQ!=E(0)for p in[2,3,5]);kk=key(xy(E,QQ));nr=rec[2][j];im=ims[nr[1]-1];assert kk==im[0]and kk[0]==im[2]==nr[3]and nr[4]*kk[0]==d and nr[2]==[nr[4]*z for z in im[1]];rows.append(kk)
 out.append([ci,d,rows])
 if ci%200==0:print('D31B_C60_GENUINE_MAP_PROGRESS',ci,time.monotonic()-T,flush=True)
assert len(out)==1674 and len(fibs)==1848 and len(cands)==3414
(P/'independent_genuine_source_maps.json').write_text(json.dumps({'status':'PASS','source_points':1674,'marked_target_maps':5022,'all_exact_orders_and_independence_verified':True,'all_intrinsic_target_keys_and_norms_match':True,'observations':out,'elapsed_seconds':time.monotonic()-T},indent=2)+'\n');print('D31B_C60_GENUINE_SOURCE_MAPS_COMPLETE',1674,5022,time.monotonic()-T,flush=True)
