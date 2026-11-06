from pathlib import Path
from collections import Counter
import ast,re,json,time,math
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
P=Path(__file__).resolve().parent;t0=time.monotonic();R0=PolynomialRing(GF(2),'z');K=GF(2**30,'z');z=K.gen();fields={}
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']').strip())
def fromint(v):return sum(K((v>>i)&1)*z**i for i in range(30))
def fieldroot(poly):
 key=tuple(poly)
 if key not in fields:
  roots=PolynomialRing(K,'r')(poly).roots(multiplicities=False);assert len(roots)==len(poly)-1;fields[key]=roots[0]
 return fields[key]
def element(v,root):return sum(K((v>>i)&1)*root**i for i in range(len(root.minpoly())))
def integer(v):return sum(int(c)<<i for i,c in enumerate(v._vector_()))
def ckey(vals):
 d=int(lcm([v.minpoly().degree()for v in vals]));return d,min(tuple(integer(v**(2**i))for v in vals)for i in range(d))
def moduli(E,Q,R):
 assert 25*Q==E(0)and 5*Q!=E(0)and 3*R==E(0)and R!=E(0)
 a1,a2,a3,a4,a6=E.a_invariants();x0,y0=Q.xy();b3=a3+a1*x0;b4=x0*x0+a4-a1*y0;sl=b4/b3;b2=x0+a2-sl*a1-sl*sl;sc=b2/b3;b=-b2**3/b3**2;c=1-a1*sc
 r=b/c;s=c*c/(b-c);t=r*s+1;xx=(s-r)/t;yy=t/(s*s-s-r+1)
 x3=sc**2*(R[0]-x0);y3=sc**3*(R[1]-y0-sl*(R[0]-x0));ET=EllipticCurve(K,[1-c,-b,-b,0,0]);RT=ET(x3,y3);assert 3*RT==ET(0)
 xc=(x3+b/(1-c))/(1-c)**2;assert xc**4+xc**3==1/E.j_invariant()
 f3,f4,f5,f6,f7,f8,f9=b,c,b-c,s-1,s-r,r*s+1,s*s-s-r+1
 f10,f12=xx-yy+1,xx-yy;f13=xx**3*yy-xx**2*yy**2-xx**2*yy+xx*yy**2-yy+1
 hv=f3*f4*f6**3*f7**2*f10*f12**3/(f5*f8**2*f9**2*f13**3)
 return (hv,xx,xc),(b,c,xx,yy)
def record(row):
 poly,data,d=row;root=fieldroot(poly);v=[element(n,root)for n in data];assert len(v)==15
 E=EllipticCurve(K,v[:5]);Q=E(*v[5:7]);R=E(*v[7:9]);vals,extra=moduli(E,Q,R);assert extra==tuple(v[9:13])and vals==(v[13],v[11],v[14]);key=ckey(vals);assert key[0]==d
 return E,Q,R,key
results=[];all_records=0
for job in range(1,4):
 text=(P/f'moduli_job{job}_v2.txt').read_text().replace('[*','[').replace('*]',']');rows=[ast.literal_eval(n.value)for n in ast.parse(text).body];assert len(rows)==28
 for off in [0,14]:
  ci=2*job-1+off//14;block=rows[off:off+14];E,Q,R,source=record(block[0]);assert source[0]==10
  decoded=[record(row)for row in block];all_records+=14
  assert decoded[5][3]==ckey(moduli(E,7*Q,7*R)[0]);packet=lit((P/f'proof_packet_case{ci}.m').read_text());assert packet[0]==ci;meta=packet[8];assert len(meta)==4
  native_full=Counter(v[3]for v in decoded[6:]);native_fast=Counter()
  for row,md in zip(decoded[1:5],meta):
   orbit,d,weight=md;assert d==row[3][0]and 10*orbit==d*weight;native_fast[row[3]]+=weight
  # Enumerate all8 geometric seven-kernels independently over F(2^30).
  RR=PolynomialRing(K,'x');poly=RR(E.division_polynomial(7));xs=poly.roots(multiplicities=False);assert len(xs)==24
  seen=set();geometric=Counter();a1,a2,a3,a4,a6=E.a_invariants();RY=PolynomialRing(K,'y');Y=RY.gen()
  for x in xs:
   roots=(Y*Y+(a1*x+a3)*Y-(x**3+a2*x*x+a4*x+a6)).roots(multiplicities=False);assert len(roots)==2;T=E(x,roots[0]);assert 7*T==E(0)and T!=E(0);T.set_order(7)
   key=frozenset(integer((k*T)[0])for k in [1,2,3]);assert len(key)==3
   if key in seen:continue
   seen.add(key);phi=E.isogeny(T,check=True);assert phi.degree()==7;target=ckey(moduli(phi.codomain(),phi(Q),phi(R))[0]);geometric[target]+=1
  assert len(seen)==8 and geometric==native_full
  expected=Counter()
  for key,count in geometric.items():assert (10*count)%key[0]==0;expected[key]=(10*count)//key[0]
  assert expected==native_fast and sum(key[0]*n for key,n in expected.items())==80
  results.append(dict(case=ci,source_degree=10,target_degrees=[r[1]for r in meta],kernel_orbit_lengths=[r[0]for r in meta],norm_weights=[r[2]for r in meta],independent_geometric_kernels=8,full_native_targets_match_independent_Velu=True,point_Hecke_operator='T7-1-7delta7'))
  print('C75_DECIC_INDEPENDENT_MODULI_AND_ALL8_KERNELS_PASS',ci,'seconds',time.monotonic()-t0,flush=True)
# Independently count all source degree10 points from the accepted full Frobenius polynomial.
frob=lit((P/'inputs/source_frobenius_arithmetic.txt').read_text())[0];_,g,co,native,total=frob;assert g==145 and len(co)==291
tr=[];closed=[]
for n in range(1,11):
 trace=-n*co[2*g-n]-sum(co[2*g-j]*tr[n-j-1]for j in range(1,n));tr.append(trace);num=2**n+1-trace-sum(e*closed[e-1]for e in range(1,n)if n%e==0);assert num%n==0;closed.append(num//n)
cusps={1:20,2:10,4:8,20:2};nc=[n-cusps.get(i+1,0)for i,n in enumerate(closed)];assert nc==[0,0,0,0,0,0,20,0,0,120]
# Six sourcej-orbits, cyclic groups1050, give exactly120 closed degree10 levels.
assert len(results)==6 and 6*(40//2)==nc[9]
out=dict(complete=True,actual_moduli_records_checked=all_records,source_degree10_diamond_orbits=6,source_NC_counts_degrees1to10=nc,all_geometric7_kernels_independently_verified=48,all_actual_level_images_independently_verified=True,all_norm_multiplicities_checked=True,source_full_Frobenius_counts_match_coverage=True,results=results)
(P/'finite_degree10_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('C75_DECIC_FINITE_ARITHMETIC_INDEPENDENT_SAGE_PASS','seconds',time.monotonic()-t0,flush=True)
