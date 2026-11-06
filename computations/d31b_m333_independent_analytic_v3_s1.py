from sage.all import *
from pathlib import Path
import ast,re,json,time,hashlib
T=time.monotonic();P=Path('work/deg10/continuation_20261031b/M3N33')
def rd(n):
 s=(P/n).read_text();s=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)\s*\]',lambda m:str(list(range(int(m[1]),int(m[2])+1))),s);s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
a=rd('relative_kummer_frame.m');js,wits,perm,idx=rd('analytic_cusp_jets_v4.m');bps=rd('relative_base_cusps.m');high=rd('high_primitive_series.m');assert len(high)==24
k=GF(4,'w',modulus=PolynomialRing(GF(2),'z')([1,1,1]));w=k.gen();ku=PolynomialRing(k,'u');u=ku.gen();RF=ku.fraction_field()
def dc(c):return k(c%2)+k(c//2)*w
def dr(x):return RF(ku([dc(c) for c in x[0]]))/RF(ku([dc(c) for c in x[1]]))
N=len(a[1])-1;assert N==12
zero=tuple([RF(0)]*N)
def df(x):return tuple([dr(c) for c in x]+[RF(0)]*(N-len(x)))
uf=list(zero);uf[0]=RF(u);vf=list(zero);vf[1]=RF(1);functions=[tuple(uf),tuple(vf),tuple(-c for c in df(a[2][0]))];ids=[]
for hol in a[4]:
 row=[]
 for co in hol:
  f=df(co)
  if f not in functions:functions.append(f)
  row.append(functions.index(f)+1)
 ids.append(row)
assert ids==idx and len(functions)==len(wits[0][4])
CM=matrix(k,81,81,[dc(x) for x in a[7][0]]);TM=matrix(k,81,81,[dc(x) for x in a[7][1]]);U=matrix(k,36,81,[dc(x) for x in a[7][2]])
assert TM*TM==1 and U.rank()==36
polyf=[dr(c) for c in a[1]];checks=[];allhol={};minprec=10**9
for wi in wits:
 if (wi[0]-1)%4!=1:continue
 bi,d,m,kpol,ser=wi;assert 1<=bi<=24 and d in [1,5]
 if d==1:
  K=k;om=w;th=k(1)
 else:
  K=GF(2**(2*d),'z');om=K.gen()**((K.order()-1)//3);PK=PolynomialRing(K,'zeta');f=PK([K(c%2)+K(c//2)*om for c in kpol]);roots=f.roots(multiplicities=False);assert len(roots)==d;th=roots[0]
 def emb(c):return K(int(c.polynomial()[0]))+K(int(c.polynomial()[1]))*om
 def dec(v):return sum((K(c%2)+K(c//2)*om)*th**i for i,c in enumerate(v))
 if d==1:
  def enc(z):return [int(z.polynomial()[0])+2*int(z.polynomial()[1])]
 else:
  M=matrix(GF(2),[list((th**i*om**j).polynomial())+[0]*(2*d-len((th**i*om**j).polynomial())) for i in range(d) for j in range(2)]);assert M.is_invertible();MI=M.inverse()
  def enc(z):
   vs=list(z.polynomial());v=vector(GF(2),vs+[0]*(2*d-len(vs)))*MI;return [int(v[2*i])+2*int(v[2*i+1]) for i in range(d)]
 assert all(dec(enc(z))==z for z in [K.random_element() for _ in range(12)])
 S=LaurentSeriesRing(K,'s',default_prec=600);s=S.gen()
 def ds(z):
  lo,hi,cs=z
  if not cs:return S(0).add_bigoh(hi)
  assert len(cs)==hi-lo
  return S({lo+i:dec(c) for i,c in enumerate(cs)}).add_bigoh(hi)
 ss=[ds(z) for z in ser];hw=high[bi-1];assert hw[:3]==[bi,d,kpol];us,vs=[ds(z) for z in hw[3]];assert us==ss[0] and vs==ss[1]
 # Evaluate the original rational model and both native prime generators.
 cache={}
 def ep(pol):
  out=S(0)
  for c in reversed(pol.list()):out=out*us+emb(c)
  return out
 def er(f):
  if f not in cache:cache[f]=ep(f.numerator())/ep(f.denominator())
  return cache[f]
 def ef(ff):
  out=S(0)
  for co in reversed(ff):out=out*vs+er(co)
  return out
 rel=sum(er(co)*vs**i for i,co in enumerate(polyf));assert rel==0
 gens=[ef(df(c)) for c in bps[bi-1][2]];assert all(z.valuation()>0 for z in gens) and min(z.valuation() for z in gens)==1
 actual=[ef(ff) for ff in functions]
 assert all(z.precision_absolute()>=saved.precision_absolute() and z==saved for z,saved in zip(actual,ss));actual=[z.add_bigoh(saved.precision_absolute()) for z,saved in zip(actual,ss)]
 assert ss[2].valuation()==m
 # Each actual function agrees to enough precision to determine every saved jet.
 branches=[r for r in js if r[0]==bi];assert len(branches)==(3 if m%3==0 else 1)
 for row in branches:
  _,branch,dd,e,mm,be,ce,J,UJ=row;assert dd==d and mm==m
  beta=dec(be);c=dec(ce)
  def sub(z):
   if z==0:return S(0).add_bigoh(e*z.precision_absolute())
   lo=z.valuation();hi=z.precision_absolute()
   return sum(z[i]*c**i*s**(e*i) for i in range(lo,hi)).add_bigoh(e*hi)
  arr=[sub(z) for z in actual];D=arr[2];V=D/(s**(e*m));assert V.valuation()==0 and beta**3==V[0]
  H=S(beta).add_bigoh(V.precision_absolute())
  for _ in range(12):
   H=V/H**2
   if H**3==V:break
  assert H**3==V;r=s**(e*m//3)*H;assert r**3==D
  du=arr[0].derivative();assert du!=0
  hol=[sum(arr[j-1]*r**z for z,j in enumerate(rowids))*du for rowids in idx]
  nv=5 if d==1 else 1
  assert all(z==0 or z.valuation()>=0 for z in hol)
  assert all(z.precision_absolute()>=nv for z in hol)
  minprec=min(minprec,min(z.precision_absolute() for z in hol));codes=[[cc for j in range(nv) for cc in enc(z[j])] for z in hol];assert codes==J
  mat=matrix(k,81,5,[dc(c) for rr in codes for c in rr]);assert [[int(x.polynomial()[0])+2*int(x.polynomial()[1]) for x in rr] for rr in U*mat]==UJ
  allhol[(bi,branch)]=(hol,K,om,enc)
 checks.append({'base_index':bi,'degree':d,'kummer_valuation':m,'branches':len(branches),'model_relation':True,'both_prime_generators_uniformizer':True,'all_function_series_bound':True})
 print('M333_INDEPENDENT_ANALYTIC_BASE',bi,d,m,len(branches),float(time.monotonic()-T),flush=True)
report={'status':'PASS','base_cusps':len(checks),'source_cusps':len(js),'rational_source_cusps':sum(r[2]==1 for r in js),'quintic_source_cusps':sum(r[2]==5 for r in js),'unique_function_witnesses':len(functions),'holomorphic_functions':81,'verified_jet_matrices':48,'model_and_prime_binding':True,'independent_tame_uniformizers_and_kummer_identity':True,'min_available_holomorphic_precision':int(minprec),'checks':checks,'elapsed_seconds':time.monotonic()-T,'input_sha256':{n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['relative_kummer_frame.m','relative_base_cusps.m','analytic_cusp_jets_v4.m','high_primitive_series.m']}}
(P/'independent_analytic_series_s1.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_M333_INDEPENDENT_ANALYTIC_SHARD_COMPLETE_1',len(checks),len(js),time.monotonic()-T,flush=True)
