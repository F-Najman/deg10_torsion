from sage.all import *
from pathlib import Path
from collections import Counter
import ast,json,re,time
H=Path('work/deg10');R=H/'continuation_20261031b';out=R/'Y228_shared';start=time.monotonic();Fp=GF(3);RP=PolynomialRing(Fp,'z')
def read(p):return ast.literal_eval(re.sub(r'\bfalse\b','False',re.sub(r'\btrue\b','True',p.read_text().replace('[*','[').replace('*]',']'))))
a=read(R/'C56/projected_Y28.m');b=read(R/'M4N28/fourmap_Y28.m');assert a[0]==b[0];fd=a[0];fields={int(d):GF(3**int(d),'a',modulus=RP(pol))for d,pol in fd}
def coeff(v,d):return tuple(int(v.polynomial()[i])for i in range(d))
def ky(v,d):return min(tuple(c for z in vv for c in coeff(z,d))for vv in [[v[0]**(3**j),v[1]**(3**j)]for j in range(d)])
def intrinsic(d,key,e):
 K=fields[d];z=K.gen();v=[RP(key[:d])(z),RP(key[d:])(z)];assert lcm([c.minpoly().degree()for c in v])==e
 if d!=e:
  L=fields[e];pol=L.modulus();rr=PolynomialRing(K,'X')(pol.list()).roots(multiplicities=False)[0];M=matrix(Fp,[vector(Fp,coeff(rr**j,d))for j in range(e)])
  vv=[M.transpose().solve_right(vector(Fp,coeff(c,d)))for c in v];v=[sum((L(cc)*L.gen()**i for i,cc in enumerate(row)),L(0))for row in vv]
 return ky(v,e)
cache={}
def orbit(e,key):
 tag=(e,tuple(key))
 if tag in cache:return cache[tag]
 K=fields[e];z=K.gen();bb=RP(key[:e])(z);cc=RP(key[e:])(z);E=EllipticCurve(K,[0,cc,0,bb*(1-bb-cc),0]);P=E(0,0);Q=E(bb,bb);assert 2*P==E(0)and P!=E(0)and P!=14*Q and 28*Q==E(0)and 14*Q!=E(0)and 4*Q!=E(0)
 vals=[]
 for B in [0,1]:
  for Cc in [0,1]:
   for u in [1,3,5,9,11,13]:
    PP=P+B*14*Q;QQ=Cc*P+u*Q;r=PP[0];xx=QQ[0]-r;yy=QQ[1];assert xx!=0 and yy!=0
    vals.append(ky([xx**3/yy**2,(cc+3*r)*xx**2/yy**2],e))
 cache[tag]=vals;return vals
cases={};bindings=[]
def add(label,idx,cusp,components):
 assert sum(e*m for e,m,key in components)+cusp==10
 candidates=[]
 for g in range(24):
  parts=Counter()
  for e,m,key in components:parts[(e,orbit(e,key)[g])]+=m
  candidates.append(tuple(sorted((e,m,k)for (e,k),m in parts.items())))
 canon=min(candidates);tag=(cusp,canon)
 if tag not in cases:cases[tag]=len(cases)+1
 bindings.append([label,int(idx),cases[tag]])
rd={int(v[0]):v for v in a[1]}
for idx,ids,cusp in a[2]:
 comps=[]
 for ri in ids:
  _,d,key,e,n=rd[int(ri)];assert d==e*n;comps.append((int(e),int(n),intrinsic(int(d),key,int(e))))
 add('C56',idx,int(cusp),comps)
rb={int(v[0]):v for v in b[1]}
for idx,ids,cusp in b[2]:
 comps=[]
 for ri in ids:
  _,d,m,key,*_=rb[int(ri)];comps.append((int(d),int(m),tuple(key)))
 add('M4N28',idx,int(cusp),comps)
rows=[[i,cusp,[[e,m,list(key)]for e,m,key in canon]]for (cusp,canon),i in sorted(cases.items(),key=lambda z:z[1])]
def magma(v):
 if isinstance(v,(int,Integer)):return str(v)
 if isinstance(v,str):return json.dumps(v)
 assert isinstance(v,(list,tuple))
 if all(isinstance(z,(int,Integer))for z in v):return '['+', '.join(map(str,v))+']'
 return '[*'+', '.join(magma(z)for z in v)+'*]'
(out/'union_cases.m').write_text('UnionFields:='+magma(fd)+';\nUnionCases:='+magma(rows)+';\n')
result={'status':'PASS','C56_original_cases':len(a[2]),'M4N28_original_cases':len(b[2]),'union_distinct_cases':len(rows),'bindings':bindings,'cases':rows,'normalization':'Simultaneous actual24automorphisms and individualclosedFrobenius; lower moduli-field norms and repeated components are retained','seconds':time.monotonic()-start}
(out/'union_normalization.json').write_text(json.dumps(result,indent=2)+'\n');print('D31B_Y228_UNION_NORMALIZED',len(a[2]),len(b[2]),len(rows),time.monotonic()-start,flush=True)
