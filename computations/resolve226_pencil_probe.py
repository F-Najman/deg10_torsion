from pathlib import Path
import ast,json
from sage.all import ZZ,QQ,matrix,vector
from fractions import Fraction
from math import gcd
P=Path('work/deg10/C226_resolution_20261030');C=Path('work/deg10/continuation_20261029b/C226')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
labels,rows=read(P/'optimized_cusp_divisors.m');B=matrix(ZZ,json.loads((P/'optimized_cusp_divisors.json').read_text())['principal_relation_basis']);print('RELATION_BASIS_POSITIVE_DEGREES',[sum(max(0,x) for x in b) for b in B.rows()])
n=26;ls={(a,r,s,t) for a in range(n) for r in range(n) if gcd(gcd(a,r),n)==1 for s in range(2) for t in range(2) if (s*r-t*a)%2==1}
def orb(v):
 a,r,s,t=v;return frozenset(((e*(a+k*r))%n,e*r%n,e*(s+k*t)%2,e*t%2) for e in [-1,1] for k in range(n))
cs=sorted(set(map(orb,ls)),key=lambda x:sorted(x));index={v:i for i,c in enumerate(cs) for v in c}
def val(v):
 a,r,s,t=v;u=min(r,n-r)
 for w in range(1,n+1):
  z=((a+w*r)%n,r,(s+w*t)%2,t)
  if any(z==(e*a%n,e*r%n,e*s%2,e*t%2) for e in [-1,1]):break
 x=Fraction(w,n)*(-2*u if t==0 else u);assert x.denominator==1;return int(x)
def gal(i,k):a,r,s,t=min(cs[i]);return index[(k*a%n,r,k*s%2,t)]
rat=[i for i in range(36) if all(gal(i,k)==i for k in range(n) if gcd(k,n)==1)];byval={val(min(cs[i])):i for i in rat};assert set(byval)==set(labels)
out=[]
for idx,desc,co,candidates in rows:
 def action(v):
  a,r,s,t=v
  return ((13*s+(18 if idx==24 else 5)*a)%26,(13*t+(18 if idx==24 else 5)*r)%26,a%2 if idx==24 else s,r%2 if idx==24 else t)
 perm=[index[action(min(c))] for c in cs];assert all(perm[perm[i]]==i for i in range(36));gp={v:val(min(cs[perm[byval[v]]])) for v in labels};A=[max(0,x) for x in co];Bv=[max(0,-x) for x in co];
 def invariant(v):return all(v[labels.index(gp[x])]==v[labels.index(x)] for x in labels)
 out.append({'id':idx,'g_on_rational_cusp_labels':gp,'g_fixed_rational_cusps':[v for v in labels if gp[v]==v],'positive_part_g_invariant':invariant(A),'negative_part_g_invariant':invariant(Bv),'positive_degree':sum(A),'negative_degree':sum(Bv)})
(P/'pencil_involutions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
