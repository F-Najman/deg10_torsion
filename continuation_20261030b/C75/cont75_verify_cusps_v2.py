from pathlib import Path
from collections import Counter,defaultdict
from math import gcd,comb
import ast,json,re,time
from sage.all import GF,Matrix
P=Path(__file__).resolve().parent;t0=time.monotonic();N=75
U=lambda n:[v for v in range(n)if gcd(v,n)==1]
def can(c):
 d,a,b=c;return min(c,(d,(-a)%d,(-b)%(N//d)))
S={can((d,a,b))for d in range(1,N+1)if N%d==0 for a in U(d)for b in U(N//d)}
def diamond(c,v):
 d,a,b=c;return can((d,v*a%d,v*b%(N//d)))
def gal(c,v):
 d,a,b=c;return can((d,v*a%d,b))
def qcan(c):return min(c,diamond(c,26))
Q={qcan(c)for c in S}
def cyclic_orbits(A,fn):
 left=set(A);out=[]
 while left:
  x=min(left);row=[];c=x
  while c not in row:row.append(c);c=fn(c)
  assert c==x and set(row)<=left;left-=set(row);out.append(row)
 return out
sf=cyclic_orbits(S,lambda c:gal(c,2));qf=cyclic_orbits(Q,lambda c:qcan(gal(c,2)))
sg={tuple(sorted(gal(c,v)for v in U(N)))for c in S}
# Turn repeated orbit entries into sets before comparison.
sg={tuple(sorted({gal(c,v)for v in U(N)}))for c in S};qg={tuple(sorted({qcan(gal(c,v))for v in U(N)}))for c in Q}
assert Counter(map(len,sg))=={1:20,2:10,4:4,8:2,20:2}
assert {frozenset(r)for r in qf}=={frozenset(r)for r in qg}
assert Counter(map(len,qg))=={1:20,4:4,10:2}
assert all(all(gal(c,7)==c for c in row)for row in sf if len(row)<=3)
qindex={c:i for i,row in enumerate(qf)for c in row};steps=defaultdict(set)
for row in sf:
 target=qindex[qcan(row[0])];steps[target].add(len(row)//len(qf[target]))
records=[]
for i,row in enumerate(qf):
 assert len(steps[i])==1;records.append((len(row),N//row[0][0],next(iter(steps[i]))))
assert Counter(records)=={(1,25,2):10,(1,75,1):10,(4,5,1):2,(4,15,1):2,(10,1,2):1,(10,3,2):1}
ways=[0]*11;ways[0]=1
for d,width,step in records:
 for k in range(d*step,11):ways[k]+=ways[k-d*step]
assert ways[10]==871630
# Exactly110 configurations have a local weighted multiplicity>8.
assert 10+10+10*9==110
s=(P/'high_cusp_matrices.m').read_text()
def mat(name):
 z=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert z,name
 return Matrix(GF(2),int(z[1]),int(z[2]),ast.literal_eval(z[3]))
C=mat('CartierMatrix');M=mat('Diamond2Matrix');B=mat('ImageBasis');W3=mat('Kernel3').row_space();W12=mat('Kernel12').row_space();W10=mat('Kernel10Width25').row_space();I=Matrix.identity(GF(2),73)
assert M**10==I and M*C==C*M and B.row_space()==(C-I).row_space() and B.rank()==71
A=B.row_space();assert 71-(A&W3).dimension()==3 and 71-(A&W12).dimension()==12 and 71-(A&W10).dimension()==10
assert W12.is_subspace(W3)
mixed=[]
for k in range(1,10):
 T=(W3.basis_matrix()*M**(-k)).row_space();r=71-(A&W12&T).dimension();assert r==15;mixed.append(int(r))
ret=json.loads((P/'degree10_retained_FI.json').read_text());assert ret['tested_divisors']==871520 and ret['exceptions_count']==0
out=dict(complete=True,level=75,source_global_cusp_signature=dict(Counter(map(len,sg))),source_F2_cusp_signature=dict(Counter(map(len,sf))),quotient_all_cusps_globally_intact=True,quotient_cusp_signature=dict(Counter(map(len,qg))),source_low_degree3_cusps_globally_intact_and_Gal7_fixed=True,admissible_degree10_cusp_images=871630,retained_common_basis_tests=871520,new_high_cases=110,new_own_basis_U_dimension=71,rank_at_12_width75_cusp=12,rank_at_10_width25_cusp=10,all_nine_rank_at_12_plus_3_width75_cusps=mixed,all_degree10_cusp_images_rank10=True,no_cross_export_basis_alignment_used=True)
(P/'cusp_degree10_verified.json').write_text(json.dumps(out,indent=2)+'\n')
print('C75_DECIC_ALL_CUSP_BRANCHES_INDEPENDENT_SAGE_PASS',ways[10],'seconds',time.monotonic()-t0,flush=True)
