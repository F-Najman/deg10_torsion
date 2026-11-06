#!/usr/bin/env python3
"""Independent binary verification of the new degree-ten C51 cusp restrictions."""
from pathlib import Path
import ast,re,json,time,hashlib
H=Path(__file__).resolve().parent;t0=time.monotonic();P=H/'new_jets.json';d=json.loads(P.read_text())
assert d['level']==51 and d['degree']==10 and d['genus']==65 and d['U_dimension']==62 and d['complete']
oldpath=Path('work/degree9/cyclic/continuation_20261026/N51_geometry/cyclic_sieves/N51_frozen_jets.m');s=oldpath.read_text();assert s.count('FI_JET_DATA_READY')==1
weights=ast.literal_eval(re.search(r'ds:=(.*?);',s)[1]);assert sorted(weights)==[1]*16+[2]*8+[8]*4

def flatmatrix(flat,r,c):
 assert len(flat)==r*c and set(flat)<={0,1}
 return tuple(sum(flat[i*c+j]<<j for j in range(c))for i in range(r))
def parse(pattern):
 out={}
 for m in re.finditer(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S):
  a=m.groups();key=tuple(map(int,a[:-3]));r,c=map(int,a[-3:-1]);out[key]=(r,c,flatmatrix(ast.literal_eval(a[-1]),r,c))
 return out

def basis(vals):
 B={}
 for v in vals:
  while v:
   p=v.bit_length()-1
   if p in B:v^=B[p]
   else:B[p]=v;break
 return B

def reduce(v,B):
 while v:
  p=v.bit_length()-1
  if p not in B:return v
  v^=B[p]
 return 0

def rowtimes(v,rows):
 out=0
 while v:
  bit=v&-v;out^=rows[bit.bit_length()-1];v^=bit
 return out

def mm(A,B):return tuple(rowtimes(v,B)for v in A)
def apply(A,v):return sum(((r&v).bit_count()&1)<<i for i,r in enumerate(A))
def columns(rows,c):return tuple(sum(((r>>j)&1)<<i for i,r in enumerate(rows))for j in range(c))
def kernel(cols,n=65):
 b=basis(cols);ans=[]
 for free in range(n):
  if free in b:continue
  v=1<<free
  for p in sorted(b):
   if(v&b[p]).bit_count()&1:v^=1<<p
  assert all((v&c).bit_count()%2==0 for c in cols);ans.append(v)
 return tuple(ans)
def samespan(A,B):
 a=basis(A);b=basis(B);return len(a)==len(b)and all(reduce(v,b)==0 for v in a.values())
Cold=parse('CartierMatrix')[()][2];Uold=parse('ImageBasis')[()][2];Aold=parse(r'basisjets\[(\d+)\]')
C=flatmatrix(d['Cartier'],65,65);U=flatmatrix(d['U_basis'],62,65);T=flatmatrix(d['T_old_to_new'],65,65)
assert len(basis(T))==65 and len(basis(U))==62 and len(basis(Uold))==62
assert samespan([v^(1<<i)for i,v in enumerate(C)],U) and samespan([v^(1<<i)for i,v in enumerate(Cold)],Uold)
assert mm(Cold,T)==mm(T,C) and samespan(mm(Uold,T),U)
AR=flatmatrix(d['rational_cap10'],65,10);AQ=flatmatrix(d['quadratic_cap5'],65,10);rcols=columns(AR,10);qcols=columns(AQ,10)
oldcols=columns(Aold[1,][2],8);Wold=kernel(oldcols);Wnew=kernel(rcols[:8]);anchor=flatmatrix(d['anchor_kernel_new'],57,65)
assert samespan(mm(Wold,T),Wnew) and samespan(anchor,Wnew)
# Check uniqueness of the basis transition, independently from the native solver.
constraints=[]
for i in range(65):
 for j in range(65):
  row=0
  for k in range(65):
   if(Cold[i]>>k)&1:row^=1<<(65*k+j)
   if(C[k]>>j)&1:row^=1<<(65*i+k)
  constraints.append(row)
for w in Wold:
 for col in rcols[:8]:
  row=0
  for i in range(65):
   if(w>>i)&1:row^=col<<(65*i)
  constraints.append(row)
rankconstraints=len(basis(constraints));assert rankconstraints==4224
flatT=sum(v<<(65*i)for i,v in enumerate(T));assert all((r&flatT).bit_count()%2==0 for r in constraints)
ratindices=[i+1 for i,w in enumerate(weights)if w==1];assert ratindices==d['rational_old_indices'] and ratindices[0]==1
ratnew=[flatmatrix(a,65,1)for a in d['rational_first_jets']];assert len(ratnew)==16
for i,A in zip(ratindices,ratnew):assert apply(T,columns(A,1)[0])==columns(Aold[i,][2],8)[0]
assert columns(ratnew[0],1)[0]==rcols[0]

def project(v):return sum(((u&v).bit_count()&1)<<i for i,u in enumerate(U))
pr=tuple(project(c)for c in rcols);pq=tuple(project(c)for c in qcols)
assert len(basis(pr))==10 and len(basis(pq))==10
ranks=[]
for k in range(1,16):
 rank=len(basis(pr[:9]+(project(columns(ratnew[k],1)[0]),)));ranks.append(rank);assert rank==10
oldcheck=json.loads((H/'old_jets_degree10_replay.json').read_text());assert oldcheck['complete_enumeration']and oldcheck['tested_divisors']==9655880 and oldcheck['exceptions_count']==0
orbits=json.loads((H/'cusp_orbits_verified.json').read_text());assert orbits['status']=='PASS'and orbits['rational_cusps_transitive']and orbits['quadratic_cusps_transitive']
# Independent full generating coefficient and disjoint high-multiplicity count.
ways=[0]*11;ways[0]=1
for w in weights:
 for e in range(w,11):ways[e]+=ways[e-w]
assert ways[10]==9655880+16*15+16+8==9656144
result={'N':51,'degree':10,'status':'PASS','rank_all':10,'old_cuspidal_divisors':9655880,'new_high_multiplicity_divisors':264,'all_cuspidal_divisors':ways[10],'new_native_reduced_cases':17,'rational_cap10_rank':10,'quadratic_cap5_rank':10,'nine_plus_one_ranks':ranks,'basis_transition_invertible':True,'basis_transition_dimension':4225-rankconstraints,'Cartier_intertwiner_verified':True,'U_transition_verified':True,'anchor_degree':8,'old_jets_sha256':hashlib.sha256(oldpath.read_bytes()).hexdigest(),'new_jets_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0}
(H/'new_jets_verified.json').write_text(json.dumps(result,indent=2)+'\n');print('DECIC51_NEW_JETS_REPLAY_COMPLETE',json.dumps(result,sort_keys=True))
