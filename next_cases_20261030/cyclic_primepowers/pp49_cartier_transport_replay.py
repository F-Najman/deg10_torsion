#!/usr/bin/env python3
"""Exact Cartier-cyclic recovery of the diamond action and remaining C49 ranks."""
from pathlib import Path
from itertools import combinations
import ast,re,json,time,hashlib,math
H=Path(__file__).resolve().parent; t0=time.monotonic()
oldpath=Path('work/degree9/prime_orders/continued_20261026/P49/cartier_mixed9_jets.m')
partialpath=H/'c49_fast_cartier_jets.m'; permpath=H/'c49_cusp_permutation.m'
old=oldpath.read_text(); partial=partialpath.read_text(); ps=permpath.read_text()
assert old.count('NONIC49_MIXED9_JET_DATA_READY')==1
assert ps.count('DECIC_C49_CUSP_PERMUTATION_DATA_READY')==1
status=json.loads(Path('work/deg10/jobs/pp49_fast_cartier_export.status.json').read_text())
assert status['state']=='EXITED' and status['actual_exit_code']==-15
# The full exporter timed out. Only its completed first block is used below,
# and this is verified independently by the native Completion control.
log=Path('work/deg10/jobs/pp49_fast_cartier_export.log').read_text()
assert 'DECIC_C49_FAST_CARTIER_JET_DONE i=1 ' in log

def parse(s,pattern):
 out={}
 for match in re.finditer(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S):
  a=match.groups();key=tuple(map(int,a[:-3]));r,c=map(int,a[-3:-1]);v=ast.literal_eval(a[-1]);assert len(v)==r*c and set(v)<={0,1} and key not in out
  rows=tuple(sum(v[i*c+j]<<j for j in range(c))for i in range(r));cols=tuple(sum(v[i*c+j]<<i for i in range(r))for j in range(c));out[key]=(r,c,rows,cols)
 return out
C=parse(old,r'CartierMatrix')[()];U=parse(old,r'ImageBasis')[()]
assert C[:2]==(69,69) and U[:2]==(68,69)
A=parse(old,r'basisjets\[(\d+)\]');newA=parse(partial,r'basisjets\[(\d+)\]');newJ=parse(partial,r'localjets\[<(\d+),(\d+)>\]')
assert len(A)==24 and parse(partial,r'CartierMatrix')[()]==C and parse(partial,r'ImageBasis')[()]==U

def apply(rows,v):return sum(((row&v).bit_count()&1)<<i for i,row in enumerate(rows))
def binary_basis(vals):
 b={}
 for v in vals:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return b

def row_times(v,rows):
 ans=0
 while v:
  low=v&-v;ans^=rows[low.bit_length()-1];v^=low
 return ans
# At an F2-rational cusp the nonzero differential-line trivialization is unique.
# The identical Cartier matrices and these common canonical first jets identify
# the two native differential bases: their Cartier-cyclic span is the full space.
frame=[]
for i in range(1,9):
 assert A[i,][3][0]==newA[i,][3][0]
 v=A[i,][3][0]
 for k in range(69):frame.append(v);v=apply(C[2],v)
assert len(binary_basis(frame))==69
for i in range(1,min(18,len(newA))+1):assert A[i,][3][0]==newA[i,][3][0]
assert newA[1,][:2]==(69,7) and newA[1,][3][:6]==A[1,][3]

def listvalue(name):return ast.literal_eval(re.search(name+r':=(.*?);',ps,re.S)[1])
perm=[i-1 for i in listvalue('diamond2_permutation')];sigs=listvalue('rational_cusp_signatures');images=listvalue('diamond2_image_signatures')
assert sorted(perm)==list(range(21)) and len(set(map(tuple,sigs)))==21
assert all(images[i]==sigs[perm[i]]for i in range(21))
cycle=[];i=0
while i not in cycle:cycle.append(i);i=perm[i]
assert i==0 and len(cycle)==21 and [i+1 for i in cycle]==listvalue('diamond2_cycle')
known=[A[i+1,][3][0] for i in range(18)]
# The actual diamond pullback commutes with Cartier and sends the canonical
# evaluator at P to the evaluator at diamond2(P). These known pairs determine
# it on the full differential dual; all dependencies must be respected.
intertwiner={};pair_count=0;relations=0
for i in range(18):
 if perm[i]>=18:continue
 pair_count+=1;v=known[i];w=known[perm[i]]
 for k in range(69):
  vv,ww=v,w
  while vv:
   p=vv.bit_length()-1
   if p in intertwiner:vv^=intertwiner[p][0];ww^=intertwiner[p][1]
   else:intertwiner[p]=(vv,ww);break
  if not vv:assert ww==0;relations+=1
  v=apply(C[2],v);w=apply(C[2],w)
assert len(intertwiner)==69

def diamond(v):
 out=0
 while v:
  p=v.bit_length()-1;vv,ww=intertwiner[p];v^=vv;out^=ww
 return out
Mcols=tuple(diamond(1<<i)for i in range(69));Mrows=tuple(sum(((v>>i)&1)<<j for j,v in enumerate(Mcols))for i in range(69))
assert len(binary_basis(Mcols))==69
for i in range(69):assert diamond(apply(C[2],1<<i))==apply(C[2],diamond(1<<i))
for n in [3,7,21]:
 cols=[]
 for i in range(69):
  v=1<<i
  for k in range(n):v=diamond(v)
  cols.append(v)
 assert (tuple(cols)==tuple(1<<i for i in range(69)))==(n==21)
# Check the whole old point set, including indices whose image was missing.
ratcols=[None]*21;v=known[0]
for i in cycle:ratcols[i]=v;v=diamond(v)
assert v==known[0] and len(set(ratcols))==21 and ratcols[:18]==known
assert all(diamond(ratcols[i])==ratcols[perm[i]]for i in range(21))
# Invariance of the full Cartier image on row vectors.
ub=binary_basis(U[2])
def reduce(v,b):
 while v:
  p=v.bit_length()-1
  if p not in b:return v
  v^=b[p]
 return 0
assert all(reduce(row_times(u,Mrows),ub)==0 for u in U[2])

def project(col):return sum(((row&col).bit_count()&1)<<i for i,row in enumerate(U[2]))
rat=[project(col)for col in ratcols]
cubic=[tuple(project(col)for col in A[i,][3][:3])for i in range(19,25)]
base7=tuple(project(col)for col in newA[1,][3]);assert len(base7)==7
assert newJ[1,7][:2]==(68,7) and newJ[1,7][3]==base7
exceptions=[];distinct_count=0
for support in combinations(range(21),7):
 rc=tuple(rat[i]for i in support)
 for j,cc in enumerate(cubic):
  rank=len(binary_basis(rc+cc));distinct_count+=1
  if rank!=10:exceptions.append({'case':'seven distinct rational cusps plus cubic','rational_support':support,'cubic':j,'rank':rank})
assert distinct_count==6*math.comb(21,7)==697680
for j,cc in enumerate(cubic):
 rank=len(binary_basis(base7+cc))
 if rank!=10:exceptions.append({'case':'sevenfold fixed rational cusp plus cubic','cubic':j,'rank':rank})
result={'level':49,'degree':10,'source_genus':69,'Cartier_image_dimension':68,'basis_alignment_canonical_rational_jets':8,'basis_alignment_Cartier_cyclic_rank':69,'known_diamond_pairs':pair_count,'known_pair_Cartier_span_rank':69,'zero_dependency_relations_checked':relations,'diamond_order':21,'diamond_commutes_with_Cartier':True,'diamond_checks_all_old_rational_evaluators':18,'recovered_rational_first_jets':3,'seven_distinct_rational_support_cases':distinct_count,'sevenfold_fixed_rational_cusp_cases':6,'exceptions':exceptions,'exceptions_count':len(exceptions),'rank_all':10 if not exceptions else None,'partial_export_actual_exit':-15,'partial_export_used_as_complete_certificate':False,'saved_completed_cap7_block_independently_validated':True,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [oldpath,partialpath,permpath]},'elapsed_seconds':time.monotonic()-t0,'status':'PASS'if not exceptions else'FAIL'}
(H/'c49_cartier_transport_verified.json').write_text(json.dumps(result,indent=2)+'\n')
(H/'c49_recovered_diamond_matrix.json').write_text(json.dumps({'order':21,'nrows':69,'rows_bitpacked':Mrows,'columns_bitpacked':Mcols,'rational_first_jets_bitpacked':ratcols},indent=2)+'\n')
print('DECIC_C49_CARTIER_TRANSPORT_REPLAY_COMPLETE',json.dumps(result,sort_keys=True));assert not exceptions
