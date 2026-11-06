#!/usr/bin/env python3
"""Independent binary replay of the two new degree-ten C49 Cartier branches."""
from pathlib import Path
from itertools import combinations
import ast,re,json,time,hashlib,math
H=Path(__file__).resolve().parent;P=H/'c49_fast_cartier_jets.m';s=P.read_text();t0=time.monotonic()
assert s.count('DECIC_C49_FAST_MIXED10_JET_DATA_READY')==1
assert re.search(r'\bN:=49;D:=10;gx:=69;udim:=68;',s)
weights=ast.literal_eval(re.search(r'ds:=(.*?);',s)[1]);caps=ast.literal_eval(re.search(r'caps:=(.*?);',s)[1])
assert weights==[1]*21+[3]*6 and caps==[7]+[1]*26

def binary_basis(vals):
 b={}
 for v in vals:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return b

def reduce(v,b):
 while v:
  p=v.bit_length()-1
  if p not in b:return v
  v^=b[p]
 return 0

def parse(pattern):
 out={}
 for match in re.finditer(pattern+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S):
  a=match.groups();key=tuple(map(int,a[:-3]));r,c=map(int,a[-3:-1]);v=ast.literal_eval(a[-1]);assert len(v)==r*c and set(v)<={0,1} and key not in out
  rows=tuple(sum(v[i*c+j]<<j for j in range(c))for i in range(r));cols=tuple(sum(v[i*c+j]<<i for i in range(r))for j in range(c));out[key]=(r,c,rows,cols)
 return out
CM=parse(r'CartierMatrix')[()];U=parse(r'ImageBasis')[()]
assert CM[:2]==(69,69)and U[:2]==(68,69)
cb=binary_basis([v^(1<<i)for i,v in enumerate(CM[2])]);ub=binary_basis(U[2]);assert len(cb)==len(ub)==68 and all(reduce(v,ub)==0 for v in cb.values())
A=parse(r'basisjets\[(\d+)\]');J=parse(r'localjets\[<(\d+),(\d+)>\]');assert len(A)==27 and len(J)==33
jets={}
for i,w in enumerate(weights):
 a=A[i+1,];assert a[:2]==(69,w*caps[i])
 expected=tuple(sum(((v&col).bit_count()&1)<<k for k,v in enumerate(U[2]))for col in a[3])
 for m in range(1,caps[i]+1):
  j=J[i+1,m];assert j[:2]==(68,w*m)and j[3]==expected[:w*m]
 jets[i]=expected
exceptions=[];distinct_count=0
for support in combinations(range(21),7):
 ratcols=tuple(jets[i][0]for i in support)
 for j in range(21,27):
  rank=len(binary_basis(ratcols+jets[j]));distinct_count+=1
  if rank!=10:exceptions.append({'case':'seven distinct rational cusps plus one cubic cusp','rational_support':support,'cubic_cusp':j,'rank':rank})
assert distinct_count==6*math.comb(21,7)==697680
high_count=0
for j in range(21,27):
 rank=len(binary_basis(jets[0]+jets[j]));high_count+=1
 if rank!=10:exceptions.append({'case':'seven times fixed rational cusp plus one cubic cusp','cubic_cusp':j,'rank':rank})
assert high_count==6
result={'level':49,'degree':10,'source_genus':69,'Cartier_image_dimension':68,'one_common_differential_basis':True,'rational_cusps':21,'cubic_cusps':6,'seven_distinct_rational_support_cases':distinct_count,'sevenfold_fixed_rational_cusp_cases':high_count,'exceptions':exceptions,'exceptions_count':len(exceptions),'rank_all':10 if not exceptions else None,'geometry_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0,'status':'PASS'if not exceptions else'FAIL'}
(H/'c49_fast_jets_verified.json').write_text(json.dumps(result,indent=2)+'\n')
print('DECIC_C49_FAST_CARTIER_REPLAY_COMPLETE',json.dumps(result,sort_keys=True));assert not exceptions
