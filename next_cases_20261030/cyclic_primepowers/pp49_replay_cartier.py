#!/usr/bin/env python3
"""Independent binary proof for one complete common-basis C49 degree10 export."""
from pathlib import Path
import ast,re,json,time,hashlib,sys,math
H=Path(__file__).resolve().parent;P=H/'c49_cartier_jets.m';s=P.read_text();t0=time.monotonic()
assert s.count('DECIC_C49_TRANSPORT_MIXED10_JET_DATA_READY')==1,'native exporter incomplete'
assert re.search(r'\bN:=49;D:=10;gx:=69;udim:=68;',s)
weights=ast.literal_eval(re.search(r'ds:=(.*?);',s)[1]);caps=ast.literal_eval(re.search(r'caps:=(.*?);',s)[1])
assert weights==[1]*19+[3]*6 and caps==[7]*19+[3]*6

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
A=parse(r'basisjets\[(\d+)\]');J=parse(r'localjets\[<(\d+),(\d+)>\]');assert len(A)==25 and len(J)==151
MD=parse(r'DiamondMatrix')[()];assert MD[:2]==(69,69)
def row_times(v,rows):
 ans=0
 while v:
  low=v&-v;ans^=rows[low.bit_length()-1];v^=low
 return ans
def mm(a,b):return tuple(row_times(v,b)for v in a)
identity=tuple(1<<i for i in range(69));powers=[identity]
for k in range(21):powers.append(mm(powers[-1],MD[2]))
assert powers[21]==identity and powers[3]!=identity and powers[7]!=identity
assert mm(MD[2],CM[2])==mm(CM[2],MD[2])
assert all(reduce(row_times(v,MD[2]),ub)==0 for v in U[2])
# Every rational restriction is the exact pullback of the same base jet.
for k in range(19):assert A[k+1,][2]==mm(powers[k],A[1,][2])
# The six cubic blocks comprise two complete three-cycles.
for k in range(3):
 for j in range(2):assert A[20+2*k+j,][2]==mm(powers[k],A[20+j,][2])
for j in range(2):
 for m in range(1,4):
  cols=A[20+j,][3][:3*m];span=binary_basis(cols)
  transported=tuple(sum(((v&c).bit_count()&1)<<i for i,v in enumerate(powers[3]))for c in cols)
  assert len(span)==3*m and all(reduce(c,span)==0 for c in transported)
jets={}
for i,w in enumerate(weights):
 a=A[i+1,];assert a[:2]==(69,w*caps[i])
 expected=tuple(sum(((v&col).bit_count()&1)<<k for k,v in enumerate(U[2]))for col in a[3])
 for m in range(1,caps[i]+1):
  j=J[i+1,m];assert j[:2]==(68,w*m)and j[3]==expected[:w*m]
 jets[i]=expected
# Put the cubic blocks first so the pure rational branch can be discarded early.
order=list(range(19,25))+list(range(19));basis=[0]*68;vector=[0]*25;count=0;bycubic={3:0,6:0,9:0};exceptions=[]
def visit(k,remaining,rank,cubicdegree):
 global count
 if remaining==0:
  assert cubicdegree in bycubic
  count+=1;bycubic[cubicdegree]+=1
  if rank<10:exceptions.append({'rank':rank,'multiplicities':vector.copy()})
  return
 if k==25 or(k==6 and cubicdegree==0):return
 i=order[k];w=weights[i];maximum=min(caps[i],remaining//w);added=[];current=rank
 for m in range(maximum+1):
  if m:
   for col in jets[i][(m-1)*w:m*w]:
    while col:
     p=col.bit_length()-1
     if basis[p]:col^=basis[p]
     else:basis[p]=col;added.append(p);current+=1;break
  vector[i]=m
  visit(k+1,remaining-m*w,current,cubicdegree+(m*w if k<6 else 0))
 vector[i]=0
 for p in added:basis[p]=0
visit(0,10,0,0)
assert bycubic=={3:6*math.comb(25,7),6:math.comb(7,2)*math.comb(22,4),9:math.comb(8,3)*19}and count==3038879
result={'level':49,'degree':10,'source_genus':69,'Cartier_image_dimension':68,'one_common_differential_basis':True,'selected_rational_cusps':19,'omitted_rational_cusps':2,'all_cubic_cusps':6,'tested_divisors':count,'by_cubic_degree':bycubic,'exceptions':exceptions,'exceptions_count':len(exceptions),'rank_all':10 if not exceptions else None,'geometry_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0,'status':'PASS'if not exceptions else'FAIL'}
if '--readonly' in sys.argv:
 saved=json.loads((H/'c49_mixed10_verified.json').read_text());assert all(saved[k]==v for k,v in result.items()if k not in ['elapsed_seconds','by_cubic_degree']);assert saved['by_cubic_degree']=={str(k):v for k,v in bycubic.items()}
else:(H/'c49_mixed10_verified.json').write_text(json.dumps(result,indent=2)+'\n')
print('DECIC_C49_MIXED10_REPLAY_COMPLETE',json.dumps(result,sort_keys=True));assert not exceptions
