"""Independent binary arithmetic for the degree-ten C91 relative fibres."""
from pathlib import Path
import re,json,hashlib
H=Path(__file__).resolve().parent

def values(s):return list(map(int,re.findall(r'\b[01]\b',s)))
def rows(s,n,m):
 a=values(s);assert len(a)==n*m
 return [a[i*m:(i+1)*m] for i in range(n)]
def rank(a):
 p={}
 for row in a:
  v=sum(b<<i for i,b in enumerate(row))
  while v:
   i=v.bit_length()-1
   if i in p:v^=p[i]
   else:p[i]=v;break
 return len(p)
s=(H/'tate7_differential_matrices.m').read_text()
ubmatches=list(re.finditer(r'UB:=Matrix\(GF\(2\),(\d+),23,\[(.*?)\]\);',s,re.S));assert ubmatches
us=[rows(x.group(2),int(x.group(1)),23)for x in ubmatches];assert all(u==us[0]for u in us)
u=us[0];assert len(u)==rank(u)==14
s=(H/'tate7_cusp_blocks.m').read_text();ms={}
for x in re.finditer(r'localjets\[<(\d+),(\d+)>\]:=Matrix\(GF\(2\),(\d+),(\d+),\[(.*?)\]\);',s,re.S):
 i,m,n,d=map(int,x.group(1,2,3,4));ms[i,m]=rows(x.group(5),n,d)
assert len(ms)==33
cusp=[]
for a in range(11):
 for b in range(11-a):
  for c in range((10-a-b)//3+1):
   d=10-a-b-3*c;v=[a,b,c,d];rr=[[]for _ in range(14)]
   for i,m in enumerate(v,1):
    if m:
     for j,row in enumerate(ms[i,m]):rr[j]+=row
   r=rank(rr);k=min(a,b,d)+c;assert r==10-k
   cusp.append(dict(vector=v,rank=r,degree_downstairs=k,fixed_cusp_degree=10-3*k))
assert len(cusp)==120 and min(c['fixed_cusp_degree']for c in cusp)==1
s=(H/'tate7_CM_blocks.m').read_text();cm=[]
for x in re.finditer(r'Append\(~CMRecords,<([12]),([123]),([123]),8,Matrix\(GF\(2\),(\d+),23,\[(.*?)\]\)>\);',s,re.S):
 j,a,b,n=map(int,x.group(1,2,3,4));w=rows(x.group(5),n,23);assert rank(w)==n
 r=rank(u+w)-rank(w);assert r==8
 cm.append(dict(CM=j,ratA=a,ratB=b,rank=r,degree_downstairs=2,fixed_cusp_degree=2))
assert len(cm)==12 and len({(x['CM'],x['ratA'],x['ratB'])for x in cm})==12
record=dict(status='PASS',method='Independent plain-Python Gaussian elimination over F2',cusp_patterns=cusp,CM_patterns=cm,
 inputs={q:hashlib.sha256((H/q).read_bytes()).hexdigest()for q in ['tate7_differential_matrices.m','tate7_cusp_blocks.m','tate7_CM_blocks.m']})
(H/'independent_arithmetic_replay.json').write_text(json.dumps(record,indent=2)+'\n')
print('C91_INDEPENDENT_ARITHMETIC_PASS cusp120 CM12 all_relative_ranks_exact=true')
