from sage.all import *
from pathlib import Path
from collections import Counter
import json,time
T=time.monotonic();P=Path('work/deg10/continuation_20261031b/C60')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
a=rd('target30_proven_picard.m');points,candidates,fibres=rd('mixed_fibres.json');records={}
for s in range(3):
 aa=rd(f'three_maps_s{s}_v2.m')
 for r in aa[0]:
  ci,d,rows=r;assert ci not in records and d==points[ci-1][0]
  for v,ix,cl,e,n in rows:
   im=aa[1][ix-1];assert im[0][0]==e==im[2]and n*e==d and cl==[n*z for z in im[1]]
  records[ci]=r
assert sorted(records)==list(range(1,len(points)+1))
n=len(a[5][0]);L=matrix(ZZ,a[5]+a[12]+[[0]*(n-1)+[1]]).row_module();B=L.basis_matrix();M=abs(B.det());assert M==4352;BI=B.inverse()
def q(v):return tuple(int(M*z)%int(M)for z in vector(ZZ,v)*BI)
allowed=[set()for _ in range(7)];allowed[0].add(q([0]*n))
for d,cl in zip(a[8],a[9]):
 for j in range(d,7):allowed[j]|={tuple((x+y)%M for x,y in zip(q(cl),v))for v in allowed[j-d]}
tests=[];left=[]
for ci,cand in enumerate(candidates,1):
 nc=sum(points[i-1][0]*mul for i,mul in cand);cd=10-nc;assert 0<=cd<=6;flags=[]
 for j in range(3):
  cl=sum((mul*vector(ZZ,records[i][2][j][2])for i,mul in cand),vector(ZZ,[0]*n));flags.append(q(-cl)in allowed[cd])
 tests.append([ci,cand,cd,flags])
 if all(flags):left.append(tests[-1])
out={'status':'PASS','source_points':len(points),'source_patterns':len(candidates),'projected_maps':3*len(points),'quotient_index':int(M),'allowed_complement_sizes':list(map(len,allowed)),'remaining':left,'remaining_count':len(left),'remaining_types':[[list(k),v]for k,v in Counter(tuple((points[i-1][0],m)for i,m in r[1])for r in left).items()],'all_tests':tests,'elapsed_seconds':time.monotonic()-T}
(P/'three_maps_summary.json').write_text(json.dumps(out,indent=2)+'\n')
def mag(v):
 if isinstance(v,bool):return'true'if v else'false'
 if isinstance(v,(int,Integer)):return str(v)
 a=all(isinstance(z,(int,Integer))and not isinstance(z,bool)for z in v);return('['if a else'[*')+','.join(mag(z)for z in v)+(']'if a else'*]')
(P/'three_maps_remaining.m').write_text(mag([points,left])+'\n');print('D31B_C60_THREE_MAP_MERGE_COMPLETE',len(points),len(candidates),len(left),out['remaining_types'],out['allowed_complement_sizes'],time.monotonic()-T,flush=True)
