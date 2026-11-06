LEVEL=91;CASEPATH='work/deg10/continuation_20261030c/C91/source_backup'
from pathlib import Path
import ast,re,json,time
from sage.all import GF,Matrix,vector,PolynomialRing,PowerSeriesRing
N=LEVEL;P=Path(CASEPATH);t0=time.monotonic();F2=GF(2)
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(F2,int(m[1]),int(m[2]),lit(m[3]))
s=(P/'inputs'/'canonical_unit_stage1_data.m').read_text();Q=mat(s,'qcoeffs');g=Q.nrows();D2=mat(s,'Diamond1Matrix');DD=mat(s,'Diamond2Matrix');I=Matrix.identity(F2,g)
def hecke(n,D):
 l=(Q.ncols()-1)//n+1;A=Q[:,:l];assert A.rank()==g
 B=Matrix(F2,[[Q[i,n*k]+(n%2)*int(D[i,:]*Q[:,k//n])[0]if False else 0 for k in range(l)]for i in []]) if False else Matrix(F2,g,l)
 DQ=D*Q
 for i in range(g):
  for k in range(l):B[i,k]=Q[i,n*k]+(DQ[i,k//n]if k%n==0 and n%2 else 0)
 T=A.solve_left(B);assert T*A==B;return T
T2=hecke(2,D2)
W=((T2-I)*(DD-I)).row_space().basis_matrix();Wq=W
group=[D2**a*DD**b for b in range(3)for a in range(12)];assert len({tuple(m.list())for m in group})==36
arcspec=[('q3',3,3,12)]
print('CONTC_CUSP_OPERATOR_RANK',N,W.nrows(),Wq.nrows(),flush=True)
def bits(v):return sum(int(x)<<i for i,x in enumerate(v))
def normarc(A):
 k=next(i for i in range(g)if A[i][0]);h=A[k];return [v/h for v in A]
def key(A):
 zz=[v[0]for v in A];h=next(v for v in zz if v);zz=[v/h for v in zz];d=zz[0].parent().degree();return min(tuple(bits(v**(2**j))for v in zz)for j in range(d))
result={'N':N,'g':g,'Wrank':W.nrows(),'Wqrank':Wq.nrows(),'rational':{},'nonrational':{}}
for label,BM in [('local',W),('global',Wq)]:
 rows=[]
 for k,M in enumerate(group):
  A=BM*M*Q[:,1:11];rows.append([[bits(A.column(j))for j in range(10)],k])
 result['rational'][label]=rows

rows=result['rational']['local'][0][0][:10]+sum((v[0][:5]for v in result['rational']['local'][1:]),[]);M=Matrix(F2,[[row>>j&1 for j in range(W.nrows())]for row in rows]);print('C91_SOURCE_PREFLIGHT_COMPLETE',W.nrows(),len(rows),M.rank(),M.left_kernel().dimension(),flush=True)
(P/'rational_preflight.json').write_text(json.dumps(result,indent=2)+'\n')
