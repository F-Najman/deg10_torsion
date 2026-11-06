LEVEL=85;CASEPATH='work/deg10/continuation_20261030c/C85'
from pathlib import Path
import ast,re,json,time
from sage.all import GF,Matrix,vector,PolynomialRing,PowerSeriesRing
N=LEVEL;P=Path(CASEPATH);t0=time.monotonic();F2=GF(2)
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(F2,int(m[1]),int(m[2]),lit(m[3]))
s=(P/'inputs'/('canonical_unit_stage1_data.m'if N==85 else'canonical_unit_data.m')).read_text();Q=mat(s,'qcoeffs');g=Q.nrows();D2=mat(s,'Diamond2Matrix'if N==85 else'Diamond1Matrix');DD=mat(s,'Diamond3Matrix'if N==85 else'Diamond2Matrix');I=Matrix.identity(F2,g)
def hecke(n,D):
 l=(Q.ncols()-1)//n+1;A=Q[:,:l];assert A.rank()==g
 B=Matrix(F2,[[Q[i,n*k]+(n%2)*int(D[i,:]*Q[:,k//n])[0]if False else 0 for k in range(l)]for i in []]) if False else Matrix(F2,g,l)
 DQ=D*Q
 for i in range(g):
  for k in range(l):B[i,k]=Q[i,n*k]+(DQ[i,k//n]if k%n==0 and n%2 else 0)
 T=A.solve_left(B);assert T*A==B;return T
T2=hecke(2,D2)
if N==85:
 T3=hecke(3,DD);K3=T3**4;K2=T2**4+I;W=Matrix(F2,list((K3*(T2-I)).rows())+list((K2*(T2-I)).rows())).row_space().basis_matrix()
 Wq=Matrix(F2,list((K3*(T3-I-DD)).rows())+list((K2*(T3-I-DD)).rows())).row_space().basis_matrix()
 group=[D2**a*DD**b for b in range(16)for a in range(2)];arcspec=[('q4',4,2,8),('q8width5',8,1,4),('q8width1',8,1,4)]
else:
 W=Matrix(F2,list(((T2-I)*(D2-I)).rows())+list(((T2-I)*(DD-I)).rows())).row_space().basis_matrix();Wq=W
 group=[D2**a*DD**b for b in range(3)for a in range(10)];arcspec=[('q5',5,2,6)]
assert (W*D2).row_space()==W.row_space() and (W*DD).row_space()==W.row_space()
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
for name,deg,cap,count in arcspec:
 a=lit((P/f'cusp_arc_{name}.m').read_text());pol=a[1];R=PolynomialRing(F2,'a');K=GF(2**deg,'a',modulus=R(pol));PS=PowerSeriesRing(K,'u',default_prec=cap);u=PS.gen();Z=[PS([sum(K(c)*K.gen()**j for j,c in enumerate(v))for v in co]).add_bigoh(cap)for co in a[6]];assert len(Z)==g
 normZ=normarc(Z)
 if N==85:
  labels=[(pow(2,-aa,85)*pow(3,-bb,85))%85 for bb in range(16)for aa in range(2)];physical={'q4':52,'q8width5':36,'q8width1':2}[name]
 else:labels=[(pow(2,-aa,93)*pow(5,-bb,93))%93 for bb in range(3)for aa in range(10)];physical=64
 mi=next(j for j,v in enumerate(labels)if v in [physical,N-physical]);GM=group[mi]
 transported=normarc([sum((normZ[j]for j in range(g)if GM[i,j]),PS(0))for i in range(g)])
 frobenius=normarc([PS([v[j]**2 for j in range(cap)]).add_bigoh(cap)for v in normZ]);assert transported==frobenius,(N,name,'actual Frobenius action')
 print('CONTC_CUSP_FROBENIUS_TRANSPORT_PASS',N,name,physical,flush=True)
 cusps={}
 for M in group:
  A=[sum((Z[j]for j in range(g)if M[i,j]),PS(0))for i in range(g)];A=normarc(A);ky=key(A)
  if ky in cusps:continue
  values=[sum((A[j]for j in range(g)if W[i,j]),PS(0))for i in range(W.nrows())];cols=[]
  for n in range(cap):
   for component in range(deg):cols.append(bits([v[n]._vector_()[component]for v in values]))
  cusps[ky]=cols
 assert len(cusps)==count,(name,len(cusps),count)
 result['nonrational'][name]=list(cusps.values());print('CONTC_CUSP_ARCS_TRANSPORTED',N,name,count,flush=True)
# Report normalized rational-prefix relation dimensions for both targets.
for label in ['local','global']:
 B=result['rational'][label];rows=B[0][0][:10]+sum((v[0][:5]for v in B[1:]),[]);ambient=W.nrows()if label=='local'else Wq.nrows();M=Matrix(F2,[[int(row>>j&1)for j in range(ambient)]for row in rows]);print('CONTC_RATIONAL_PREFIX_RANK',N,label,len(rows),M.rank(),M.left_kernel().dimension(),flush=True)
assert result==json.loads((P/'cusp_blocks.json').read_text())
ops={'level':N,'mu_local_killer':'T2-I','W':[hex(bits(row))for row in W.rows()],'Wq':[hex(bits(row))for row in Wq.rows()],'T2':[hex(bits(row))for row in T2.rows()],'D2':[hex(bits(row))for row in D2.rows()],'Dother':[hex(bits(row))for row in DD.rows()],'canonical_cusp_frobenius_verified':True,'diamond_stability_verified':True,'all_blocks_recomputed_equal':True}
(P/'operator_and_frobenius_control.json').write_text(json.dumps(ops,indent=2)+'\n');print('CONTC_CUSP_BLOCKS_COMPLETE',N,time.monotonic()-t0,flush=True)
