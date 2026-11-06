"""Rebuild C65 jets and check all prefixes using independent full-rank elimination."""
from pathlib import Path
import ast,re,json,time,subprocess,hashlib
from sage.all import GF,PolynomialRing,Matrix,vector
H=Path(__file__).resolve().parents[1];P=H/'continuation_20261030c/C65';t0=time.time();k=GF(2);g=121
parse=lambda s:ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def read_matrix(txt,name):
 m=re.search(r'\b'+name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',txt,re.S);assert m,name
 return Matrix(k,int(m[1]),int(m[2]),parse(m[3]))
s=(P/'inputs/canonical_source_data.m').read_text();q=read_matrix(s,'qcoeffs');A=read_matrix(s,'Diamond1Matrix');B=read_matrix(s,'Diamond2Matrix')
w=read_matrix((P/'inputs/annihilator_data.m').read_text(),'AnnihilatingDifferentials');assert w.nrows()==106 and w.rank()==106
maps=[A**a*B**b for a in range(2)for b in range(12)];assert len({tuple(m.list())for m in maps})==24
width,pol,prec,pivot,inds,sidx,co,xc,yc=parse((P/'inputs/quartic_cusp_arc_data.txt').read_text());assert width==5 and prec==800
F=GF(16,'u',modulus=PolynomialRing(k,'z')(pol));u=F.gen()
base=Matrix(F,[[sum(F(c)*u**i for i,c in enumerate(row[j]))for j in range(2)]for row in co]);assert base.nrows()==g
assert base[pivot-1,0]==1 and base[pivot-1,1]==0 and any(base[:,1].list())
def projective_orbit_key(v):
 v=list(v);c=next(x for x in v if x);v=[x/c for x in v]
 return min(tuple(int(a)for x in v for a in (x**(2**i))._vector_())for i in range(4))
representatives={}
for m in maps:
 v=m.change_ring(F)*base;representatives.setdefault(projective_orbit_key(v.column(0)),v)
assert len(representatives)==6
packed=[]
def encode(v):
 b=sum(int(c)<<i for i,c in enumerate(v));return [b&((1<<64)-1),b>>64]
for i,m in enumerate(maps):
 rr=w*m*q.matrix_from_columns(range(1,11 if i==0 else 6));packed.append([1,10 if i==0 else 5,[encode(rr.column(j))for j in range(rr.ncols())]])
# Coordinate functionals form another basis of Hom_F2(F16,F2), independent
# of the primary implementation's four absolute-trace functionals.
for arc in representatives.values():
 ev=w.change_ring(F)*arc;vec=[]
 for j in range(2):
  for component in range(4):vec.append(encode([v._vector_()[component]for v in ev.column(j)]))
 packed.append([4,2,vec])
assert len(packed)==30
lines=['30 106']
for cost,cap,vs in packed:
 lines.append(f'{cost} {cap} {len(vs)}');lines.extend(f'{lo} {hi}'for lo,hi in vs)
vecfile=P/'parent_prefix_vectors.txt';vecfile.write_text('\n'.join(lines)+'\n');cpp=H/'computations/contc65_parent_prefix_native.cpp';exe=P/'parent_prefix_native'
subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
r=subprocess.run([str(exe),str(vecfile)],capture_output=True,text=True,check=True);print(r.stdout,flush=True)
counts=[list(map(int,v.split()))for v in r.stdout.splitlines()if re.match(r'^\d+ ',v)];assert len(counts)==10
for n0,ct in counts:
 dp=[1]+[0]*10
 for cost,cap,vs in packed[1:]:
  lim=min(n0,5)if cost==1 else cap;nd=[0]*11
  for j,v in enumerate(dp):
   for m in range(lim+1):
    if j+cost*m<11:nd[j+cost*m]+=v
  dp=nd
 assert ct==dp[10-n0],(n0,ct,dp[10-n0])
assert sum(v for _,v in counts)==7170048
out={'status':'PASS','case':'C65','method':'Independent Sage jet vectors using F16 coordinate functionals and compiled rank-at-leaf Gaussian elimination; no incremental matroid state reused','total_configurations':7170048,'counts_by_max_rational_multiplicity':counts,'all_prefix_matrices_have_rank10':True,'cpp_sha256':hashlib.sha256(cpp.read_bytes()).hexdigest(),'vector_data_sha256':hashlib.sha256(vecfile.read_bytes()).hexdigest(),'seconds':time.time()-t0}
(P/'parent_prefix_replay.json').write_text(json.dumps(out,indent=2)+'\n');print('C65_PARENT_PREFIX_REPLAY_COMPLETE',7170048,round(time.time()-t0,2),flush=True)
