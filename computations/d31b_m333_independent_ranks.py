from sage.all import *
from pathlib import Path
from collections import Counter
import ast,re,json,itertools,time,random
P=Path('work/deg10/continuation_20261031b/M3N33');T=time.monotonic()
def rd(n):
 s=(P/n).read_text();s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
a=rd('analytic_cusp_jets_v4.m');rows,wits,perm,idx=a;frame=rd('relative_kummer_frame.m');k=GF(4,'w',modulus=PolynomialRing(GF(2),'x')([1,1,1]));w=k.gen();dc=[k(0),k(1),w,w+1]
CM=matrix(k,81,81,[dc[c]for c in frame[7][0]]);TM=matrix(k,81,81,[dc[c]for c in frame[7][1]]);U=matrix(k,36,81,[dc[c]for c in frame[7][2]]);I=identity_matrix(k,81);C2=CM.apply_map(lambda c:c*c)*CM
assert TM*TM==I and C2*TM==TM*C2 and (TM-I).rank()==40 and U.rank()==36
assert U.row_space()==((C2-I)*(TM-I)).row_space()
Js=[matrix(k,r[7]) for r in []] # explicit encoded conversion below
Js=[matrix(k,81,5,[dc[c]for v in r[7]for c in v])for r in rows]
Us=[matrix(k,36,5,[dc[c]for v in r[8]for c in v])for r in rows]
assert all(U*J==R and J.rank()==5 for J,R in zip(Js,Us))
spaces=[J[:,:r[2]].transpose().row_space()for J,r in zip(Js,rows)]
assert all(spaces[i]!=spaces[j]for i in range(48)for j in range(i))
for i,(J,r)in enumerate(zip(Js,rows)):
 target=(TM*J[:,:r[2]]).transpose().row_space();hits=[j for j,S in enumerate(spaces)if S==target]
 assert hits==[perm[i]-1]and hits[0]!=i and perm[hits[0]]==i+1
rat=[i for i,r in enumerate(rows)if r[2]==1];quin=[i for i,r in enumerate(rows)if r[2]==5]
assert len(rat)==40 and len(quin)==8 and all(Us[i].rank()==5 for i in quin)
pos={j:i for i,j in enumerate(rat)};hp=[pos[perm[i]-1]for i in rat]
def bitcol(codes):
 lo=hi=0
 for i,c in enumerate(codes):
  if c&1:lo|=1<<i
  if c&2:hi|=1<<i
 return lo,hi
def mul(v,c):
 a,b=v
 return (a,b)if c==1 else ((b,a^b)if c==2 else (a^b,a))
def grank(vs):
 piv={}
 for aa,bb in vs:
  while aa|bb:
   z=(aa|bb)&-(aa|bb);p=z.bit_length()-1;c=((aa>>p)&1)+2*((bb>>p)&1)
   if p in piv:
    a,b=mul(piv[p],c);aa^=a;bb^=b
   else:
    piv[p]=mul((aa,bb),[0,1,3,2][c]);break
 return len(piv)
cols=[[bitcol([rows[i][8][r][j]for r in range(36)])for j in range(5)]for i in rat]
rng=random.Random(33320261031)
for _ in range(2000):
 ids=sorted(rng.randrange(40)for z in range(5));co=Counter(ids);vs=[];sm=[]
 for i,n in co.items():vs+=cols[i][:n];sm += [Us[rat[i]].column(j)for j in range(n)]
 assert grank(vs)==matrix(k,sm).rank()
print('M333_INDEPENDENT_FRAME_H_AND_GF4_ARITHMETIC_PASS',flush=True)
hist=[0,0,0];tested=0
for ids in itertools.combinations_with_replacement(range(40),5):
 co=Counter(ids);kk=sum(min(n,co.get(hp[i],0))for i,n in co.items()if i<hp[i]);assert 0<=kk<=2
 vs=[]
 for i,n in co.items():vs.extend(cols[i][:n])
 rank=grank(vs);assert rank==5-kk,(ids,rank,kk);hist[kk]+=1;tested+=1
assert tested==1086008 and hist==[864008,213600,8400]
report={'status':'PASS','frame_Cartier_square_and_h_matrices_independently_recomputed':True,'h_cotangent_rank':40,'U_dimension':36,'all48_full_canonical_first_column_spaces_distinct':True,'all48_h_partners_uniquely_reconstructed':True,'h_free_involution':True,'all48_saved_UJ_products_verified':True,'Sage_matrix_rank_crosschecks':2000,'independent_GF4_integer_gaussian_ranks':tested,'rational_cusp_divisors':tested,'rank_distribution_by_pulledback_degree':hist,'quintic_cusp_divisors':8,'quintic_ranks':[5]*8,'bad':0,'elapsed_seconds':time.monotonic()-T}
(P/'independent_ranks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print('D31B_M333_INDEPENDENT_ALL_RANKS_COMPLETE',time.monotonic()-T,flush=True)
