from pathlib import Path
import json,ast,re,math,collections,hashlib
from sage.all import ZZ,QQ,GF,matrix,vector,PolynomialRing,NumberField,EllipticCurve,companion_matrix
H=Path(__file__).resolve().parents[1];O=H/'exclusions/315d5'
read=lambda p:ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
a=read(O/'diagonal_data.txt');d4=read(O/'M3N15_D4_exact_cusp_lifts.txt');gg=read(O/'galois3_data.txt')
n4,n7=a[2][:-1],a[4][:-1];k=6
for inv,rel in [(a[2],a[3]),(a[4],a[5])]:assert rel==[[inv[i] if i==j else 0 for j in range(7)] for i in range(6)]
assert math.prod(n4)==5971968 and math.prod(n7)==149299200
# Galois saturation of the rational 3-torsion, independently over F3.
assert gg[0]==a[2] and gg[1]==a[3]
F3=GF(3);V=matrix(F3,gg[4]);W=matrix(F3,gg[5]);M=matrix(F3,6,6,gg[6]);perm=gg[3]
assert V.rank()==6 and W.rank()==5 and M*M==matrix.identity(F3,6)
assert V*M==matrix(F3,[gg[4][i-1] for i in perm])
assert (M-matrix.identity(F3,6)).rank()==1 and W*(M-matrix.identity(F3,6))==0
assert collections.Counter(gg[2])=={1:16,2:8} and sorted(perm)==list(range(1,25))
print('DECIC315_GALOIS_SATURATION_REPLAY_PASS rational_J_equals_C=true',flush=True)
# Check the independent modular-symbol census using Newton/Frobenius traces.
raw=(O/'frobenius_data.txt').read_text();ff=[ast.literal_eval(s.replace('[*','[').replace('*]',']')) for s in re.findall(r'\[\*.*?\*\]',raw,re.S)];R=PolynomialRing(ZZ,'x');frob=dict((p,R(cs)) for p,cs in ff)
Jorders={p:int(frob[p](1)) for p in [7,13,19]};assert math.gcd(*Jorders.values())==5971968
assert frob[2](1)*frob[2](-1)==math.prod(n4)
expected={4:[0,0,0,102,192],7:[0,24,144,396,3312]}
for q,p,e,cuspds in [(4,2,2,[1]*16+[2]*8),(7,7,1,[1]*16+[4]*4)]:
 C=companion_matrix(frob[p]);closed=[]
 for d in range(1,6):
  count=int(q**d+1-(C**(e*d)).trace())-sum(t for t in cuspds if d%t==0)
  new=(count-sum(j*closed[j-1] for j in range(1,d) if d%j==0))//d
  assert count-sum(j*closed[j-1] for j in range(1,d) if d%j==0)==d*new
  closed.append(new)
 assert closed==expected[q];print('DECIC315_INDEPENDENT_CENSUS_PASS',q,closed,flush=True)
# Integer Hermite normal form gives a coordinate-free diagonal subgroup test.
h4=[v[:-1] for v in a[6]];h7=[v[:-1] for v in a[7]];mods=n4+n7
rows=[v+w for v,w in zip(h4,h7)]+[[mods[i] if i==j else 0 for j in range(12)] for i in range(12)]
B=matrix(ZZ,rows).row_module().basis_matrix();assert B.nrows()==12;BI=B.inverse();den=int(BI.denominator());N=matrix(ZZ,BI*den);nr=[[int(x) for x in row] for row in N.rows()]
assert math.prod(mods)//abs(B.det())==1990656
key4=lambda v,scale:tuple(sum(scale*v[i]*nr[i][j] for i in range(6))%den for j in range(12))
key7=lambda v,scale:tuple(sum(-scale*v[i]*nr[6+i][j] for i in range(6))%den for j in range(12))
assert len(a[8])==24432 and len(a[9])==57072
assert collections.Counter(st[0] for st in a[8])=={0:22608,4:1632,5:192}
assert collections.Counter(st[0] for st in a[9])=={2:19584,3:19584,4:6336,5:3312,22:4800,23:3456}
for ss in [a[8],a[9]]:assert all(st[2][-1]==0 for st in ss)
found={}
for scale in [3,1]:
 S=collections.defaultdict(list)
 for i,st in enumerate(a[8],1):S[key4(st[2],scale)].append(i)
 pairs=[]
 for j,st in enumerate(a[9],1):
  for i in S.get(key7(st[2],scale),[]):pairs.append((i,j))
 found[scale]=pairs
assert set(found[3])==set(map(tuple,a[10])) and len(found[3])==4080
assert len(found[1])==240
# Generate all global degree-five cusp divisors and the 12 degree-four
# sporadic points plus every rational cusp, independently of native groups.
def weighted(ds,d):
 def it(i,r,v):
  if i==len(ds):
   if r==0:yield v
   return
  for k in range(r//ds[i]+1):yield from it(i+1,r-k*ds[i],v+[k])
 return it(0,d,[])
norm4=lambda v:tuple(int(v[i])%n4[i] for i in range(6));norm7=lambda v:tuple(int(v[i])%n7[i] for i in range(6))
image=lambda v,h,n:tuple(sum(v[i]*h[i][j] for i in range(20))%n[j] for j in range(6))
pure=[image(v,h4,n4) for v in weighted(a[1],5)];assert len(pure)==15568
known=set(pure)
for v in d4[6]:
 for i,d in enumerate(a[1]):
  if d==1:known.add(norm4([image(v,h4,n4)[j]+h4[i][j] for j in range(6)]))
assert len(known)==15760 and known=={norm4(v) for v in a[11]}
survivors={norm4(a[8][i-1][2]) for i,j in found[1]};unknown=survivors-known
assert len(survivors)==240 and len(unknown)==48 and len(survivors & known)==192
assert survivors.isdisjoint(set(pure))
# All 72640 effective degree-five divisors over F7 have distinct classes.
all7=[norm7(st[2]) for st in a[9]]+[image(v,h7,n7) for v in weighted(a[1],5)]
assert len(all7)==len(set(all7))==72640
print('DECIC315_INTEGER_DIAGONAL_REPLAY_PASS threefold=4080 exact=240 known_reducible=192 remaining=48 gonality_F7_gt5=true',flush=True)
# Check the exact sporadic decic curve and an unramified good reduction.
d=json.loads((H/'occurrences/new315_exact_candidate.json').read_text());S=PolynomialRing(QQ,'z');F=NumberField(S.gen()**2+S.gen()+1,'z');z=F.gen();T=PolynomialRing(F,'t');el=lambda c:QQ(c[0])+QQ(c[1])*z
f=T([el(c) for c in d['relative_polynomial']]);qpoly=T([el(c) for c in d['q_coefficients']]);assert f.is_irreducible();K=F.extension(f,'w');w=K.gen();q=qpoly(w)
E=EllipticCurve(K,[(1+q)*w+(2-q),0,q*w*w+(1-q)*w,0,0]);P=E([0,0]);Q=E([-w,w*w]);assert 3*P==E(0) and 15*Q==E(0) and 3*Q!=E(0) and 5*Q!=E(0) and 5*Q not in [P,-P]
jp=E.j_invariant().absolute_minpoly();assert jp.degree()==5 and any(c not in ZZ for c in jp.list())
Fp=GF(43);red=lambda a:Fp(QQ(F(a).list()[0]))+(Fp(QQ(F(a).list()[1])) if len(F(a).list())>1 else Fp(0))*Fp(6)
Rp=PolynomialRing(Fp,'t');fp=Rp([red(c) for c in f.list()]);qp=Rp([red(c) for c in qpoly.list()]);assert fp(Fp(28))==0 and fp.derivative()(Fp(28))!=0
qr=qp(Fp(28));wr=Fp(28);Er=EllipticCurve(Fp,[(1+qr)*wr+(2-qr),0,qr*wr*wr+(1-qr)*wr,0,0]);assert Er.cardinality()==45
print('DECIC315_EXACT_WITNESS_REPLAY_PASS full_torsion=[3,15] j_degree=5 nonCM=true prime=43 order=45',flush=True)
report={'complete':True,'degree':10,'rational_Jacobian':'[3,3,24,96,96]','full_torsion':[3,15],'all_containing_curves_same_full_torsion':True,'local_F4_divisors':24432,'local_F7_noncuspidal_divisors':57072,'gonality_bound_F7':6,'exact_diagonal_classes':240,'known_reducible_classes':192,'remaining_classes':48,'new_curve_automorphism_orbit':48,'excluded_strict_supergroup_examples':[[3,30],[3,45]]}
(O/'independent_replay.json').write_text(json.dumps(report,indent=2)+'\n')
print('DECIC315_COMPLETE_REPLAY_PASS',flush=True)
