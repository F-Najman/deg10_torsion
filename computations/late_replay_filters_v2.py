from pathlib import Path
import ast,json,time
from collections import Counter
from sage.all import GF,PolynomialRing,EllipticCurve,divisors,ZZ,matrix
P=Path('work/deg10/continuation_20261029b');t0=time.monotonic()
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quot:
 def __init__(self,rels,gens=[]):
  self.n=len(rels[0])-1;self.H=matrix(ZZ,[v[:-1] for v in rels+gens]).row_module().basis_matrix();assert self.H.nrows()==self.n;self.H=self.H.hermite_form();self.index=abs(self.H.det());self.cols=self.H.transpose()
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0
  return tuple(self.H.row_module().coordinate_vector(ZZ(v[0])*self.H.row(0))) if False else tuple(self.cols.solve_right(matrix(ZZ,self.n,1,list(v[:-1]))).list())
 def contains(self,v):return all(x.denominator()==1 for x in self.key(v))
 def residue(self,v):return tuple(x-x.floor() for x in self.key(v))
def sums(rows):return [sum(x) for x in zip(*rows)]
def allowed(rel,hg,ds,cc,D=6):
 Q=Quot(rel,hg);out=[{Q.residue([0]*(Q.n+1))}]+[set() for _ in range(D)]
 for w,v in zip(ds,cc):
  for d in range(w,D+1):out[d].update(Q.residue([x+y for x,y in zip(list(z)+[0],v)]) for z in []) if False else out[d].update(tuple(((x+y)-(x+y).floor()) for x,y in zip(z,Q.residue(v))) for z in list(out[d-w]))
 return Q,out
R=PolynomialRing(GF(3),'z');fields={}
def field(d,pol):
 k=(d,tuple(pol))
 if k not in fields:fields[k]=GF(3**d,'z',modulus=R(pol))
 return fields[k]
def src(d,ky,K):
 z=K.gen();b=sum(K(v)*z**i for i,v in enumerate(ky[:d]));c=sum(K(v)*z**i for i,v in enumerate(ky[d:]));E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);return E,E(0,0),E(b,b)
def iso2(E,T,points):
 r=T[0];A=E.a2()+3*r;B=3*r*r+2*E.a2()*r+E.a4();F=EllipticCurve(E.base_field(),[0,-2*A,0,A*A-4*B,0]);out=[]
 for U in points:x=U[0]-r;out.append(F(x+A+B/x,U[1]*(1-B/x**2)))
 return F,out
def bc(E,T,Q):x=Q[0]-T[0];y=Q[1];return [x**3/y**2,(E.a2()+3*T[0])*x*x/y**2]
def tate(E,Q):
 x,y=Q[0],Q[1];a3=2*y;a2=3*x+E.a2();a4=3*x*x+2*x*E.a2()+E.a4();a1=2*a4/a3;a2-=(a4/a3)**2;b=-a2**3/a3**2;c=-(a1*a2-a3)/a3;r=b/c;s=c*c/(b-c);u=r*s-2*r+1;return [(s-r)/u,u/(s*s-s-r+1)]
def deg(v,d):return min(int(e) for e in divisors(d) if all(x**(3**int(e))==x for x in v))
# The two low-genus filters at 40.
f=read(P/'C240/filter_quotient20_linear_arithmetic.txt');g,fs,single,small,pairs,triples,first,cov=f;iv,rel,ds,cc,hg,hinv=g;G=Quot(rel);H,A=allowed(rel,hg,ds,cc);assert H.index==4 and all(H.contains(v) for v in cc)
F={d:field(d,p) for d,p in fs};own={(r[0],r[1]):r for r in single};sm={(r[0],r[1]):r for r in small};assert len(single)==301 and len(pairs)==480 and not triples
for row in single+small:
 d,i,ky,e,w,cl=row[:6];E,T,Q=src(d,ky,F[d]);assert deg(bc(E,T,2*Q),d)==e and w==d//e
expected=[[1,r[0],r[1]] for r in single if H.contains(r[5])]
assert all(bool(r[6])==H.contains(r[5]) for r in single)
for D,V in [(4,4),(4,5),(4,6),(5,5)]:
 rows=[r for r in pairs if r[0]==D and r[2]==V];assert [(r[1],r[3]) for r in rows]==[(i,j) for d,i in own if d==D for v,j in sm if v==V]
 for d,i,v,j,cl,yes in rows:
  actual=sums([own[d,i][5],sm[v,j][5]]);assert G.contains([x-y for x,y in zip(cl,actual)]) and bool(yes)==H.contains(cl)
  if yes:expected.append([2,d,i,v,j])
assert expected==first and len(first)==288
ig,ims,tests,second=read(P/'C240/filter_quotient20_isogeny_arithmetic.txt');iiv,ir,ids,icc,ihg,ihin=ig;IG=Quot(ir);IH,IA=allowed(ir,ihg,ids,icc);cache={}
for kk,e,w,cl in ims:
 d=kk[0];E,T,Q=src(d,kk[1:],F[d]);EP,(TP,U)=iso2(E,20*Q,[T,Q]);assert 20*U==EP(0) and 10*U!=EP(0) and 4*U!=EP(0) and TP!=10*U;assert deg(bc(EP,TP,U),d)==e and w==d//e;cache[tuple(kk)]=cl
for k,desc,cl,yes in tests:
 assert desc==first[k-1];rows=[own[desc[1],desc[2]]]+([sm[desc[3],desc[4]]] if desc[0]==2 else []);actual=sums([cache[tuple([r[0]]+r[2])] for r in rows]);assert IG.contains([x-y for x,y in zip(cl,actual)]) and bool(yes)==IH.contains(cl)
assert second==[r[1] for r in tests if r[3]] and len(second)==118
# All eight degree-independent Picard presentations and every new class sum.
seen=set();picard_images=0;witnesses=0;picgroup=None
for shard in range(8):
 group,images,tests,left=read(P/f'C240/portable/shard{shard}_arithmetic.txt');assert not left
 if picgroup is None:picgroup=group
 else:assert group==picgroup
 iv,rel,hg,ds,cc,native=group;G40=Quot(rel);H40,A40=allowed(rel,hg,ds,cc);assert G40.index==31579088486400 and H40.index==4592 and list(map(len,A40))==list(map(len,native))
 c={}
 for kk,e,w,cl in images:
  m,d=kk[:2];E,T,Q=src(d,kk[2:],F[d]);assert m==1;assert deg(tate(E,Q),d)==e and w==d//e;c[tuple(kk)]=cl
 for k,desc,cd,checks,passed in tests:
  assert k not in seen and (k-1)%8==shard;seen.add(k);assert desc==second[k-1] and not passed
  rows=[own[desc[1],desc[2]]]+([sm[desc[3],desc[4]]] if desc[0]==2 else []);assert cd==10-sum(r[0] for r in rows)
  for m,cl,yes in checks:
   actual=sums([c[tuple([m,r[0]]+r[2])] for r in rows]);assert G40.contains([x-y for x,y in zip(cl,actual)]) and bool(yes)==(H40.residue([-x for x in cl]) in A40[cd])
  assert not checks[-1][2]
 picard_images+=len(images)
 geom={n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse((P/'C240/portable/imported/factorbase_geometry.py').read_text()).body};fb=geom['factorbase_classes'];fd=geom['factorbase_degrees']
 ww=ast.literal_eval(ast.parse((P/f'C240/portable/shard{shard}_RR_identity.log.witnesses.py').read_text()).body[0].value)
 for trial,d,result,ev,rv,va,b in ww:
  assert sum(e*fd[j-1] for j,e in ev)==sum(e*fd[j-1] for j,e in rv)==25
  assert sum(e*fd[j-1] for j,e in va)==1
  actual=[0]*9
  for sign,sp in [(1,ev),(-1,rv),(d,va)]:
   for j,e in sp:actual=[a+sign*e*x for a,x in zip(actual,fb[j-1])]
  assert G40.contains([a-b for a,b in zip(actual,result)]);witnesses+=1
assert seen==set(range(1,119))
# Low-genus C2 x C38 filter, all mixed types, and every final Hecke decision.
group,source,small,cover,images,tests,remaining=read(P/'C238/filter19_v2.m');iv,rel,hg,ds,cc,native=group;G19=Quot(rel);H19,A19=allowed(rel,hg,ds,cc);assert H19.index==2
F38={d:field(d,p) for d,p,c in cover};o={(r[0],r[1]):r for r in source};s={(r[0],r[1]):r for r in small};cache={}
for kk,e,w,cl in images:
 m,d=kk[:2];E,T,Q=src(d,kk[2:],F38[d]);V=2*Q
 if m==1:EP=E;U=V
 else:EP,(U,)=iso2(E,[T,19*Q,T+19*Q][m-2],[V])
 assert 19*U==EP(0) and U!=EP(0) and deg(tate(EP,U),d)==e and w==d//e;cache[tuple(kk)]=cl
support=[[1,d,i] for d,i in o]+[[2,5,i,5,j] for d,i in o if d==5 for v,j in s if v==5];assert len(support)==233
for k,desc,cd,checks,passed in tests:
 assert desc==support[k-1];rs=[o[desc[1],desc[2]]]+([s[desc[3],desc[4]]] if desc[0]==2 else []);assert cd==10-sum(r[0] for r in rs)
 for m,cl,yes in checks:
  actual=sums([cache[tuple([m,r[0]]+r[2])] for r in rs]);assert G19.contains([a-b for a,b in zip(actual,cl)]) and bool(yes)==(H19.residue([-x for x in cl]) in A19[cd])
 assert bool(passed)==all(x[2] for x in checks)
assert remaining==[r[1] for r in tests if r[4]] and len(remaining)==35
seen=set()
for shard in range(4):
 images,tests,left=read(P/f'C238/hecke_shard{shard}.m');assert not left
 for kk,e,w,xy in images:
  m,d=kk[:2];E,T,Q=src(d,kk[2:],F38[d]);assert m==1;pair=tate(E,Q);assert [[int(x) for x in z._vector_()] for z in pair]==xy and deg(pair,d)==e and w==d//e
 for k,desc,cd,checks,passed in tests:
  assert k not in seen and (k-1)%4==shard and desc==remaining[k-1] and not passed;seen.add(k);assert all(positive==(q+1)*(10-cd) for m,q,positive,yes in checks);assert not checks[-1][3]
assert seen==set(range(1,36))
report={'status':'PASS','independently_replayed_groups':[[2,38],[2,40]],'cyclic_isogeny_consequences':[[1,76],[1,80]],'C240_initial_supports':781,'C240_first_survivors':288,'C240_second_survivors':118,'C240_final_survivors':0,'C240_picard_images':picard_images,'C240_exact_RR_witness_class_sums':witnesses,'C238_initial_supports':233,'C238_low_genus_survivors':35,'C238_final_Hecke_survivors':0,'elapsed_seconds':time.monotonic()-t0}
(P/'four_exclusions_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
