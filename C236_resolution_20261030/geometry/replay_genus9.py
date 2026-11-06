from pathlib import Path
import ast,json
from sage.all import ZZ,matrix,vector,FreeModule,GF,PolynomialRing,EllipticCurve,divisors
P=Path(__file__).resolve().parent

def rd(p):
 return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
class Quotient:
 def __init__(self,R,hs):
  self.n=len(R[0])-1
  assert all(r[-1]==0 for r in R+hs)
  self.H=matrix(ZZ,[r[:-1] for r in R+hs]).row_module()
  self.index=self.H.index_in(FreeModule(ZZ,self.n))
 def zero(self,v):
  assert len(v)==self.n+1 and v[-1]==0
  return vector(ZZ,v[:-1]) in self.H
iv,rels,hc,acs,a,aa,possible=rd(P/'quotient9_image_fullcontrol_arithmetic.m')
assert iv==[3,3,3,24,24,192,0] and all(not any(v) for v in acs)
G=Quotient(rels,[]);H=Quotient(rels,hc)
assert G.index==2985984 and H.index==8
assert all(not H.zero([n*v for v in a]) for n in range(1,8)) and H.zero([8*v for v in a])
assert G.zero([8*v for v in aa]) and all(not G.zero([n*v for v in aa]) for n in [1,2,4])
r=G.n;rr=[v[:-1] for v in rels];gens=hc+[a];outs=acs+[aa]
rows=[v+[0]*r for v in rr]+[[0]*r+v for v in rr]+[u[:-1]+v[:-1] for u,v in zip(gens,outs)]
graph=matrix(ZZ,rows).row_module();assert graph.index_in(FreeModule(ZZ,2*r))==G.index
assert Quotient(rels,gens).index==1 and Quotient(rels,outs).index==G.index//8
out=dict(complete=True,finite_picard_order=int(G.index),rational_image_upper_index=8,A7_image_order=8,A7_kernel_is_all_finite_cuspidal_subgroup=True,graph_relations_checked_over_Z=True)
if (P/'genus9_stage2_arithmetic.m').exists():
 data=rd(P/'genus9_stage2_arithmetic.m');fields,images,tests,remaining=data;pic=rd(P/'quotient9_proven_picard_v2.m');G2=Quotient(pic[5],[]);H2=Quotient(pic[5],pic[9]);assert G2.index==G.index and H2.index==8
 R=PolynomialRing(GF(5),'z');fs={d:GF(5**d,'z',modulus=R(pol)) for d,pol in fields};cache={}
 def fkey(vals,d):
  ans=[]
  for _ in range(d):
   ans.append(tuple(int(x) for v in vals for x in v._vector_()));vals=[v**5 for v in vals]
  return min(ans)
 for d,key,e,w,cl,projkey in images:
  F=fs[d];z=F.gen();b=R(key[:d])(z);c=R(key[d:])(z);E=EllipticCurve(F,[0,c,0,b*(1-b-c),0]);Q=E(b,b);assert 36*Q==E(0) and 18*Q!=E(0) and 12*Q!=E(0)
  Q12=3*Q;T=12*Q;R9=4*Q;lam=Q12[1]/Q12[0];xp=R9[0]+(R9+T)[0]+(R9-T)[0]-2*T[0];vals=[Q12[0]/lam**2,c/lam**2,xp/lam**2]
  assert fkey(vals,d)==tuple(projkey)
  de=min(int(j) for j in divisors(d) if all(v**(5**int(j))==v for v in vals));assert e==de and w==d//e and w in [1,3]
  cache[(d,tuple(key))]=cl
 original=rd(P.parent/'enumeration/first_survivors_unified.m')[1];orig={row[0]:row for row in original};expected=[]
 for case_id,cl,yes in tests:
  summed=[0]*7
  for d,key in orig[case_id][1]:summed=[x+y for x,y in zip(summed,cache[d,tuple(key)])]
  assert G2.zero([x-y for x,y in zip(summed,cl)]) and bool(yes)==H2.zero(cl)
  if yes:expected.append(orig[case_id])
 assert expected==remaining
 out.update(stage2_supports_checked=len(tests),stage2_remaining=len(remaining),stage2_unique_projected_places_checked=len(images),all_input_checked=len(tests)==len(original))
(P/'genus9_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_D10_GENUS9_A7_LATTICE_AND_FILTER_REPLAY_PASS',json.dumps(out),flush=True)
