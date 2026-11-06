from pathlib import Path
import ast,json
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,EllipticCurve,prime_range
P=Path('work/deg10/C226_resolution_20261030');C=Path('work/deg10/continuation_20261029b/C226')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def q(v):return QQ(v[0])/QQ(v[1])
data=read(C/'filter26.m');own={(r[0],r[1]):r for r in data[1]};small={(r[0],r[1]):r for r in data[2]};reports=[];jpolys=[]
def key(bc,d):
 out=[]
 for _ in range(d):out.append(tuple(int(v) for x in bc for v in x._vector_()));bc=[x**3 for x in bc]
 return min(out)
for idx,desc in [(30,[2,5,1,5,19]),(24,[1,10,54])]:
 fc,bc,cc,m=read(P/f'exact_points_{idx}.m')[0];R=PolynomialRing(QQ,'x');f=R([q(v) for v in fc]);S=PolynomialRing(GF(3),'z');fp=S(f);assert fp.gcd(fp.derivative()).degree()==0;fac=fp.factor();ds=sorted([int(h.degree()) for h,e in fac for _ in range(e)]);assert ds==([5,5] if idx==30 else [10]);d=ds[0];modulus=next(r[1] for r in data[3] if r[0]==d);L=GF(3**d,'z',modulus=S(modulus));T=PolynomialRing(L,'x');xs=T(fp).roots(multiplicities=False);assert len(xs)==10;keys=set()
 for x in xs:
  b=sum(L(q(v))*x**i for i,v in enumerate(bc));c=sum(L(q(v))*x**i for i,v in enumerate(cc));keys.add(key([b,c],d))
 expected={tuple(own[desc[1],desc[2]][2])}
 if desc[0]==2:expected.add(tuple(small[desc[3],desc[4]][2]))
 assert keys==expected,(idx,keys,expected)
 K=NumberField(f,'a');b=K([q(v) for v in bc]);c=K([q(v) for v in cc]);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);jpoly=E.j_invariant().minpoly();jpolys.append(jpoly);certs=[];roots=[K(0),(13*E(b,b))[0]];roots.append(-c-roots[1]);differences=[(i,j,roots[j]-roots[i]) for i in range(3) for j in range(3) if i!=j];done=set()
 for pp in prime_range(17,3000):
  pp=int(pp);k=GF(pp);U=PolynomialRing(k,'x')
  try:ff=U(f);vs=[U([k(v) for v in a.list()]) for i,j,a in differences]
  except (ZeroDivisionError,ValueError):continue
  if ff.gcd(ff.derivative()).degree()!=0:continue
  for rr,_ in ff.roots():
   for no,((i,j,a),pol) in enumerate(zip(differences,vs)):
    if no in done:continue
    v=pol(rr)
    if v!=0 and v**((pp-1)//2)==k(-1):certs.append({'kernel_root_index':i,'other_root_index':j,'prime':pp,'field_generator_residue':int(rr),'difference_residue':int(v),'euler_symbol':-1});done.add(no)
   if len(done)==6:break
  if len(done)==6:break
 assert len(done)==6
 reports.append({'class':idx,'descriptor':desc,'field_mod3_factor_degrees':ds,'exact_Frobenius_keys_match_saved_class':True,'j_polynomial':[str(v) for v in jpoly.list()],'six_nonsquare_certificates':certs,'all_three_two_isogenous_curves_have_no_K_rational_point_of_order4':True})
assert jpolys[0]!=jpolys[1]
report={'status':'PASS','both_classes_effective_over_Q':True,'each_unique_effective_divisor_irreducible_degree10':True,'distinct_Galois_orbits_of_j_invariants':True,'class_records':reports,'cyclic52_excluded_by_complete_C226_classification_and_two_isogeny_tests':True}
(P/'independent_class_binding.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
