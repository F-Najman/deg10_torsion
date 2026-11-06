from pathlib import Path
import ast,json
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,EllipticCurve
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
c=json.loads((D/'M6N12_absolute_polynomial_coefficients.json').read_text())[::-1];R=PolynomialRing(QQ,'a');f=R(c);assert f.degree()==10 and f.is_monic() and f.is_irreducible();K=NumberField(f,'a');a=K.gen();w=1/a
A=3*w**5+4*w**4+7*w**3-7*w**2-6*w-1;B=5*w**4+8*w**3+7*w**2-6*w-2;z=A/B;assert z**2-z+1==0
raw=ast.literal_eval((D/'M6N12_exact_witness.m').read_text().replace('[*','[').replace('*]',']'))
df=lambda row:K(QQ(row[0][0])/QQ(row[0][1]))+K(QQ(row[1][0])/QQ(row[1][1]))*z
dk=lambda row:sum((df(c)*w**i for i,c in enumerate(row)),K(0))
assert sum(df(c)*w**i for i,c in enumerate(raw[0]))==0
r=(w*w+1)/(w*w-w);s=(w*w+w+1)/(w*w);E=EllipticCurve(K,[1+s-r*s,r*s*(1-r),r*s*(1-r),0,0]);assert E.a_invariants()==tuple(dk(c) for c in raw[1])
P=E([dk(c) for c in raw[2]]);Q=E([dk(c) for c in raw[3]]);assert 6*P==E(0) and 2*P!=E(0) and 3*P!=E(0) and 12*Q==E(0) and 6*Q!=E(0) and 4*Q!=E(0);assert P.weil_pairing(2*Q,ZZ(6))==z
p=61;F=GF(p);aa=F(23);ff=f.change_ring(F);assert ff(aa)==0 and ff.derivative()(aa)!=0
red=lambda x:sum((F(v)*aa**i for i,v in enumerate(K(x).list())),F(0));assert red(w)==8 and red(z)==14
Er=EllipticCurve(F,[red(c) for c in E.a_invariants()]);assert Er.discriminant()!=0 and Er.cardinality()==72
enc=lambda x:[str(c) for c in K(x).list()]
report={'group':[6,12],'absolute_degree':10,'absolute_polynomial_coefficients':c,'absolute_field_polynomial':str(f),'curve_a_invariants':[enc(x) for x in E.a_invariants()],'P6':[enc(x) for x in P],'Q12':[enc(x) for x in Q],'zeta6':enc(z),'u':enc(w),'good_reduction':{'p':61,'absolute_generator_root':23,'count':72,'unramified':True},'independent_absolute_field_isomorphism_points_pairing_and_full_torsion':'PASS'}
(H/'occurrences/new612_absolute.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC612_ABSOLUTE_WITNESS_REPLAY_COMPLETE',report['absolute_field_polynomial'],flush=True)
