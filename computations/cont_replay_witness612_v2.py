from pathlib import Path
import ast,json,math
from sage.all import QQ,ZZ,GF,PolynomialRing,NumberField,EllipticCurve
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
a=ast.literal_eval((D/'M6N12_exact_witness.m').read_text().replace('[*','[').replace('*]',']'))
R=PolynomialRing(QQ,'Z');F=NumberField(R([1,-1,1]),'z');z=F.gen();S=PolynomialRing(F,'t')
def decF(row):return F(QQ(row[0][0])/QQ(row[0][1]))+F(QQ(row[1][0])/QQ(row[1][1]))*z
g=S([decF(c) for c in a[0]]);assert g.degree()==5 and g.is_irreducible();K=F.extension(g,'w');w=K.gen();assert K.absolute_degree()==10
def decK(row):return sum((K(decF(c))*w**i for i,c in enumerate(row)),K(0))
ai=[decK(c) for c in a[1]];E=EllipticCurve(K,ai);assert E.discriminant()!=0
P=E([decK(c) for c in a[2]]);Q=E([decK(c) for c in a[3]]);assert 6*P==E(0) and 3*P!=E(0) and 2*P!=E(0);assert 12*Q==E(0) and 6*Q!=E(0) and 4*Q!=E(0)
pairing=P.weil_pairing(2*Q,ZZ(6));assert pairing==K(z) and pairing.multiplicative_order()==6
assert E.j_invariant()==decK(a[4]);bound=0;reductions=[]
for p,zr,wr,count,invs in a[5]:
 p=int(p);fp=GF(p);zz=fp(zr);ww=fp(wr);assert zz**2-zz+1==0 and 2*zz-1!=0
 redF=lambda x:sum((fp(c)*zz**j for j,c in enumerate(F(x).list())),fp(0));Sp=PolynomialRing(fp,'x');gp=Sp([redF(c) for c in g.list()]);assert gp(ww)==0 and gp.derivative()(ww)!=0
 redK=lambda x:sum((redF(c)*ww**j for j,c in enumerate(K(x).list())),fp(0));Er=EllipticCurve(fp,[redK(c) for c in ai]);assert Er.discriminant()!=0 and Er.cardinality()==count
 # Count every Weierstrass fiber using only modular arithmetic.
 aa=[int(v) for v in Er.a_invariants()];total=1
 for xx in range(p):
  disc=((aa[0]*xx+aa[2])**2+4*(xx**3+aa[1]*xx**2+aa[3]*xx+aa[4]))%p;leg=pow(disc,(p-1)//2,p);total+=1+(0 if leg==0 else 1 if leg==1 else -1)
 assert total==count;bound=math.gcd(bound,count);reductions.append({'p':p,'z_root':int(zr),'w_root':int(wr),'count':int(count),'independent_fiber_count':total,'unramified_simple_root':True})
assert bound==72 and a[6]==72 and a[7]==[6,12]
report={'group':[6,12],'absolute_degree':10,'relative_degree':5,'relative_field_polynomial':str(g),'exact_point_orders':[6,12],'primitive_Weil_pairing_order':6,'known_subgroup_order':72,'good_reductions':reductions,'full_torsion_order_upper_bound':72,'full_torsion_invariants':[6,12],'independent_field_point_pairing_and_full_torsion_replay':'PASS','sporadic_reason':'rank J1(6,12)(Q(zeta6))=0 and relative degree 5 < gonality at least 6'}
(H/'occurrences/new612_verified.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC612_EXACT_WITNESS_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
