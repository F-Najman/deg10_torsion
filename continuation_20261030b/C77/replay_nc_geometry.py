"""Independent source, marked3/11-isogeny and quadratic lift geometry."""
from sage.all import GF,PolynomialRing,EllipticCurve,lcm
from pathlib import Path
from collections import Counter
import ast,json,time,sys
H=Path(__file__).resolve().parent;k=GF(2);R=PolynomialRing(k,'x');x=R.gen();t0=time.monotonic()
def read(p):return ast.literal_eval(Path(p).read_text().replace('[*','[').replace('*]',']'))
proj=read(H/'projection_model_data.m')
def bc(E,P):
 X,Y=P.xy();a1,a2,a3,a4,a6=E.a_invariants();u3=2*Y+a3+a1*X;u2=3*X+a2;u4=3*X**2+2*X*a2+a4-a1*Y
 v1=2*u4/u3+a1;v2=u2-(u4/u3)**2-a1*u4/u3
 return(-v2**3/u3**2,-(v1*v2-u3)/u3)
def t7(E,P):
 assert 7*P==E(0)and P!=E(0);b,c=bc(E,P);t=b/c;assert b==t**3-t**2 and c==t**2-t;return t
def ord77(P):return 77*P==P.curve()(0)and 7*P!=P.curve()(0)and 11*P!=P.curve()(0)
def j11(E,P):
 P11=7*P;assert 11*P11==E(0)and P11!=E(0);T=PolynomialRing(E.base_field(),'w');w=T.gen();f=T(1)
 for i in range(1,6):f*=w-(i*P11)[0]
 return E.isogeny(f).codomain().j_invariant()
def key(vs):
 F=vs[0].parent();n=F.degree()
 return min(tuple(int((v**(2**i)).polynomial()[j])for v in vs for j in range(n))for i in range(n))
def deg(vs):
 for d in range(1,vs[0].parent().degree()+1):
  if all(v**(2**d)==v for v in vs):return d
 raise AssertionError
def rf(code,z):
 T=PolynomialRing(z.parent(),'w');return T(code[0])(z)/T(code[1])(z)
def ef3(code,tv,xv):return sum(rf(c,tv)*xv**i for i,c in enumerate(code))
def x11_check(E,PN,xv):
 F=E.base_field();T=PolynomialRing(F,'w');w=T.gen();jj=E.j_invariant();jj11=j11(E,PN)
 jp=[rf(c,xv)for c in proj[1]];kp=[rf(c,xv)for c in proj[0]];dp=[rf(c,xv)for c in proj[2]]
 ys=[(jj-jp[0])/jp[1]]if jp[1]else T(kp).roots(multiplicities=False)
 ys=[y for y in ys if T(kp)(y)==0 and jp[0]+jp[1]*y==jj and dp[0]+dp[1]*y==jj11];assert len(ys)==1
 return jj11
CASES=[int(sys.argv[1])]if len(sys.argv)>1 else [1,2,3]
stats=[];kernels=0;sourceprojections=[]
for case in CASES:
 packet=read(H/f'nc_case{case}_v3_input_packet.m');row=packet[7];obs=packet[8];n=int(lcm([len(r[4][0])-1 for r in obs]));F=GF(2**n,'a');T=PolynomialRing(F,'w');w=T.gen()
 def root(defn):
  p=T(defn);assert R(defn).is_irreducible();rs=p.roots(multiplicities=False);assert len(rs)==len(defn)-1;return rs[0]
 beta=root(row[1]);fv=lambda v:T(v)(beta)
 E=EllipticCurve(F,[fv(v)for v in row[2]]);PN=E([fv(v)for v in row[3]]);assert ord77(PN)
 tv=t7(E,11*PN);assert deg([E.j_invariant()])==10 and t7(E,22*PN)==(tv+1)/tv and t7(E,33*PN)==1/(tv+1)
 # Independently solve the exact X0(11) j/j11 fibre for the input.
 PZ=PolynomialRing(F,'z');Q=PZ.fraction_field();z=Q.gen()
 def qr(c):return Q(PZ(c[0]))/Q(PZ(c[1]))
 j0,j1=[qr(c)for c in proj[1]];aa=qr(proj[0][1]);bb=qr(proj[0][0]);N=((E.j_invariant()+j0)**2+aa*j1*(E.j_invariant()+j0)+bb*j1*j1).numerator()
 valid=[]
 for xv in N.roots(multiplicities=False):
  try:x11_check(E,PN,xv)
  except (AssertionError,ZeroDivisionError):continue
  valid.append(xv)
 assert len(valid)==1;sval=valid[0];assert sum(rf(c,tv)*sval**i for i,c in enumerate(packet[0]))==0
 assert len(T([ef3(c,tv,sval)for c in packet[1]]).roots(multiplicities=False))==2
 sourceprojections.append(dict(case=case,field_polynomial=[int(c)for c in F.modulus()],coordinates=[[int(v.polynomial()[i])for i in range(n)]for v in [tv,sval]]))
 roots=E.division_polynomial(3).roots(multiplicities=False);assert len(roots)==4;expect=Counter();product=T(1)
 for r in roots:
  product*=w-r;phi=E.isogeny(w-r);PP=phi(PN);assert ord77(PP);expect[key([t7(phi.codomain(),11*PP),j11(phi.codomain(),PP)])]+=1
 assert T(E.division_polynomial(3))==E.division_polynomial(3).leading_coefficient()*product
 actual=Counter();details=[]
 for idx,wgt,d,mult,ff in obs:
  gamma=root(ff[0]);vals=[T(v)(gamma)for v in ff[1]];EE=EllipticCurve(F,vals[:5]);PP=EE(vals[5:7]);assert ord77(PP)
  assert t7(EE,11*PP)==vals[7]and x11_check(EE,PP,vals[8])==vals[9]
  assert sum(rf(c,vals[7])*vals[8]**i for i,c in enumerate(packet[0]))==0
  assert len(T([ef3(c,vals[7],vals[8])for c in packet[1]]).roots(multiplicities=False))==2
  assert deg(vals[7:9])==d and 10*wgt==d*mult
  actual[key([vals[7],vals[9]])]+=wgt;details.append([int(wgt),int(d),int(mult)])
 assert expect==actual and sum(actual.values())==4
 kernels+=4;stats.append(dict(case=case,isogeny_orbits=details,lift_allocations=len(packet[9])))
assert pow(36,5,77)==1 and 64==pow(36,2,77)and 71==pow(36,3,77)
assert (2*3)%77==(-pow(36,3,77))%77
out=dict(status='PASS',cases=stats,independent_marked_Velu3_kernels=kernels,independent_marked_Velu11_quotients=kernels+len(CASES),all_division3_polynomials_complete=True,all_equal_degree_AS_lifts_split=True,delta2_cubic_exponent=2,delta3_is_inverse_delta2=True,source_projection_values=sourceprojections,seconds=time.monotonic()-t0)
(H/('independent_nc_geometry_case%d.json'%CASES[0]if len(CASES)==1 else 'independent_nc_geometry.json')).write_text(json.dumps(out,indent=2)+'\n');print('C77_INDEPENDENT_NC_GEOMETRY_PASS',json.dumps(out,sort_keys=True))
