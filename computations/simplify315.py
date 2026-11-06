from sage.all import *
from pathlib import Path
import json
H=Path(__file__).resolve().parents[1];d=json.loads((H/'occurrences/new315_verified.json').read_text());R=PolynomialRing(QQ,'z');F=NumberField(R.gen()**2+R.gen()+1,'z');z=F.gen();T=PolynomialRing(F,'t');t=T.gen();el=lambda a:QQ(a[0])+QQ(a[1])*z;f=T([el(a) for a in d['relative_polynomial']]);q=T([el(a) for a in d['q_coefficients']]);cols=[list((q*t**i)%f)+[F(0)]*(5-len(list((q*t**i)%f))) for i in range(3)]+[list(-t**i)+[F(0)]*(5-len(list(-t**i))) for i in range(3)];V=matrix(F,cols).left_kernel();assert V.dimension()==1;v=V.basis()[0];den=T(list(v[:3]));num=T(list(v[3:]));lc=den.leading_coefficient();num=T(num/lc);den=T(den/lc);assert (q*den-num)%f==0
print('Simplified q:',num,' / ',den)
S=PolynomialRing(QQ,'w');w=S.gen();aa=w**5+13*w**4-22*w**3+4*w**2+4*w+1;bb=13*w**4-24*w**3+4*w**2+w+1;g=aa**2-aa*bb+bb**2;assert g.is_irreducible();zz=(-aa*bb.inverse_mod(g))%g;assert (zz*zz+zz+1)%g==0
out={'absolute_field_polynomial':list(map(str,g.list())),'zeta3_polynomial':list(map(str,zz.list())),'q_rational_function_numerator':[[str(c) for c in F(a).list()]+['0']*(2-len(F(a).list())) for a in num.list()],'q_rational_function_denominator':[[str(c) for c in F(a).list()]+['0']*(2-len(F(a).list())) for a in den.list()]};(H/'occurrences/new315_simplified.json').write_text(json.dumps(out,indent=2)+'\n');print('Absolute polynomial:',g)
