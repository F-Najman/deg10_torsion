# Written for the repository: polredabs of the decic fields of Sections 3 and 11 and their LMFDB labels.
from sage.all import QQ, PolynomialRing, NumberField, pari
import psycopg2
R = PolynomialRing(QQ, 'x'); x = R.gen()
def absolute(rel_base, rel):
    F = NumberField(rel_base, 'z'); z = F.gen()
    S = PolynomialRing(F, 'W'); W = S.gen()
    L = F.extension(eval(rel), 'w')
    return L.absolute_field('a').defining_polynomial()
fields = {
 '6x12': absolute(x**2-x+1, "3*W**5+(4-5*z)*W**4+(7-8*z)*W**3-7*(z+1)*W**2+6*(z-1)*W+2*z-1"),
 '3x15': absolute(x**2+x+1, "W**5+13*(z+1)*W**4-(24*z+22)*W**3+4*(z+1)*W**2+(z+4)*W+z+1"),
 '3x18': absolute(x**2+x+1, "3*W**5+(1+5*z)*W**4+(6*z-3)*W**3-(6*z+9)*W**2-(5*z+4)*W+3"),
 '2x26a': 2*x**10+11*x**9+17*x**8-7*x**7-42*x**6-20*x**5+9*x**4+9*x**3+10*x**2+6*x+1,
 '2x26b': x**10+3*x**9-3*x**8-23*x**7-26*x**6+3*x**5+13*x**4-2*x**3-9*x**2-5*x-1,
}
# Z/2 x Z/28 and Z/2 x Z/30: quadratic extension of the quintic field by the quadratic factor of the 2-division polynomial
def full2(fq, bco, cco):
    F = NumberField(fq, 'a'); a = F.gen()
    b = sum(c*a**i for i, c in enumerate(bco)); c = sum(cc*a**i for i, cc in enumerate(cco))
    S = PolynomialRing(F, 'X'); X = S.gen()
    # E(b,c): y^2 + (1-c)xy - b y = x^3 - b x^2 ; 2-division polynomial 4x^3 + b2 x^2 + 2 b4 x + b6
    a1, a2, a3, a4, a6 = 1-c, -b, -b, 0, 0
    b2 = a1**2+4*a2; b4 = 2*a4+a1*a3; b6 = a3**2+4*a6
    psi2 = 4*X**3+b2*X**2+2*b4*X+b6
    h = [g for g, e in psi2.factor() if g.degree() == 2][0]
    return F.extension(h, 'w').absolute_field('t').defining_polynomial()
fields['2x28'] = full2(x**5-x**4-2*x**3-x**2+2*x+2, [1032,-239,-222,-759,419], [42,-11,-8,-31,17])
fields['2x30'] = full2(x**5+x**4-3*x**3+3*x+1, [-9118,-14321,20475,-1917,-6378], [-162,-257,366,-34,-114])
conn = psycopg2.connect(host='devmirror.lmfdb.xyz', dbname='lmfdb', user='lmfdb', password='lmfdb', connect_timeout=20)
cur = conn.cursor()
for k, f in fields.items():
    g = R(pari(f).polredabs())
    K = NumberField(g, 'a')
    coeffs = [int(c) for c in g.list()]
    cur.execute("select label from nf_fields where coeffs = %s::numeric[]", (coeffs,))
    print('FIELD', k, 'polredabs', g, 'disc', K.discriminant(), 'signature', K.signature(), 'LMFDB', cur.fetchall())
