// Lemma lem:sw-60-quintic. For the two pairs (K_30,E_30) and (K_30,E_30') of [N5, Proposition 9.1], and M = K_30(E[2]),
// each of the 3 curves E/<T>, T in E[2] \ {0}, which are the curves 2-isogenous to E over M, has torsion subgroup
// Z/30Z over M (Magma's TorsionSubgroup).
// The curves are E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2, with (b,c) as in [N5, (7) and (8)], the same values as in
// continuation_20261031b/C60/d31b_C60_quintic_2tower.m.
// Run from the repository root: magma -b additions/C60/tors60.m
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
R<t>:=PolynomialRing(Rationals());F<a>:=NumberField(t^5+t^4-3*t^3+3*t+1);
// coefficients of 1, a, ..., a^4 of b and c
rows:=[*[*[-9118,-14321,20475,-1917,-6378],[-162,-257,366,-34,-114]*],[*[295,238,-382,-38,190],[28,9,-36,13,7]*]*];
S<x>:=PolynomialRing(F);
for i in [1,2] do
  b:=Evaluate(R!rows[i][1],a);c:=Evaluate(R!rows[i][2],a);
  E:=EllipticCurve([1-c,-b,-b,0,0]);P:=E![0,0];assert Order(P) eq 30;
  fac:=Factorization(DivisionPolynomial(E,2));assert Sort([Degree(f[1]):f in fac]) eq [1,2];
  h:=[f[1]:f in fac|Degree(f[1]) eq 2][1];
  M<w>:=ext<F|h>;assert AbsoluteDegree(M) eq 10;EM:=BaseChange(E,M);
  printf "TORS60 curve %o: E(M)_tors has invariants %o\n",i,Invariants(TorsionSubgroup(EM));
  SM<y>:=PolynomialRing(M);
  rts:=[r[1]:r in Roots(DivisionPolynomial(EM,2))];assert #rts eq 3;
  for r in rts do
    E2,phi:=IsogenyFromKernel(EM,y-r);assert Degree(phi) eq 2;
    inv:=Invariants(TorsionSubgroup(E2));
    printf "TORS60 curve %o: E/<T> for x(T) = %o has torsion invariants %o over M\n",i,r,inv;
    assert inv eq [30];
  end for;
end for;
printf "TORS60_DONE: the 6 curves 2-isogenous over M to E_30 or E_30' all have torsion Z/30Z over M (%o s)\n",Cputime(tt);
quit;
