SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);tt:=Cputime();
F<z>:=CyclotomicField(3);E:=EllipticCurve([F|0,0,0,0,-1]);f:=DivisionPolynomial(E,13);fac:=Factorization(f);ds:=[Degree(r[1]):r in fac];assert not 1 in ds and not 5 in ds;
function En(v) return [[Numerator(c),Denominator(c)]:c in Eltseq(F!v)];end function;
Write("work/deg10/continuation_20261029b/C226/j0_twists.m",Sprint([*Coefficients(DefiningPolynomial(F)),[En(c):c in Coefficients(f)],[*[*Degree(r[1]),r[2],[En(c):c in Coefficients(r[1])]*]:r in fac*]*]));
print "J0_QUADRATIC_TWISTS_NO_QUINTIC_13_TORSION",ds,Cputime(tt);quit;
