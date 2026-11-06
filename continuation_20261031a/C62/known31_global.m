SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);tt:=Cputime();
R<t>:=PolynomialRing(Rationals());f:=t^9-2*t^8+t^6-t^5+14*t^4-28*t^3+19*t^2-2*t-1;assert IsIrreducible(f) and IsIrreducible(Polynomial(GF(3),f));
K<a>:=NumberField(f);x:=a;y:=(9*a^8-28*a^7+16*a^6+46*a^5-11*a^4+125*a^3-440*a^2+386*a+150)/119;
r:=(x^2*y-x*y+y-1)/(x^2*y-x);s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];assert Discriminant(E) ne 0 and 31*P eq E!0 and P ne E!0;
printf "D31A_C62_GLOBAL_NONIC31_COMPLETE fielddegree9=true exactpoint31=true irreduciblemod3=true time=%o\n",Cputime(tt);quit;
