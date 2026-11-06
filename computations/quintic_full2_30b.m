SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
R<t>:=PolynomialRing(Rationals());f:=R![1, 3, 0, -3, 1, 1];assert IsIrreducible(f);K<a>:=NumberField(f);S<x>:=PolynomialRing(K);
b:=Evaluate(R![295, 238, -382, -38, 190],a);c:=Evaluate(R![28, 9, -36, 13, 7],a);E:=EllipticCurve([1-c,-b,-b,0,0]);P:=E![0,0];N:=30;assert Order(P) eq N;
function ReductionBound(E,h,N)
 K:=BaseRing(E);ff:=DefiningPolynomial(K);coeffs:=&cat[Eltseq(z):z in Coefficients(h)] cat &cat[Eltseq(z):z in aInvariants(E)] cat Coefficients(ff);bound:=0;
 for p in PrimesUpTo(401) do
  if p eq 2 or not &and[Denominator(z) mod p ne 0:z in coeffs] then continue;end if;
  F:=GF(p);fff:=Polynomial(F,ff);
  for rec in Factorization(fff) do
   if Degree(rec[1]) ne 1 or rec[2] ne 1 then continue;end if;r:=-Coefficient(rec[1],0)/Coefficient(rec[1],1);
   red:=func<z|Evaluate(Polynomial(F,Eltseq(K!z)),r)>;
   hh:=Polynomial(F,[red(z):z in Coefficients(h)]);
   if Degree(hh) ne Degree(h) or IsZero(Discriminant(hh)) or #Roots(hh) eq 0 then continue;end if;
   aa:=[red(z):z in aInvariants(E)];a1,a2,a3,a4,a6:=Explode(aa);b2:=a1^2+4*a2;b4:=a1*a3+2*a4;b6:=a3^2+4*a6;b8:=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
   if -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6 eq 0 then continue;end if;
   Eq:=EllipticCurve(aa);z:=#Eq;bound:=GCD(bound,z);printf "DECIC_TOWER_REDUCTION order=%o p=%o quintic_root=%o relative_root=%o count=%o gcd=%o\n",N,p,r,Roots(hh)[1][1],z,bound;
   if bound eq N then return bound;end if;
  end for;
 end for;
 return bound;
end function;
hs:=[r[1]:r in Factorization(DivisionPolynomial(E,2))|Degree(r[1]) eq 2];assert #hs eq 1;h:=hs[1];assert IsIrreducible(h);
L<w>:=NumberField(h);assert AbsoluteDegree(L) eq 10;EL:=ChangeRing(E,L);PL:=EL!P;
twos:=DivisionPoints(EL!0,2);assert #twos eq 4 and Order(PL) eq N;
bound:=ReductionBound(E,h,2*N);assert bound eq 2*N;
printf "DECIC_FULL2_WITNESS N=%o full_group=[2,%o] exact=true quintic_polynomial=%o relative_polynomial=%o\n",N,N,f,h;
print "DECIC_FULL2_WITNESS_COMPLETE";quit;
