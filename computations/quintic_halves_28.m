SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
R<t>:=PolynomialRing(Rationals());f:=R![2, 2, -1, -2, -1, 1];assert IsIrreducible(f);K<a>:=NumberField(f);S<x>:=PolynomialRing(K);
b:=Evaluate(R![1032, -239, -222, -759, 419],a);c:=Evaluate(R![42, -11, -8, -31, 17],a);E:=EllipticCurve([1-c,-b,-b,0,0]);P:=E![0,0];N:=28;assert Order(P) eq N;
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
a1,a2,a3,a4,a6:=Explode(aInvariants(E));b2:=a1^2+4*a2;b4:=a1*a3+2*a4;b6:=a3^2+4*a6;b8:=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
halfpoly:=x^4-b4*x^2-2*b6*x-b8-P[1]*(4*x^3+b2*x^2+2*b4*x+b6);facts:=Factorization(halfpoly);
printf "DECIC_HALVING_FACTOR_DEGREES N=%o degrees=%o\n",N,[Degree(r[1]):r in facts];
for rec in facts do
 h:=rec[1];if Degree(h) ne 2 then continue;end if;L<w>:=NumberField(h);EL:=ChangeRing(E,L);PL:=EL!P;
 aa1,aa2,aa3,aa4,aa6:=Explode(aInvariants(EL));iss,ss:=IsSquare((aa1*w+aa3)^2+4*(w^3+aa2*w^2+aa4*w+aa6));
 if iss then
  Q:=EL![w,(-aa1*w-aa3+ss)/2];assert 2*Q in {PL,-PL} and Order(Q) eq 2*N;
  bound:=ReductionBound(E,h,2*N);printf "DECIC_HALVING_WITNESS N=%o point_order=%o torsion_order_bound=%o exact=%o quintic_polynomial=%o relative_polynomial=%o\n",N,2*N,bound,bound eq 2*N,f,h;
 end if;
end for;
print "DECIC_HALVING_SEARCH_COMPLETE";quit;
