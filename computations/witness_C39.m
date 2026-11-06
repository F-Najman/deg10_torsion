SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);t0:=Cputime();
R<t>:=PolynomialRing(Rationals());N:=39;f:=R![1,-4,5,-1,-9,13,-6,-4,7,-4,1];assert Degree(f) eq 10 and IsIrreducible(f);
K<a>:=NumberField(f);bb:=R![248,12,-2110,2122,-284,-5359,4932,801,-3393,1574];cc:=R![11,-51,73,-7,-126,193,-60,-93,100,-32];b:=Evaluate(bb,a);c:=Evaluate(cc,a);
E:=EllipticCurve([1-c,-b,-b,0,0]);assert Discriminant(E) ne 0;
P:=E![0,0,1];assert N*P eq E!0 and &and[(N div ell)*P ne E!0:ell in PrimeDivisors(N)];
bound:=0;checks:=0;
for q in PrimesUpTo(251) do
 if q eq 2 or N mod q eq 0 or not &and[Denominator(x) mod q ne 0:x in Coefficients(f) cat Coefficients(bb) cat Coefficients(cc)] then continue;end if;
 fq:=Polynomial(GF(q),f);
 for rec in Factorization(fq) do
  if Degree(rec[1]) ne 1 or rec[2] ne 1 then continue;end if;
  x:=-Coefficient(rec[1],0)/Coefficient(rec[1],1);
  bq:=Evaluate(Polynomial(GF(q),bb),x);cq:=Evaluate(Polynomial(GF(q),cc),x);
  Eq:=EllipticCurve([1-cq,-bq,-bq,0,0]);if Discriminant(Eq) eq 0 then continue;end if;
  count:=#Eq;assert count mod N eq 0;bound:=GCD(bound,count);checks+:=1;
  printf "DECIC_WITNESS_REDUCTION N=%o p=%o root=%o residue_degree=1 count=%o gcd=%o unramified_simple_root=true\n",N,q,x,count,bound;
  if bound eq N then break;end if;
 end for;
 if bound eq N then break;end if;
end for;
printf "DECIC_WITNESS_RESULT N=%o field_degree=10 point_order=%o torsion_order_bound=%o reductions=%o exact=%o time=%o\n",N,N,bound,checks,bound eq N,Cputime(t0);
assert bound eq N;print "DECIC_WITNESS_COMPLETE";quit;
