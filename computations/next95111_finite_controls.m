SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);tt:=Cputime();all:=[* *];PP:=PolynomialRing(GF(2));
for d in [1..10] do
 K:=GF(2^d);rows:=[* *];orders:={Integers()|};
 for j in K do
  es:=Twists(EllipticCurveWithjInvariant(j));
  for E in es do
   a:=aInvariants(E);n:=#E;assert n mod 95 ne 0 and n mod 111 ne 0;
   for v in a cat [j] do assert Evaluate(PP!Eltseq(v,GF(2)),K.1) eq v;end for;
   Append(~rows,[* [Integers()!x:x in Eltseq(j,GF(2))],[[Integers()!x:x in Eltseq(v,GF(2))]:v in a],n*]);Include(~orders,n);
  end for;
 end for;
 assert Degree(DefiningPolynomial(K,GF(2))) eq d;
 Append(~all,[*d,[Integers()!x:x in Coefficients(DefiningPolynomial(K,GF(2)))],rows*]);
 Write("work/deg10/next_cases_20261030/dyadic_local/all_curves.m",Sprint(all):Overwrite:=true);
 print "NEXT95111_FINITE_DEGREE_COMPLETE",d,#rows,Sort(Setseq(orders)),Cputime(tt);
end for;
print "NEXT95111_FINITE_EXHAUSTION_COMPLETE",Cputime(tt);quit;
