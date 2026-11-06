SetColumns(0);SetQuitOnError(true);SetNthreads(1);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
J:=JH(72,[13]);assert Dimension(J) eq 9;rows:=[* *];
for p in [7,13,19] do
 f:=FrobeniusPolynomial(J,p);Append(~rows,[*p,Coefficients(f),Integers()!Evaluate(f,1)*]);
 print "DECIC612_FROBENIUS",p,Coefficients(f),Evaluate(f,1),Cputime(tt);
end for;
Write("work/deg10/continuation_20261029/frob612_arithmetic.txt",Sprint(rows));print "DECIC612_FROBENIUS_COMPLETE";quit;
