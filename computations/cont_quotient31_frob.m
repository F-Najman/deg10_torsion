SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
J:=JH(31,[2]);assert Dimension(J) eq 6;
rows:=[* *];
for p in [2,3,5,7,11,13] do
 f:=FrobeniusPolynomial(J,p);print "DECIC_Y31_FROBENIUS",p,Evaluate(f,1),Coefficients(f);Append(~rows,[*p,Coefficients(f)*]);
end for;
Write("work/deg10/continuation_20261029/quotient31_frob.txt",Sprint(rows));
print "DECIC_Y31_FROBENIUS_COMPLETE";quit;
