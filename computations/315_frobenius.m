SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
AttachSpec("work/deg10/shared/mdmagma/Magma/magma.spec");J:=JH(45,[16]);assert Dimension(J) eq 9;
for p in [2,7,13,19] do f:=FrobeniusPolynomial(J,p);Write("work/deg10/exclusions/315d5/frobenius_data.txt",Sprint([*p,[Integers()!c:c in Coefficients(f)]*]));end for;
print "DECIC315_FROBENIUS_POLYNOMIALS_COMPLETE";quit;
