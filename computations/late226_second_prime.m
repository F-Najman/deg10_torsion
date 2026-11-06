SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/continuation_20261029b/C226/small_cusp_divisors.m");out:=[* *];
for p in [5,7,11] do
 X:=MDX11(2,26,GF(p));C:=Curve(X);assert Genus(C) eq 25;cs:=[c:c in Cusps(X)|Degree(c) eq 1];b:=X`_coordinates[1];labels:=[Valuation(b,c):c in cs];assert #cs eq 18 and Seqset(labels) eq Seqset(a[1]);
 for row in a[2] do
  DD:=&+[row[3][i]*Divisor(cs[Index(labels,a[1][i])]):i in [1..18]];assert Degree(DD) eq 10;V,mu:=RiemannRochSpace(DD);dim:=Dimension(V);Append(~out,[*p,row[1],row[2],dim*]);print "C226_SECOND_PRIME_RR",p,row[1],dim,Cputime(tt);
 end for;
end for;
Write("work/deg10/continuation_20261029b/C226/second_prime_RR.m",Sprint(out));print "C226_SECOND_PRIME_RR_COMPLETE",Cputime(tt);quit;
