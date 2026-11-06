SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/C226_resolution_20261030/optimized_cusp_divisors.m");out:=[* *];
for p in [29] do
 X:=MDX11(2,26,GF(p));C:=Curve(X);cs:=[c:c in Cusps(X)|Degree(c) eq 1];b:=X`_coordinates[1];labels:=[Valuation(b,c):c in cs];assert Seqset(labels) eq Seqset(a[1]);
 for row in a[2] do
  v:=row[3];AA:=&+[Maximum(0,v[i])*Divisor(cs[Index(labels,a[1][i])]):i in [1..18]];BB:=&+[Maximum(0,-v[i])*Divisor(cs[Index(labels,a[1][i])]):i in [1..18]];L:=AA-BB;
  hA:=Dimension(RiemannRochSpace(AA));hL:=Dimension(RiemannRochSpace(L));print "Q226_FINITE_PENCIL",p,row[1],Degree(AA),Degree(BB),hA,hL,Cputime(tt);Append(~out,[p,row[1],Degree(AA),Degree(BB),hA,hL]);
 end for;
end for;
Write("work/deg10/C226_resolution_20261030/finite_pencils_29.m",Sprint(out));print "Q226_FINITE_PENCILS_COMPLETE",Cputime(tt);quit;
