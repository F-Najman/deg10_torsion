SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);records:=[* *];left:=[];
for d in [2,3,4] do
 xs:=NoncuspidalPlaces(X,d);rr,labels,keys:=Decic420Orbits(X,xs);print "DECIC420_ORBITS",d,#xs,#rr,Cputime(tt);tests:=[];
 for i in [1..#rr] do
  x:=rr[i];A:=HeckeOperatorFast(X,29,x)-Divisor(x)-29*Divisor(DiamondOperator(X,29,x));assert Degree(A) eq 0;pr:=IsPrincipal(A);
  Append(~tests,[i,29,pr select 1 else 0]);if pr then Append(~left,[d,i]);end if;
  print "DECIC420_A29",d,i,pr,Cputime(tt);
  Write("work/deg10/continuation_20261029/M4N20_lower_partial.m",Sprint([*records,[*d,keys,labels,[Index(xs,x):x in rr],tests*],left*]):Overwrite:=true);
 end for;
 Append(~records,[*d,keys,labels,[Index(xs,x):x in rr],tests*]);
end for;
Write("work/deg10/continuation_20261029/M4N20_lower_arithmetic.txt",Sprint([*records,left*]));
print "DECIC420_LOWER_HECkE_COMPLETE",left,Cputime(tt);quit;
