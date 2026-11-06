SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_6_12.m");
X:=MDX612(GF(7));C:=Curve(X);exports:=[* *];left:=[];
for d in [3..5] do
 xs:=NoncuspidalPlaces(X,d);rs,labels,keys:=Decic612Orbits(X,xs);print "DECIC612_ORBITS",d,#xs,#rs,Cputime(tt);tests:=[];
 for i in [1..#rs] do
  killed:=false;x:=rs[i];
  for q in [13] do
   A:=HeckeOperator612(X,q,x)-Divisor(x)-q*Divisor(Diamond612(X,q,x));assert Degree(A) eq 0;
   pr:=IsPrincipal(A);Append(~tests,[i,q,pr select 1 else 0]);print "DECIC612_HECKE_TEST",d,i,q,pr,Cputime(tt);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~left,[d,i]);end if;
 end for;
 Append(~exports,[*d,keys,labels,[Index(xs,x):x in rs],tests*]);
end for;
Write("work/deg10/continuation_20261029/sieve612_D5_v2_arithmetic.txt",Sprint([*exports,left*]));
print "DECIC612_SINGLE_SUPPORT_COMPLETE",left,Cputime(tt);quit;
