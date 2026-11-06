SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_4_20.m");
X:=MDX11(4,20,GF(9));xs:=NoncuspidalPlaces(X,4);rr,labels,keys:=Decic420Orbits(X,xs);assert #rr eq 12;tests:=[];left:=[];
for i in [1,2,4,12] do
 killed:=false;
 for q in [13,17,41] do
  x:=rr[i];A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;pr:=IsPrincipal(A);Append(~tests,[i,q,pr select 1 else 0]);print "DECIC420_QUARTIC_FOLLOWUP",i,q,pr,Cputime(tt);
  if not pr then killed:=true;break;end if;
 end for;if not killed then Append(~left,i);end if;
end for;
Write("work/deg10/continuation_20261029/M4N20_quartic_followup_arithmetic.txt",Sprint([*4,keys,labels,[Index(xs,x):x in rr],tests,left*]));print "DECIC420_QUARTIC_FOLLOWUP_COMPLETE",left,Cputime(tt);quit;
