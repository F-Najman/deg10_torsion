SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_4_20.m");
X:=MDX11(4,20,GF(9));C:=Curve(X);xs:=NoncuspidalPlaces(X,2);ys:=NoncuspidalPlaces(X,3);assert #xs eq 96 and #ys eq 320;
rr,labels,keys:=Decic420Orbits(X,xs);assert #rr eq 2;F:=BaseRing(X);K:=ext<F|3>;ykeys:=[M4SixteenRelativeKey([Evaluate(f,y):f in X`_coordinates],K,F,3):y in ys];assert #Seqset(ykeys) eq #ys;
function AQ(X,q,x)
 return HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));
end function;
rep29:=[AQ(X,29,x):x in rr];rep13:=[AQ(X,13,x):x in rr];tests:=[];left:=[];
for typ in [22,23] do
 zs:=typ eq 22 select xs else ys;q:=typ eq 22 select 29 else 13;aa:=typ eq 22 select rep29 else rep13;
 for j in [1..#zs] do
  if (j-1) mod 8 ne 1 then continue;end if;
  B:=AQ(X,q,zs[j]);
  for i in [1..#rr] do
   A:=aa[i]+B;assert Degree(A) eq 0;pr:=IsPrincipal(A);Append(~tests,[typ,i,j,q,pr select 1 else 0]);if pr then Append(~left,[typ,i,j]);end if;
   print "DECIC420_MIXED_TEST",typ,i,j,q,pr,Cputime(tt);
  end for;
 end for;
end for;
Write("work/deg10/continuation_20261029/cont_sieve_4_20_mixed_shard1_arithmetic.txt",Sprint([*keys,labels,[Index(xs,x):x in rr],ykeys,tests,left*]));print "DECIC420_MIXED_SHARD_COMPLETE",1,left,Cputime(tt);quit;
