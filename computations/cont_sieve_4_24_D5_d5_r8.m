SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_24_root.m");
X:=MDX11(4,24,GF(5));C:=Curve(X);exports:=[* *];remaining:=[];
for dd in [5] do
 xs:=NoncuspidalPlaces(X,dd);rs,labels,keys:=PlacesUpToFullAutM4TwentyFourRoot(X,xs);
 print "DECIC_HECKE_ORBIT_COUNTS",4,24,dd,#xs,#rs,Cputime(t0);
 tests:=[];
 assert #rs eq 11;
 for i in [8] do
  x:=rs[i];killed:=false;
  for q in [13,17,19] do
   if q mod 4 ne 1 or 24 mod q eq 0 or 5 mod q eq 0 then continue;end if;
   A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));assert Degree(A) eq 0;
   pr:=IsPrincipal(A);Append(~tests,[i,q,pr select 1 else 0]);
   print "DECIC_HECKE_SINGLE_SUPPORT",4,24,dd,i,q,pr,Cputime(t0);
   if not pr then killed:=true;break;end if;
  end for;
  if not killed then Append(~remaining,[dd,i]);end if;
 end for;
 Append(~exports,[*dd,keys,labels,[Index(xs,x):x in rs],tests*]);
end for;
Write("work/deg10/continuation_20261029/cont_sieve_4_24_D5_d5_r8_arithmetic.txt",Sprint([*exports,remaining*]));
print "DECIC_SINGLE_SUPPORT_SIEVE_COMPLETE",4,24,remaining,Cputime(t0);quit;
