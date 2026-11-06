SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
Attach("work/degree9/shared/fast_hecke.m");Attach("work/degree9/shared/fast_hecke_11.m");Attach("work/degree9/shared/fast_hecke_cyclic_commonzero.m");
for t in [<46,3,1958,5>,<50,3,2040,7>,<51,2,80,5>] do
 n,p,total,q:=Explode(t);X:=MDX1(n,GF(p));S:=NoncuspidalPlaces(X,9);assert #S eq total;P:=S[1];D:=Divisor(P);
 printf "NONIC_NATIVE_CONTROL_INPUT N=%o p=%o degree=%o totalNC9=%o time=%o\n",n,p,Degree(P),#S,Cputime(t0);
 F:=HeckeOperatorCommonZeros(X,q,D);T:=HeckeOperator(X,q,D);assert F eq T;
 B:=T-D-q*DiamondOperator(X,q,D);flag:=IsPrincipal(B);
 printf "NONIC_NATIVE_HECKE_EQUAL N=%o q=%o degree=%o principal=%o time=%o\n",n,q,Degree(T),flag,Cputime(t0);
end for;
print "NONIC_NATIVE_CONTROLS_COMPLETE levels=46,50,51";quit;
