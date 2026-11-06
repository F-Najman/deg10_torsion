SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(4*10^9);Alarm(120);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";a:=eval Read(root cat"global_function_field_char0.m");k<z>:=CyclotomicField(3);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK![k!c:c in x[1]]))/(B!(PK![k!c:c in x[2]]));end function;
F<v>:=FunctionField(R![DR(x):x in a[1]]);print "M318_DECOMP_RESTORE",Cputime(tt);
p:=Zeros(B!u)[1];print "M318_DECOMP_LOCAL_BEGIN",Cputime(tt);dec:=Decomposition(F,p);print "D31B_M318_LOCAL_DECOMP_COMPLETE",dec,Cputime(tt);quit;
