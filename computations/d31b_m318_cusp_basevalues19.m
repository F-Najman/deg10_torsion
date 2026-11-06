SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(180);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";
a:=eval Read(root cat"global_cusp_model_p19.m");rep:=eval Read(root cat"endomorphism_gluing/short_cusp_divisors.m");assert a[1]eq 19 and &and[a[4][i]eq rep[1][i]:i in[1..24]] and a[5]eq rep[2];k:=GF(19);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK!x[1]))/(B!(PK!x[2]));end function;
f:=R![DR(x):x in a[2]];assert IsIrreducible(f);F<v>:=FunctionField(f);assert Genus(F)eq 10;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(x)return F![DR(z):z in x];end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
out:=[* *];
for i in[1..24]do if a[5][i]ne 1 then continue;end if;P:=Prime(a[3][i][1]);v:=Valuation(F!u,P);if v lt 0 then val:=[*"infinity",v*];else val:=[*Integers()!Evaluate(F!u,P),v*];end if;Append(~out,[*i,val*]);end for;
Write(root cat"rational_cusp_basevalues19.m",Sprint(out):Overwrite:=true);print out;print "D31B_M318_CUSP_BASEVALUES_COMPLETE";quit;
