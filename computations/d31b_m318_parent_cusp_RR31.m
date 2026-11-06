SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(180);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";
a:=eval Read(root cat"global_cusp_model_p31.m");rep:=eval Read(root cat"endomorphism_gluing/short_cusp_divisors.m");assert a[1]eq 31 and &and[a[4][i]eq rep[1][i]:i in[1..24]] and a[5]eq rep[2];k:=GF(31);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK!x[1]))/(B!(PK!x[2]));end function;
f:=R![DR(x):x in a[2]];assert IsIrreducible(f);F<v>:=FunctionField(f);assert Genus(F)eq 10;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(x)return F![DR(z):z in x];end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
blocks:=[&+[Divisor(Prime(p)):p in ps]:ps in a[3]];assert [Degree(D):D in blocks]eq a[5];out:=[* *];
for row in rep[3]do id,coeff:=Explode(row);D:=&+[coeff[i]*blocks[i]:i in[1..24]];assert Degree(D)eq 5;V,m:=RiemannRochSpace(D);h:=Dimension(V);Append(~out,[*id,coeff,h*]);print "M318_PARENT_CUSP_RR31",id,h,Cputime(tt);end for;
Write(root cat"parent_cusp_RR31.m",Sprint(out):Overwrite:=true);print "D31B_M318_PARENT_CUSP_RR31_COMPLETE",#out,Cputime(tt);quit;
