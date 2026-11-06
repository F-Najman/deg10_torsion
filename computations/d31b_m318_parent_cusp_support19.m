SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(180);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";
a:=eval Read(root cat"global_cusp_model_p19.m");rep:=eval Read(root cat"endomorphism_gluing/short_cusp_divisors.m");assert a[1]eq 19 and &and[a[4][i]eq rep[1][i]:i in[1..24]] and a[5]eq rep[2];k:=GF(19);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK!x[1]))/(B!(PK!x[2]));end function;
f:=R![DR(x):x in a[2]];assert IsIrreducible(f);F<v>:=FunctionField(f);assert Genus(F)eq 10;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(x)return F![DR(z):z in x];end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
blocks:=[&+[Divisor(Prime(p)):p in ps]:ps in a[3]];assert [Degree(D):D in blocks]eq a[5];qt:=[DF(z):z in a[6]];q,t:=Explode(qt);E:=EllipticCurve([q*t-q+t+2,F!0,q*t^2-q*t+t,F!0,F!0]);jf:=jInvariant(E);out:=[* *];
for row in rep[3]do id,coeff:=Explode(row);D:=&+[coeff[i]*blocks[i]:i in[1..24]];assert Degree(D)eq 5;V,m:=RiemannRochSpace(D);h:=Dimension(V);assert h eq 1;Z:=Divisor(m(V.1))+D;ps,es:=Support(Z);assert Degree(Z)eq 5 and &and[e gt 0:e in es];checks:=[* [*Degree(ps[i]),es[i],Valuation(jf,ps[i])*]:i in[1..#ps]*];Append(~out,[*id,coeff,h,checks*]);print "M318_PARENT_CUSP_SUPPORT19",id,checks,Cputime(tt);end for;
Write(root cat"parent_cusp_support19.m",Sprint(out):Overwrite:=true);print "D31B_M318_PARENT_CUSP_SUPPORT19_COMPLETE",#out,Cputime(tt);quit;
