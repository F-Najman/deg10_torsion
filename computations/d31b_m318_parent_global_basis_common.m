SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(240);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";a:=eval Read(root cat"global_cusp_model_char0_rational.m");rr:=eval Read(root cat SNAPSHOT);k<z>:=CyclotomicField(3);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK![k!c:c in x[1]]))/(B!(PK![k!c:c in x[2]]));end function;
F<v>:=FunctionField(R![DR(x):x in a[2]]);function DF(x)return F![DR(y):y in x];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(x):x in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
ps:=[* #p eq 0 select 0 else Prime(p[1]):p in a[3]*];out:=[];
for row in rr do id,co,h,bs:=Explode(row);if id ne WHICH then continue;end if;assert h eq 1 and #bs eq 1;A:=&+[co[i]*Divisor(ps[i]):i in[1..24]|co[i]ne 0];assert Degree(A)eq 5;f:=DF(bs[1]);assert f ne 0;I,J:=Ideals(A);assert IsIntegral(I*(f*OF))and IsIntegral(J*(f*OI));Append(~out,[id,h,1,1]);print "M318_PARENT_GLOBAL_BASIS",id,1,1,Cputime(tt);end for;
assert #out eq 1;Write(root cat Sprintf("parent_global_basis_%o.m",WHICH),Sprint(out):Overwrite:=true);print "D31B_M318_PARENT_GLOBAL_BASIS_COMPLETE",WHICH,Cputime(tt);quit;
