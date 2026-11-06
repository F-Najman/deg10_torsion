SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(7*10^9);Alarm(400);tt:=Cputime();k:=GF(5);R<x>:=FunctionField(k);S<y>:=PolynomialRing(R);root:="work/deg10/continuation_20261031b/C48/";m:=eval Read(root cat "hecke7_s0_model.m");a:=eval Read(root cat "cusp_corrections_v2.m");
function DR(v)return (PolynomialRing(k)!v[1])(x)/(PolynomialRing(k)!v[2])(x);end function;
F<t>:=FunctionField(S![DR(v):v in m[1]]);assert Genus(F)eq 37;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(z):z in v];end function;
b,c:=Explode([DF(v):v in m[3]]);E:=EllipticCurve([F|1-c,-b,-b,0,0]);j:=jInvariant(E);pds:=[* *];
for i in[1,2]do ps:=[];for pd in a[1][i]do O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);p:=Place(I);assert Degree(p)eq pd[1]and pd[1]eq 2*i;assert -Valuation(j,p)eq (i eq 1 select 2 else 1);Append(~ps,p);end for;assert ps[1]ne ps[2];B:=12*(Divisor(ps[1])-Divisor(ps[2]));assert a[2][i][1][1]eq 1 and a[2][i][1][2]eq 1;v:=DF(a[2][i][1][3]);IB,JB:=Ideals(B);assert IB*(v*OF)eq ideal<OF|1> and JB*(v*OI)eq ideal<OI|1>;assert Divisor(v)eq -B;Append(~pds,[i,2*i,12,1]);print "D31B_C48_NATIVE_CUSP_PRINCIPAL_PASS",i,Degree(B),Cputime(tt);end for;
Write(root cat "independent_cusp_principal.m",Sprint(pds):Overwrite:=true);print "D31B_C48_NATIVE_CUSP_PRINCIPAL_COMPLETE",2,Cputime(tt);quit;
