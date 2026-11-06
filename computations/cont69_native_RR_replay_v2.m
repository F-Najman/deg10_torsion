SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
m:=eval Read("work/deg10/next_cases_20261030/C69/cover_model_v2.m");k:=GF(2);PP:=PolynomialRing(k);U<x>:=FunctionField(k);R<Y>:=PolynomialRing(U);
function DR(v) return U!(PP!v[1])/U!(PP!v[2]);end function;
K<y>:=FunctionField(R![DR(v):v in m[1]]);assert Genus(K) eq 2;
function DK(v) return K![DR(a):a in v];end function;
RA<T>:=PolynomialRing(K);AF<t>:=FunctionField(RA![DK(v):v in m[2]]);assert Genus(AF) eq 7;
function DA(v) return AF![DK(a):a in v];end function;
RY<W>:=PolynomialRing(AF);CY<z>:=FunctionField(RY![DA(v):v in m[3]]);assert Genus(CY) eq 13;
function DY(v) return CY![DA(a):a in v];end function;
OAF:=MaximalOrderFinite(AF);OAI:=MaximalOrderInfinite(AF);OY:=MaximalOrderFinite(CY);OYI:=MaximalOrderInfinite(CY);
function BaseDiv(rows)
 D:=DivisorGroup(AF)!0;
 for row in rows do p,e:=Explode(row);O:=p[1] eq 1 select OAF else OAI;J:=ideal<O|[DA(v):v in p[3]]>;assert IsPrime(J);Q:=Place(J);assert Degree(Q) eq p[2];D+:=e*Divisor(Q);end for;assert Degree(D) eq 0;return D;
end function;
function CoverDiv(rows)
 D:=DivisorGroup(CY)!0;
 for row in rows do p,e:=Explode(row);O:=p[1] eq 1 select OY else OYI;J:=ideal<O|[DY(v):v in p[3]]>;assert IsPrime(J);Q:=Place(J);assert Degree(Q) eq p[2];D+:=e*Divisor(Q);end for;assert Degree(D) eq 0;return D;
end function;
a:=eval Read("work/deg10/next_cases_20261030/C69/cover_results_v2.m");results:=[* *];
for row in a[4] do D:=BaseDiv(row[3]);h:=#Basis(D:Reduction:=false,Simplification:="None");assert h eq(row[4] select 1 else 0);Append(~results,[*5,7,row[1],row[2],0,h*]);end for;
for row in a[5] do D:=CoverDiv(row[4]);h:=#Basis(D:Reduction:=false,Simplification:="None");assert h eq(row[5] select 1 else 0);Append(~results,[*5,13,row[1],row[2],row[3],h*]);print "CONT69_NATIVE_RR_REPLAY_PROGRESS",#results,h,Cputime(tt);end for;
Write("work/deg10/next_cases_20261030/C69/native_RR_replay_v2.m",Sprint(results):Overwrite:=true);print "CONT69_NATIVE_RR_REPLAY_V2_COMPLETE",#results,Cputime(tt);quit;
