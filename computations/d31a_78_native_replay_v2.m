SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);Alarm(850);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031a/C78/";
modelDat:=eval Read(root cat "model39.m");raw:=eval Read(root cat "A5_RR.m");src:=eval Read(root cat "source39.m");k:=GF(2);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<Z>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;ff:=R![DR(q):q in modelDat[1]];assert Degree(ff) eq 17 and IsIrreducible(ff);F<u>:=FunctionField(ff);assert Genus(F) eq 33;
function DF(v)return F![DR(q):q in v];end function;xy:=[DF(v):v in modelDat[2]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);zero:=DivisorGroup(F)!0;
function Prime(pd)O:=pd[2] eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I) and IsPrime(I);P:=Place(I);assert Degree(P) eq pd[1] and Ideal(P) eq I;a,b:=TwoGenerators(P);assert ideal<O|a,b> eq I;return P;end function;
function Key2(vals)
 L:=Parent(vals[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);found:=false;
 for bits in[[1,0],[0,1],[1,1]]do g:=&+[L!bits[i]*vals[i]:i in[1..2]];mp:=MinimalPolynomial(g,k);if Degree(mp) eq d then chosen:=bits;found:=true;break;end if;end for;assert found;
 M:=Matrix(k,[Eltseq(g^i,k):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vals];return[*d,chosen,[Integers()!q:q in Coefficients(mp)],[[Integers()!q:q in row]:row in co]*];
end function;
out:=[* *];npr:=0;
for row in raw do
 ci:=row[1];d:=row[2];assert src[3][ci][1] eq d;poly:=[r[2]:r in src[1]|r[1] eq d][1];L<z>:=ext<k|PP!poly>;ky:=src[3][ci][2];b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([L|1-c,-b,-b,0,0]);Q:=E![0,0];assert 39*Q eq E!0 and 13*Q ne E!0 and 3*Q ne E!0;
 P:=Prime(row[3]);Da:=Prime(row[4]);npr+:=2;assert Degree(P) eq d and Degree(Da) eq d;kp:=Key2([Evaluate(v,P):v in xy]);kd:=Key2([Evaluate(v,Da):v in xy]);assert kp eq Key2(TateNormalForm_xy(E,Q)) and kd eq Key2(TateNormalForm_xy(E,5*Q));T:=zero;hec:=[* *];assert &+[r[1]:r in row[5]] eq 6;
 for h in row[5] do Z:=Prime(h[3]);npr+:=1;assert d*h[1] eq h[2]*Degree(Z);T+:=h[2]*Divisor(Z);Append(~hec,[*h[1],h[2],Key2([Evaluate(v,Z):v in xy])*]);end for;
 assert Degree(T) eq 6*d;D:=T-Divisor(P)-5*Divisor(Da);assert Degree(D) eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");assert h0 eq row[6];Append(~out,[*ci,d,kp,kd,hec,h0*]);printf "D31A_C78_NATIVE_RR_CASE ci=%o degree=%o h0=%o time=%o\n",ci,d,h0,Cputime(tt);
end for;
assert #out eq 10 and [v[6]:v in out] eq [1,0,0,0,0,0,0,0,0,1];Write(root cat "native_replay.m",Sprint(out):Overwrite:=true);printf "D31A_C78_NATIVE_RR_COMPLETE tests=10 primechecks=%o time=%o\n",npr,Cputime(tt);quit;
