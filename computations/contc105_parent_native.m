SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(500);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261030c/C105/";p:=eval Read(root cat "quotient35_RR.m");src:=eval Read(root cat "source105.m");
k:=GF(2);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
F<y>:=FunctionField(R![DR(v):v in p[1]]);assert Degree(F) eq 15 and Genus(F) eq 25;
function DF(v)return F![DR(c):c in v];end function;
xy:=[DF(v):v in p[2]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function Prime(row)
 O:=row[2] eq 1 select OF else OI;I:=ideal<O|[DF(v):v in row[3]]>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq row[1];return Q;
end function;
function Key2(v)
 L:=Parent(v[1]);d:=LCM([Degree(MinimalPolynomial(s,k)):s in v]);found:=false;
 for b in [[1,0],[0,1],[1,1]] do
  g:=&+[L!b[i]*v[i]:i in[1..2]];mp:=MinimalPolynomial(g,k);
  if Degree(mp) eq d then bits:=b;found:=true;break;end if;
 end for;assert found;B:=Matrix(k,[Eltseq(g^i,k):i in[0..d-1]]);assert Rank(B) eq d;
 co:=[Eltseq(Solution(B,Vector(k,Eltseq(s,k)))):s in v];assert &and[&+[L!co[j][i+1]*g^i:i in[0..d-1]] eq v[j]:j in[1..2]];
 return [*d,bits,[Integers()!c:c in Coefficients(mp)],[[Integers()!c:c in r]:r in co]*];
end function;
pol:=[v[2]:v in src[1]|v[1] eq 10][1];L<z>:=ext<k|PP!pol>;out:=[* *];npr:=0;
for row in p[3] do
 ci:=row[1];d,ky:=Explode(src[3][ci]);assert d eq 10;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([L|1-c,-b,-b,0,0]);P:=E![0,0];assert Order(P) eq 105;Q:=3*P;assert Order(Q) eq 35;
 ps:=Prime(row[2]);pd:=Prime(row[4]);npr+:=2;assert row[3] eq 1 and row[5] eq 3 and Degree(ps) eq 10 and Degree(pd) eq 10;
 ks:=Key2([Evaluate(v,ps):v in xy]);kd:=Key2([Evaluate(v,pd):v in xy]);assert ks eq Key2(TateNormalForm_xy(E,Q));assert kd eq Key2(TateNormalForm_xy(E,3*Q));
 T:=DivisorGroup(F)!0;obs:=[* *];assert &+[r[1]:r in row[6]] eq 4;
 for h in row[6] do
  qp:=Prime(h[3]);npr+:=1;assert 10*h[1] eq h[2]*Degree(qp);T+:=h[2]*Divisor(qp);Append(~obs,[*h[1],h[2],Degree(qp),Key2([Evaluate(v,qp):v in xy])*]);
 end for;
 assert Degree(T) eq 40;B:=T-row[3]*Divisor(ps)-row[5]*Divisor(pd);assert Degree(B) eq 0;h0:=#Basis(B:Reduction:=false,Simplification:="None");assert h0 eq 0 and row[7] eq 0;
 Append(~out,[*ci,ks,kd,obs,h0*]);print "C105_PARENT_NATIVE_CASE",ci,h0,Cputime(tt);
end for;
assert #out eq 6 and npr eq 24;Write(root cat "parent_native_RR.m",Sprint(out):Overwrite:=true);print "C105_PARENT_NATIVE_COMPLETE",6,npr,Cputime(tt);quit;
