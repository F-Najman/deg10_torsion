SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(600);clockStart:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";outdir:=root cat "endomorphism_gluing/";
a:=eval Read(root cat"global_function_field_char0.m");rr:=eval Read(outdir cat"global_RR_section13900.snapshot.m");row:=[v:v in rr|v[1]eq 13900][1];assert row[3]eq 1 and #row[4]eq 1;k<zet>:=CyclotomicField(3);PK<U>:=PolynomialRing(k);B<u>:=FunctionField(k);R<V>:=PolynomialRing(B);
function DR(x)return (B!(PK![k!c:c in x[1]]))/(B!(PK![k!c:c in x[2]]));end function;
F<v>:=FunctionField(R![DR(x):x in a[1]]);function DF(x)return F![DR(y):y in x];end function;
f:=DF(row[4][1]);assert f ne 0;nf:=Norm(f);num:=PK!Numerator(nf);den:=PK!Denominator(nf);removed:=[];
for c in[k!0,k!-1,zet,zet^2]do ell:=U-c;e:=0;while IsDivisibleBy(num,ell)do num div:=ell;e+:=1;end while;Append(~removed,e);end for;
assert Degree(num)eq 5;g:=num/LeadingCoefficient(num);assert IsIrreducible(g);print "M318_GLOBAL_QUINTIC_EXTRACTED",Degree(F),g,removed,Cputime(clockStart);L<u0>:=ext<k|g>;assert Degree(L)eq 5;
PL<W>:=PolynomialRing(L);
function SpecialPolynomial(cs)
 dd:=LCM([Denominator(c):c in cs]);pol:=PL![Evaluate(PK!Numerator(dd*c),u0):c in cs];return pol,Evaluate(PK!dd,u0);
end function;
ff,fd:=SpecialPolynomial(Coefficients(DefiningPolynomial(F)));nf0,fd0:=SpecialPolynomial(Eltseq(f));assert ff ne 0 and nf0 ne 0 and fd ne 0 and fd0 ne 0;gg:=GCD(ff,nf0);assert Degree(gg)eq 1;v0:=-Coefficient(gg,0)/Coefficient(gg,1);assert Evaluate(ff,v0)eq 0 and Evaluate(nf0,v0)eq 0;
function EvalF(z)
 cs:=Eltseq(F!z);out:=L!0;for i in[1..#cs]do out+:=Evaluate(PK!Numerator(cs[i]),u0)/Evaluate(PK!Denominator(cs[i]),u0)*v0^(i-1);end for;return out;
end function;
qt:=[EvalF(DF(x)):x in a[2]];q,t:=Explode(qt);EE:=EllipticCurve([q*t-q+t+2,L!0,q*t^2-q*t+t,L!0,L!0]);PP:=EE![0,0];QQ:=EE![-t,t^2];assert 3*PP eq EE!0 and PP ne EE!0 and 18*QQ eq EE!0 and 9*QQ ne EE!0 and 6*QQ ne EE!0 and PP ne 6*QQ and PP ne -6*QQ;wp:=WeilPairing(PP,6*QQ,3);assert wp ne 1 and wp^3 eq 1;assert EvalF(DF(a[3]))eq jInvariant(EE);
function Conj(c)z:=Eltseq(k!c);return z[1]+z[2]*zet^2;end function;
found:=false;PQ<X>:=PolynomialRing(Rationals());
for c in[1..5]do h:=Evaluate(g,U-c*zet);hc:=PK![Conj(z):z in Coefficients(h)];pabs:=PQ![Rationals()!z:z in Coefficients(h*hc)];assert Degree(pabs)eq 10;theta:=u0+c*L!zet;assert Evaluate(pabs,theta)eq 0;if IsIrreducible(pabs)then found:=true;chosen:=c;break;end if;end for;assert found;
function Flat(z)return &cat[Eltseq(k!c):c in Eltseq(L!z)];end function;
M:=Matrix(Rationals(),[Flat(theta^i):i in[0..9]]);assert Rank(M)eq 10;
function Rat(z)return[Numerator(Rationals()!z),Denominator(Rationals()!z)];end function;
function Abs(z)c:=Eltseq(Solution(M,Vector(Rationals(),Flat(z))));assert &+[L!c[i]*theta^(i-1):i in[1..10]]eq z;return[Rat(w):w in c];end function;
function Rel(z)return[[Rat(w):w in Eltseq(k!c)]:c in Eltseq(L!z)];end function;
ptpayload:=[*13900,[Rat(c):c in Coefficients(pabs)],Abs(L!zet),[*Abs(q),Abs(t)*],[Abs(c):c in aInvariants(EE)],[Abs(PP[1]),Abs(PP[2])],[Abs(QQ[1]),Abs(QQ[2])],Abs(jInvariant(EE)),[*Abs(u0),Abs(v0)*],[[Rat(w):w in Eltseq(k!c)]:c in Coefficients(g)],Abs(wp),chosen*];
Write(outdir cat"global_point13900_absolute.m",Sprint(ptpayload):Overwrite:=true);Write(outdir cat"global_point13900_extraction.m",Sprint([*row[1],removed,[[Eltseq(k!c):c in Coefficients(Numerator(nf))],[Eltseq(k!c):c in Coefficients(Denominator(nf))]],Rel(u0),Rel(v0),Degree(gg),Rel(q),Rel(t)*]):Overwrite:=true);print "D31B_M318_GLOBAL_POINT13900_COMPLETE",Degree(g),Degree(pabs),chosen,[Integers()!Order(PP),Integers()!Order(QQ)],Cputime(clockStart);quit;
