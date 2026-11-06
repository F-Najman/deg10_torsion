SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
R<X>:=PolynomialRing(Rationals());F<z>:=NumberField(X^2-X+1);S<T>:=PolynomialRing(F);
a:=eval Read("work/deg10/continuation_20261029/M6N12_char0_norm_certificate.m");
function DecodeF(row)
 return F!(Rationals()!(row[1][1]/row[1][2]))+F!(Rationals()!(row[2][1]/row[2][2]))*z;
end function;
rows:=[r:r in a[10]|r[3] and r[4] eq [5] and r[5] eq [2]];assert #rows eq 1;poly:=S![DecodeF(c):c in rows[1][2]];ok,g:=IsSquare(poly);assert ok;g/:=LeadingCoefficient(g);assert Degree(g) eq 5 and IsIrreducible(g);
K<w>:=NumberField(g);assert Degree(K) eq 5 and AbsoluteDegree(K) eq 10;
u:=w;v:=(4*z+1)*w^4+(-z+4)*w^3+(-2*z+10)*w^2+(-13*z+2)*w+2*z-5;
raw:=Split(Read("work/deg10/shared/mdmagma/models_X1_m_n/X1_6_12.txt"));
for line in raw do
 value:=Split(Split(line,"=")[2],";")[1];
 if line[1] eq "X" then equation:=eval(value);end if;
 if line[1] eq "r" then r:=eval(value);end if;
 if line[1] eq "s" then s:=eval(value);end if;
 if line[1] eq "q" then q:=eval(value);end if;
 if line[1] eq "E" then ai:=eval(value);end if;
end for;assert equation eq 0;E:=EllipticCurve([K!c:c in ai]);assert Discriminant(E) ne 0;Q:=E![0,0];assert 12*Q eq E!0 and 6*Q ne E!0 and 4*Q ne E!0;
PP<Y>:=PolynomialRing(K);aa:=aInvariants(E);yy:=Roots(Y^2+(aa[1]*q+aa[3])*Y-q^3-aa[2]*q^2-aa[4]*q-aa[5]);assert #yy eq 2;P:=E![q,yy[1][1]];
assert 6*P eq E!0 and 3*P ne E!0 and 2*P ne E!0;if WeilPairing(P,2*Q,6) ne K!z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq K!z;
print "DECIC612_EXACT_LEVEL",g,MinimalPolynomial(jInvariant(E)),Cputime(tt);
function EncF(x)
 cc:=Eltseq(F!x);return [[Numerator(r),Denominator(r)]:r in cc];
end function;
function EncK(x)
 return [EncF(c):c in Eltseq(K!x)];
end function;
reductions:=[* *];bound:=0;g1:=0;g2:=0;
for p in PrimesUpTo(2000) do
 if p mod 6 ne 1 then continue;end if;
 flat:=&cat[Eltseq(c):c in Coefficients(g)] cat &cat[&cat[Eltseq(c):c in Eltseq(x)]:x in aInvariants(E)];
 if not &and[Denominator(c) mod p ne 0:c in flat] then continue;end if;
 FF:=GF(p);RR<Z>:=PolynomialRing(FF);
 for zr in Roots(Z^2-Z+1) do
  zz:=zr[1];red:=func<x|Evaluate(Polynomial(FF,Eltseq(F!x)),zz)>;gp:=RR![red(c):c in Coefficients(g)];
  for root in Roots(gp) do
   ww:=root[1];if Evaluate(Derivative(gp),ww) eq 0 then continue;end if;
   redK:=func<x|Evaluate(RR![red(c):c in Eltseq(K!x)],ww)>;Er:=EllipticCurve([redK(c):c in aInvariants(E)]);if Discriminant(Er) eq 0 then continue;end if;
   N:=#Er;iv:=Invariants(AbelianGroup(Er));assert #iv eq 2 and iv[1] mod 6 eq 0 and iv[2] mod 12 eq 0;bound:=GCD(bound,N);g1:=GCD(g1,iv[1]);g2:=GCD(g2,iv[2]);Append(~reductions,[*p,Integers()!zz,Integers()!ww,N,iv*]);
   print "DECIC612_GOOD_REDUCTION",reductions[#reductions],bound,[g1,g2],Cputime(tt);
   if bound eq 72 or (g1 eq 6 and g2 eq 12) then break;end if;
  end for;if bound eq 72 or (g1 eq 6 and g2 eq 12) then break;end if;
 end for;if bound eq 72 or (g1 eq 6 and g2 eq 12) then break;end if;
end for;
assert bound eq 72 or (g1 eq 6 and g2 eq 12);
Write("work/deg10/continuation_20261029/M6N12_exact_witness.m",Sprint([*[EncF(c):c in Coefficients(g)],[EncK(c):c in aInvariants(E)],[EncK(c):c in Eltseq(P)],[EncK(c):c in Eltseq(Q)],EncK(jInvariant(E)),reductions,bound,[g1,g2]*]));
print "DECIC612_FULL_TORSION_VERIFIED",[6,12],AbsoluteDegree(K),Cputime(tt);quit;
