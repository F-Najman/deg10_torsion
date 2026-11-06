SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/C232_resolution_20261030/hecke11_geometry_v2.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 37;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
function EP(pl) I:=Ideal(pl);return [*Order(I) eq OF select 1 else 0,Degree(pl),[EF(v):v in Basis(I)]*];end function;
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
geom:=eval Read("work/deg10/C232_resolution_20261030/exponent_cusp_geometry.m");A:=Divisor(DP(geom[1]));N:=Divisor(DP(a[5]));
function ED(D) I,J:=Ideals(D);return [*[EF(v):v in Basis(I)],[EF(v):v in Basis(J)]*];end function;
certs:=[* *];
procedure Red(~D,~positive,~certs)
 if positive le 37 then return;end if;T:=D+37*A;basis:=Basis(T:Reduction:=false,Simplification:="None");assert #basis gt 0;h:=basis[1];IT,JT:=Ideals(T);IE:=IT*(h*OF);JE:=JT*(h*OI);assert IsIntegral(IE) and IsIntegral(JE);DD:=Divisor(IE,JE)-37*A;I,J:=Ideals(D);ID,JD:=Ideals(DD);assert ID/(h*OF) eq I and JD/(h*OI) eq J;assert Degree(DD) eq 0;
 Append(~certs,[*ED(D),ED(DD),EF(h)*]);D:=DD;positive:=37;
end procedure;
D:=N-6*A;u:=D;ud:=6;v:=DivisorGroup(F)!0;vd:=0;m:=23280;
while m gt 0 do
 if IsOdd(m) then v+:=u;vd+:=ud;Red(~v,~vd,~certs);end if;m div:=2;if m gt 0 then u:=2*u;ud*:=2;Red(~u,~ud,~certs);end if;
 print "C232_EXPONENT_BINARY_STEP",m,#certs,Cputime(tt);
end while;
basis:=Basis(v:Reduction:=false,Simplification:="None");Write("work/deg10/C232_resolution_20261030/exponent_source.m",Sprint([*ED(v),#basis,certs*]):Overwrite:=true);print "C232_EXPONENT_SOURCE_COMPLETE",#basis,Cputime(tt);quit;
