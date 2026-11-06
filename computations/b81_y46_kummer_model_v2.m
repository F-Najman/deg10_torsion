SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261030b/C81/";k<eta>:=GF(4);Y:=MDX1(27,k);C:=Curve(Y);F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 13;E:=EllipticCurve([F!mp(v):v in Y`_E]);P:=E![0,0];S:=9*P;aa:=aInvariants(E);xs:=S[1];ys:=S[2];m:=(3*xs^2+2*aa[2]*xs+aa[4]-aa[1]*ys)/(2*ys+aa[1]*xs+aa[3]);assert 3*xs+aa[2]-m^2-aa[1]*m eq 0;fk:=P[2]-ys-m*(P[1]-xs);assert fk ne 0;
RX<X>:=PolynomialRing(F);EG,phi:=IsogenyFromKernel(E,X-xs);dual:=DualIsogeny(phi);dp:=DefiningPolynomials(dual);nx:=Evaluate(dp[1],[X,RX!0,RX!1]);dx:=Evaluate(dp[3],[X,RX!0,RX!1]);g:=GCD(nx,dx);nx div:=g;dx div:=g;nx/:=LeadingCoefficient(nx);assert Degree(nx) eq 3 and IsIrreducible(nx);aa2:=Coefficient(nx,2);pp:=aa2^2+Coefficient(nx,1);qq:=aa2*Coefficient(nx,1)+Coefficient(nx,0);assert qq ne 0;
print "C81_KUMMER_BASE_AND_DUAL_READY",Cputime(tt);rr:=Roots(X^2+qq*X+pp^3);assert #rr eq 2;found:=false;chosenPower:=0;chosenA:=F!0;chosenB:=F!0;trials:=[* *];
for r in rr do
 U:=r[1];if U eq 0 then continue;end if;
 for power in [1,2] do
  roots:=Roots(X^3-U/fk^power);Append(~trials,[power,#roots]);if #roots eq 0 then continue;end if;A:=roots[1][1];B:=pp/(A*fk^power);assert A^3*fk^power+B^3*fk^(2*power) eq qq;chosenPower:=power;chosenA:=A;chosenB:=B;found:=true;break;
 end for;if found then break;end if;
end for;assert found;A:=chosenA;B:=chosenB;print "C81_KUMMER_CUBE_CLASS_VERIFIED",chosenPower,trials,Cputime(tt);
RL<W>:=PolynomialRing(F);L<w>:=FunctionField(W^3-fk);uz:=chosenPower eq 1 select w else w^2;xx:=L!aa2+(L!A)*uz+(L!B)*uz^2;assert Evaluate(nx,xx) eq 0;assert Degree(MinimalPolynomial(xx,F)) eq 3;
Df:=Divisor(fk);ps,es:=Support(Df);branch:=[p:p in ps|Valuation(fk,p) mod 3 ne 0];bdeg:=&+[Degree(p):p in branch];assert 1+3*(Genus(F)-1)+bdeg eq 46;jf:=F!mp(jInvariantMap(Y));assert &and[Valuation(jf,p) lt 0:p in branch];
function EC(c) return [Integers()!a:a in Eltseq(k!c,GF(2))];end function;
function ER(v) return [[EC(c):c in Coefficients(Numerator(v))],[EC(c):c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);out:=[*[ER(c):c in Coefficients(DefiningPolynomial(F))],[EF(F!mp(v)):v in Y`_coordinates_xy],EF(fk),[EF(v):v in Coefficients(nx)],chosenPower,EF(A),EF(B),EF(aa2),[[*Degree(p),Valuation(fk,p),Order(Ideal(p)) eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*]:p in ps],[EF(v):v in aInvariants(E)]*];Write(root cat "y46_kummer_model.m",Sprint(out):Overwrite:=true);print "C81_Y46_KUMMER_MODEL_COMPLETE",Degree(F),Genus(F),bdeg,[Degree(p):p in branch],chosenPower,Cputime(tt);quit;
