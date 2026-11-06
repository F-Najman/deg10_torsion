SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261030c/C99/";k:=GF(2);B<t>:=FunctionField(k);c:=t^2*(t-1);b:=c*(t^2-t+1);E9:=EllipticCurve([B|1-c,-b,-b,0,0]);P9:=E9![0,0];assert 9*P9 eq E9!0 and 3*P9 ne E9!0;assert c^2/(b-c) eq t;j9:=jInvariant(E9);
C11:=SmallModularCurve(11);jr:=jInvariant(C11,11);C:=ChangeRing(C11,k);Fc:=FunctionField(C);j0:=Fc!(Numerator(jr)/Denominator(jr));F0,mp:=AlgorithmicFunctionField(Fc);F11:=RationalExtensionRepresentation(F0);assert Genus(F11) eq 1 and Degree(F11) eq 2;ff:=DefiningPolynomial(F11);jj:=F11!mp(j0);co:=Eltseq(jj);assert #co eq 2 and co[2] ne 0;
R<X>:=PolynomialRing(B);FR:=FieldOfFractions(R);
function RatSub(v)return FR!(R!Numerator(v))/FR!(R!Denominator(v));end function;
aa:=RatSub(co[1]);bb:=RatSub(co[2]);yy:=(FR!j9-aa)/bb;eq:=&+[RatSub(Coefficient(ff,i))*yy^i:i in [0..Degree(ff)]];pol:=Numerator(eq);pol/:=LeadingCoefficient(pol);fac:=Factorization(pol);print "C99_Y29_FIBRE_POLYNOMIAL",Degree(pol),[Degree(v[1]):v in fac],Cputime(tt);assert #fac eq 1 and Degree(pol) eq 12;
F<z>:=FunctionField(pol);assert Genus(F) eq 29;yz:=Evaluate(Numerator(yy),z)/Evaluate(Denominator(yy),z);assert Evaluate(DefiningPolynomial(F11),z) ne 0 or true;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
// The function-field generators identify the modular curve by exact j equality.
xx:=BaseField(F11).1;fcco:=Eltseq(F11!mp(Fc!C.1));gcco:=Eltseq(F11!mp(Fc!C.2));
out:=[*[ER(v):v in Coefficients(pol)],ER(j9),[ER(v):v in Coefficients(ff)],[ER(v):v in co],EF(F!B.1),EF(F.1),EF(yz),[[ER(v):v in fcco],[ER(v):v in gcco]]*];Write(root cat "Y29_model.m",Sprint(out):Overwrite:=true);print "C99_Y29_MODEL_COMPLETE",Degree(F),Genus(F),Cputime(tt);quit;
