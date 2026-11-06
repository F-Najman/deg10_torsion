SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(700);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");k<w>:=GF(4);X:=MDX1(33,k);FF:=FunctionField(Curve(X));F,mp:=AlgorithmicFunctionField(FF);B:=BaseField(F);E:=EllipticCurve([mp(a):a in X`_E]);Q:=E![mp(a):a in X`_P];assert 33*Q eq E!0;printf "M333_COVER_BASE_READY g=%o degree=%o cpu=%o\n",Genus(F),Degree(F),Cputime(tt);
PZ<Z>:=PolynomialRing(F);A3pt:=11*Q;a1,a2,a3,a4,a6:=Explode(aInvariants(E));x3:=A3pt[1];y3:=A3pt[2];s3:=(3*x3^2+2*a2*x3+a4-a1*y3)/(2*y3+a1*x3+a3);A1:=a1+2*s3;A3:=a3+a1*x3+2*y3;dd:=A1^3-27*A3;f3:=Z^3-dd;print "M333_KUMMER_POLYNOMIAL_READY",Cputime(tt);L<r>:=FunctionField(f3:Check:=false);rx:=L!x3+r^2+(L!A1)*r+L!A1^2;assert Evaluate(PolynomialRing(L)!DivisionPolynomial(E,3),rx) eq 0;print "M333_KUMMER_IDENTITY_VERIFIED",Cputime(tt);assert Genus(L) eq 81;print "M333_KUMMER_COVER_GENUS_READY",Genus(L),Cputime(tt);
XY:=TateNormalForm_xy(E,10*Q);
function Sub(f)
 ff:=f@@mp;nn:=Numerator(ff);dd:=Denominator(ff);args:=Rank(Parent(nn)) eq 3 select XY cat [F!1] else XY;return Evaluate(nn,args)/Evaluate(dd,args);
end function;
up:=Sub(F!B.1);vp:=Sub(F.1);bm:=hom<B->L|L!up>;fm:=hom<F->L|bm,L!vp>;
A:=10*Q;a1,a2,a3,a4,a6:=Explode(aInvariants(E));xr:=A[1];yr:=A[2];ss:=(3*xr^2+2*a2*xr+a4-a1*yr)/(2*yr+a1*xr+a3);a2p:=a2+3*xr-a1*ss-ss^2;a3p:=a3+a1*xr+2*yr;sc:=a2p/a3p;rp:=(L!sc)*r;tau:=hom<L->L|fm,rp>;
for a in [L.1,L!F.1,L!(F!B.1)] do assert tau(tau(a)) eq a;end for;assert tau(L.1) ne L.1 or tau(L!F.1) ne L!F.1;print "M333_DELTA10_INVOLUTION_READY",Cputime(tt);
j:=L!jInvariant(E);unitfs:=[mp(v):v in X`_coordinates_xy];unitfs cat:=[1-v:v in unitfs];basecs:=Setseq(&join[Seqset(Support(Divisor(v))):v in unitfs]);basecs:=[p:p in basecs|Valuation(jInvariant(E),p) lt 0];assert &+[Degree(p):p in basecs] eq 40;print "M333_BASE_CUSPS_FROM_SMALL_FUNCTIONS",[Degree(p):p in basecs],Cputime(tt);cs:=[];for bP in basecs do
 aps:=Decomposition(L,bP);a,b:=TwoGenerators(bP);above:=[p:p in aps|Valuation(L!a,p) gt 0 and Valuation(L!b,p) gt 0];assert &+[Valuation(jInvariant(E),bP)*0+Valuation(L!UniformizingElement(bP),p)*Degree(p):p in above] eq 3*Degree(bP);cs cat:=above;print "M333_KUMMER_CUSP_FIBRE",Degree(bP),[Degree(p):p in above],Cputime(tt);end for;assert {*Degree(P):P in cs*} eq {*1^^40,5^^8*};print "M333_COVER_CUSPS_READY",[Degree(P):P in cs],Cputime(tt);
VO,om:=SpaceOfHolomorphicDifferentials(L);BB:=[om(v):v in Basis(VO)];assert #BB eq 81;CM:=Matrix(k,[Eltseq(Cartier(b)@@om):b in BB]);C2:=Matrix(k,81,81,[a^2:a in Eltseq(CM)])*CM;assert C2 eq Matrix(k,[Eltseq(Cartier(Cartier(b))@@om):b in BB]);
u:=L!(F!B.1);du:=Differential(u);bb:=[b/du:b in BB];dt:=Differential(tau(u))/du;TM:=Matrix(k,[Eltseq((tau(b)*dt*du)@@om):b in bb]);assert TM^2 eq 1;assert C2*TM eq TM*C2;U:=RowSpace((C2-IdentityMatrix(k,81))*(TM-IdentityMatrix(k,81)));print "M333_CARTIER2_RELATIVE_FORMS_READY",Dimension(U),Rank(TM-1),Cputime(tt);
function EC(a)v:=Eltseq(k!a,GF(2));return Integers()!v[1]+2*Integers()!v[2];end function;
function ER(v)return [[EC(c):c in Coefficients(Numerator(v))],[EC(c):c in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(c):c in Eltseq(F!v)];end function;
function EL(v)return [EF(c):c in Eltseq(L!v)];end function;
function EP(P)a,b:=TwoGenerators(P);assert ideal<Order(Ideal(P))|a,b> eq Ideal(P);return [*Degree(P),IsFinite(P) select 1 else 0,[EL(a),EL(b)]*];end function;
packet:=[*Coefficients(DefiningPolynomial(k)),[ER(c):c in Coefficients(DefiningPolynomial(F))],[EF(c):c in Coefficients(f3)],[EL(tau(a)):a in[L.1,L!F.1,u]],[EL(b):b in bb],EL(j),[EP(P):P in cs],[[EC(a):a in Eltseq(CM)],[EC(a):a in Eltseq(TM)],[EC(a):a in Eltseq(BasisMatrix(U))]]*];Write("work/deg10/continuation_20261031b/M3N33/relative_kummer_cover_model_v3.m",Sprint(packet):Overwrite:=true);print "D31B_M333_KUMMER_COVER_V3_COMPLETE",Dimension(U),Cputime(tt);quit;
