SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(700);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");k<w>:=GF(4);X:=MDX1(33,k);FF:=FunctionField(Curve(X));F,mp:=AlgorithmicFunctionField(FF);B:=BaseField(F);E:=EllipticCurve([mp(a):a in X`_E]);Q:=E![mp(a):a in X`_P];assert 33*Q eq E!0;printf "M333_COVER_BASE_READY g=%o degree=%o cpu=%o\n",Genus(F),Degree(F),Cputime(tt);
PZ<Z>:=PolynomialRing(F);f:=DivisionPolynomial(E,3);f3,r0:=Quotrem(f,Z-(11*Q)[1]);assert r0 eq 0 and Degree(f3) eq 3 ;print "M333_CUBIC_POLYNOMIAL_READY",Cputime(tt);L<r>:=FunctionField(f3:Check:=false);print "M333_CUBIC_FIELD_CREATED",Cputime(tt);assert Genus(L) eq 81;print "M333_CUBIC_COVER_READY",Degree(L),Genus(L),Cputime(tt);
XY:=TateNormalForm_xy(E,10*Q);
function Sub(f)
 ff:=f@@mp;nn:=Numerator(ff);dd:=Denominator(ff);args:=Rank(Parent(nn)) eq 3 select XY cat [F!1] else XY;return Evaluate(nn,args)/Evaluate(dd,args);
end function;
up:=Sub(F!B.1);vp:=Sub(F.1);bm:=hom<B->L|L!up>;fm:=hom<F->L|bm,L!vp>;
A:=10*Q;a1,a2,a3,a4,a6:=Explode(aInvariants(E));xr:=A[1];yr:=A[2];ss:=(3*xr^2+2*a2*xr+a4-a1*yr)/(2*yr+a1*xr+a3);a2p:=a2+3*xr-a1*ss-ss^2;a3p:=a3+a1*xr+2*yr;sc:=a2p/a3p;rp:=(L!sc)^2*(r-L!xr);tau:=hom<L->L|fm,rp>;
for a in [L.1,L!F.1,L!(F!B.1)] do assert tau(tau(a)) eq a;end for;assert tau(L.1) ne L.1 or tau(L!F.1) ne L!F.1;print "M333_DELTA10_INVOLUTION_READY",Cputime(tt);
j:=L!jInvariant(E);cs:=Poles(j);assert {*Degree(P):P in cs*} eq {*1^^40,5^^8*};print "M333_COVER_CUSPS_READY",[Degree(P):P in cs],Cputime(tt);
VO,om:=SpaceOfHolomorphicDifferentials(L);BB:=[om(v):v in Basis(VO)];assert #BB eq 81;CM:=Matrix(k,[Eltseq(Cartier(b)@@om):b in BB]);C2:=Matrix(k,81,81,[a^2:a in Eltseq(CM)])*CM;assert C2 eq Matrix(k,[Eltseq(Cartier(Cartier(b))@@om):b in BB]);
u:=L!(F!B.1);du:=Differential(u);bb:=[b/du:b in BB];dt:=Differential(tau(u))/du;TM:=Matrix(k,[Eltseq((tau(b)*dt*du)@@om):b in bb]);assert TM^2 eq 1;assert C2*TM eq TM*C2;U:=RowSpace((C2-IdentityMatrix(k,81))*(TM-IdentityMatrix(k,81)));print "M333_CARTIER2_RELATIVE_FORMS_READY",Dimension(U),Rank(TM-1),Cputime(tt);
function EC(a)v:=Eltseq(k!a,GF(2));return Integers()!v[1]+2*Integers()!v[2];end function;
function ER(v)return [[EC(c):c in Coefficients(Numerator(v))],[EC(c):c in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(c):c in Eltseq(F!v)];end function;
function EL(v)return [EF(c):c in Eltseq(L!v)];end function;
function EP(P)a,b:=TwoGenerators(P);assert ideal<Order(Ideal(P))|a,b> eq Ideal(P);return [*Degree(P),IsFinite(P) select 1 else 0,[EL(a),EL(b)]*];end function;
packet:=[*Coefficients(DefiningPolynomial(k)),[ER(c):c in Coefficients(DefiningPolynomial(F))],[EF(c):c in Coefficients(f3)],[EL(tau(a)):a in[L.1,L!F.1,u]],[EL(b):b in bb],EL(j),[EP(P):P in cs],[[EC(a):a in Eltseq(CM)],[EC(a):a in Eltseq(TM)],[EC(a):a in Eltseq(BasisMatrix(U))]]*];Write("work/deg10/continuation_20261031b/M3N33/relative_cover_model.m",Sprint(packet):Overwrite:=true);print "D31B_M333_COVER_PREFLIGHT_COMPLETE",Dimension(U),Cputime(tt);quit;
