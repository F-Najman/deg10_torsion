SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(7*10^9);tt:=Cputime();
SetUseMontes(FldFun,true);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX11(2,26,Rationals());C:=Curve(X);FF:=FunctionField(C);AF,mp:=AlgorithmicFunctionField(FF);F0:=RationalExtensionRepresentation(AF);B:=BaseField(F0);R:=PolynomialRing(B);f:=R!DefiningPolynomial(F0);f/:=LeadingCoefficient(f);n:=Degree(f);assert n eq 12;
print "Q226_FUNCTION_FIELD_READY",n,Cputime(tt);
fac:=Factorization(LCM([Denominator(a):a in Coefficients(f)]));D:=B!1;
for pe in fac do
 pp:=pe[1];ee:=Maximum([Ceiling(Valuation(Denominator(Coefficient(f,i)),pp)/(n-i)):i in [0..n-1]|Coefficient(f,i) ne 0]);if ee gt 0 then D*:=pp^ee;end if;
end for;
g:=R![Coefficient(f,i)*D^(n-i):i in [0..n]];assert LeadingCoefficient(g) eq 1 and &and[Denominator(v) eq 1:v in Coefficients(g)];F<z>:=FunctionField(g);hm:=hom<F0->F|z/D>;
b:=hm(F0!mp(FF!X`_coordinates[1]));c:=hm(F0!mp(FF!X`_coordinates[2]));
print "Q226_MONIC_INTEGRAL_MODEL",D,Cputime(tt);
OF:=MaximalOrderFinite(F:Al:="Montes");print "Q226_FINITE_MAXIMAL_ORDER_READY",Cputime(tt);
OI:=MaximalOrderInfinite(F);print "Q226_INFINITE_MAXIMAL_ORDER_READY",Cputime(tt);print "Q226_KNOWN_GENUS_25_FROM_MODULAR_MODEL",Cputime(tt);
bb:=Divisor(b);pp,ee:=Support(bb);cs:=[pp[i]:i in [1..#pp]|Degree(pp[i]) eq 1];labels:=[ee[i]:i in [1..#pp]|Degree(pp[i]) eq 1];assert #cs eq 18 and #Seqset(labels) eq 18;print "Q226_RATIONAL_CUSPS_READY",labels,Cputime(tt);
a:=eval Read("work/deg10/continuation_20261029b/C226/small_cusp_divisors_v2.m");assert Seqset(labels) eq {v:v in a[1]};out:=[* *];
function EQ(v) return [Numerator(v),Denominator(v)];end function;
function ER(v) return [[EQ(x):x in Coefficients(Numerator(v))],[EQ(x):x in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(x):x in Eltseq(F!v)];end function;
model:=[*[ER(x):x in Coefficients(g)],EF(b),EF(c),[EF(v):v in Basis(OF)],[EF(v):v in Basis(OI)]*];
Write("work/deg10/C226_resolution_20261030/Q_model_nogenus.m",Sprint(model));
for row in a[2] do
 L:=&+[row[3][i]*Divisor(cs[Index(labels,a[1][i])]):i in [1..18]];assert Degree(L) eq 10;print "Q226_RR_START",row[1],Cputime(tt);
 bs:=Basis(L:Reduction:=false,Simplification:="None");print "Q226_RR_DIMENSION",row[1],#bs,Cputime(tt);Append(~out,[*row[1],row[2],#bs*]);
 Write("work/deg10/C226_resolution_20261030/Q_RR_nogenus.m",Sprint(out):Overwrite:=true);
 if #bs eq 1 then
  h:=bs[1];EE:=L+Divisor(h);assert IsEffective(EE) and Degree(EE) eq 10;pls,mults:=Support(EE);print "Q226_EFFECTIVE_SUPPORT",row[1],[Degree(p):p in pls],mults,Cputime(tt);
  Write(Sprintf("work/deg10/C226_resolution_20261030/Q_function_nogenus_%o.m",row[1]),Sprint(EF(h)));
  records:=[* *];
  for i in [1..#pls] do
   p:=pls[i];K:=ResidueClassField(p);bk:=Evaluate(b,p);ck:=Evaluate(c,p);E:=EllipticCurve([K|0,ck,0,bk*(1-bk-ck),0]);P:=E![0,0];Q:=E![bk,bk];assert Discriminant(E) ne 0 and Order(P) eq 2 and Order(Q) eq 26 and P ne 13*Q;
   print "Q226_NONCUSPIDAL_POINT",row[1],i,Degree(p),mults[i],Cputime(tt);
   if Type(K) eq FldRat then kp:=[Rationals()|0,1];else kp:=Coefficients(DefiningPolynomial(K));end if;
   Append(~records,[*[EQ(v):v in kp],[EQ(v):v in Eltseq(K!bk)],[EQ(v):v in Eltseq(K!ck)],mults[i]*]);
  end for;
  Write(Sprintf("work/deg10/C226_resolution_20261030/Q_points_nogenus_%o.m",row[1]),Sprint(records));
 end if;
end for;
print "Q226_CHARACTERISTIC_ZERO_RR_COMPLETE",Cputime(tt);quit;
