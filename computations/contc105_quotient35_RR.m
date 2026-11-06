SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030c/C105/inputs/fast_hecke.m");
load "work/deg10/continuation_20261030c/C105/inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261030c/C105/";k:=GF(2);PP:=PolynomialRing(k);X:=MDX1(35,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 25;coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
print "C105_TARGET35_NATIVE_MODEL_READY",Degree(F),Genus(F),Cputime(tt);
function ER(v)return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(c):c in Eltseq(F!v)];end function;
function PD(p) return [*Degree(p),Order(Ideal(p)) eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*];end function;
function Get(E,P)
 vals:=TateNormalForm_xy(E,P);p:=C48PairAFPlaceLinear(coords,vals);assert Degree(p) eq LCM([Degree(MinimalPolynomial(v,k)):v in vals]);return p;
end function;
data:=eval Read(root cat "source105.m");pol:=[r[2]:r in data[1]|r[1] eq 10][1];K<z>:=ext<k|PP!pol>;tests:=[* *];
for i in [1..#data[3]] do
 d,ky:=Explode(data[3][i]);assert d eq 10;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];assert Order(P) eq 105;Q:=3*P;assert Order(Q) eq 35;ps:=Get(E,Q);pd:=Get(E,3*Q);ns:=d div Degree(ps);nd:=3*d div Degree(pd);assert ns*Degree(ps) eq d and nd*Degree(pd) eq 3*d;
 D:=-ns*Divisor(ps)-nd*Divisor(pd);hec:=[* *];fast:=MDIsogeniesFast(E,Q,3);assert &+[v[3]:v in fast] eq 4;
 for row in fast do ET,QT,e:=Explode(row);assert Order(QT) eq 35;pt:=Get(ET,QT);n:=d*e div Degree(pt);assert n*Degree(pt) eq d*e;D+:=n*Divisor(pt);Append(~hec,[*e,n,PD(pt),[[Integers()!v:v in Eltseq(t,k)]:t in TateNormalForm_xy(ET,QT)]*]);end for;
 assert Degree(D) eq 0;h0:=Dimension(D);Append(~tests,[*i,PD(ps),ns,PD(pd),nd,hec,h0*]);Write(root cat Sprintf("quotient35_progress%o.m",i),Sprint(tests):Overwrite:=true);print "C105_TARGET35_A3_RR_CASE",i,Degree(ps),Degree(pd),h0,Cputime(tt);
end for;
out:=[*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords],tests*];Write(root cat "quotient35_RR.m",Sprint(out):Overwrite:=true);print "C105_TARGET35_A3_RR_COMPLETE",#tests,[v[7]:v in tests],Cputime(tt);quit;
