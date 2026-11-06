SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(1100);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031a/C78/inputs/fast_hecke.m");
load "work/deg10/continuation_20261031a/C78/inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031a/C78/";k:=GF(2);PP:=PolynomialRing(k);X:=MDX1(39,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 33;coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
print "D31A_C78_X139_NATIVE_MODEL_READY",Degree(F),Genus(F),Cputime(tt);
function ER(v)return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(c):c in Eltseq(F!v)];end function;
function PD(p) a,b:=TwoGenerators(p);O:=IsFinite(p) select OF else OI;assert ideal<O|a,b> eq Ideal(p);return [*Degree(p),IsFinite(p) select 1 else 0,[EF(a),EF(b)]*];end function;
function Get(E,P)
 vals:=TateNormalForm_xy(E,P);p:=C48PairAFPlaceLinear(coords,vals);assert Degree(p) eq LCM([Degree(MinimalPolynomial(v,k)):v in vals]);return p;
end function;

data:=eval Read(root cat "source39.m");tests:=[* *];Write(root cat "model39.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords]*]):Overwrite:=true);
for i in[1..#data[3]] do
 d,ky:=Explode(data[3][i]);assert d in[9,10];pol:=[r[2]:r in data[1]|r[1] eq d][1];K<z>:=ext<k|PP!pol>;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);Q:=E![0,0];assert Order(Q) eq 39;
 ps:=Get(E,Q);pd:=Get(E,5*Q);assert Degree(ps) eq d and Degree(pd) eq d;D:=-Divisor(ps)-5*Divisor(pd);hec:=[* *];fast:=MDIsogeniesFast(E,Q,5);assert &+[v[3]:v in fast] eq 6;
 for row in fast do ET,QT,e:=Explode(row);assert Order(QT) eq 39;pt:=Get(ET,QT);n:=d*e div Degree(pt);assert n*Degree(pt) eq d*e;D+:=n*Divisor(pt);Append(~hec,[*e,n,PD(pt),[[Integers()!v:v in Eltseq(t,k)]:t in TateNormalForm_xy(ET,QT)],[Integers()!v:v in Coefficients(DefiningPolynomial(BaseRing(ET),k))]*]);end for;
 assert Degree(D) eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");Append(~tests,[*i,d,PD(ps),PD(pd),hec,h0*]);Write(root cat "A5_progress.m",Sprint(tests):Overwrite:=true);printf "D31A_C78_X139_A5_CASE ci=%o degree=%o h0=%o time=%o\n",i,d,h0,Cputime(tt);
end for;
Write(root cat "A5_RR.m",Sprint(tests):Overwrite:=true);print "D31A_C78_X139_A5_COMPLETE",#tests,[v[6]:v in tests],Cputime(tt);quit;
