SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(550);tt:=Cputime();root:="work/deg10/continuation_20261031b/C66/";AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/continuation_20261031b/C66/inputs/pair_af_place_linear.m";
X:=MDX1(33,GF(2));C:=Curve(X);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);assert Genus(AF)eq21;cs:=Cusps(X);os:=CuspOrbitsQ(X);assert Sort([Degree(o):o in os])eq[1^^10,2^^5,10,10];c0:=Support([o:o in os|Degree(o)eq1][1])[1];A:=FunctionFieldDivisor(Divisor(c0));
G,_,phi:=ClassGroup(AF:ReductionDivisor:=A,Proof:=true);assert Invariants(G)eq[15,42373650,0];base:=phi(A);qc:=[phi(FunctionFieldDivisor(o))-Degree(o)*base:o in os];H:=sub<G|qc>;assert Invariants(H)eq[5,42373650];Q,qmap:=quo<G|sub<G|qc cat[base]>>;assert #Q eq3;coords:=[am(FF.i):i in[1,2]];
printf "D31B_C66_F2_PICARD_READY G=%o H=%o quotient3=true time=%o\n",Invariants(G),Invariants(H),Cputime(tt);
function ER(f)return [[Integers()!q:q in Coefficients(Numerator(f))],[Integers()!q:q in Coefficients(Denominator(f))]];end function;
function EF(f)return [ER(q):q in Eltseq(AF!f)];end function;
function PD(p)g,h:=TwoGenerators(p);OO:=IsFinite(p)select MaximalOrderFinite(AF)else MaximalOrderInfinite(AF);assert ideal<OO|g,h>eq Ideal(p);return[*Degree(p),IsFinite(p)select 1 else 0,[EF(g),EF(h)]*];end function;
function SD(D)ps,es:=Support(D);return [*[*PD(ps[i]),es[i]*]:i in[1..#ps]*];end function;
function Cl(D)return [Integers()!q:q in Eltseq(phi(D)-Degree(D)*base)];end function;
css:=[* *];for p in cs do D:=FunctionFieldDivisor(Divisor(p));Append(~css,[*Degree(p),SD(D),Cl(D),Eltseq(qmap(phi(D)-Degree(D)*base))*]);end for;
qcss:=[*[*Degree(o),SD(FunctionFieldDivisor(o)),Cl(FunctionFieldDivisor(o)),CuspSignature(X,Support(o)[1])*]:o in os*];
frame:=[*[ER(q):q in Coefficients(DefiningPolynomial(AF))],[EF(q):q in coords],SD(A),[Eltseq(v):v in Rows(RelationMatrix(G))],[Eltseq(g):g in qc],[[Integers()!v:v in Eltseq(qmap(G.i))]:i in[1..Ngens(G)]],css,qcss*];Write(root cat "F2_proven_frame.m",Sprint(frame):Overwrite:=true);
nc:=[* *];PF:=PolynomialRing(GF(2));counts:=[];
for d in[1..10] do
 ps:=NoncuspidalPlaces(X,d);Append(~counts,#ps);rows:=[* *];
 for p in ps do vals:=[Evaluate(FF.i,p):i in[1,2]];q:=C48PairAFPlaceLinear(coords,vals);assert Degree(q)eq d;D:=Divisor(q);L:=Parent(vals[1]);pol:=Degree(L)eq 1 select[0,1]else[Integers()!v:v in Coefficients(DefiningPolynomial(L,GF(2)))];xy:=[[Integers()!v:v in Eltseq(w,GF(2))]:w in vals];cl:=phi(D)-d*base;Append(~rows,[*PD(q),[Integers()!v:v in Eltseq(cl)],[Integers()!v:v in Eltseq(qmap(cl))],pol,xy*]);end for;
 Append(~nc,[*d,rows*]);Write(root cat "F2_NC_classes_progress.m",Sprint(nc):Overwrite:=true);printf "D31B_C66_F2_NC_CLASS_DEGREE d=%o count=%o time=%o\n",d,#ps,Cputime(tt);
end for;
assert counts[1..9]eq[0,0,0,0,2,20,30,25,50];Write(root cat "F2_NC_classes.m",Sprint(nc):Overwrite:=true);printf "D31B_C66_F2_FRAME_COMPLETE counts=%o all_TwoGenerators=true time=%o\n",counts,Cputime(tt);quit;
