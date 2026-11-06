SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(9*10^9);Alarm(850);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261031a/C54/";X:=MDX1(54,GF(5));C:=Curve(X);assert Genus(C)eq 52;F,mp:=AlgorithmicFunctionField(FunctionField(C));assert Genus(F)eq 52;print "D31A_C54_CUSP_MODEL_READY",Degree(F),Genus(F),Cputime(tt);
cs:=Cusps(X);ps:=[Support(DivisorGroup(F)!FunctionFieldDivisor(Divisor(P)))[1]:P in cs];assert {*Degree(P):P in ps*}eq {*1^^18,2^^6,6^^2,9^^2*};
V,om:=SpaceOfHolomorphicDifferentials(F);assert Dimension(V)eq 52;ws:=[om(v):v in Basis(V)];data:=[];
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function PD(p)return [*Degree(p),Order(Ideal(p))eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*];end function;
Write(root cat "cusp_completion_model.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[PD(P):P in ps]*]):Overwrite:=true);
for i in[1..#ps]do
 p:=ps[i];d:=Degree(p);lim:=10 div d;u:=UniformizingElement(p);du:=Differential(u);L,ev:=Completion(F,p:Precision:=160);s:=[ev(w/du):w in ws];assert &and[Valuation(w)ge 0 and AbsolutePrecision(w)ge lim:w in s];K:=BaseRing(L);assert #K eq 5^d;
 jets:=[[[Integers()!a:a in Eltseq(K!Coefficient(w,n))]:w in s]:n in[0..lim-1]];
 if d eq 1 then def:=[0,1];else def:=[Integers()!a:a in Coefficients(DefiningPolynomial(K))];end if;
 Append(~data,[*i,d,def,jets,EF(u)*]);Write(root cat "cusp_completion_progress.m",Sprint(data):Overwrite:=true);print "D31A_C54_CUSP_COMPLETION_PLACE",i,d,lim,Cputime(tt);
end for;
Write(root cat "cusp_completion_flags.m",Sprint(data):Overwrite:=true);print "D31A_C54_CUSP_COMPLETION_FLAGS_COMPLETE",#ps,Cputime(tt);quit;
