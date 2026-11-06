SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(10*10^9);Alarm(2800);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031a/C46/parent_inputs/fast_hecke.m");load "work/deg10/continuation_20261031a/C46/parent_inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031a/C46/";load "work/deg10/continuation_20261031a/C46/exceptional_cases.m";
k:=GF(3);PP:=PolynomialRing(k);X:=MDX1(46,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 45;coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);print "C46_PARENT_AF_MODEL_READY",Degree(F),Genus(F),Cputime(tt);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)a,b:=TwoGenerators(p);O:=IsFinite(p)select OF else OI;assert ideal<O|a,b>eq Ideal(p);return[*Degree(p),IsFinite(p)select 1 else 0,[EF(a),EF(b)]*];end function;
function Get(E,Q)return C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,Q));end function;
function Vec(a)return[Integers()!v:v in Eltseq(a,k)];end function;
Write(root cat "parent_model46.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords]*]):Overwrite:=true);
alltests:=[* *];
for caseRow in ExceptionalCases do
 typ,idx,parts:=Explode(caseRow);ds:=[];models:=[* *];
 for part in parts do
  d,mul,pol,co:=Explode(part);K<z>:=ext<k|PP!pol>;pt:=Curve(X)(K)![Evaluate(PP!row,z):row in co];ps:=Places(pt);assert #ps eq 1 and Degree(ps[1])eq d;E:=EllipticCurve(X,ps[1]);Q:=LevelStructure(X,ps[1])`P;assert Order(Q)eq 46;Append(~models,[*d,mul,E,Q*]);Append(~ds,d*mul);
 end for;
 assert &+ds eq 10;attempts:=[* *];
 for q in[5,7,11,13,17,19]do
  D:=DivisorGroup(F)!0;rows:=[* *];
  for modrow in models do
   d,mul,E,Q:=Explode(modrow);ps:=Get(E,Q);pd:=Get(E,q*Q);assert Degree(ps)eq d and Degree(pd)eq d;D-:=mul*(Divisor(ps)+q*Divisor(pd));fast:=MDIsogeniesFast(E,Q,q);assert &+[v[3]:v in fast]eq q+1;hec:=[* *];
   for h in fast do ET,QT,e:=Explode(h);assert Order(QT)eq 46;V:=Get(ET,QT);n:=d*e div Degree(V);assert n*Degree(V)eq d*e;D+:=mul*n*Divisor(V);Append(~hec,[*e,n,PD(V)*]);end for;
   L:=BaseRing(E);Append(~rows,[*d,mul,[Integers()!a:a in Coefficients(DefiningPolynomial(L,k))],[Vec(a):a in aInvariants(E)],[Vec(Q[1]),Vec(Q[2])],PD(ps),PD(pd),hec*]);
  end for;
  assert Degree(D)eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");Append(~attempts,[*q,rows,h0*]);print "C46_PARENT_EXCEPTIONAL_ATTEMPT",typ,idx,q,h0,Cputime(tt);if h0 eq 0 then break;end if;
 end for;
 assert attempts[#attempts][3]eq 0;Append(~alltests,[*typ,idx,attempts*]);Write(root cat "parent_exceptional_hecke_progress.m",Sprint(alltests):Overwrite:=true);
end for;
Write(root cat "parent_exceptional_hecke.m",Sprint(alltests):Overwrite:=true);print "C46_PARENT_EXCEPTIONAL_HECKE_COMPLETE",#alltests,Cputime(tt);quit;
