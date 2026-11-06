SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
N:=69;
assert N eq 69;m:=N eq 55 select 5 else 3;ell:=N div m;h:=(ell-1) div 2;tstart:=Cputime();
CQ:=SmallModularCurve(ell);jraw:=jInvariant(CQ,ell);wraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,ell,ell));
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq (ell eq 23 select 2 else 1);
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);ww:=[RC!e:e in wraw];
jell:=FunctionField(CC)!(jn/jd);jNell:=FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww));
K,mK:=AlgorithmicFunctionField(FunctionField(CC));jK:=mK(jell);
jdegree:=-&+[Valuation(jK,P)*Degree(P):P in Poles(jK)];assert jdegree eq ell+1;
printf "DATABASE_ELL_MODEL N=%o ell=%o genus=%o good_reduction=true j_degree=%o time=%o\n",N,ell,Genus(K),jdegree,Cputime(tstart);
R<v>:=PolynomialRing(K);f:=N eq 55 select v^6+jK*v+jK else v^4+v^3+1/jK;assert IsIrreducible(f);
AF<t3>:=FunctionField(f);gg:=Genus(AF);assert gg eq (N eq 69 select 7 else 5);
O:=MaximalOrderFinite(AF);OK:=MaximalOrderFinite(K);
numY,denY:=IntegralSplit(K.1,OK);numJ,denJ:=IntegralSplit(jK,OK);
IX:=AF!BaseField(K).1;IY:=AF!K!numY;multT:=N eq 55 select K!denJ else K!numJ;IT:=(AF!multT)*t3;
integral_coordinates:=[O!IX,O!IY,O!IT];
fx:=(K!BaseField(K).1)@@mK;fy:=(K!numY)@@mK;fJ:=multT@@mK;
printf "X0_3ELL_MODEL N=%o genus=%o degree_over_X0ell=4 integral_coordinates_ready=true time=%o\n",N,gg,Cputime(tstart);
function ThreeParameter(E,R3)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x:=R3[1];y:=R3[2];
 aa3:=2*y+a3+a1*x;aa2:=3*x+a2;aa4:=3*x^2+2*x*a2+a4-a1*y;
 s:=aa4/aa3;aone:=a1+2*s;assert aa2-s^2-a1*s eq 0;
 tt:=aa3/aone^3;assert 1/(tt^3*(tt+1)) eq jInvariant(E);return tt;
end function;
function SmallParameter(E,R)
 if N eq 55 then
  b,c:=Explode(TateNormalForm_bc(E,R));assert b eq c;u:=b+1/b;assert u^6/(u+1) eq jInvariant(E);return u;
 else return ThreeParameter(E,R);end if;
end function;
function EllPoint(E,Rl)
 L:=BaseRing(E);PL<X>:=PolynomialRing(L);polyker:=&*[X-(k*Rl)[1]:k in [1..(ell-1) div 2]];
 Et,phi:=IsogenyFromKernel(E,polyker);jj:=jInvariant(E);jjN:=jInvariant(Et);
 A:=Evaluate(jn,[X,PL!0,PL!1]);BB:=Evaluate(jn,[X,PL!1,PL!1])-A;DD:=Evaluate(jd,[X,PL!0,PL!1]);
 assert Evaluate(jd,[X,PL!1,PL!1]) eq DD;
 ce:=DefiningEquation(CC);c0:=Evaluate(ce,[X,PL!0,PL!1]);c2:=PL!1;
 if ell eq 19 then c1:=PL!1;else c1:=X^3+X+1;end if;
 assert Evaluate(ce,[X,PL!1,PL!1]) eq c0+c1+c2;
 YY:=jj*DD-A;res:=c2*YY^2+c1*YY*BB+c0*BB^2;
 out:=[];
 for rt in Roots(res) do
  xx:=rt[1];bb:=Evaluate(BB,xx);
  if bb ne 0 then ys:=[Evaluate(YY,xx)/bb];else ys:=[zz[1]:zz in Roots(c2*X^2+Evaluate(c1,xx)*X+Evaluate(c0,xx))];end if;
  for yy in ys do
   pt:=CC(L)![xx,yy,L!1];
   if Evaluate(jd,Eltseq(pt)) ne 0 and Evaluate(jell,pt) eq jj and Evaluate(jNell,pt) eq jjN then Append(~out,pt);end if;
  end for;
 end for;
 assert #out eq 1;return out[1];
end function;
function Coordinates0(E,R3,Rl)
 pt:=EllPoint(E,Rl);return <Evaluate(fx,pt),Evaluate(fy,pt),Evaluate(fJ,pt)*SmallParameter(E,R3)>;
end function;
function Key0(r)
 return Minimum([&cat[Eltseq(s^(2^i),GF(2)):s in r]:i in [0..Degree(Parent(r[1]))-1]]);
end function;
function CoordPlace(r)
 L:=Parent(r[1]);degrees:=[Degree(MinimalPolynomial(s,GF(2))):s in r];d:=LCM(degrees);
 candidates:=[<IX,r[1]>,<IY,r[2]>,<IT,r[3]>,<IX+IY,r[1]+r[2]>,<IX+IT,r[1]+r[3]>,<IY+IT,r[2]+r[3]>,<IX+IY+IT,r[1]+r[2]+r[3]>];
 found:=false;
 for q in candidates do
  mp:=MinimalPolynomial(q[2],GF(2));if Degree(mp) eq d then gen:=q[1];val:=q[2];found:=true;break;end if;
 end for;
 assert found;
 MB:=Matrix(GF(2),[Eltseq(val^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 gens:=[O!Evaluate(mp,gen)];
 for i in [1..3] do
  xv:=Vector(GF(2),Eltseq(r[i],GF(2)));assert Rank(VerticalJoin(MB,Matrix(xv))) eq d;
  co:=Solution(MB,xv);poly:=&+[AF!co[k+1]*gen^k:k in [0..d-1]];
  Append(~gens,O!([IX,IY,IT][i]-poly));
 end for;
 II:=ideal<O|gens>;assert IsPrime(II);P:=Place(II);assert Degree(P) eq d;return P;
end function;

function SR(a)
 return [[Integers()!v:v in Coefficients(Numerator(a))],[Integers()!v:v in Coefficients(Denominator(a))]];
end function;
function SK(a) return [SR(v):v in Eltseq(K!a)];end function;
function SA(a) return [SK(v):v in Eltseq(AF!a)];end function;

HP:=7;EPS:=KroneckerSymbol(69,HP);assert EPS eq -1;
model:=eval Read("work/deg10/next_cases_20261030/C69/cover_model.m");
assert [SR(v):v in Coefficients(DefiningPolynomial(K))] eq model[1];assert [SK(v):v in Coefficients(DefiningPolynomial(AF))] eq model[2];
PR:=PolynomialRing(GF(2));
function DR(v) return BaseField(K)!(PR!v[1])/BaseField(K)!(PR!v[2]);end function;
function DK(v) return K![DR(a):a in v];end function;
function DA(v) return AF![DK(a):a in v];end function;
CY<z>:=FunctionField(PolynomialRing(AF)![DA(v):v in model[3]]);assert Genus(CY) eq 13;
function SY(v) return [SA(a):a in Eltseq(CY!v)];end function;
function PY(P) return [*IsFinite(P) select 1 else 0,Degree(P),[SY(v):v in Basis(Ideal(P))]*];end function;
function DY(D) ps,es:=Support(D);return [*[*PY(ps[i]),es[i]*]:i in [1..#ps]*];end function;
function RelativeDecomposition(F,P)
 assert BaseField(F) cmpeq FunctionField(P);
 f,g:=TwoGenerators(P);u:=UniformizingElement(P);
 assert Valuation(f,P) gt 0 and Valuation(g,P) gt 0 and Valuation(u,P) eq 1;
 raw:=Decomposition(F,P);
 out:=[Q:Q in raw|Valuation(F!f,Q) gt 0 and Valuation(F!g,Q) gt 0];
 assert #out gt 0 and &+[Degree(Q):Q in out] eq 2*Degree(P);
 assert &and[Valuation(F!u,Q) eq 1:Q in out];
 return out;
end function;

D:=10;FD:=GF(2^D);lookupD:=AssociativeArray();modelsD:=AssociativeArray();countD:=0;
for j in FD do for E in Twists(EllipticCurveFromjInvariant(j)) do
 if #E mod 69 ne 0 then continue;end if;countD+:=1;G,mp:=AbelianGroup(E);assert Invariants(G) eq [966];P69:=14*mp(G.1);assert Order(P69) eq 69;P3:=23*P69;P23:=3*P69;r:=Coordinates0(E,P3,P23);ky:=Key0(r);lookupD[ky]:=r;modelsD[ky]:=<E,P3,P23>;
end for;end for;assert countD eq 20 and #Keys(lookupD) eq 2;
results:=[* *];rr:=[* *];idx:=0;
for ky in Sort(Setseq(Keys(lookupD))) do
 idx+:=1;E,P3,P23:=Explode(modelsD[ky]);P:=CoordPlace(lookupD[ky]);assert Degree(P) eq 10;out:=[];
 for ph in MDIsogenies(E,HP) do Append(~out,CoordPlace(Coordinates0(Codomain(ph),ph(P3),ph(P23))));end for;assert #out eq HP+1;
 T:=DivisorGroup(AF)!0;for Q in Seqset(out) do n:=#[a:a in out|a eq Q];assert(D*n)mod Degree(Q) eq 0;T+:=((D*n)div Degree(Q))*Divisor(Q);end for;assert Degree(T) eq (HP+1)*D;
 basepr:=IsPrincipal(T-(HP+1)*Divisor(P));print "CONT69_H7_BASE",idx,basepr,Cputime(tstart);if not basepr then Append(~results,[*idx,ky,false,0,0*]);continue;end if;
 dp:=RelativeDecomposition(CY,P);assert #dp eq 2 and &and[Degree(a) eq 10:a in dp];neg:=Divisor(dp[1])+HP*Divisor(EPS eq 1 select dp[1] else dp[2]);choices:=[DivisorGroup(CY)!0];
 for Q in Seqset(out) do
  n:=(D*#[a:a in out|a eq Q])div Degree(Q);dd:=RelativeDecomposition(CY,Q);
  if #dd eq 1 then options:=n mod 2 eq 0 select [(n div 2)*Divisor(dd[1])] else [];
  else assert #dd eq 2;options:=[a*Divisor(dd[1])+(n-a)*Divisor(dd[2]):a in [0..n]];end if;
  choices:=[A+B:A in choices,B in options];
 end for;
 bad:=0;for i in [1..#choices] do A:=choices[i]-neg;assert Degree(A) eq 0;pr:=IsPrincipal(A);if pr then bad+:=1;end if;Append(~rr,[*idx,i,DY(A),pr*]);print "CONT69_H7_COVER",idx,i,#choices,pr,Cputime(tstart);end for;
 Append(~results,[*idx,ky,true,#choices,bad*]);Write("work/deg10/next_cases_20261030/C69/hecke7_results.m",Sprint([*results,rr*]):Overwrite:=true);
end for;
Write("work/deg10/next_cases_20261030/C69/hecke7_results.m",Sprint([*results,rr*]):Overwrite:=true);print "CONT69_H7_COMPLETE",results,Cputime(tstart);quit;
