SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030b/C87/deps/fast_hecke.m");Alarm(600);
N:=87;
assert N eq 87;m:=N eq 55 select 5 else 3;ell:=N div m;h:=(ell-1) div 2;tstart:=Cputime();
CQ:=SmallModularCurve(ell);jraw:=jInvariant(CQ,ell);wraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,ell,ell));
CC:=ChangeRing(CQ,GF(2));assert IsNonsingular(CC) and Genus(CC) eq 2;
RC:=CoordinateRing(Ambient(CC));jn:=RC!Numerator(jraw);jd:=RC!Denominator(jraw);ww:=[RC!e:e in wraw];
jell:=FunctionField(CC)!(jn/jd);jNell:=FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww));
K,mK:=AlgorithmicFunctionField(FunctionField(CC));jK:=mK(jell);
jdegree:=-&+[Valuation(jK,P)*Degree(P):P in Poles(jK)];assert jdegree eq ell+1;
printf "DATABASE_ELL_MODEL N=%o ell=%o genus=%o good_reduction=true j_degree=%o time=%o\n",N,ell,Genus(K),jdegree,Cputime(tstart);
R<v>:=PolynomialRing(K);f:=N eq 55 select v^6+jK*v+jK else v^4+v^3+1/jK;assert IsIrreducible(f);
AF<t3>:=FunctionField(f);gg:=Genus(AF);assert gg eq 9;
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
 QY<Y0>:=PolynomialRing(PL);cp:=Evaluate(ce,[QY!X,Y0,QY!1]);assert Degree(cp) eq 2;c2:=Coefficient(cp,2);c1:=Coefficient(cp,1);
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


function Key0(r)
 return Minimum([&cat[Eltseq(s^(2^i),GF(2)):s in r]:i in [0..Degree(Parent(r[1]))-1]]);
end function;
function SR(v)return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function SK(v)return [SR(c):c in Eltseq(K!v)];end function;
function SA(v)return [SK(c):c in Eltseq(AF!v)];end function;
function PB(P)
 f,g:=TwoGenerators(P);assert ideal<Order(Ideal(P))|f,g> eq Ideal(P);return [*IsFinite(P) select 1 else 0,Degree(P),[SA(f),SA(g)]*];
end function;
function DB(D)
 ps,es:=Support(D);return [*[*PB(ps[i]),es[i]*]:i in [1..#ps]*];
end function;
function EC(v)return [Integers()!c:c in Eltseq(v,GF(2))];end function;
function CurveData(E,P)
 L:=BaseRing(E);return [*Coefficients(DefiningPolynomial(L,GF(2))),[EC(v):v in aInvariants(E)],[EC(v):v in Eltseq(P)]*];
end function;
root:="work/deg10/continuation_20261030b/C87/";
Write(root cat "base_model.m",Sprint([*[SR(v):v in Coefficients(DefiningPolynomial(K))],[SK(v):v in Coefficients(DefiningPolynomial(AF))],[SA(v):v in [IX,IY,IT]],SK(jK)*]):Overwrite:=true);
results:=[* *];rr:=[* *];inputdata:=[* *];census:=[* *];isogdata:=[* *];
for d in [9,10] do
 F<a>:=GF(2^d);lookup:=AssociativeArray();models:=AssociativeArray();counts:=AssociativeArray();nc:=0;nl:=0;
 for j in F do for E in Twists(EllipticCurveFromjInvariant(j)) do
  n:=#E;if n mod 87 ne 0 then continue;end if;G,mp:=AbelianGroup(E);iv:=Invariants(G);if iv[#iv] mod87 ne0 then continue;end if;
  nc+:=1;kyi:=iv;if not IsDefined(counts,kyi) then counts[kyi]:=0;end if;counts[kyi]+:=1;
  R29:=(iv[#iv] div29)*mp(G.#iv);assert Order(R29) eq29;gs:=[(iv[i] div3)*mp(G.i):i in[1..#iv]|iv[i] mod3 eq0];
  tors:=#gs eq1 select [a*gs[1]:a in[0..2]] else [a*gs[1]+b*gs[2]:a,b in[0..2]];
  seen:={};for R3 in tors do if R3 eq E!0 then continue;end if;kk:=EC(R3[1]);if kk in seen then continue;end if;Include(~seen,kk);nl+:=1;
   assert Order(R3) eq3;PN:=R3+R29;assert Order(PN) eq87;r:=Coordinates0(E,R3,R29);ky:=Key0(r);lookup[ky]:=r;models[ky]:=<E,PN>;
  end for;
 end for;end for;
 rows:=[*[*iv,counts[iv]*]:iv in Keys(counts)*];Append(~census,[*d,nc,nl,#Keys(lookup),rows*]);print "B87_BASE_CENSUS",d,nc,nl,#Keys(lookup),rows,Cputime(tstart);
 idx:=0;for ky in Sort(Setseq(Keys(lookup))) do
  idx+:=1;E,PN:=Explode(models[ky]);P:=CoordPlace(lookup[ky]);assert Degree(P) eq d;Append(~inputdata,[*d,idx,ky,CurveData(E,PN),PB(P)*]);
  ins:=MDIsogeniesFast(E,PN,5);assert &+[v[3]:v in ins] eq6;T:=DivisorGroup(AF)!0;idata:=[* *];
  for row in ins do EG,PG,w:=Explode(row);r:=Coordinates0(EG,29*PG,3*PG);Q:=CoordPlace(r);de:=Degree(Q);assert (d*w) mod de eq0;mult:=d*w div de;T+:=mult*Divisor(Q);Append(~idata,[*w,de,mult,CurveData(EG,PG),PB(Q)*]);end for;
  assert Degree(T) eq6*d;D:=T-6*Divisor(P);assert Degree(D) eq0;hh:=#Basis(D:Reduction:=false,Simplification:="None");assert hh in[0,1];Append(~rr,[*d,idx,DB(D),hh*]);Append(~isogdata,[*d,idx,idata*]);Append(~results,[d,idx,hh]);print "B87_BASE_RR",d,idx,hh,Cputime(tstart);
  Write(root cat "base_results.m",Sprint([*census,inputdata,isogdata,rr,results*]):Overwrite:=true);
 end for;
end for;
print "B87_BASE_COMPLETE",#results,#[v:v in results|v[3] eq1],Cputime(tstart);quit;
