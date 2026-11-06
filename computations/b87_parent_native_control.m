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

function IntrinsicKey(r)
 L:=Parent(r[1]);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in r]);found:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]] do
  g:=&+[L!bits[i]*r[i]:i in[1..3]];mp:=MinimalPolynomial(g,GF(2));if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(GF(2),[Eltseq(g^i,GF(2)):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(GF(2),Eltseq(v,GF(2))))):v in r];assert &and[&+[L!co[j][i+1]*g^i:i in[0..d-1]] eq r[j]:j in[1..3]];
 return [d] cat bb cat[Integers()!v:v in Coefficients(mp)] cat[Integers()!v:v in &cat co];
end function;
root:="work/deg10/continuation_20261030b/C87/";data:=eval Read(root cat "base_results.m");raw:=eval Read(root cat "raw_hecke_isogenies.m");m:=eval Read(root cat "base_model.m");PR:=PolynomialRing(GF(2));
function DR(v)return BaseField(K)!(PR!v[1])/BaseField(K)!(PR!v[2]);end function;
function DK(v)return K![DR(t):t in v];end function;
function DA(v)return AF![DK(t):t in v];end function;
assert Coefficients(DefiningPolynomial(K)) eq[DR(v):v in m[1]] and Coefficients(DefiningPolynomial(AF)) eq[DK(v):v in m[2]];assert [IX,IY,IT] eq[DA(v):v in m[3]];
OI:=MaximalOrderInfinite(AF);
function ParsePrime(row)
 OO:=row[1] eq 1 select O else OI;I:=ideal<OO|[DA(v):v in row[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq row[2];return P;
end function;
function DecodeCurve(cd)
 F<a>:=ext<GF(2)|PR!cd[1]>;fv:=func<v|&+[F!v[i]*a^(i-1):i in[1..#v]]>;E:=EllipticCurve([fv(v):v in cd[2]]);P:=E![fv(v):v in cd[3]];return E,P;
end function;
function DecodeDiv(rows)
 D:=DivisorGroup(AF)!0;for row in rows do D+:=row[2]*Divisor(ParsePrime(row[1]));end for;return D;
end function;
out:=[* *];nk:=0;
for ii in[1..#data[2]] do
 rec:=data[2][ii];d,idx,ky,cd,pp:=Explode(rec);E,PN:=DecodeCurve(cd);assert Order(PN) eq 87;P:=ParsePrime(pp);assert CoordPlace(Coordinates0(E,29*PN,3*PN)) eq P;nat:=AssociativeArray();fast:=AssociativeArray();rrow:=[r:r in raw|r[1] eq d and r[2] eq idx][1];assert #rrow[4] eq 6;
 for pack in rrow[4] do
  fd:=pack[1];L<z>:=ext<GF(2)|PR!fd[1]>;vals:=[Evaluate(PR!v,z):v in fd[2]];r:=vals[22..24];kk:=IntrinsicKey(r);assert kk[1] eq pack[2];if not IsDefined(nat,kk) then nat[kk]:=0;end if;nat[kk]+:=1;nk+:=1;
 end for;
 T:=DivisorGroup(AF)!0;rr:=[r:r in data[3]|r[1] eq d and r[2] eq idx][1];assert &+[r[1]:r in rr[3]] eq 6;
 for im in rr[3] do w,de,mult,ec,qp:=Explode(im);EG,PG:=DecodeCurve(ec);assert Order(PG) eq 87;r:=Coordinates0(EG,29*PG,3*PG);kk:=IntrinsicKey(r);assert kk[1] eq de and d*w eq de*mult;Q:=ParsePrime(qp);assert CoordPlace(r) eq Q and Degree(Q) eq de;T+:=mult*Divisor(Q);if not IsDefined(fast,kk) then fast[kk]:=0;end if;fast[kk]+:=w;end for;
 assert Keys(nat) eq Keys(fast);for ky in Keys(nat) do assert nat[ky] eq fast[ky];end for;
 row:=[r:r in data[4]|r[1] eq d and r[2] eq idx][1];D:=DecodeDiv(row[3]);assert D eq T-6*Divisor(P) and Degree(D) eq 0;hh:=#Basis(D:Reduction:=false,Simplification:="None");assert hh eq 0 and row[4] eq 0;Append(~out,[d,idx,#Keys(nat),hh]);print "B87_PARENT_NATIVE_NORM_RR",d,idx,#Keys(nat),hh,Cputime(tstart);
end for;assert #out eq 12 and nk eq 72;Write(root cat "parent_native_norm_RR.m",Sprint(out):Overwrite:=true);print "B87_PARENT_NATIVE_COMPLETE",12,72,Cputime(tstart);quit;
