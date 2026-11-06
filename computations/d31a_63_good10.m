SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(1140);
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
e:=10;
assert e in [7,10];
t0:=Cputime();
function E9at(t)
 b:=t^2*(t-1)*(t^2-t+1);c:=t^2*(t-1);K:=Parent(t);
 return EllipticCurve([1-c,-b,-b,K!0,K!0]);
end function;
function NormalizePair(E,R,S)
 assert R[3] eq 1 and S[3] eq 1;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x0:=R[1];y0:=R[2];
 aa3:=2*y0+a3+a1*x0;aa2:=3*x0+a2;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;
 slope:=aa4/aa3;atwo:=aa2-slope^2-a1*slope;scale:=aa3/atwo;
 b:=-atwo^3/aa3^2;c:=-((2*slope+a1)*atwo-aa3)/aa3;
 t:=c^2/(b-c);x:=(S[1]-x0)/scale^2;y:=(S[2]-y0-slope*(S[1]-x0))/scale^3;
 En:=E9at(t);assert aInvariants(En) eq [1-c,-b,-b,Parent(t)!0,Parent(t)!0];
 Q:=En![x,y];assert 7*Q eq En!0 and Q ne En!0;
 return t,x;
end function;
F<z>:=GF(2^e);pairs:=[];
for tt in F do
 if tt in [F!0,F!1] then continue;end if;
 try E:=E9at(tt);catch e;continue;end try;
 if #E mod 63 ne 0 then continue;end if;
 assert #E eq (e eq 7 select 126 else 1008);
 for rr in Roots(DivisionPolynomial(E,7)) do
  assert #Points(E,rr[1]) eq 2;
  Append(~pairs,<tt,rr[1]>);
 end for;
end for;
assert #pairs eq (e eq 7 select 189 else 270);
function PairKey(t,x)
 return Minimum([Eltseq(t^(2^j),GF(2)) cat Eltseq(x^(2^j),GF(2)):j in [0..e-1]]);
end function;
lookup:=AssociativeArray();for r in pairs do lookup[PairKey(r[1],r[2])]:=r;end for;
assert #Keys(lookup) eq 27;
assert &and[#{<r[1]^(2^j),r[2]^(2^j)>:j in [0..e-1]} eq e:r in pairs];
left:=Keys(lookup);reps:=[];
while #left gt 0 do
 key:=Minimum(Setseq(left));Append(~reps,lookup[key]);todo:=[key];orb:={key};
 while #todo gt 0 do
  k:=todo[1];Remove(~todo,1);r:=lookup[k];E:=E9at(r[1]);P:=E![0,0];Q:=Points(E,r[2])[1];
  for a in [2,5] do
   ta,xa:=NormalizePair(E,a*P,a*Q);ka:=PairKey(ta,xa);assert IsDefined(lookup,ka);
   if ka notin orb then Include(~orb,ka);Append(~todo,ka);end if;
  end for;
 end while;
 assert #orb eq 9;left diff:=orb;
end while;
assert #reps eq 3;
printf "Y63_NC_EXHAUSTIVE e=%o geometric_pairs=%o closed_places=%o diamond_orbits=%o time=%o\n",e,#pairs,#Keys(lookup),#reps,Cputime(t0);
for i in [1..#reps] do printf "Y63_NC_REP index=%o fieldpoly=%o t=%o x7=%o\n",i,Eltseq(DefiningPolynomial(F,GF(2))),Eltseq(reps[i][1],GF(2)),Eltseq(reps[i][2],GF(2));end for;
F0<t>:=FunctionField(GF(2));E9:=E9at(t);f:=DivisionPolynomial(E9,7);assert IsIrreducible(f) and Degree(f) eq 24;
AF<u>:=FunctionField(f);assert Genus(AF) eq 49;
O:=MaximalOrderFinite(AF);
printf "Y63_MODEL genus=49 degree=24 maximal_order_ready=true time=%o\n",Cputime(t0);
function PairPlace(t0,x0)
 K:=Parent(t0);ambient_d:=Degree(K,GF(2));mp:=MinimalPolynomial(t0,GF(2));d:=Degree(mp);
 MB:=Matrix(GF(2),[Eltseq(t0^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 xv:=Vector(GF(2),Eltseq(x0,GF(2)));
 assert Rank(VerticalJoin(MB,Matrix(xv))) eq d;
 co:=Solution(MB,xv);
 g:=&+[AF!co[i+1]*(AF!t)^i:i in [0..d-1]];
 I:=ideal<O|O!Evaluate(mp,AF!t),O!(u-g)>;
 assert IsPrime(I);
 P:=Place(I);assert Degree(P) eq d;
 return P;
end function;

function ER(v)return [[Integers()!q:q in Coefficients(Numerator(v))],[Integers()!q:q in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(q):q in Eltseq(AF!v)];end function;
function PD(p)
 g,h:=TwoGenerators(p);OO:=IsFinite(p) select MaximalOrderFinite(AF) else MaximalOrderInfinite(AF);assert ideal<OO|g,h> eq Ideal(p);return [*Degree(p),IsFinite(p) select 1 else 0,[EF(g),EF(h)]*];
end function;
function SD(D)pp,ee:=Support(D);return [*[*PD(pp[i]),ee[i]*]:i in[1..#pp]*];end function;
model:=[* [ER(q):q in Coefficients(f)],EF(AF!jInvariant(E9)),e,[*[*[Integers()!q:q in Eltseq(r[1],GF(2))],[Integers()!q:q in Eltseq(r[2],GF(2))]*]:r in reps*],[Integers()!q:q in Coefficients(DefiningPolynomial(F,GF(2)))]*];
Write("work/deg10/continuation_20261031a/C63/good10_model.m",Sprint(model):Overwrite:=true);
cs:=[];jaf:=AF!jInvariant(E9);
for p0 in Poles(jInvariant(E9)) do for p in Decomposition(AF,p0) do if p notin cs and Valuation(jaf,p) lt 0 then Append(~cs,p);end if;end for;end for;
badcs:=[p:p in cs|Degree(p) eq 3 and Valuation(jaf,p) eq -1];assert #badcs eq 3;
printf "D31A_C63_NC_SETUP e=%o reps3=true bad_cubics3=true time=%o\n",e,Cputime(t0);
answers:=[* *];
for i in[1..#reps] do
 r:=reps[i];E:=E9at(r[1]);P9:=E![0,0];P7:=Points(E,r[2])[1];P:=PairPlace(r[1],r[2]);assert Degree(P) eq e;
 ta,xa:=NormalizePair(E,5*P9,5*P7);Da:=PairPlace(ta,xa);assert Degree(Da) eq e;
 isogs:=MDIsogenies(E,5);assert #isogs eq 6;out:=[];
 for phi in isogs do
  tp,xp:=NormalizePair(Codomain(phi),phi(P9),phi(P7));Append(~out,PairPlace(tp,xp));
 end for;
 uq:=Setseq(Seqset(out));T:=DivisorGroup(AF)!0;targetrows:=[* *];
 for Q in uq do
  mult:=#[z:z in out|z eq Q];assert (e*mult) mod Degree(Q) eq 0;n:=(e*mult) div Degree(Q);T+:=n*Divisor(Q);Append(~targetrows,[*mult,n,PD(Q)*]);
 end for;
 assert Degree(T) eq 6*e;A:=T-Divisor(P)-5*Divisor(Da);assert Degree(A) eq 0;
 if e eq 10 then
  h0:=#Basis(A:Reduction:=false,Simplification:="None");assert h0 eq 0;
  Append(~answers,[*i,PD(P),PD(Da),targetrows,SD(A),h0*]);
  printf "D31A_C63_NC10_NONPRINCIPAL index=%o h0=%o hecke_degree=%o time=%o\n",i,h0,Degree(T),Cputime(t0);
 else
  tests:=[* *];
  for ci,cj in[1..3] do
   if ci eq cj then continue;end if;
   B:=A+4*Divisor(badcs[ci])-4*Divisor(badcs[cj]);h0:=#Basis(B:Reduction:=false,Simplification:="None");assert h0 eq 0;
   Append(~tests,[*ci,cj,PD(badcs[ci]),PD(badcs[cj]),SD(B),h0*]);
   printf "D31A_C63_NC7_CUBIC_NONPRINCIPAL index=%o ci=%o cj=%o h0=%o time=%o\n",i,ci,cj,h0,Cputime(t0);
  end for;
  Append(~answers,[*i,PD(P),PD(Da),targetrows,tests*]);
 end if;
 Write("work/deg10/continuation_20261031a/C63/good10_RR.m",Sprint(answers):Overwrite:=true);
end for;
printf "D31A_C63_NC_COMPLETE e=%o reps=%o h0_all_zero=true time=%o\n",e,#reps,Cputime(t0);quit;
