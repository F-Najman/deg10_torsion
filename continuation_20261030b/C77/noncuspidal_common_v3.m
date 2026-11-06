SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030b/C77/deps/fast_hecke.m");
load "work/deg10/continuation_20261030b/C77/geometry_fast_nc.m";
assert KroneckerSymbol(77,2) eq -1 and KroneckerSymbol(77,3) eq -1;
HH:={Integers()|36^i mod 77:i in [0..4]};assert 6 in HH or 71 in HH;
bb,cc:=Explode(TateNormalForm_bc(E7u,2*P7u));td2:=bb/cc;
if td2 eq 1/(t+1) then de:=1;else assert td2 eq (t+1)/t;de:=2;end if;
deck2:=hom<CY->CY|map<E3->CY|a:->CY!a>,CY.1+1>;
deck3:=hom<CY->CY|map<E3->CY|a:->CY!tau(a)>,CY.1>;
function Sig(a,n)
 b:=a;for i in [1..n mod 6] do for j in [1..de] do b:=deck3(b);end for;b:=deck2(b);end for;return b;
end function;
for a in [CY.1,CY!E3.1,CY!FB.1] do assert Sig(a,6) eq a;end for;
print "C77_DIAMONDS_VERIFIED delta2_cubic_exponent",de,"delta3=delta2inverse";
wwraw:=DefiningPolynomials(AtkinLehnerInvolution(CQ,11,11));ww:=[RC!a:a in wwraw];jdual:=kmap(FunctionField(CC)!(Evaluate(jn,ww)/Evaluate(jd,ww)));jdp:=Eltseq(jdual);assert #jdp eq 2;
function EV(q,xx)
 L:=Parent(xx);P:=PolynomialRing(L);return Evaluate(P!Numerator(q),xx)/Evaluate(P!Denominator(q),xx);
end function;
function X11Coordinates(E,P11)
 assert Order(P11) eq 11;L:=BaseRing(E);PL<X>:=PolynomialRing(L);ker:=&*[X-(i*P11)[1]:i in [1..5]];EE,ph:=IsogenyFromKernel(E,ker);jj:=jInvariant(E);jj11:=jInvariant(EE);
 AX<X>:=FunctionField(L);rx:=hom<BaseField(K)->AX|X>;j0:=rx(jparts[1]);j1:=rx(jparts[2]);a:=rx(Coefficient(kpoly,1));b:=rx(Coefficient(kpoly,0));poly:=Numerator((jj+j0)^2+a*j1*(jj+j0)+b*j1^2);out:=[];
 for z in Roots(poly) do
  xx:=z[1];if &or[Evaluate(PL!Denominator(q),xx) eq 0:q in jparts cat jdp cat Coefficients(kpoly)] then continue;end if;
  jf:=EV(jparts[2],xx);if jf ne 0 then ys:=[(jj+EV(jparts[1],xx))/jf];else ys:=[r[1]:r in Roots(PL![EV(q,xx):q in Coefficients(kpoly)])];end if;
  for yy in ys do if EV(jparts[1],xx)+EV(jparts[2],xx)*yy eq jj and EV(jdp[1],xx)+EV(jdp[2],xx)*yy eq jj11 then Append(~out,<xx,yy>);end if;end for;
 end for;assert #out eq 1;return out[1][1],out[1][2],jj11;
end function;
function Coordinates(E,PN)
 assert Order(PN) eq 77;P7:=11*PN;P11:=7*PN;b,c:=Explode(TateNormalForm_bc(E,P7));tv:=b/c;assert b eq tv^3-tv^2 and c eq tv^2-tv;
 xx,yy,jj11:=X11Coordinates(E,P11);assert EV(jparts[1],xx)+EV(jparts[2],xx)*yy eq jInvariant(E);
 return [tv,xx],jj11;
end function;
function PairKey(vs)
 L:=Parent(vs[1]);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in vs]);ok:=false;
 for bs in [[1,0],[0,1],[1,1]] do v:=bs[1]*vs[1]+bs[2]*vs[2];mp:=MinimalPolynomial(v,GF(2));if Degree(mp) eq d then bits:=bs;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(GF(2),[Eltseq(v^i,GF(2)):i in [0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(GF(2),Eltseq(a,GF(2))))):a in vs];
 return [d] cat bits cat [Integers()!a:a in Coefficients(mp)] cat [Integers()!a:a in &cat co];
end function;
ycoords:=[CY!FB.1,CY!E3.1];
function PlacesAt(vs)
 key:=PairKey(vs);d:=key[1];mp:=MinimalPolynomial(vs[1],GF(2));raw:=Zeros(Evaluate(mp,ycoords[1]));out:=[];
 for q in raw do if Degree(q) ne d or &or[Valuation(v,q) lt 0:v in ycoords] then continue;end if;if PairKey([Evaluate(v,q):v in ycoords]) eq key then Append(~out,q);end if;end for;
 assert #out eq 2 and &and[Degree(q) eq d:q in out];return out;
end function;
known:=[];images:=[];
function Image(q,n)
 for i in [1..#known] do if known[i] eq <q,n mod 6> then return images[i];end if;end for;
 target:=PairKey([Evaluate(Sig(v,n),q):v in ycoords]);a,b:=TwoGenerators(q);aa:=Sig(a,-n);bb:=Sig(b,-n);raw:=Zeros(aa);out:=[];
 for z in raw do if Degree(z) ne Degree(q) or Valuation(bb,z) le 0 or &or[Valuation(v,z) lt 0:v in ycoords] then continue;end if;if PairKey([Evaluate(v,z):v in ycoords]) eq target then Append(~out,z);end if;end for;
 assert #out eq 1;return out[1];
end function;
input:=eval Read("work/deg10/continuation_20261030b/C77/census.m");rows:=input[3];assert #rows eq 3;row:=rows[CASE];R0<X>:=PolynomialRing(GF(2));L<a>:=ext<GF(2)|R0!row[2]>;vv:=func<v|Evaluate(PolynomialRing(L)!v,a)>;
E:=EllipticCurve([vv(v):v in row[3]]);PN:=E![vv(v):v in row[4]];assert Order(PN) eq 77 and #E eq 1078;
vals,jj11:=Coordinates(E,PN);start:=PlacesAt(vals);P:=start[1];assert Degree(P) eq 10;
Pdelta:=Image(P,1);Pminus:=Image(P,5);assert Image(Pdelta,5) eq P and Image(Pminus,1) eq P;
print "C77_NC_INPUT_READY",CASE,Cputime(tt);
choices:=[DivisorGroup(CY)!0];observations:=[* *];targetPrimes:=[];isogs:=MDIsogeniesFast(E,PN,3);assert &+[r[3]:r in isogs] eq 4;
function FF(vals)
 F:=Parent(vals[1]);theta:=PrimitiveElement(F);d:=Degree(F,GF(2));M:=Matrix(GF(2),[Eltseq(theta^i,GF(2)):i in [0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(GF(2),Eltseq(v,GF(2))))):v in vals];return[*[Integers()!a:a in Coefficients(MinimalPolynomial(theta,GF(2)))],[[Integers()!a:a in c]:c in co]*];
end function;
for i in [1..#isogs] do
 EE,PP,w:=Explode(isogs[i]);vs,jj11:=Coordinates(EE,PP);opts:=PlacesAt(vs);d:=Degree(opts[1]);assert (10*w) mod d eq 0;mult:=10*w div d;targetPrimes cat:=opts;choices:=[T+mult*Divisor(Q):T in choices,Q in opts];
 Append(~observations,[*i,w,d,mult,FF(aInvariants(EE) cat [PP[1],PP[2]] cat vs cat [jj11])*]);
 print "C77_NC_HEKCE_ORBIT",CASE,i,w,d,mult,Cputime(tt);
end for;
function ER(v)return[[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function E3S(v)return[ER(c):c in Eltseq(E3!v)];end function;
function EYS(v)return[E3S(c):c in Eltseq(CY!v)];end function;
function DR(v)R:=PolynomialRing(GF(2));return FB!(R!v[1])/FB!(R!v[2]);end function;
function D3(v)return E3![DR(c):c in v];end function;
function DY(v)return CY![D3(c):c in v];end function;
function EP(q)
 I:=Ideal(q);a,b:=TwoGenerators(q);assert ideal<Order(I)|a,b> eq I;vs:=[EYS(a),EYS(b)];assert DY(vs[1]) eq a and DY(vs[2]) eq b;return[*IsFinite(q) select 1 else 0,Degree(q),vs*];
end function;
function ED(D)pp,ee:=Support(D);return[*[*EP(pp[i]),ee[i]*]:i in [1..#pp]*];end function;
targetPrimes:=Setseq(Seqset(targetPrimes));targetImages:=[Image(Q,1):Q in targetPrimes];
print "C77_NC_DELTA_TARGETS_CACHED",CASE,#targetPrimes,Cputime(tt);

tests:=[* *];Bdivs:=[];
for i in [1..#choices] do
 T:=choices[i];assert Degree(T) eq 40;ps,es:=Support(T);assert &and[Q in targetPrimes:Q in ps];DT:=&+[es[j]*Divisor(targetImages[Index(targetPrimes,ps[j])]):j in [1..#ps]];
 B:=T+2*DT-7*Divisor(P)-2*Divisor(Pdelta)-3*Divisor(Pminus);assert Degree(B) eq 0;Append(~Bdivs,B);Append(~tests,[*i,ED(T),ED(DT),ED(B)*]);
end for;
core:=[*[ER(c):c in Coefficients(DefiningPolynomial(E3))],[E3S(c):c in Coefficients(DefiningPolynomial(CY))],de,[EYS(v):v in ycoords],EP(P),EP(Pdelta),EP(Pminus),row,observations*];
Write(OUTPUT cat "_input_packet.m",Sprint(core cat [*tests*]):Overwrite:=true);printf "C77_NC_B_INPUT_PACKET_READY case=%o choices=%o cpu=%o\n",CASE,#tests,Cputime(tt);
results:=[* *];remaining:=0;
for i in [1..#tests] do
 rr:=#Basis(Bdivs[i]:Reduction:=false,Simplification:="None");assert rr in [0,1];if rr gt 0 then remaining+:=1;end if;
 Append(~results,tests[i] cat [*rr*]);Write(OUTPUT cat "_partial_RR.m",Sprint(results):Overwrite:=true);printf "C77_NC_B_RR case=%o choice=%o/%o h0=%o cpu=%o\n",CASE,i,#tests,rr,Cputime(tt);
end for;
Write(OUTPUT cat "_packet.m",Sprint(core cat [*results*]):Overwrite:=true);printf "C77_NC_B_COMPLETE case=%o choices=%o remaining=%o cpu=%o\n",CASE,#tests,remaining,Cputime(tt);quit;
