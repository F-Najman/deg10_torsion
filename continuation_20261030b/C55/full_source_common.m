SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030b/C55/deps/fast_hecke.m");Attach("work/deg10/continuation_20261030b/C55/deps/fast_hecke_11.m");Attach("work/deg10/continuation_20261030b/C55/deps/fast_hecke_cyclic_commonzero.m");
input:=eval Read("work/deg10/continuation_20261030b/C55/census.m");rows:=[r:r in input[3]|r[1] eq 10];assert #rows eq 5;row:=rows[CASE];
k:=GF(2);R<x>:=PolynomialRing(k);L<a>:=ext<k|R!row[2]>;V:=func<seq|Evaluate(PolynomialRing(L)!seq,a)>;
E:=EllipticCurve([V(v):v in row[3]]);PN:=E![V(v):v in row[4]];assert Order(PN) eq 55 and #E eq 990;
X:=MDX1(55,k);assert Genus(X) eq 81;p:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=PN>);assert Degree(p) eq 10;
printf "C55_SOURCE_A3_SETUP case=%o degree=%o cpu=%o\n",CASE,Degree(p),Cputime(tt);
T:=DivisorGroup(Curve(X))!0;observations:=[* *];
function FF(vals)
 F:=Parent(vals[1]);z:=PrimitiveElement(F);n:=Degree(F,k);MB:=Matrix(k,[Eltseq(z^i,k):i in [0..n-1]]);assert Rank(MB) eq n;
 vecs:=[Eltseq(Solution(MB,Vector(k,Eltseq(v,k)))):v in vals];assert &and[Evaluate(PolynomialRing(F)!vecs[i],z) eq vals[i]:i in [1..#vals]];
 return[*[Integers()!v:v in Coefficients(MinimalPolynomial(z,k))],[[Integers()!v:v in c]:c in vecs]*];
end function;
isogs:=MDIsogeniesFast(E,PN,3);assert &+[r[3]:r in isogs] eq 4;
for i in [1..#isogs] do
 EE,PP,w:=Explode(isogs[i]);assert Order(PP) eq 55;q:=ModuliPointCyclicCommonZeros(X,EE,rec<recformat<P:PtEll>|P:=PP>);assert (10*w) mod Degree(q) eq 0;mult:=10*w div Degree(q);T+:=mult*Divisor(q);
 vals:=aInvariants(EE) cat [PP[1],PP[2]] cat TateNormalForm_xy(EE,PP);
 Append(~observations,[*i,w,Degree(q),mult,FF(vals)*]);
 printf "C55_SOURCE_A3_IMAGE case=%o image=%o obs=%o degree=%o norm=%o cpu=%o\n",CASE,i,w,Degree(q),mult,Cputime(tt);
end for;
dp:=ModuliPointCyclicCommonZeros(X,E,rec<recformat<P:PtEll>|P:=3*PN>);assert Degree(dp) eq 10;assert Degree(T) eq 40;
A:=T-Divisor(p)-3*Divisor(dp);assert Degree(A) eq 0;
F,am:=AlgorithmicFunctionField(FunctionField(Curve(X)));assert Type(BaseField(F)) eq FldFunRat;DA:=DivisorGroup(F)!FunctionFieldDivisor(A);assert Degree(DA) eq 0;
function ER(v)return[[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v)return[ER(c):c in Eltseq(F!v)];end function;
function DR(v)R:=PolynomialRing(k);return BaseField(F)!(R!v[1])/BaseField(F)!(R!v[2]);end function;
function DF(v)return F![DR(c):c in v];end function;
function EP(q)
 I:=Ideal(q);f,g:=TwoGenerators(q);assert ideal<Order(I)|f,g> eq I;vs:=[EF(f),EF(g)];assert DF(vs[1]) eq f and DF(vs[2]) eq g;
 return[*IsFinite(q) select 1 else 0,Degree(q),vs*];
end function;
function ED(D)
 pp,ee:=Support(D);return[*[*EP(pp[i]),ee[i]*]:i in [1..#pp]*];
end function;
function FD(D)return DivisorGroup(F)!FunctionFieldDivisor(D);end function;
packet:=[*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(am(FunctionField(Curve(X)).i)):i in [1..2]],ED(FD(Divisor(p))),ED(FD(Divisor(dp))),ED(FD(T)),ED(DA),row,observations*];
Write(OUTPUT cat "_packet.m",Sprint(packet):Overwrite:=true);
printf "C55_SOURCE_A3_PACKET case=%o fielddegree=%o cpu=%o\n",CASE,Degree(F),Cputime(tt);
rr:=Dimension(RiemannRochSpace(DA));assert rr in [0,1];Write(OUTPUT cat "_result.m",Sprint([*CASE,rr,Cputime(tt)*]):Overwrite:=true);
printf "C55_SOURCE_A3_COMPLETE case=%o h0=%o cpu=%o\n",CASE,rr,Cputime(tt);quit;
