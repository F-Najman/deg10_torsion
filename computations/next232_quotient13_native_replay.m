SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/C232_resolution_20261030/deps/portable_forward.m";
a:=eval Read("work/deg10/C232_resolution_20261030/quotient13_proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 13;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
fb:=[];for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A0:=Divisor(fb[a[8]]);A:=Divisor(DP(a[12]));cs:=[DP(v):v in a[11]];cc:=[q(Z!v):v in a[10]];coords:=[GFUN(v):v in a[7]];b,c,w:=Explode(coords);
load "work/deg10/C232_resolution_20261030/deps/native_group_generators.m";
ivs:=ClassGroupAbelianInvariants(F:ReductionDivisor:=A,Proof:=true);assert F`cld`proven;MDProvenClassGroupGenerators(F);GN:=F`cld`Cl;fn:=F`cld`ZnToCl;ZN:=Domain(fn);fbn:=Reverse([x`place:x in F`cld`fb]);AB:=FreeAbelianGroup(#fb);cn:=[];
for p in fb do j:=Index(fbn,p);assert j gt 0;Append(~cn,fn(ZN.j));end for;
oldmap:=hom<AB->G|fc>;newmap:=hom<AB->GN|cn>;assert Kernel(oldmap) eq Kernel(newmap);assert Invariants(G) eq Invariants(GN);print "C232_QUOTIENT13_PICARD_FRAME_REPROVED",#fb,Invariants(G),Cputime(tt);
basepl:=Support(A)[1];ai:=Index(fb,basepl);assert ai gt 0;ca:=fc[ai];
for i in [1..#cs] do j:=Index(fb,cs[i]);assert j gt 0;assert cc[i] eq fc[j]-Degree(cs[i])*ca;end for;
qd:=eval Read("work/deg10/C232_resolution_20261030/quotient13_filter.m");units:=[GFUN(v):v in qd[7]];assert [[Valuation(u,p):u in units]:p in cs] eq qd[1];
hg:=[&+[cc[i]:i in bl]:bl in qd[2]];assert [Eltseq(v):v in hg] eq qd[3];H:=sub<G|hg>;assert #H eq 512000;
function Canonical(vals)
 kk:=PrimeField(Parent(vals[1]));e:=LCM([Degree(MinimalPolynomial(v,kk)):v in vals]);j:=[i:i in [1..#vals]|Degree(MinimalPolynomial(vals[i],kk)) eq e][1];gen:=vals[j];mat:=Matrix(kk,[Eltseq(gen^i,kk):i in [0..e-1]]);pol:=MinimalPolynomial(gen,kk);repr:=[Eltseq(Solution(mat,Vector(kk,Eltseq(v,kk)))):v in vals];return [*j,Coefficients(pol),repr*];
end function;
src:=eval Read("work/deg10/C232_resolution_20261030/source6.m");K<aa>:=ext<k|PP!src[2]>;rows:=(eval Read("work/deg10/C232_resolution_20261030/quotient13_orbit_filter.m"))[1];
for row in rows do
 pl:=fb[row[4]];assert Degree(pl) eq 3;assert Canonical([Evaluate(v,pl):v in coords]) eq Canonical([Evaluate(PP!v,aa):v in row[2]]);assert q(Z!row[3]) eq fc[row[4]]-3*ca;
end for;
Write("work/deg10/C232_resolution_20261030/quotient13_native_replay.m",Sprint([*#fb,#cs,#rows,true*]):Overwrite:=true);print "C232_QUOTIENT13_NATIVE_REPLAY_COMPLETE",#fb,#cs,#rows,Cputime(tt);quit;
