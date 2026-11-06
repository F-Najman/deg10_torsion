SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261029b/deps/native_group_generators.m";
X:=MDX11(2,26,GF(3));C:=Curve(X);assert Genus(C) eq 25;cs:=Cusps(X);rat:=[c:c in cs|Degree(c) eq 1];cub:=[c:c in cs|Degree(c) eq 3];assert #rat eq 18 and #cub eq 6;c0:=rat[1];F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);base:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c0));
print "C226_PROVEN_PICARD_SETUP_START",Cputime(tt);
ivs:=ClassGroupAbelianInvariants(F:ReductionDivisor:=base,Proof:=true);assert assigned F`cld and F`cld`proven;print "C226_CLASSGROUP_PROVEN",ivs,Cputime(tt);
MDProvenClassGroupGenerators(F);G:=F`cld`Cl;ff:=F`cld`ZnToCl;Z:=Domain(ff);fb:=Reverse([a`place:a in F`cld`fb]);assert #fb eq Ngens(Z);classes:=[ff(Z.j):j in [1..#fb]];order:=&*[e:e in Invariants(G)|e ne 0];assert order eq 124501550040000;
AA:=F`cld`A1;pp,ee:=Support(AA);assert #pp eq 1 and ee[1] eq 1;ai:=Index(fb,pp[1]);assert ai gt 0 and Degree(pp[1]) eq 1;
function FBIndex(c)
 DD:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c));pp,ee:=Support(DD);assert #pp eq 1 and ee[1] eq 1;j:=Index(fb,pp[1]);assert j gt 0;return j;
end function;
ri:=[FBIndex(c):c in rat];ci:=[FBIndex(c):c in cub];H:=sub<G|[classes[j]-classes[ai]:j in ri]>;assert #H eq 355718714400;
function EncRat(v) return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EncFun(v) return [EncRat(c):c in Eltseq(F!v)];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);flags:=[Order(Ideal(p)) eq OF select 1 else 0:p in fb];degs:=[Degree(p):p in fb];ideals:=[*[EncFun(v):v in Basis(Ideal(p))]:p in fb*];
coords:=[EncFun(F!mp(v)):v in X`_coordinates];fpoly:=[EncRat(v):v in Coefficients(DefiningPolynomial(F))];fp:=FrobeniusPolynomial(JH(52,[27]),3);assert Integers()!Evaluate(fp,1) eq order;
Write("work/deg10/continuation_20261029b/C226/proven_picard.m",Sprint([*fpoly,flags,degs,ideals,[Eltseq(v):v in classes],[Eltseq(r):r in Rows(RelationMatrix(G))],coords,ai,ri,ci,Coefficients(fp),#H,Invariants(H)*]));
print "C226_PROVEN_PORTABLE_PICARD_COMPLETE",#fb,order,#H,Cputime(tt);quit;
