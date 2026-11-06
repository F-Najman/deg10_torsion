SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);tt:=Cputime();
load "work/deg10/continuation_20261031a/C62/inputs/C62_X131/inputs.m";
assert N eq 31 and p eq 3 and ExpectedGenus eq 26 and KnownOrder eq 134097116106800;
assert ClassGroupGenerationBound(p,ExpectedGenus) le 9;
k:=GF(p);K<t>:=RationalFunctionField(k);P<z>:=PolynomialRing(K);Pr:=PolynomialRing(k);
function DecRat35(v)return K!(Pr!v[1])/(K!(Pr!v[2]));end function;
f:=P![DecRat35(a):a in FieldPolynomialData];assert IsMonic(f);
F<u>:=FunctionField(f);assert Genus(F) eq ExpectedGenus;
function DecElt35(v)return F![DecRat35(a):a in v];end function;
coords:=[DecElt35(a):a in ModularCoordinatesData];

qq:=Places(F,10);
printf "D31A_TARGET_PLACES10_COMPLETE level=%o count=%o time=%o\n",N,#qq,Cputime(tt);
Write("work/deg10/continuation_20261031a/C62/places10_preflight.m", Sprintf("Degrees10:=%o;",[Degree(Q):Q in qq]));
quit;
