SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/C232_resolution_20261030/hecke11_geometry_v2.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 37;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
function EP(pl) I:=Ideal(pl);return [*Order(I) eq OF select 1 else 0,Degree(pl),[EF(v):v in Basis(I)]*];end function;
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
cg:=eval Read("work/deg10/C232_resolution_20261030/hecke5_cusp_geometry.m");ps:=[DP(s):s in cg[2]];vectors:=eval Read("work/deg10/C232_resolution_20261030/hecke7_cusp_candidates.m");hg:=eval Read("work/deg10/C232_resolution_20261030/hecke7_geometry.m");assert hg[1] eq a[1];N:=Divisor(DP(hg[5]));DN:=Divisor(DP(hg[6]));T:=&+[r[10]*Divisor(DP(r[11])):r in hg[7]];assert Degree(T) eq 48;AA:=T-N-7*DN;assert Degree(AA) eq 0;tests:=[* *];remaining:=[];
for i in [1..#vectors] do
 v:=vectors[i];C:=&+[DivisorGroup(F)|v[j]*Divisor(ps[j]):j in [1..#ps]];assert Degree(C) eq 0;D:=AA+C;basis:=Basis(D:Reduction:=false,Simplification:="None");yes:=#basis gt 0;Append(~tests,[*i,#basis,yes*]);if yes then Append(~remaining,i);end if;
 print "C232_HECKE7_CUSP_TEST",i,#vectors,yes,Cputime(tt);Write("work/deg10/C232_resolution_20261030/hecke7_cusp_filter.m",Sprint([*tests,remaining*]):Overwrite:=true);
end for;
print "C232_HECKE7_CUSP_FILTER_COMPLETE",#tests,remaining,Cputime(tt);quit;
