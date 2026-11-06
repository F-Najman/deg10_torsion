SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);SetNthreads(1);Alarm(650);tt:=Cputime();root:="work/deg10/continuation_20261031b/C66/";
load "work/deg10/continuation_20261031b/C66/inputs/F5_rr/input33.m";
k:=GF(5);K<t>:=FunctionField(k);PK:=PolynomialRing(k);RY<y>:=PolynomialRing(K);
function Rat(v)return (K!PK!v[1])/(K!PK!v[2]);end function;
F<u>:=FunctionField(RY![Rat(v):v in FieldEquation]);assert Degree(F)eq 12 and Genus(F)eq 21;
function Fun(v)return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);blocks:=[* *];
for z in GlobalCuspBlocks do
 sig,deg,pp:=Explode(z);DD:=DivisorGroup(F)!0;
 for v in pp do dd,m,flag,ib:=Explode(v);O:=flag eq 1 select OF else OI;II:=ideal<O|[Fun(b):b in ib]>;assert IsPrime(II);P:=Place(II);assert Degree(P)eq dd;DD+:=m*Divisor(P);end for;
 assert Degree(DD)eq deg;Append(~blocks,DD);
end for;

load "work/deg10/continuation_20261031b/C66/gonal_cusp_input.m";
D0:=&+[GonalCuspVectors[1][i]*blocks[i]:i in[1..17]];V,vm:=RiemannRochSpace(D0);assert Dimension(V) eq 2;fs:=[vm(w):w in Basis(V)];g:=[f:f in fs|#Poles(f) gt 0][1];assert &+[-Valuation(g,p)*Degree(p):p in Poles(g)] eq 10;
function ER(f)return [[Integers()!q:q in Coefficients(Numerator(f))],[Integers()!q:q in Coefficients(Denominator(f))]];end function;
function EF(f)return [ER(q):q in Eltseq(F!f)];end function;
Write(root cat "gonal_first_function_oldfield.m",Sprint(EF(g)):Overwrite:=true);
for a in[0..4] do pp:=Zeros(g-GF(5)!a);vv:=[Valuation(g-GF(5)!a,p):p in pp];assert &+[vv[i]*Degree(pp[i]):i in[1..#pp]] eq 10;printf "D31B_C66_CLASSICAL_FIBRE a=%o degrees=%o multiplicities=%o time=%o\n",a,[Degree(p):p in pp],vv,Cputime(tt);end for;
printf "D31B_C66_CLASSICAL_PENCIL_PROFILE_COMPLETE time=%o\n",Cputime(tt);quit;
