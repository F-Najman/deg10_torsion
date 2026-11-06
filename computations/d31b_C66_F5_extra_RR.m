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

rats:=[j:j in[1..17]|Degree(blocks[j]) eq 1];assert #rats eq 10;out:=[* *];
for i in[1..10] do
 C:=&+[CuspRepresentatives[i][j]*blocks[j]:j in[1..17]];assert Degree(C) eq 9;
 for j in rats do V,mp:=RiemannRochSpace(C+blocks[j]);Append(~out,[i,j,Dimension(V)]);printf "D31B_C66_F5_EXTRA index=%o cusp=%o h0=%o time=%o\n",i,j,Dimension(V),Cputime(tt);end for;
end for;
Write(root cat "F5_degree10_signed_RR.m",Sprint(out):Overwrite:=true);printf "D31B_C66_F5_SIGNED_RR_COMPLETE tests=100 dims=%o no_classgroup=true time=%o\n",[o[3]:o in out],Cputime(tt);quit;
