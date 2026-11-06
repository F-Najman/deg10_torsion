SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);SetNthreads(1);Alarm(400);tt:=Cputime();root:="work/deg10/continuation_20261031b/C66/";
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

rats:=[j:j in[1..17]|Degree(blocks[j]) eq 1];targets:=[];
for i in[1..10] do for j in rats do vv:=CuspRepresentatives[i];vv[j]+:=1;Append(~targets,vv);end for;end for;assert #targets eq 100;
ac:=[Max([0]cat[targets[i][j]:i in[1..100]]):j in[1..17]];A:=&+[ac[j]*blocks[j]:j in[1..17]];assert Degree(A)eq 58;
V,mp:=RiemannRochSpace(A);assert Dimension(V)eq 38;fs:=[mp(v):v in Basis(V)];
function CodeRat(a)return [[Integers()!c:c in Coefficients(Numerator(a))],[Integers()!c:c in Coefficients(Denominator(a))]];end function;
function CodeF(a)return [CodeRat(c):c in Eltseq(F!a)];end function;
PrintFile("work/deg10/continuation_20261031b/C66/F5_common_basis58.m",Sprintf("ACoefficients:=%o;BasisFunctions:=%o;",ac,[CodeF(f):f in fs]):Overwrite:=true);
for f in fs do for P in Poles(f)do assert Valuation(f,P)+Valuation(A,P)ge 0;end for;end for;
out:=Open("work/deg10/continuation_20261031b/C66/F5_common_basis58_jets.jsonl","w");tables:=[* *];
for j in[1..17]do
 used:=[ac[j]-targets[i][j]:i in[1..100]];assert &and[v ge 0:v in used];if Max(used)eq 0 then continue;end if;
 ps,es:=Support(blocks[j]);
 for k0 in[1..#ps]do
  P:=ps[k0];e:=es[k0];dd:=Degree(P);ord:=ac[j]*e;maxb:=Max(used)*e;LS,comp:=Completion(F,P:Precision:=100);zz:=LS.1;
  vv:=[zz^ord*comp(f):f in fs];assert &and[Valuation(v)ge 0 and AbsolutePrecision(v)ge maxb:v in vv];
  rows:=[&cat[[Integers()!c:c in Eltseq(Coefficient(v,n),k)]:n in[0..maxb-1]]:v in vv];M:=Matrix(k,rows);assert Nrows(M)eq 38 and Ncols(M)eq dd*maxb;
  orders:=[e*v:v in used];Append(~tables,<j,k0,dd,orders,M>);
  fprintf out,"{\"block_index\":%o,\"place_index\":%o,\"degree\":%o,\"A_order\":%o,\"orders\":%o,\"rows\":%o}\n",j,k0,dd,ord,orders,rows;Flush(out);
 end for;
end for;delete out;rr:=[];
for i in[1..100]do
 J:=ZeroMatrix(k,38,0);for z in tables do if z[4][i]gt 0 then J:=HorizontalJoin(J,Submatrix(z[5],1,1,38,z[3]*z[4][i]));end if;end for;
 assert Ncols(J)eq 48 and Rank(J)eq 38;Append(~rr,Rank(J));
end for;
PrintFile("work/deg10/continuation_20261031b/C66/F5_common_basis58_native.json",Sprintf("{\"complete\":true,\"A_degree\":58,\"basis_dimension\":38,\"B_degrees\":%o,\"ranks\":%o}",[48:i in[1..100]],rr):Overwrite:=true);
assert not assigned F`cld;printf "D31B_C66_COMMON_BASIS58_COMPLETE Adegree58_dim38=true all100rank38_on48=true exact_poles_checked=true no_classgroup=true time=%o\n",Cputime(tt);quit;
