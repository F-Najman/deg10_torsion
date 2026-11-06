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
data:=eval Read("work/deg10/C232_resolution_20261030/quotient13_filter.m");hg:=[q(Z!v):v in data[3]];H:=sub<G|hg>;assert #H eq 512000;QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 32;
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A0,"work/deg10/C232_resolution_20261030/quotient13_orbits_RR.log");src:=eval Read("work/deg10/C232_resolution_20261030/source6.m");ky:=(eval Read("work/deg10/C232_resolution_20261030/filterprobe6.m"))[3][1][3];K<aa>:=ext<k|PP!src[2]>;b32:=Evaluate(PP!ky[1..6],aa);c32:=Evaluate(PP!ky[7..12],aa);EX:=EllipticCurve([K|0,c32,0,b32*(1-b32-c32),0]);P:=EX![0,0];S:=EX![b32,b32];
function Place3(vals)
 e:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);assert e eq 3;chosen:=[j:j in [1..3]|Degree(MinimalPolynomial(vals[j],k)) eq 3][1];gen:=vals[chosen];mat:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);vs:=[Eltseq(Solution(mat,Vector(k,Eltseq(v,k)))):v in vals];pol:=MinimalPolynomial(gen,k);gg:=coords[chosen];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!vs[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];I:=ideal<OF|fs>;fa:=Factorization(I);pp:=[Place(r[1]):r in fa|r[2] gt 0];assert #pp eq 1 and Degree(pp[1]) eq e;return pp[1];
end function;
rows:=[* *];common:={QH|};first:=true;
for bb,dd in [0,1],u in [1,3,5,7,9,11,13,15] do
 PS:=P+bb*(16*S);SS:=dd*P+u*S;Q16:=2*SS;T:=16*SS;vals:=DerickxNormalForm_bc(EX,PS,Q16);alpha:=(Q16[1]-PS[1])/Q16[2];bt:=3*T[1]^2+2*c32*T[1]+b32*(1-b32-c32);ww:=alpha*((SS[1]-T[1])^2-bt)/(2*SS[2]);Append(~vals,ww);assert &and[v^27 eq v:v in vals];npl:=Place3(vals);v:=phi(Divisor(npl)-3*A);j:=Index(fb,npl);assert j gt 0;assert v eq fc[j]-3*phi(A);
 poss:={QH|};checks:=[* *];
 for i in [1..#cs] do if Degree(cs[i]) ne 2 then continue;end if;z:=qh(v+cc[i]);good:=2*z eq QH!0;if good then Include(~poss,z);end if;Append(~checks,[*i,Eltseq(z),good*]);end for;
 if first then common:=poss;first:=false;else common meet:=poss;end if;
 Append(~rows,[*[bb,dd,u],[[Integers()!x:x in Eltseq(a,k)]:a in vals],Eltseq(v),j,checks,[Eltseq(z):z in poss],[Eltseq(z):z in common]*]);print "C232_QUOTIENT13_ORBIT_TEST",[bb,dd,u],#poss,#common,Cputime(tt);
end for;
Write("work/deg10/C232_resolution_20261030/quotient13_orbit_filter.m",Sprint([*rows,[Eltseq(z):z in common]*]):Overwrite:=true);fprintf witness,"]\n";Flush(witness);print "C232_QUOTIENT13_ORBIT_FILTER_COMPLETE",#rows,#common,Cputime(tt);quit;
