SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(180);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";
a:=eval Read(root cat"global_cusp_model_p19.m");rep:=eval Read(root cat"endomorphism_gluing/short_cusp_divisors.m");assert a[1]eq 19 and a[4]eq rep[1] and a[5]eq rep[2];k:=GF(19);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(x)return (B!(PK!x[1]))/(B!(PK!x[2]));end function;
f:=R![DR(x):x in a[2]];assert IsIrreducible(f);F<v>:=FunctionField(f);assert Genus(F)eq 10;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(x)return F![DR(z):z in x];end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;

function ER(z)return [[Integers()!w:w in Coefficients(Numerator(z))],[Integers()!w:w in Coefficients(Denominator(z))]];end function;
function EF(z)return [ER(c):c in Eltseq(F!z)];end function;
function PD(P)g,h:=TwoGenerators(P);O:=IsFinite(P)select OF else OI;assert ideal<O|g,h>eq Ideal(P);return[*Degree(P),IsFinite(P)select 1 else 0,[EF(g),EF(h)]*];end function;
function FX(vals)
 K:=Parent(vals[1]);theta:=PrimitiveElement(K);d:=Degree(K,k);M:=Matrix(k,[Eltseq(theta^i,k):i in[0..d-1]]);assert Rank(M)eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(z,k)))):z in vals];return[*[Integers()!x:x in Coefficients(MinimalPolynomial(theta,k))],[[Integers()!x:x in row]:row in co]*];
end function;
blocks:=[&+[Divisor(Prime(p)):p in ps]:ps in a[3]];assert [Degree(D):D in blocks]eq a[5];allc:=[Prime(pd):ps in a[3],pd in ps];out:=[* *];qt:=[DF(v):v in a[6]];
for row in rep[3]do
 id,coeff:=Explode(row);A:=&+[coeff[i]*blocks[i]:i in[1..24]];assert Degree(A)eq 5;vv,mp:=RiemannRochSpace(A);assert Dimension(vv)eq 1;f:=mp(vv.1);D:=A+Divisor(f);assert Degree(D)eq 5;ps,es:=Support(D);assert &and[n gt 0:n in es];IA,JA:=Ideals(A);ID,JD:=Ideals(D);assert ID eq IA*(f*OF)and JD eq JA*(f*OI);obs:=[* *];
 for j in[1..#ps]do
  P:=ps[j];isc:=P in allc;assert not isc;vals:=[Evaluate(v,P):v in qt];qq,tt:=Explode(vals);EE:=EllipticCurve([qq*tt-qq+tt+2,Parent(qq)!0,qq*tt^2-qq*tt+tt,Parent(qq)!0,Parent(qq)!0]);PP:=EE![0,0];QQ:=EE![-tt,tt^2];assert Order(PP)eq 3 and Order(QQ)eq 18 and WeilPairing(PP,6*QQ,3)ne 1;gg,hh:=AbelianGroup(EE);Append(~obs,[*Degree(P),es[j],PD(P),FX(vals),[Integers()!z:z in Coefficients(MinimalPolynomial(jInvariant(EE),k))],Invariants(gg)*]);
 end for;
 Append(~out,[*id,coeff,EF(f),obs*]);print "M318_EFFECTIVE_F19_DIVISOR",id,[(<v[1],v[2],v[5],v[6]>):v in obs],Cputime(tt);
end for;
Write(root cat"endomorphism_gluing/effective_cusp_divisors_p19.m",Sprint(out):Overwrite:=true);print "D31B_M318_EFFECTIVE_F19_COMPLETE",#out,Cputime(tt);quit;
