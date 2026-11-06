SetColumns(0);SetNthreads(1);SetQuitOnError(true);Alarm(500);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N27/";md:=eval Read(root cat "quotient_model.m");rows:=eval Read(root cat "quotient_A7_RR.m");src:=eval Read(root cat "physical_source5.m");F2:=GF(2);P2:=PolynomialRing(F2);k<w>:=ext<F2|P2!md[1]>;assert #k eq 4;
function DC(n)return k!(n mod 2)+k!(n div 2)*w;end function;
function EC(a)v:=Eltseq(k!a,F2);return Integers()!v[1]+2*Integers()!v[2];end function;
B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);PK:=PolynomialRing(k);
function DR(a)return B!(PK![DC(n):n in a[1]])/(B!(PK![DC(n):n in a[2]]));end function;
f:=R![DR(a):a in md[2]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 13;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;coords:=[DF(v):v in md[3]];
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
function Key4(v)
 L:=Parent(v[1]);Embed(k,L);d:=LCM([Degree(MinimalPolynomial(a,k)):a in v]);ok:=false;
 for bits in[[k!1,k!0],[k!0,k!1]]cat[[k!1,a]:a in k|a ne 0]do z:=&+[L!bits[i]*v[i]:i in[1,2]];mp:=MinimalPolynomial(z,k);if Degree(mp)eq d then chosen:=bits;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(k,[Eltseq(z^i,k):i in[0..d-1]]);assert Rank(M)eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(a,k)))):a in v];return[*d,[EC(a):a in chosen],[EC(a):a in Coefficients(mp)],[[EC(a):a in c]:c in co]*];
end function;
K<z>:=ext<F2|P2!src[1]>;Embed(k,K);out:=[* *];np:=0;
for row in rows do
 i:=row[1];P:=Prime(row[3]);V:=Prime(row[5]);np+:=2;assert row[4]*Degree(P)eq 5 and row[6]*Degree(V)eq 35;
 kp:=Key4([Evaluate(v,P):v in coords]);kd:=Key4([Evaluate(v,V):v in coords]);assert kp eq Key4([Evaluate(P2!v,z):v in row[7]])and kd eq Key4([Evaluate(P2!v,z):v in row[8]]);
 D:=-row[4]*Divisor(P)-row[6]*Divisor(V);hec:=[* *];assert &+[h[1]:h in row[9]]eq 8;
 for h in row[9]do T:=Prime(h[3]);np+:=1;assert 5*h[1]eq h[2]*Degree(T);L<a>:=ext<F2|P2!h[5]>;Embed(k,L);kt:=Key4([Evaluate(v,T):v in coords]);assert kt eq Key4([Evaluate(P2!v,a):v in h[4]]);D+:=h[2]*Divisor(T);Append(~hec,[*h[1],h[2],kt*]);end for;
 assert Degree(D)eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");assert h0 eq row[10]and h0 eq 0;Append(~out,[*i,kp,kd,hec,h0*]);print "M3N27_PARENT_NATIVE_CASE",i,h0,Cputime(tt);
end for;
assert #out eq 6;Write(root cat "parent_native_replay.m",Sprint(out):Overwrite:=true);print "D31B_M3N27_PARENT_NATIVE_COMPLETE",#out,np,Cputime(tt);quit;
