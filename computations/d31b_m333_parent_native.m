SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(800);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N33/";md:=eval Read(root cat "rational_trace_model.m");initial:=eval Read(root cat "rational_trace_A5_RR.m");three:=eval Read(root cat "rational_trace_3isog_A5_RR.m");src:=eval Read(root cat "physical_source5.m");k:=GF(2);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(a)return B!(PK!a[1])/(B!(PK!a[2]));end function;
f:=R![DR(a):a in md[2]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 21;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;coords:=[DF(v):v in md[3]];
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
function Key(v)
 L:=Parent(v[1]);dd:=LCM([Degree(MinimalPolynomial(a,k)):a in v]);ok:=false;
 for bits in[[1,0],[0,1],[1,1]]do gen:=&+[L!bits[j]*v[j]:j in[1,2]];pol:=MinimalPolynomial(gen,k);if Degree(pol)eq dd then chosen:=bits;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(k,[Eltseq(gen^i,k):i in[0..dd-1]]);assert Rank(M)eq dd;co:=[Eltseq(Solution(M,Vector(k,Eltseq(a,k)))):a in v];return[*dd,chosen,[Integers()!a:a in Coefficients(pol)],[[Integers()!a:a in r]:r in co]*];
end function;
selected:=[*r:r in initial|r[10]eq 0*];assert #selected eq 9;
for i in[1,6,10]do rr:=[*r:r in three|r[1]eq i and r[13]eq 0*];assert #rr gt 0;r:=rr[1];Append(~selected,[*r[1],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12],r[13],r[14]*]);end for;
assert #selected eq 12 and Seqset([r[1]:r in selected])eq {1..12};K<z>:=ext<k|PK!src[1]>;np:=0;out:=[* *];
for row in selected do
 P:=Prime(row[3]);V:=Prime(row[5]);np+:=2;assert row[4]*Degree(P)eq 10 and row[6]*Degree(V)eq 50;
 kp:=Key([Evaluate(v,P):v in coords]);kd:=Key([Evaluate(v,V):v in coords]);assert kp eq Key([Evaluate(PK!v,z):v in row[7]])and kd eq Key([Evaluate(PK!v,z):v in row[8]]);
 D:=-row[4]*Divisor(P)-row[6]*Divisor(V);assert &+[h[1]:h in row[9]]eq 6;hec:=[* *];
 for h in row[9]do T:=Prime(h[3]);np+:=1;assert h[2]*Degree(T)eq 10*h[1];L<a>:=ext<k|PK!h[5]>;kt:=Key([Evaluate(v,T):v in coords]);assert kt eq Key([Evaluate(PK!v,a):v in h[4]]);D+:=h[2]*Divisor(T);Append(~hec,[*h[1],h[2],kt*]);end for;
 assert Degree(D)eq 0;h0:=#Basis(D);assert h0 eq row[10]and h0 eq 0;Append(~out,[*row[1],kp,kd,hec,h0*]);print "M333_PARENT_NATIVE_CASE",row[1],h0,Cputime(tt);
end for;
Write(root cat "parent_native_selected12.m",Sprint(out):Overwrite:=true);print "D31B_M333_PARENT_NATIVE_COMPLETE",#out,np,Cputime(tt);quit;
