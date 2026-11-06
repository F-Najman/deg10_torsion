SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(600);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261030b/C55/";packets:=[eval Read(root cat Sprintf("source_case%o_packet.m",i)):i in[1..5]];m:=packets[1];assert &and[p[1] eq m[1]:p in packets];
k:=GF(2);PP:=PolynomialRing(k);K<x>:=FunctionField(k);R<Y>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
F<y>:=FunctionField(R![DR(v):v in m[1]]);assert Degree(F) eq 37 and Genus(F) eq 81;
function DF(v)return F![DR(c):c in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);coords:=[DF(v):v in m[2]];assert #coords eq 2 and &and[p[2] eq m[2]:p in packets];
print "C55_PARENT_NATIVE_MODEL_READY",Degree(F),Genus(F),Cputime(tt);Alarm(600);
cache:=AssociativeArray();
for p in packets do for j in[3..6] do for row in p[j] do
 ep:=row[1];ky:=Sprint(ep);if IsDefined(cache,ky) then continue;end if;O:=ep[1] eq 1 select OF else OI;J:=ideal<O|[DF(v):v in ep[3]]>;assert IsPrime(J);q:=Place(J);assert Degree(q) eq ep[2];cache[ky]:=q;
end for;end for;end for;
function DD(rows)
 D:=DivisorGroup(F)!0;for row in rows do D+:=row[2]*Divisor(cache[Sprint(row[1])]);end for;return D;
end function;
function PairKey(vs)
 L:=Parent(vs[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vs]);ok:=false;
 for bits in[[1,0],[0,1],[1,1]] do u:=L!bits[1]*vs[1]+L!bits[2]*vs[2];mp:=MinimalPolynomial(u,k);if Degree(mp) eq d then bb:=bits;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(k,[Eltseq(u^i,k):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vs];assert &and[&+[L!co[j][i+1]*u^i:i in[0..d-1]] eq vs[j]:j in[1..2]];
 return [d] cat bb cat[Integers()!v:v in Coefficients(mp)] cat[Integers()!v:v in &cat co];
end function;
function PointKey(D)
 ps,es:=Support(D);assert #ps eq 1 and es eq[1];return PairKey([Evaluate(c,ps[1]):c in coords]);
end function;
out:=[* *];
for i in[1..5] do
 p:=packets[i];Ns:=DD(p[3]);Nd:=DD(p[4]);Th:=DD(p[5]);A3:=DD(p[6]);assert Degree(Ns) eq 10 and Degree(Nd) eq 10 and Degree(Th) eq 40 and Degree(A3) eq 0 and A3 eq Th-Ns-3*Nd;
 row:=p[7];L<a>:=ext<k|PP!row[2]>;fv:=func<v|Evaluate(PolynomialRing(L)!v,a)>;E:=EllipticCurve([fv(v):v in row[3]]);PN:=E![fv(v):v in row[4]];assert Order(PN) eq 55;assert PointKey(Ns) eq PairKey(TateNormalForm_xy(E,PN));assert PointKey(Nd) eq PairKey(TateNormalForm_xy(E,3*PN));
 recorded:=AssociativeArray();actual:=AssociativeArray();ps,es:=Support(Th);for j in[1..#ps] do ky:=PairKey([Evaluate(c,ps[j]):c in coords]);if not IsDefined(actual,ky) then actual[ky]:=0;end if;actual[ky]+:=es[j];end for;
 assert &+[r[2]:r in p[8]] eq 4;
 for r in p[8] do ix,w,d,mult,ff:=Explode(r);T<b>:=ext<k|PP!ff[1]>;vv:=[Evaluate(PolynomialRing(T)!v,b):v in ff[2]];EE:=EllipticCurve(vv[1..5]);PG:=EE![vv[6],vv[7]];assert Order(PG) eq 55;xy:=TateNormalForm_xy(EE,PG);assert xy eq vv[8..9];ky:=PairKey(xy);assert ky[1] eq d and 10*w eq d*mult;if not IsDefined(recorded,ky) then recorded[ky]:=0;end if;recorded[ky]+:=mult;end for;
 assert Keys(actual) eq Keys(recorded);for ky in Keys(actual) do assert actual[ky] eq recorded[ky];end for;
 hh:=#Basis(A3:Reduction:=false,Simplification:="None");assert hh eq 0;Append(~out,[*i,Degree(Ns),Degree(Nd),Degree(Th),hh*]);print "C55_PARENT_NATIVE_CASE",i,hh,Cputime(tt);
end for;
Write(root cat "parent_native_replay.m",Sprint(out):Overwrite:=true);print "C55_PARENT_NATIVE_COMPLETE",5,#Keys(cache),Cputime(tt);quit;
