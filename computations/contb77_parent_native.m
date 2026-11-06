SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);Alarm(600);tt:=Cputime();
root:="work/deg10/continuation_20261030b/C77/";packets:=[eval Read(root cat Sprintf("nc_case%o_v3_input_packet.m",i)):i in[1..3]];p:=packets[1];assert &and[q[1] eq p[1] and q[2] eq p[2] and q[3] eq 2 and q[4] eq p[4]:q in packets];
k:=GF(2);PP:=PolynomialRing(k);B<t>:=FunctionField(k);R<S>:=PolynomialRing(B);
function DR(v)return B!(PP!v[1])/B!(PP!v[2]);end function;
E3<s>:=FunctionField(R![DR(v):v in p[1]]);assert Degree(E3) eq 12 and Genus(E3) eq 19;
function DE(v)return E3![DR(c):c in v];end function;
RY<Y>:=PolynomialRing(E3);CY<z>:=FunctionField(RY![DE(v):v in p[2]]);assert Degree(CY) eq 2 and Genus(CY) eq 37;
function DY(v)return CY![DE(c):c in v];end function;
coords:=[DY(v):v in p[4]];assert coords eq [CY!t,CY!s];coords cat:=[z];OF:=MaximalOrderFinite(CY);OI:=MaximalOrderInfinite(CY);
function RationalAt(v,x)return Evaluate(PP!Numerator(v),x)/Evaluate(PP!Denominator(v),x);end function;
function SigE(v,x)return &+[RationalAt(c,x)*(CY!s)^(i-1):i->c in Eltseq(E3!v)];end function;
for x in [(CY!t+1)/(CY!t),1/(CY!t+1)] do
 assert &+[RationalAt(c,x)*(CY!s)^(i-1):i->c in Coefficients(DefiningPolynomial(E3))] eq 0;
 assert &+[SigE(c,x)*(z+1)^(i-1):i->c in Coefficients(DefiningPolynomial(CY))] eq 0;
end for;
print "C77_PARENT_NATIVE_MODEL_READY",Genus(E3),Genus(CY),Cputime(tt);Alarm(1200);
cache:=AssociativeArray();
for p in packets do
 rows:=[*p[5],p[6],p[7]*];for test in p[10] do for j in[2,3,4] do for r in test[j] do Append(~rows,r[1]);end for;end for;end for;
 for ep in rows do ky:=Sprint(ep);if IsDefined(cache,ky) then continue;end if;OO:=ep[1] eq 1 select OF else OI;I:=ideal<OO|[DY(v):v in ep[3]]>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq ep[2];cache[ky]:=Q;end for;
end for;
function PlaceOf(ep)return cache[Sprint(ep)];end function;
function DD(rows)
 D:=DivisorGroup(CY)!0;for r in rows do D+:=r[2]*Divisor(PlaceOf(r[1]));end for;return D;
end function;
function Key(vs)
 K:=Parent(vs[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vs]);ok:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]] do gen:=&+[K!bits[i]*vs[i]:i in[1..#vs]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq d then bb:=bits;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(k,[Eltseq(gen^i,k):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vs];assert &and[&+[K!co[j][i+1]*gen^i:i in[0..d-1]] eq vs[j]:j in[1..#vs]];return [d] cat bb cat[Integers()!c:c in Coefficients(pol)]cat[Integers()!v:v in &cat co];
end function;
function PointValues(P)
 assert &and[Valuation(f,P) ge 0:f in coords];return[Evaluate(f,P):f in coords];
end function;
function ClassKeys(D,shift)
 ans:=AssociativeArray();ps,ms:=Support(D);
 for i in[1..#ps] do v:=PointValues(ps[i]);if shift eq 1 then v:=[(v[1]+1)/v[1],v[2],v[3]+1];elif shift eq -1 then v:=[1/(v[1]+1),v[2],v[3]+1];end if;ky:=Key(v);if not IsDefined(ans,ky) then ans[ky]:=0;end if;ans[ky]+:=ms[i];end for;return ans;
end function;
procedure SameKeys(A,C)
 assert Keys(A) eq Keys(C);for ky in Keys(A) do assert A[ky] eq C[ky];end for;
end procedure;
oracles:=eval Read(root cat "parent_source_oracles.m");res:=[* *];n:=0;
for ci in[1..3] do
 p:=packets[ci];P:=PlaceOf(p[5]);DP:=PlaceOf(p[6]);IP:=PlaceOf(p[7]);assert Degree(P) eq 10 and Degree(DP) eq 10 and Degree(IP) eq 10;SameKeys(ClassKeys(Divisor(P),1),ClassKeys(Divisor(DP),0));SameKeys(ClassKeys(Divisor(P),-1),ClassKeys(Divisor(IP),0));
 oc:=oracles[ci];assert oc[1] eq ci;K<a>:=ext<k|PP!oc[2]>;vv:=[Evaluate(PP!v,a):v in oc[3]];native:=PointValues(P);assert Key(vv) eq Key(native[1..2]);
 obs:=p[9];assert &+[r[2]:r in obs] eq 4 and &and[10*r[2] eq r[3]*r[4] and r[4] eq 1:r in obs];assert #p[10] eq (ci eq 2 select 16 else 4);
 for tr in p[10] do
  ix:=tr[1];T:=DD(tr[2]);DT:=DD(tr[3]);D:=DD(tr[4]);assert Degree(T) eq 40 and Degree(DT) eq 40 and Degree(D) eq 0 and D eq T+2*DT-7*Divisor(P)-2*Divisor(DP)-3*Divisor(IP);SameKeys(ClassKeys(T,1),ClassKeys(DT,0));ps,ms:=Support(T);
  for dg in Seqset([r[3]:r in obs] cat[Degree(Q):Q in ps]) do assert &+[Integers()|r[4]:r in obs|r[3] eq dg] eq &+[Integers()|ms[i]:i in[1..#ps]|Degree(ps[i]) eq dg];end for;
  hh:=#Basis(D:Reduction:=false,Simplification:="None");assert hh eq 0;n+:=1;Append(~res,[ci,ix,hh]);print "C77_PARENT_NATIVE_CASE",ci,ix,hh,Cputime(tt);
 end for;
end for;
assert n eq 24;Write(root cat "parent_native_RR.m",Sprint([*res,#Keys(cache)*]):Overwrite:=true);print "C77_PARENT_NATIVE_COMPLETE",24,#Keys(cache),Cputime(tt);quit;
