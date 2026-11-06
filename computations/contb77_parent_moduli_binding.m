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
for p in packets do for tr in p[10] do for row in tr[2] do ep:=row[1];ky:=Sprint(ep);if IsDefined(cache,ky) then continue;end if;OO:=ep[1] eq 1 select OF else OI;I:=ideal<OO|[DY(v):v in ep[3]]>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq ep[2];cache[ky]:=Q;print "C77_PARENT_T3_BINDING_PRIME",#Keys(cache),Cputime(tt);end for;end for;end for;
function Key(vs)
 K:=Parent(vs[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vs]);ok:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]] do gen:=&+[K!bits[i]*vs[i]:i in[1..#vs]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq d then bb:=bits;ok:=true;break;end if;end for;assert ok;
 M:=Matrix(k,[Eltseq(gen^i,k):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vs];assert &and[&+[K!co[j][i+1]*gen^i:i in[0..d-1]] eq vs[j]:j in[1..#vs]];return [d] cat bb cat[Integers()!c:c in Coefficients(pol)]cat[Integers()!v:v in &cat co];
end function;

procedure SameKeys(A,C)
 assert Keys(A) eq Keys(C);for ky in Keys(A) do assert A[ky] eq C[ky];end for;
end procedure;
res:=[];
for ci in[1..3] do
 p:=packets[ci];wanted:=AssociativeArray();
 for r in p[9] do
  ff:=r[5];K<a>:=ext<k|PP!ff[1]>;vals:=[Evaluate(PP!v,a):v in ff[2]];assert #vals eq 10;ky:=Key(vals[8..9]);assert ky[1] eq r[3] and r[4] eq 1 and 10*r[2] eq r[3]*r[4];if not IsDefined(wanted,ky) then wanted[ky]:=0;end if;wanted[ky]+:=r[4];
 end for;
 for tr in p[10] do
  actual:=AssociativeArray();for row in tr[2] do Q:=cache[Sprint(row[1])];assert &and[Valuation(v,Q) ge 0:v in coords[1..2]];ky:=Key([Evaluate(v,Q):v in coords[1..2]]);assert ky[1] eq Degree(Q);if not IsDefined(actual,ky) then actual[ky]:=0;end if;actual[ky]+:=row[2];end for;SameKeys(actual,wanted);Append(~res,[ci,tr[1],#Keys(actual)]);
 end for;print "C77_PARENT_T3_MODULI_CASE",ci,#p[10],Cputime(tt);
end for;
assert #res eq 24;Write(root cat "parent_T3_moduli_binding.m",Sprint(res):Overwrite:=true);print "C77_PARENT_T3_MODULI_BINDING_COMPLETE",24,#Keys(cache),Cputime(tt);quit;
