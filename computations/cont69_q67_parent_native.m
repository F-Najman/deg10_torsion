SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(16*10^9);tt:=Cputime();Alarm(300);
root:="work/deg10/next_cases_20261030/C69/quotient67/";p:=eval Read(root cat "proof_packet_v2.m");assert p[1] eq "Q67_PACKET_V2";
k:=GF(2);PP:=PolynomialRing(k);B<x>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(v)return B!(PP!v[1])/B!(PP!v[2]);end function;
K<y>:=FunctionField(R![DR(v):v in p[2]]);assert Genus(K) eq 12;
function DK(v)return K![DR(t):t in v];end function;
RZ<W>:=PolynomialRing(K);Z<z>:=FunctionField(RZ![DK(v):v in p[3]]);assert Degree(Z) eq 4 and Genus(Z) eq 67;
function DZ(v)return Z![DK(t):t in v];end function;
jK:=DK(p[4]);assert Degree(jK) eq 264 and DefiningPolynomial(Z) eq RZ![1/jK,0,0,1,1];
OK:=MaximalOrderFinite(K);OZ:=MaximalOrderFinite(Z);
function ParseK(row)
 O:=row[1] eq 1 select OK else MaximalOrderInfinite(K);I:=ideal<O|[DK(v):v in row[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq row[2];return P;
end function;
function ParseZ(row)
 O:=row[1] eq 1 select OZ else MaximalOrderInfinite(Z);I:=ideal<O|[DZ(v):v in row[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq row[2];return P;
end function;
rc:=[P:P in Poles(jK)|Degree(P) eq 1];assert #rc eq 11;f:=DefiningPolynomial(Z);assert Valuation(Evaluate(f,K!1),rc[1]) gt 0 and Valuation(Evaluate(Derivative(f),K!1),rc[1]) eq 0;
print "C69_Q67_PARENT_MODEL_READY",Cputime(tt);Alarm(1200);
zc:=AssociativeArray();kc:=AssociativeArray();zr:=[* *];for ii in [5..8] do for r in p[ii] do Append(~zr,r[1]);end for;end for;for r in p[10] do Append(~zr,r[2]);end for;
for row in zr do ky:=Sprint(row);if not IsDefined(zc,ky) then zc[ky]:=ParseZ(row);print "C69_Q67_PARENT_PRIME_LOADED",#Keys(zc),Cputime(tt);end if;end for;
for row in p[10] do ky:=Sprint(row[1]);if not IsDefined(kc,ky) then kc[ky]:=ParseK(row[1]);end if;end for;
function PK(row)return kc[Sprint(row)];end function;
function PZ(row)return zc[Sprint(row)];end function;
function DD(rows)
 D:=DivisorGroup(Z)!0;for r in rows do P:=PZ(r[1]);D+:=r[2]*Divisor(P);end for;return D;
end function;
Ns:=DD(p[5]);Nd:=DD(p[6]);Th:=DD(p[7]);A5:=DD(p[8]);assert Degree(Ns) eq 10 and Degree(Nd) eq 10 and Degree(Th) eq 60 and Degree(A5) eq 0;assert A5 eq Th-Ns-5*Nd;
pairs:=p[10];assert #pairs eq #p[9]+2;wantT:=DivisorGroup(Z)!0;checks:=[* *];nw:=0;
for i in [1..#pairs] do
 r:=pairs[i];P:=PK(r[1]);Q:=PZ(r[2]);d:=r[3];w:=r[4];mult:=r[5];assert Degree(P) eq d and Degree(Q) eq d;
 g,h:=TwoGenerators(P);un:=UniformizingElement(P);assert Valuation(Z!g,Q) gt 0 and Valuation(Z!h,Q) gt 0 and Valuation(Z!un,Q) eq 1;
 assert Valuation(jK,P) eq 0;L:=ResidueClassField(P);RL<V>:=PolynomialRing(L);ff:=RL![Evaluate(c,P):c in Coefficients(f)];fac:=Factorization(ff);assert Sort([Degree(t[1]):t in fac]) eq [1,3] and &and[t[2] eq 1:t in fac];
 if i eq 1 then assert Divisor(Q) eq Ns and d eq 10;elif i eq 2 then assert Divisor(Q) eq Nd and d eq 10;else
  ob:=p[9][i-2];assert ob[1] eq i-2 and ob[2] eq d and ob[3] eq d and ob[4] eq d and ob[5] eq w and ob[6] eq mult;
  assert d in [10,50] and mult*d eq 10*w;wantT+:=mult*Divisor(Q);nw+:=w;
 end if;
 Append(~checks,[*i,d,w,mult*]);print "C69_Q67_PARENT_PRIME_RESTRICTION",i,d,Cputime(tt);
end for;assert nw eq 6 and wantT eq Th;
print "C69_Q67_PARENT_OPERATOR_VERIFIED",Degree(Th),Degree(A5),Cputime(tt);
rr:=Cputime();bb:=Basis(A5:Reduction:=false,Simplification:="None");assert #bb eq 0;
Write(root cat "parent_native_replay.m",Sprint([*checks,Degree(Th),Degree(A5),#bb,Cputime(rr)*]):Overwrite:=true);
print "C69_Q67_PARENT_NATIVE_RR_COMPLETE",#pairs,#bb,Cputime(tt);quit;
