SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);Alarm(300);tt:=Cputime();
root:="work/deg10/continuation_20261030b/C75/";packets:=[eval Read(root cat Sprintf("proof_packet_case%o.m",i)):i in[1..6]];p:=packets[1];assert &and[q[2] eq p[2] and q[3] eq p[3] and q[4] eq p[4]:q in packets];
k:=GF(2);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
L<a>:=FunctionField(R![DR(v):v in p[2]]);assert Degree(L) eq 5 and Genus(L) eq 12;
function DL(v)return L![DR(c):c in v];end function;
RM<W>:=PolynomialRing(L);M<w>:=FunctionField(RM![DL(v):v in p[3]]);j25:=DL(p[4]);assert DefiningPolynomial(M) eq RM![1/j25,0,0,1,1];assert Degree(j25) eq 300 and Genus(M) eq 73;
function DM(v)return M![DL(c):c in v];end function;
OF:=MaximalOrderFinite(M);OI:=MaximalOrderInfinite(M);print "C75_PARENT_NATIVE_MODEL_READY",Genus(L),Genus(M),Cputime(tt);Alarm(600);
cache:=AssociativeArray();
for p in packets do for j in[5..8] do for row in p[j] do
 ep:=row[1];ky:=Sprint(ep);if IsDefined(cache,ky) then continue;end if;O:=ep[1] eq 1 select OF else OI;I:=ideal<O|[DM(v):v in ep[3]]>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq ep[2];cache[ky]:=Q;
end for;end for;end for;
function DD(rows)
 D:=DivisorGroup(M)!0;for row in rows do D+:=row[2]*Divisor(cache[Sprint(row[1])]);end for;return D;
end function;
res:=[* *];
for ci in[1..6] do
 p:=packets[ci];assert p[1] eq ci;Ns:=DD(p[5]);Nd:=DD(p[6]);Th:=DD(p[7]);A7:=DD(p[8]);assert Degree(Ns) eq 10 and Degree(Nd) eq 10 and Degree(Th) eq 80 and Degree(A7) eq 0 and A7 eq Th-Ns-7*Nd;
 ps,ms:=Support(Ns);assert #ps eq 1 and ms eq [1];jp:=MinimalPolynomial(Evaluate(M!j25,ps[1]),k);assert jp eq PP!p[10][6] and Degree(jp) eq 10;
 qs,vs:=Support(Nd);assert #qs eq 1 and vs eq [1] and MinimalPolynomial(Evaluate(M!j25,qs[1]),k) eq jp;
 obs:=p[9];assert &+[v[1]:v in obs] eq 8 and &and[10*v[1] eq v[2]*v[3]:v in obs];ps,ms:=Support(Th);
 for dg in Seqset([v[2]:v in obs] cat[Degree(Q):Q in ps]) do assert &+[Integers()|v[3]:v in obs|v[2] eq dg] eq &+[Integers()|ms[i]:i in[1..#ps]|Degree(ps[i]) eq dg];end for;
 hh:=#Basis(A7:Reduction:=false,Simplification:="None");assert hh eq 0;Append(~res,[*ci,Degree(Ns),Degree(Nd),Degree(Th),hh*]);print "C75_PARENT_NATIVE_CASE",ci,hh,Cputime(tt);
end for;
Write(root cat "parent_native_replay.m",Sprint(res):Overwrite:=true);print "C75_PARENT_NATIVE_COMPLETE",6,#Keys(cache),Cputime(tt);quit;
