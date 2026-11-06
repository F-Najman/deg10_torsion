SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
R<x>:=PolynomialRing(Rationals());K<a>:=NumberField(R!(x^10 - 2*x^7 - 3*x^6 - 3*x^5 - 3*x^4 - 2*x^3 + 1));E:=EllipticCurve([K|25*a^9 + 14*a^8 + 7*a^7 - 47*a^6 - 103*a^5 - 134*a^4 - 151*a^3 - 136*a^2 - 77*a - 43,1376*a^9 + 765*a^8 + 403*a^7 - 2579*a^6 - 5658*a^5 - 7402*a^4 - 8381*a^3 - 7540*a^2 - 4288*a - 2435,1376*a^9 + 765*a^8 + 403*a^7 - 2579*a^6 - 5658*a^5 - 7402*a^4 - 8381*a^3 - 7540*a^2 - 4288*a - 2435,0,0]);
P:=E![5*a^9 + 19/4*a^8 + 25/4*a^7 - 1/2*a^6 - 9*a^5 - 61/4*a^4 - 41/2*a^3 - 81/4*a^2 - 12*a - 29/4,-1731/8*a^9 + 111/2*a^8 + 1383/4*a^7 + 9449/8*a^6 + 7731/4*a^5 + 18295/8*a^4 + 9439/4*a^3 + 7839/4*a^2 + 8659/8*a + 4467/8];Q:=E![0,0];assert Order(P) eq 2 and Order(Q) eq 30;
ks:=[P,(15)*Q,P+(15)*Q];T:=ks[1];assert Order(T) eq 2;
S<t>:=PolynomialRing(K);EE,phi:=IsogenyFromKernel(E,t-T[1]);
A,mp:=TorsionSubgroup(EE);invs:=Invariants(A);
print "DECIC_ISOGENY_FULL_TORSION",30,0,invs,Cputime(tt);
Write("work/deg10/continuation_20261029/cont_isogeny_30_a_0_torsion.txt",Sprint([*invs,[Eltseq(z):z in aInvariants(EE)],[[Eltseq(z):z in Eltseq(mp(A.j))]:j in [1..Ngens(A)]]*]));
print "DECIC_ISOGENY_EXACT_COMPLETE";quit;
