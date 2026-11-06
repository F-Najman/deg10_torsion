SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
R<x>:=PolynomialRing(Rationals());K<a>:=NumberField(R!(x^10 + x^8 - 2*x^6 + x^4 + 2*x^2 - 2));E:=EllipticCurve([K|-17*a^8 - 31*a^6 + 8*a^4 - 11*a^2 - 41,-419*a^8 - 759*a^6 + 222*a^4 - 239*a^2 - 1032,-419*a^8 - 759*a^6 + 222*a^4 - 239*a^2 - 1032,0,0]);
P:=E![-143/4*a^9 + 117/4*a^8 - 517/8*a^7 + 423/8*a^6 + 77/4*a^5 - 125/8*a^4 - 81/4*a^3 + 17*a^2 - 88*a + 289/4,-2445/2*a^9 + 9725/8*a^8 - 35441/16*a^7 + 35243/16*a^6 + 1293/2*a^5 - 10283/16*a^4 - 5583/8*a^3 + 2775/4*a^2 - 12047/4*a + 23955/8];Q:=E![0,0];assert Order(P) eq 2 and Order(Q) eq 28;
ks:=[P,(14)*Q,P+(14)*Q];T:=ks[3];assert Order(T) eq 2;
S<t>:=PolynomialRing(K);EE,phi:=IsogenyFromKernel(E,t-T[1]);
A,mp:=TorsionSubgroup(EE);invs:=Invariants(A);
print "DECIC_ISOGENY_FULL_TORSION",28,2,invs,Cputime(tt);
Write("work/deg10/continuation_20261029/cont_isogeny_28_a_2_torsion.txt",Sprint([*invs,[Eltseq(z):z in aInvariants(EE)],[[Eltseq(z):z in Eltseq(mp(A.j))]:j in [1..Ngens(A)]]*]));
print "DECIC_ISOGENY_EXACT_COMPLETE";quit;
