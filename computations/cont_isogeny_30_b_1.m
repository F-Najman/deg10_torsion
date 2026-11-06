SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
R<x>:=PolynomialRing(Rationals());K<a>:=NumberField(R!(x^10 + x^8 + 3*x^6 - 3*x^2 + 1));E:=EllipticCurve([K|-173*a^8 - 249*a^6 - 629*a^4 - 277*a^2 + 396,2592*a^8 + 3680*a^6 + 9318*a^4 + 3906*a^2 - 6149,2592*a^8 + 3680*a^6 + 9318*a^4 + 3906*a^2 - 6149,0,0]);
P:=E![-2174*a^9 + 1353*a^8 - 3106*a^7 + 1933*a^6 - 7854*a^5 + 4888*a^4 - 3368*a^3 + 2096*a^2 + 5076*a - 3159,147755*a^9 - 93249*a^8 + 210816*a^7 - 133038*a^6 + 533236*a^5 - 336510*a^4 + 227574*a^3 - 143580*a^2 - 346157*a + 218500];Q:=E![0,0];assert Order(P) eq 2 and Order(Q) eq 30;
ks:=[P,(15)*Q,P+(15)*Q];T:=ks[2];assert Order(T) eq 2;
S<t>:=PolynomialRing(K);EE,phi:=IsogenyFromKernel(E,t-T[1]);
A,mp:=TorsionSubgroup(EE);invs:=Invariants(A);
print "DECIC_ISOGENY_FULL_TORSION",30,1,invs,Cputime(tt);
Write("work/deg10/continuation_20261029/cont_isogeny_30_b_1_torsion.txt",Sprint([*invs,[Eltseq(z):z in aInvariants(EE)],[[Eltseq(z):z in Eltseq(mp(A.j))]:j in [1..Ngens(A)]]*]));
print "DECIC_ISOGENY_EXACT_COMPLETE";quit;
