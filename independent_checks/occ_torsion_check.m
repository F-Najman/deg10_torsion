// Written for this repository: torsion subgroups, computed with TorsionSubgroup, of the explicit curves of
// Propositions 3.4 and 3.5, from the models used in computations/witness_C*.m and computations/quintic_full2_*_v2.m.
SetColumns(0);SetMemoryLimit(8*10^9);
R<t>:=PolynomialRing(Rationals());

t0:=Cputime(); N:=29; f:=R![1,-7,20,-37,51,-56,47,-31,15,-5,1]; K<a>:=NumberField(f); b:=Evaluate(R![44,-435,1759,-3899,5512,-5368,3787,-1868,616,-112],a); c:=Evaluate(R![-3,34,-135,291,-405,390,-273,134,-44,8],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=31; f:=R![1,0,2,3,3,7,8,7,7,4,1]; K<a>:=NumberField(f); b:=Evaluate(R![104,290,157,513,1041,1155,1670,1825,956,207],a); c:=Evaluate(R![-25,26,-12,-46,24,-14,-6,102,90,26],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=35; f:=R![-1,1,3,3,3,3,1,4,8,5,1]; K<a>:=NumberField(f); b:=Evaluate(R![-109,165,241,201,222,211,-1,435,645,209],a); c:=Evaluate(R![10,-17,-23,-20,-22,-21,-1,-44,-61,-19],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=37; f:=R![1,-3,0,9,-7,-8,10,3,-5,0,1]; K<a>:=NumberField(f); b:=Evaluate(R![126,-555,741,201,-1205,537,640,-453,-75,102],a); c:=Evaluate(R![8,-41,59,14,-99,44,56,-39,-7,9],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=38; f:=R![-1,-2,0,7,8,-1,-9,-6,1,3,1]; K<a>:=NumberField(f); b:=Evaluate(R![-8,-18,8,69,44,-44,-68,-17,27,13],a); c:=Evaluate(R![-1,-4,1,14,7,-8,-13,-2,5,2],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=39; f:=R![1,-4,5,-1,-9,13,-6,-4,7,-4,1]; K<a>:=NumberField(f); b:=Evaluate(R![248,12,-2110,2122,-284,-5359,4932,801,-3393,1574],a); c:=Evaluate(R![11,-51,73,-7,-126,193,-60,-93,100,-32],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=40; f:=R![-1,0,-2,0,0,0,0,0,1,0,1]; K<a>:=NumberField(f); b:=Evaluate(R![16,25/2,34,33/2,-9,-51/2,-25,-9,-11/2,17/2],a); c:=Evaluate(R![-3,-3/2,-5,-1/2,4,11/2,7/2,-1/2,-1/2,-3],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=42; f:=R![-1,0,2,2,2,-4,-5,-3,3,4,1]; K<a>:=NumberField(f); b:=Evaluate(R![-9,6,29,-6,-30,-25,22,19,-2,-2],a); c:=Evaluate(R![6,7,-11,-19,-7,13,13,3,0,0],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=42; f:=R![-1,0,2,2,2,-4,-5,-3,3,4,1]; K<a>:=NumberField(f); b:=Evaluate(R![-9,6,29,-6,-30,-25,22,19,-2,-2],a); c:=Evaluate(R![6,7,-11,-19,-7,13,13,3,0,0],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=45; f:=R![-1,3,-1,1,-12,10,6,-14,11,-5,1]; K<a>:=NumberField(f); b:=Evaluate(R![-694/27,3205/27,-6203/27,1786/9,122/3,-6782/27,6995/27,-4072/27,481/9,-236/27],a); c:=Evaluate(R![86/27,-368/27,1006/27,-66,244/9,1255/27,-1912/27,1325/27,-61/3,103/27],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); N:=50; f:=R![1,0,-2,2,3,-10,15,-14,9,-4,1]; K<a>:=NumberField(f); b:=Evaluate(R![96/5,117/5,-104/5,-17/5,299/5,-118,584/5,-77,184/5,-36/5],a); c:=Evaluate(R![-18/5,-44/5,14/5,1/5,-56/5,92/5,-69/5,7,-12/5,-4/5],a);
E:=EllipticCurve([1-c,-b,-b,0,0]); T,m:=TorsionSubgroup(E);
printf "CHECK cyclic N=%o irreducible=%o deg=%o torsion=%o order_of_00=%o signature=%o disc=%o time=%o\n", N, IsIrreducible(f), Degree(K), Invariants(T), Order(E![0,0]), Signature(K), Discriminant(MaximalOrder(K)), Cputime(t0);

t0:=Cputime(); f:=R![2, 2, -1, -2, -1, 1]; K5<a>:=NumberField(f); b:=Evaluate(R![1032, -239, -222, -759, 419],a); c:=Evaluate(R![42, -11, -8, -31, 17],a);
E5:=EllipticCurve([1-c,-b,-b,0,0]); T5:=TorsionSubgroup(E5);
fac:=[g[1]:g in Factorization(DivisionPolynomial(E5,2))]; h:=[g:g in fac|Degree(g) eq 2][1];
L:=AbsoluteField(NumberField(h)); EL:=ChangeRing(E5,L); TL:=TorsionSubgroup(EL);
printf "CHECK full2 N=28 quintic_torsion=%o two_div_degrees=%o decic_torsion=%o absdeg=%o signature=%o disc=%o time=%o\n", Invariants(T5), [Degree(g):g in fac], Invariants(TL), Degree(L), Signature(L), Discriminant(MaximalOrder(L)), Cputime(t0);

t0:=Cputime(); f:=R![1, 3, 0, -3, 1, 1]; K5<a>:=NumberField(f); b:=Evaluate(R![-9118, -14321, 20475, -1917, -6378],a); c:=Evaluate(R![-162, -257, 366, -34, -114],a);
E5:=EllipticCurve([1-c,-b,-b,0,0]); T5:=TorsionSubgroup(E5);
fac:=[g[1]:g in Factorization(DivisionPolynomial(E5,2))]; h:=[g:g in fac|Degree(g) eq 2][1];
L:=AbsoluteField(NumberField(h)); EL:=ChangeRing(E5,L); TL:=TorsionSubgroup(EL);
printf "CHECK full2 N=30 quintic_torsion=%o two_div_degrees=%o decic_torsion=%o absdeg=%o signature=%o disc=%o time=%o\n", Invariants(T5), [Degree(g):g in fac], Invariants(TL), Degree(L), Signature(L), Discriminant(MaximalOrder(L)), Cputime(t0);
quit;