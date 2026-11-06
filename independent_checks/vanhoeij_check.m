// Written for this repository: the curves of Proposition 3.4 built directly from van Hoeij's coordinates (as in the
// census records inputs/1.N.10.*.json), their torsion subgroups, and the isomorphism of Q(x) with the polredabs field.
SetColumns(0);SetMemoryLimit(8*10^9);
R<t>:=PolynomialRing(Rationals());
S<X,Y>:=PolynomialRing(Rationals(),2);

t0:=Cputime(); N:=29; f:=t^10 - 3*t^9 + 8*t^7 - 5*t^6 - 6*t^5 + 5*t^4 - 2*t^2 + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=29*Y-2*X^9-X^8+11*X^7+8*X^6-49*X^5+77*X^3-35*X^2-46*X+13;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 2*t^7 + t^5 + 2*t^4 + 3*t^3 - 3*t^2 - 2*t + 1);
printf "CHECK vanHoeij id=1.29.10.b N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=31; f:=t^10 - 4*t^9 + 3*t^8 + 6*t^7 - 2*t^6 - 8*t^5 - 8*t^4 + 11*t^3 + 6*t^2 - 5*t + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=215*Y-109*X^9+266*X^8+207*X^7-801*X^6-606*X^5+304*X^4+1470*X^3+410*X^2-1000*X-71;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 + 2*t^8 - 3*t^7 + 3*t^6 - 7*t^5 + 8*t^4 - 7*t^3 + 7*t^2 - 4*t + 1);
printf "CHECK vanHoeij id=1.31.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=35; f:=t^10 + 2*t^8 + t^7 - 9*t^6 - 6*t^5 + 7*t^4 + 8*t^3 - t^2 - 3*t - 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=191*Y-39*X^9+192*X^8-127*X^7+454*X^6+452*X^5-1286*X^4-862*X^3+479*X^2+443*X-7;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 2*t^9 + 4*t^7 - 4*t^6 + 3*t^4 - t^3 - t^2 + 1);
printf "CHECK vanHoeij id=1.35.10.b N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=37; f:=t^10 - 5*t^8 + 3*t^7 + 10*t^6 - 8*t^5 - 7*t^4 + 9*t^3 - 3*t + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=Y-2*X^9-X^8+10*X^7-X^6-23*X^5+6*X^4+22*X^3-11*X^2-9*X+5;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 3*t^9 + 7*t^8 - 11*t^7 + 12*t^6 - 8*t^5 + 2*t^4 + 4*t^3 - 6*t^2 + 4*t - 1);
printf "CHECK vanHoeij id=1.37.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=38; f:=t^10 + 3*t^9 + t^8 - 6*t^7 - 9*t^6 - t^5 + 8*t^4 + 7*t^3 - 2*t - 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=Y-2*X^9-5*X^8+10*X^6+11*X^5-2*X^4-10*X^3-5*X^2+X;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 3*t^9 + 3*t^8 - t^7 + t^5 - 3*t^4 + 4*t^3 - t^2 - t + 1);
printf "CHECK vanHoeij id=1.38.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=39; f:=t^10 + 2*t^8 + 14*t^7 - 6*t^6 - 20*t^5 + 4*t^4 + 7*t^3 - 2*t^2 + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=449*Y-748*X^9-643*X^8-2185*X^7-12166*X^6-6406*X^5+8021*X^4+6917*X^3+192*X^2-243*X+889;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 4*t^9 + 7*t^8 - 4*t^7 - 6*t^6 + 13*t^5 - 9*t^4 - t^3 + 5*t^2 - 4*t + 1);
printf "CHECK vanHoeij id=1.39.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=40; f:=t^10 - 7*t^8 - 5*t^7 + 17*t^6 + t^5 - 11*t^4 + t^3 + t^2 + 2*t - 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=137*Y-712*X^9-515*X^8+4624*X^7+6920*X^6-7151*X^5-6022*X^4+3395*X^3+1608*X^2+565*X-912;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 + t^8 - 2*t^2 - 1);
printf "CHECK vanHoeij id=1.40.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=42; f:=t^10 - 5*t^9 + 4*t^8 + 14*t^7 - 24*t^6 + 6*t^5 + 5*t^4 - 6*t^3 + 8*t^2 - 5*t + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=3*Y-16*X^9+68*X^8-12*X^7-237*X^6+206*X^5+74*X^4-35*X^3+61*X^2-80*X+16;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 2*t^8 - 2*t^7 - 2*t^6 + 4*t^5 + 5*t^4 + 3*t^3 - 3*t^2 - 4*t - 1);
printf "CHECK vanHoeij id=1.42.10.b N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=45; f:=t^10 - 5*t^9 + 11*t^8 - 14*t^7 + 6*t^6 + 10*t^5 - 12*t^4 + t^3 - t^2 + 3*t - 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=135*Y-19*X^9+123*X^8-335*X^7+496*X^6-331*X^5-243*X^4+609*X^3-217*X^2-157*X-158;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - t^9 + 2*t^8 - 2*t^6 + 5*t^5 - 2*t^4 + 2*t^2 - t + 1);
printf "CHECK vanHoeij id=1.45.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);

t0:=Cputime(); N:=50; f:=t^10 - 2*t^9 + t^8 + 7*t^6 - 10*t^5 + 4*t^4 - 2*t^3 + 5*t^2 - 4*t + 1; K<a>:=NumberField(f); P<u>:=PolynomialRing(K);
G:=49*Y+255*X^9-435*X^8+104*X^7+45*X^6+1804*X^5-2005*X^4+266*X^3-426*X^2+1118*X-622;
yy:=Roots(Evaluate(G,[P!a,u]))[1][1]; xx:=a;
r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx); s:=(xx*yy-yy+1)/(xx*yy); b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); T:=TorsionSubgroup(E);
L:=NumberField(t^10 - 4*t^9 + 9*t^8 - 14*t^7 + 15*t^6 - 10*t^5 + 3*t^4 + 2*t^3 - 2*t^2 + 1);
printf "CHECK vanHoeij id=1.50.10.a N=%o order_of_00=%o torsion=%o iso_to_polredabs=%o time=%o\n", N, Order(E![0,0]), Invariants(T), IsIsomorphic(K,L), Cputime(t0);
quit;