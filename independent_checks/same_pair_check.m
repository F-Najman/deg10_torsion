// Written for this repository: for each N, the curve of computations/witness_CN.m (record model) and the curve
// built from van Hoeij's coordinates (census record inputs/1.N.10.*.json) are defined over isomorphic fields and
// have the same j-invariant, and the curves are isomorphic over K (after an isomorphism of the fields), so they
// give the same pair (K,E).
SetColumns(0);
R<t>:=PolynomialRing(Rationals());
S<X,Y>:=PolynomialRing(Rationals(),2);

K1<a>:=NumberField(R![1,-7,20,-37,51,-56,47,-31,15,-5,1]); E1:=EllipticCurve([1-Evaluate(R![-3,34,-135,291,-405,390,-273,134,-44,8],a),-Evaluate(R![44,-435,1759,-3899,5512,-5368,3787,-1868,616,-112],a),-Evaluate(R![44,-435,1759,-3899,5512,-5368,3787,-1868,616,-112],a),0,0]);
K2<x>:=NumberField(t^10 - 3*t^9 + 8*t^7 - 5*t^6 - 6*t^5 + 5*t^4 - 2*t^2 + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(29*Y-2*X^9-X^8+11*X^7+8*X^6-49*X^5+77*X^3-35*X^2-46*X+13,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=29 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![1,0,2,3,3,7,8,7,7,4,1]); E1:=EllipticCurve([1-Evaluate(R![-25,26,-12,-46,24,-14,-6,102,90,26],a),-Evaluate(R![104,290,157,513,1041,1155,1670,1825,956,207],a),-Evaluate(R![104,290,157,513,1041,1155,1670,1825,956,207],a),0,0]);
K2<x>:=NumberField(t^10 - 4*t^9 + 3*t^8 + 6*t^7 - 2*t^6 - 8*t^5 - 8*t^4 + 11*t^3 + 6*t^2 - 5*t + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(215*Y-109*X^9+266*X^8+207*X^7-801*X^6-606*X^5+304*X^4+1470*X^3+410*X^2-1000*X-71,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=31 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![-1,1,3,3,3,3,1,4,8,5,1]); E1:=EllipticCurve([1-Evaluate(R![10,-17,-23,-20,-22,-21,-1,-44,-61,-19],a),-Evaluate(R![-109,165,241,201,222,211,-1,435,645,209],a),-Evaluate(R![-109,165,241,201,222,211,-1,435,645,209],a),0,0]);
K2<x>:=NumberField(t^10 + 2*t^8 + t^7 - 9*t^6 - 6*t^5 + 7*t^4 + 8*t^3 - t^2 - 3*t - 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(191*Y-39*X^9+192*X^8-127*X^7+454*X^6+452*X^5-1286*X^4-862*X^3+479*X^2+443*X-7,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=35 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![1,-3,0,9,-7,-8,10,3,-5,0,1]); E1:=EllipticCurve([1-Evaluate(R![8,-41,59,14,-99,44,56,-39,-7,9],a),-Evaluate(R![126,-555,741,201,-1205,537,640,-453,-75,102],a),-Evaluate(R![126,-555,741,201,-1205,537,640,-453,-75,102],a),0,0]);
K2<x>:=NumberField(t^10 - 5*t^8 + 3*t^7 + 10*t^6 - 8*t^5 - 7*t^4 + 9*t^3 - 3*t + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(Y-2*X^9-X^8+10*X^7-X^6-23*X^5+6*X^4+22*X^3-11*X^2-9*X+5,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=37 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![-1,-2,0,7,8,-1,-9,-6,1,3,1]); E1:=EllipticCurve([1-Evaluate(R![-1,-4,1,14,7,-8,-13,-2,5,2],a),-Evaluate(R![-8,-18,8,69,44,-44,-68,-17,27,13],a),-Evaluate(R![-8,-18,8,69,44,-44,-68,-17,27,13],a),0,0]);
K2<x>:=NumberField(t^10 + 3*t^9 + t^8 - 6*t^7 - 9*t^6 - t^5 + 8*t^4 + 7*t^3 - 2*t - 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(Y-2*X^9-5*X^8+10*X^6+11*X^5-2*X^4-10*X^3-5*X^2+X,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=38 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![1,-4,5,-1,-9,13,-6,-4,7,-4,1]); E1:=EllipticCurve([1-Evaluate(R![11,-51,73,-7,-126,193,-60,-93,100,-32],a),-Evaluate(R![248,12,-2110,2122,-284,-5359,4932,801,-3393,1574],a),-Evaluate(R![248,12,-2110,2122,-284,-5359,4932,801,-3393,1574],a),0,0]);
K2<x>:=NumberField(t^10 + 2*t^8 + 14*t^7 - 6*t^6 - 20*t^5 + 4*t^4 + 7*t^3 - 2*t^2 + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(449*Y-748*X^9-643*X^8-2185*X^7-12166*X^6-6406*X^5+8021*X^4+6917*X^3+192*X^2-243*X+889,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=39 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![-1,0,-2,0,0,0,0,0,1,0,1]); E1:=EllipticCurve([1-Evaluate(R![-3,-3/2,-5,-1/2,4,11/2,7/2,-1/2,-1/2,-3],a),-Evaluate(R![16,25/2,34,33/2,-9,-51/2,-25,-9,-11/2,17/2],a),-Evaluate(R![16,25/2,34,33/2,-9,-51/2,-25,-9,-11/2,17/2],a),0,0]);
K2<x>:=NumberField(t^10 - 7*t^8 - 5*t^7 + 17*t^6 + t^5 - 11*t^4 + t^3 + t^2 + 2*t - 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(137*Y-712*X^9-515*X^8+4624*X^7+6920*X^6-7151*X^5-6022*X^4+3395*X^3+1608*X^2+565*X-912,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=40 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![-1,0,2,2,2,-4,-5,-3,3,4,1]); E1:=EllipticCurve([1-Evaluate(R![6,7,-11,-19,-7,13,13,3,0,0],a),-Evaluate(R![-9,6,29,-6,-30,-25,22,19,-2,-2],a),-Evaluate(R![-9,6,29,-6,-30,-25,22,19,-2,-2],a),0,0]);
K2<x>:=NumberField(t^10 - 5*t^9 + 4*t^8 + 14*t^7 - 24*t^6 + 6*t^5 + 5*t^4 - 6*t^3 + 8*t^2 - 5*t + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(3*Y-16*X^9+68*X^8-12*X^7-237*X^6+206*X^5+74*X^4-35*X^3+61*X^2-80*X+16,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=42 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![-1,3,-1,1,-12,10,6,-14,11,-5,1]); E1:=EllipticCurve([1-Evaluate(R![86/27,-368/27,1006/27,-66,244/9,1255/27,-1912/27,1325/27,-61/3,103/27],a),-Evaluate(R![-694/27,3205/27,-6203/27,1786/9,122/3,-6782/27,6995/27,-4072/27,481/9,-236/27],a),-Evaluate(R![-694/27,3205/27,-6203/27,1786/9,122/3,-6782/27,6995/27,-4072/27,481/9,-236/27],a),0,0]);
K2<x>:=NumberField(t^10 - 5*t^9 + 11*t^8 - 14*t^7 + 6*t^6 + 10*t^5 - 12*t^4 + t^3 - t^2 + 3*t - 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(135*Y-19*X^9+123*X^8-335*X^7+496*X^6-331*X^5-243*X^4+609*X^3-217*X^2-157*X-158,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=45 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;

K1<a>:=NumberField(R![1,0,-2,2,3,-10,15,-14,9,-4,1]); E1:=EllipticCurve([1-Evaluate(R![-18/5,-44/5,14/5,1/5,-56/5,92/5,-69/5,7,-12/5,-4/5],a),-Evaluate(R![96/5,117/5,-104/5,-17/5,299/5,-118,584/5,-77,184/5,-36/5],a),-Evaluate(R![96/5,117/5,-104/5,-17/5,299/5,-118,584/5,-77,184/5,-36/5],a),0,0]);
K2<x>:=NumberField(t^10 - 2*t^9 + t^8 + 7*t^6 - 10*t^5 + 4*t^4 - 2*t^3 + 5*t^2 - 4*t + 1); P<u>:=PolynomialRing(K2); yy:=Roots(Evaluate(49*Y+255*X^9-435*X^8+104*X^7+45*X^6+1804*X^5-2005*X^4+266*X^3-426*X^2+1118*X-622,[P!x,u]))[1][1];
r:=(x^2*yy-x*yy+yy-1)/(x^2*yy-x); s:=(x*yy-yy+1)/(x*yy); E2:=EllipticCurve([1-s*(r-1),-r*s*(r-1),-r*s*(r-1),0,0]);
ok,phi:=IsIsomorphic(K1,K2); samej:= ok and exists{psi : psi in Automorphisms(K2) | psi(phi(jInvariant(E1))) eq jInvariant(E2)};
isocurve:= ok and exists{psi : psi in Automorphisms(K2) | IsIsomorphic(EllipticCurve([psi(phi(c)) : c in aInvariants(E1)]), E2)};
printf "SAMEPAIR N=50 fields_isomorphic=%o same_j=%o curves_isomorphic=%o\n", ok, samej, isocurve;
print "SAMEPAIR_DONE";
quit;
