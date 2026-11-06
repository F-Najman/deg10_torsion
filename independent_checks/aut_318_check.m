// Written for this repository: the automorphism group of the field K_{3,18} of Proposition 11.12 has order 2, and its
// non-trivial element does not fix zeta_3 (used in Remark after Proposition 11.12).
SetColumns(0);
Qx<x>:=PolynomialRing(Rationals());
F<z>:=NumberField(x^2+x+1); FW<W>:=PolynomialRing(F);
L:=NumberField(3*W^5+(1+5*z)*W^4+(6*z-3)*W^3-(6*z+9)*W^2-(5*z+4)*W+3);
A:=AbsoluteField(L); zA:=A!z;
auts:=Automorphisms(A);
printf "AUT318 order=%o moves_zeta3=%o\n", #auts, [a(zA) ne zA : a in auts];
quit;
