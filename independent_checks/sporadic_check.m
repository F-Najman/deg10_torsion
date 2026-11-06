// Independent check (written for the repository) of the explicit curves of Section 11:
// torsion subgroups with TorsionSubgroup, j-invariants (degree over Q, integrality, minimal polynomials),
// and the absolute fields (discriminant, signature) for comparison with the LMFDB labels.
SetColumns(0); SetMemoryLimit(16*10^9);
Qx<x>:=PolynomialRing(Rationals());
procedure Report(name, E)
  K:=BaseRing(E); t0:=Cputime();
  T:=TorsionSubgroup(E); j:=jInvariant(E); mj:=MinimalPolynomial(j);
  printf "CHECK %o torsion=%o deg_j=%o j_integral=%o minpoly_j=%o disc=%o signature=%o time=%o\n",
    name, Invariants(T), Degree(mj), IsIntegral(j), mj, Discriminant(MaximalOrder(K)), Signature(K), Cputime(t0);
end procedure;

// Z/6 x Z/12, F = Q(z), z^2-z+1=0
F<z>:=NumberField(x^2-x+1); FW<W>:=PolynomialRing(F);
g:=3*W^5+(4-5*z)*W^4+(7-8*z)*W^3-7*(z+1)*W^2+6*(z-1)*W+2*z-1;
L<w>:=NumberField(g); A:=AbsoluteField(L); w:=A!w; zA:=A!z;
r:=(w^2+1)/(w^2-w); s:=(w^2+w+1)/w^2;
E612:=EllipticCurve([1+s-r*s, r*s*(1-r), r*s*(1-r), 0, 0]);
Report("6x12", E612);

// Z/3 x Z/15, F = Q(z), z^2+z+1=0
F<z>:=NumberField(x^2+x+1); FW<W>:=PolynomialRing(F);
g:=W^5+13*(z+1)*W^4-(24*z+22)*W^3+4*(z+1)*W^2+(z+4)*W+z+1;
L<w>:=NumberField(g); A:=AbsoluteField(L); w:=A!w; z:=A!z;
q:=((293*z+23)*w^2+(105-238*z)*w+357-28*z)/(217*w^2+(603*z+271)*w+283+86*z);
E315:=EllipticCurve([(1+q)*w+2-q, 0, q*w^2+(1-q)*w, 0, 0]);
Report("3x15", E315);

// Z/3 x Z/18 and its 2-isogenous curve
F<z>:=NumberField(x^2+x+1); FW<W>:=PolynomialRing(F);
g:=3*W^5+(1+5*z)*W^4+(6*z-3)*W^3-(6*z+9)*W^2-(5*z+4)*W+3;
L<u>:=NumberField(g); A:=AbsoluteField(L); u:=A!u; z:=A!z;
t:=((4*z-25)*u^4+(12-65*z)*u^3+(20+31*z)*u^2+(69+30*z)*u+3*z-33)/19;
E318:=EllipticCurve([(u+2)*t+1-u, 0, (u+1)*t^2-u*t, 0, 0]);
Q:=E318![-t,t^2]; printf "CHECK 3x18 order(0,0)=%o order(Q)=%o\n", Order(E318![0,0]), Order(Q);
T2:=9*Q; R<X>:=PolynomialRing(A);
E318b:=IsogenyFromKernel(E318, X-T2[1]);
Report("3x18", E318); Report("3x18_isogenous", E318b);
printf "CHECK 3x18 distinct_j_minpolys=%o\n", MinimalPolynomial(jInvariant(E318)) ne MinimalPolynomial(jInvariant(E318b));

// Z/2 x Z/26
for data in [* <[-621631,-2586377,-2555940,-4632176,2207386,13819447,4462921,-4689287,-3688261,-721870],
                [129155,551533,1269576,2528052,1757370,-2078127,-1932589,-133261,347217,92982],
                406832, 203416, 2*x^10+11*x^9+17*x^8-7*x^7-42*x^6-20*x^5+9*x^4+9*x^3+10*x^2+6*x+1, "2x26_a">,
             <[-27507,-95105,-100181,101187,190643,-236919,-349828,-76640,41312,17265],
                [90619,307177,325777,-321491,-622235,755143,1128928,249316,-133040,-55877],
                784, 1960, x^10+3*x^9-3*x^8-23*x^7-26*x^6+3*x^5+13*x^4-2*x^3-9*x^2-5*x-1, "2x26_b"> *] do
  K<a>:=NumberField(data[5]);
  b:=Evaluate(Qx!data[1],a)/data[3]; c:=Evaluate(Qx!data[2],a)/data[4];
  E:=EllipticCurve([0,c,0,b*(1-b-c),0]);
  printf "CHECK %o order(P)=%o order(Q)=%o\n", data[6], Order(E![0,0]), Order(E![b,b]);
  Report(data[6], E);
end for;
quit;
