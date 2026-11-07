/* Z/6 x Z/12 (Proposition prop:6x12 b)): the pencil |L| modulo a prime above 13.
   Run from the repository root:  magma -b additions/C6x12/square13.m
   L is the effective F-rational cuspidal divisor of degree 10 of Lemma lem:spor-612-double, given by its coefficients
   on the 24 Galois orbits of cusps in jobs/cont_reduce612_quintic_lattice_v2.log, which is in the class number 2791 of
   continuation_20261029/M6N12_D5_all_source_lifts.m. We check that the two divisors are linearly equivalent over F_13.
   We reduce L modulo the prime of F=Q(zeta_6) above 13 at which zeta_6 = 4 (the Weil pairing e_6(P,2Q) is normalized
   to 4 as in computations/cont_RR612_norm13_shard*.m), compute l(L) on X_1(6,12) over F_13, and go through all the
   (13^l-1)/12 effective divisors T in |L| over F_13. T = 2A forces the polynomial of the push-forward of T by u
   (outside u = infinity) to be a square, and for these T we compute the support and multiplicities of T.
   This is computations/cont_RR612_norm13_shard6.m restricted to the class 2791, without the early stop at the first
   square member, and with the model of continuation_20261029/deps/model_6_12.m read from a repository-relative path. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("shared/mdmagma/v2/mdmagma.spec");
Attach("continuation_20261029/deps/fast_hecke.m");
Attach("continuation_20261029/deps/fast_hecke_11.m");
// the model of MDX612 in continuation_20261029/deps/model_6_12.m
function MDX612Repo(F)
 A<u,v>:=AffineSpace(F,2);
 data:=Split(Read("shared/mdmagma/models_X1_m_n/X1_6_12.txt"));
 for line in data do
  val:=Split(Split(line,"=")[2],";")[1];
  if line[1] eq "X" then equation:=eval(val);end if;
  if line[1] eq "r" then r:=eval(val);end if;
  if line[1] eq "s" then s:=eval(val);end if;
  if line[1] eq "q" then q:=eval(val);end if;
  if line[1] eq "E" then es:=eval(val);end if;
 end for;
 C:=ProjectiveClosure(Curve(A,equation));FF:=FunctionField(C);
 X:=New(MDCrvMod11);X`M:=6;X`N:=12;
 X`_P:=[FF!q];X`_Q:=[FF!0,FF!0];X`_coordinates:=[FF.1,FF.2];
 _InitMDCrvMod(X,12,C,[FF!a:a in es]);
 return X;
end function;
function CuspBlocks(X)
C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);Q:=E![0,0];px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);assert #ys eq 2;P:=E![px,ys[1][1]];
z:=PrimitiveElement(BaseRing(X))^((#BaseRing(X)-1) div 6);if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq z;
print "SQUARE13_ZETA6",z;
fs:=[-ai[2],1-ai[1],Discriminant(E)];anchor:=P[1];
for i in [0..5] do for j in [0..11] do
 A:=i*P+j*Q;if A eq E!0 then continue;end if;
 if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
 if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
end for;end for;
cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
return blocks,sigs,ds;
end function;
a:=eval Read("continuation_20261029/M6N12_D5_all_source_lifts.m");lifts:=a[3];
X:=MDX612Repo(GF(13));C:=Curve(X);FF:=FunctionField(C);AF,af:=AlgorithmicFunctionField(FF);B:=BaseField(AF);OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);
assert Genus(C) eq 9;
blocks,sigs,ds:=CuspBlocks(X);assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];cs:=Seqset(Cusps(X));
print "SQUARE13_CUSP_BLOCK_DEGREES",ds;
function PushPolynomial(D)
 result:=B!1;pp,mm:=Support(D);
 for i in [1..#pp] do
  I:=Ideal(pp[i]);
  if Order(I) eq OF then f:=Minimum(I);assert Degree(f) gt 0 and Degree(pp[i]) mod Degree(f) eq 0;result*:=(B!f)^(mm[i]*(Degree(pp[i]) div Degree(f)));
  else assert Order(I) eq OI;end if;
 end for;return result;
end function;
h0:=PushPolynomial(FunctionFieldDivisor(Divisor(c0)));hh:=[PushPolynomial(FunctionFieldDivisor(b)):b in blocks];
function Tuples(F,k)
 if k eq 0 then return [[]];end if;return [[a] cat v:a in F,v in Tuples(F,k-1)];
end function;
k:=2791;
v:=lifts[k][2];L:=10*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];
// the effective representative of the class of L given in jobs/cont_reduce612_quintic_lattice_v2.log
e:=[0,0,0,1,0,0,2,0,0,0,2,0,1,0,0,0,0,0,0,0,1,1,0,0];Leff:=&+[e[i]*blocks[i]:i in [1..24]];
assert Degree(Leff) eq 10 and IsPrincipal(Leff-L);
print "SQUARE13_L_EFFECTIVE",k,e,"degree",Degree(Leff);
V,mp:=RiemannRochSpace(L);dim:=Dimension(V);print "SQUARE13_ELL",k,dim;assert dim eq 2;
hl:=h0^10*&*[(hh[i]/h0^ds[i])^v[i]:i in [1..24]];
members:=0;squarepol:=0;positive:=[];
for j in [1..dim] do for tail in Tuples(GF(13),dim-j) do
 u:=[GF(13)!0:i in [1..j-1]] cat [GF(13)!1] cat tail;f:=mp(V!u);members+:=1;
 zz:=hl*(B!Norm(af(f)));assert Degree(Denominator(zz)) eq 0;
 pol:=Numerator(zz);pol/:=LeadingCoefficient(pol);assert Degree(pol) le 10;
 if IsSquare(pol) then
  squarepol+:=1;
  T:=L+Divisor(f);assert IsEffective(T) and Degree(T) eq 10;sp,mults:=Support(T);degs:=[Degree(P):P in sp];
  even:=&and[IsEven(m):m in mults];cusp:=&and[P in cs:P in sp];
  normT:=PushPolynomial(FunctionFieldDivisor(T));assert Degree(Denominator(normT)) eq 0;pp:=Numerator(normT);pp/:=LeadingCoefficient(pp);assert pp eq pol;
  printf "SQUARE13_SQUARE_NORM member=%o degrees=%o multiplicities=%o all_even=%o all_cuspidal=%o\n",u,degs,mults,even,cusp;
  if even and not cusp then Append(~positive,<u,degs,mults>);end if;
 end if;
end for;end for;
assert members eq (13^dim-1) div 12;
printf "SQUARE13_COMPLETE class=%o ell=%o members=%o square_norms=%o twice_noncuspidal=%o time=%o\n",k,dim,members,squarepol,#positive,Cputime(tt);
quit;
