SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
F<z>:=CyclotomicField(3);R<t>:=PolynomialRing(F);f:=R![(1)+(1)*z,(4)+(1)*z,(4)+(4)*z,(-22)+(-24)*z,(13)+(13)*z,(1)+(0)*z];assert Degree(f) eq 5 and IsIrreducible(f);K<w>:=NumberField(f);q:=Evaluate(R![(1969716/1316941)+(211356/1316941)*z,(1098413/1316941)+(-3691023/1316941)*z,(-8879232/1316941)+(65987/1316941)*z,(5175218/1316941)+(994317/1316941)*z,(64606/1316941)+(-320275/1316941)*z],w);
E:=EllipticCurve([(1+q)*w+(2-q),K!0,q*w^2+(1-q)*w,K!0,K!0]);P:=E![0,0];Q:=E![-w,w^2];
assert Discriminant(E) ne 0 and 3*P eq E!0 and P ne E!0 and 15*Q eq E!0 and 3*Q ne E!0 and 5*Q ne E!0;assert WeilPairing(P,5*Q,3) ne 1;
jp:=MinimalPolynomial(jInvariant(E));assert Degree(jp) eq 5;
printf "DECIC315_EXACT_LEVEL field_degree=10 relative_polynomial=%o j_minpoly=%o pairing=%o\n",f,jp,WeilPairing(P,5*Q,3);
bound:=0;g1:=0;g2:=0;
for p in PrimesUpTo(251) do
 if p mod 3 ne 1 then continue;end if;
 coeffs:=&cat[Eltseq(a):a in Coefficients(f)] cat &cat[Eltseq(a):a in Coefficients(R![(1969716/1316941)+(211356/1316941)*z,(1098413/1316941)+(-3691023/1316941)*z,(-8879232/1316941)+(65987/1316941)*z,(5175218/1316941)+(994317/1316941)*z,(64606/1316941)+(-320275/1316941)*z])];
 if not &and[Denominator(a) mod p ne 0:a in coeffs] then continue;end if;
 FF:=GF(p);SS<Z>:=PolynomialRing(FF);
 for zr in Roots(Z^2+Z+1) do
  zz:=zr[1];red:=func<a|Evaluate(Polynomial(FF,Eltseq(F!a)),zz)>;
  fp:=Polynomial(FF,[red(a):a in Coefficients(f)]);qp:=Polynomial(FF,[red(a):a in Coefficients(R![(1969716/1316941)+(211356/1316941)*z,(1098413/1316941)+(-3691023/1316941)*z,(-8879232/1316941)+(65987/1316941)*z,(5175218/1316941)+(994317/1316941)*z,(64606/1316941)+(-320275/1316941)*z])]);
  for root in Roots(fp) do
   a:=root[1];if Evaluate(Derivative(fp),a) eq 0 then continue;end if;qr:=Evaluate(qp,a);a1:=(1+qr)*a+(2-qr);a3:=qr*a^2+(1-qr)*a;
   if a3^3*(a1^3-27*a3) eq 0 then continue;end if;
   Er:=EllipticCurve([a1,FF!0,a3,FF!0,FF!0]);cnt:=#Er;inv:=Invariants(AbelianGroup(Er));assert #inv eq 2 and inv[1] mod 3 eq 0 and inv[2] mod 15 eq 0;
   bound:=GCD(bound,cnt);g1:=GCD(g1,inv[1]);g2:=GCD(g2,inv[2]);
   printf "DECIC315_GOOD_REDUCTION p=%o zeta_root=%o t_root=%o count=%o invariants=%o gcd_order=%o gcd_invariants=[%o,%o]\n",p,zz,a,cnt,inv,bound,g1,g2;
   if bound eq 45 or (g1 eq 3 and g2 eq 15) then break;end if;
  end for;
  if bound eq 45 or (g1 eq 3 and g2 eq 15) then break;end if;
 end for;
 if bound eq 45 or (g1 eq 3 and g2 eq 15) then break;end if;
end for;
printf "DECIC315_FULL_TORSION_RESULT invariants_bound=[%o,%o] order_bound=%o exact_3x15=%o time=%o\n",g1,g2,bound,bound eq 45 or (g1 eq 3 and g2 eq 15),Cputime(tt);
assert bound eq 45 or (g1 eq 3 and g2 eq 15);print "DECIC315_FULL_TORSION_VERIFIED";quit;
