SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
K<z>:=CyclotomicField(11);A:=EllipticCurve([K|0,-1,1,0,0]);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Af:=EllipticCurve([GF(23)|0,-1,1,0,0]);assert #Af eq 25;assert (GF(23)!2)^11 eq 1 and GF(23)!2 ne 1;
printf "DECIC_FULL11_X111_TORSION_BOUND p=23 residue_degree=1 count=%o\n",#Af;
X:=MDX1(11,K);xx,yy:=Explode(X`_coordinates_xy);uu:=1/yy;vv:=xx/yy-1;cs:=Cusps(X);assert #cs eq 10 and &and[Degree(c) eq 1:c in cs];found:=false;
for cc in cs do
 if Valuation(uu,cc) lt 0 or Valuation(vv,cc) lt 0 then continue;end if;
 Q:=A![Evaluate(uu,cc),Evaluate(vv,cc)];
 if 25*Q eq A!0 and 5*Q ne A!0 then found:=true;break;end if;
end for;assert found;
coeff:=func<el|[[Numerator(v),Denominator(v)]:v in Eltseq(K!el)]>;
Write("work/deg10/exclusions/full11_generator.json",Sprintf("{\"field_polynomial\":%o,\"generator\":%o}",[[Numerator(v),Denominator(v)]:v in Coefficients(DefiningPolynomial(K))],[coeff(Q[1]),coeff(Q[2])]));
printf "DECIC_FULL11_X111_GENERATOR order=25 coordinates=%o time=%o\n",Eltseq(Q),Cputime(tt);
noncusps:=0;cusps:=0;
for i in [0..24] do
 q:=i*Q;
 if q eq A!0 or q[1] eq 0 then cusps+:=1;continue;end if;
 x:=(q[2]+1)/q[1];y:=1/q[1];assert x^2*y-x*y^2+y-1 eq 0;
 if x eq 0 or x*y eq 1 then cusps+:=1;continue;end if;
 r:=(x^2*y-x*y+y-1)/(x*(x*y-1));s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);
 // Check the discriminant polynomial before constructing a possibly singular curve.
 aa:=[1-c,-b,-b,K!0,K!0];a1,a2,a3,a4,a6:=Explode(aa);b2:=a1^2+4*a2;b4:=a1*a3+2*a4;b6:=a3^2+4*a6;b8:=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
 if -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6 eq 0 then cusps+:=1;continue;end if;
 E:=EllipticCurve(aa);P:=E![0,0];assert 11*P eq E!0 and P ne E!0;noncusps+:=1;killed:=false;
 for p in [23,67,89] do
  if not &and[Denominator(v) mod p ne 0:v in &cat[Eltseq(a):a in aa]] then continue;end if;
  F:=GF(p);f:=Polynomial(F,DefiningPolynomial(K));
  for rec in Roots(f) do
   zr:=rec[1];red:=func<a|Evaluate(Polynomial(F,Eltseq(a)),zr)>;aar:=[red(a):a in aa];
   a1,a2,a3,a4,a6:=Explode(aar);b2:=a1^2+4*a2;b4:=a1*a3+2*a4;b6:=a3^2+4*a6;b8:=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;
   if -b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6 eq 0 then continue;end if;
   Er:=EllipticCurve(aar);assert #Er lt 121;
   printf "DECIC_FULL11_MODULI_CASE i=%o j=%o good_prime=%o cyclotomic_root=%o reduction_order=%o excludes_full11=true\n",i,jInvariant(E),p,zr,#Er;killed:=true;break;
  end for;
  if killed then break;end if;
 end for;
 assert killed;
end for;
assert cusps eq 10 and noncusps eq 15;
printf "DECIC_FULL11_ALL_CASES_COMPLETE total=25 cusps=%o noncusps=%o all_excluded=true time=%o\n",cusps,noncusps,Cputime(tt);quit;
