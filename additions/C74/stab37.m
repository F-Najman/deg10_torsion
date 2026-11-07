// Lemma lem:sw-74-J c), e). Orbit lengths under Delta(37) of the non-cuspidal closed points of X_1(37) over F_2 of
// degrees 6 and 9, and the diamond operator <10>, of order 3, on them.
// The points are the representatives of the diamond orbits, in the order of
// continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/prime_good_branches_64GB.m,
// which wrote them to P37/input_records.m and whose log P37/good_branches_64GB.log lists the surviving indices
// (degree 6: [1], degree 9: [3,4]). The key of a closed point is the one of P37/finite_nc9_survivors.m.
// Run from the repository root: magma -b additions/C74/stab37.m
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("shared/mdmagma/v2/mdmagma.spec");
tt:=Cputime();N:=37;
function XYKey(E,P)
 xx,yy:=Explode(TateNormalForm_xy(E,P));
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));
  if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d] cat bb cat [Integers()!c:c in Eltseq(mp)] cat [Integers()!c:c in &cat cs];
end function;
assert Modexp(10,3,37) eq 1;
load "continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/input_records.m";
lengths:=[];
for rec in Inputs do
 PF:=PolynomialRing(GF(2));F<z>:=ext<GF(2)|PF!rec[3]>;E:=EllipticCurve([F!PF!v:v in rec[4]]);P:=E![F!PF!v:v in rec[5]];
 assert Order(P) eq N and XYKey(E,P)[1] eq rec[1];
 orb:=#{XYKey(E,k*P):k in [1..36]};fixed:=XYKey(E,10*P) eq XYKey(E,P);
 assert fixed eq (orb in [1,2,3,6]);
 Append(~lengths,<rec[1],rec[2],orb>);
 printf "STAB37 degree=%o index=%o orbit_length=%o fixed_by_<10>=%o\n",rec[1],rec[2],orb,fixed;
end for;
assert lengths eq [<6,1,3>,<6,2,18>,<7,1,18>,<9,1,18>,<9,2,18>,<9,3,2>,<9,4,1>];
printf "STAB37_DONE: orbit lengths 3,18 (degree 6), 18 (degree 7), 18,18,2,1 (degree 9), <10> fixes the orbits of length 3, 2, 1 (%o s)\n",Cputime(tt);
quit;
