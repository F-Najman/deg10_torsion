/* Z/90 (Lemma lem:fb-45 c), used in Proposition prop:90): l(G)=1 for every non-cuspidal closed point G of degree 10 of X_1(45) over F_2.
   Run from the repository root:  magma -b additions/C90/ell45.m
   Input: the census continuation_20261031a/C90/source45.m written by computations/d31a_90_source.m (the 108 closed
   points of degree 10, given by one representative (b,c) of each of the 10 orbits under the diamond operators).
   For each representative G we check that its diamond orbit consists of closed points of degree 10, that the 10
   orbits are disjoint and contain 108 points in total, and we compute l(G) on the model MDX1(45) of mdmagma.
   As the diamond operators are automorphisms of X_1(45)_{F_2}, l is constant on the orbits. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("shared/mdmagma/v2/mdmagma.spec");
load "continuation_20261031a/C90/inputs/pair_af_place_linear.m";
k:=GF(2);PP:=PolynomialRing(k);
X:=MDX1(45,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);
assert Genus(F) eq 41;coords:=[F!mp(v):v in X`_coordinates_xy];
print "C90_ELL45_MODEL genus",Genus(F),Cputime(tt);
data:=eval Read("continuation_20261031a/C90/source45.m");
cov:=data[2][1];assert cov[1] eq 10 and cov[6] eq 108 and cov[7] eq 10 and #data[3] eq 10;
aa:=[a:a in [1..22]|GCD(a,45) eq 1];assert #aa eq EulerPhi(45) div 2;
// Frobenius-invariant key of a closed point, as in computations/d31a_90_source.m
function Key(v,D) keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(z,GF(2))]:z in v]);v:=[z^2:z in v];end for;return Sort(keys)[1];end function;
allkeys:={};ells:=[];
for i in [1..#data[3]] do
 d,ky:=Explode(data[3][i]);assert d eq 10;pol:=[r[2]:r in data[1]|r[1] eq d][1];K<z>:=ext<k|PP!pol>;
 b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);Q:=E![0,0];assert Order(Q) eq 45;
 orb:={Key(TateNormalForm_bc(E,a*Q),d):a in aa};assert #(orb meet allkeys) eq 0;allkeys join:=orb;
 vals:=TateNormalForm_xy(E,Q);G:=C48PairAFPlaceLinear(coords,vals);assert Degree(G) eq d;
 l:=Dimension(RiemannRochSpace(Divisor(G)));Append(~ells,l);
 printf "C90_ELL45 orbit=%o orbit_length=%o deg(G)=%o ell(G)=%o time=%o\n",i,#orb,Degree(G),l,Cputime(tt);
end for;
assert #allkeys eq 108;
printf "C90_ELL45_COMPLETE orbits=%o points=%o all_ell_one=%o time=%o\n",#ells,#allkeys,&and[l eq 1:l in ells],Cputime(tt);
quit;
