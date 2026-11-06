// Independent check of the class data used in the proof of lem:spor-318-classes (X_1(3,18), types [2,2], [2,3]).
//
// The files
//   continuation_20261031b/M3N18/inputs/M3N18_9multiple_fresh_p7_{group,degree2,degree3}_data.m
//   continuation_20261031b/M3N18/inputs/M3N18_9multiple_aut27_p13_{group,degree2,degree3}_data.m
// store, at the primes above 7 and 13, the class group J_X(F_p) x Z (as invariants and a relation matrix), the
// classes of the 24 cuspidal blocks, and the classes [x - d c0] of ALL the non-cuspidal closed points x of degree
// d = 2, 3 (27 + 144 at 7, 54 + 846 at 13), together with the coordinates of x.
//
// This script recomputes everything from the model of X_1(3,18) (mdmagma) with Magma's proven class group, using
// Magma's own divisor-to-class map (not the map of the original computation), identifies every stored point with a
// freshly enumerated place through its coordinates (up to Frobenius), and then checks
//   (a) whether the frame is reproduced exactly; if so, that every stored class equals the recomputed one;
//   (b) frame-independently: that the kernel of Z^N -> J(F_p), e_k -> (stored class k), equals the kernel of
//       Z^N -> J(F_p), e_k -> (recomputed class k), where k runs over the 24 cusp blocks and all the stored points.
//       Equivalently, some isomorphism of groups maps every stored class to the recomputed class of the same divisor.
//
// Run from the repository root:   magma -b independent_checks/recompute_318_classes.m
SetColumns(0); SetQuitOnError(true); SetMemoryLimit(16*10^9);
AttachSpec("shared/mdmagma/v2/mdmagma.spec");
inp := "continuation_20261031b/M3N18/inputs/";
allok := true;

function KeyInto(bc, K, p)
  // all Frobenius conjugates of the coordinate vector bc (elements of one finite field), as integer sequences over K
  L := Parent(bc[1]); d := Degree(K);
  if #L eq p then
    v := [K!(Integers()!z) : z in bc];
  else
    f := DefiningPolynomial(L, PrimeField(L));
    fk := PolynomialRing(K)![K!(Integers()!a) : a in Coefficients(f)];
    r := Roots(fk); assert #r ge 1;
    mp := hom<L -> K | r[1][1]>; v := [mp(z) : z in bc];
  end if;
  keys := {@ @};
  for i in [1..d] do
    Include(~keys, [[Integers()!c : c in Eltseq(z)] : z in v]); v := [z^p : z in v];
  end for;
  return keys;
end function;

for pc in [<7, "fresh_p7">, <13, "aut27_p13">] do
  p := pc[1]; tag := pc[2]; t0 := Cputime(); ok := true;
  gd := eval Read(inp cat "M3N18_9multiple_" cat tag cat "_group_data.m");
  assert gd[1] eq p;

  // model and cusp blocks, grouped by the valuations of the units of [N69, Section 12.4]
  X := MDX11(3, 18, GF(p)); C := Curve(X); FF := FunctionField(C); F0 := AlgorithmicFunctionField(FF);
  assert Genus(X) eq 10;
  E := EllipticCurve([FF!f : f in X`_E]); P := E![FF!f : f in X`_P]; Q := E![FF!f : f in X`_Q];
  ai := aInvariants(E); fs := [ai[3], ai[1]^3 - 27*ai[3]];
  for k in [1..9] do
    for j in [0..2] do A0 := k*Q + j*P; fs cat:= [A0[1], A0[2]]; end for;
    A0 := k*Q; B0 := A0 + P; Append(~fs, A0[1] - B0[1]);
  end for;
  cs := Cusps(X); groups := AssociativeArray();
  for c in cs do
    sig := [Valuation(f, c) : f in fs];
    if IsDefined(groups, sig) then Append(~groups[sig], c); else groups[sig] := [c]; end if;
  end for;
  sigs := Sort(Setseq(Keys(groups)));
  blocks := [&+[Divisor(c) : c in groups[s]] : s in sigs]; ds := [Degree(b) : b in blocks];
  same_sigs := sigs eq gd[2] and ds eq gd[3];
  printf "R318_CUSP_BLOCKS p=%o blocks=%o degrees=%o same_signatures_as_stored=%o\n", p, #sigs, ds, same_sigs;
  ok and:= same_sigs;
  c0 := Support(blocks[Index(ds, 1)])[1];
  D0 := FunctionFieldDivisor(Divisor(c0));

  G, _, phi := ClassGroup(F0 : ReductionDivisor := D0, Proof := true);
  printf "R318_CLASS_GROUP p=%o invariants=%o stored=%o time=%o\n", p, Invariants(G), gd[4], Cputime(t0);
  ok and:= Invariants(G) eq gd[4];
  cl := func<D | phi(FunctionFieldDivisor(D))>;
  newgs := [cl(blocks[i] - ds[i]*Divisor(c0)) : i in [1..#blocks]];

  // stored group and stored classes
  Fn := FreeAbelianGroup(#gd[4]); Gold, qold := quo<Fn | [Fn!r : r in gd[5]]>;
  oldelt := func<v | qold(Fn!v)>;
  oldlist := [oldelt(v) : v in gd[6]]; newlist := newgs; oldvec := gd[6];
  newvec := [Eltseq(g) : g in newgs];
  exact_frame := [Eltseq(r) : r in Rows(RelationMatrix(G))] eq gd[5] and newvec eq gd[6];

  for d in [2, 3] do
    dd := eval Read(inp cat "M3N18_9multiple_" cat tag cat Sprintf("_degree%o_data.m", d));
    assert dd[1] eq p and dd[2] eq d;
    K := GF(p^d); assert Coefficients(DefiningPolynomial(K)) eq dd[3];
    xs := NoncuspidalPlaces(X, d);
    lookup := AssociativeArray();
    for j in [1..#xs] do
      for key in KeyInto([Evaluate(f, xs[j]) : f in X`_coordinates], K, p) do
        assert not IsDefined(lookup, key); lookup[key] := j;
      end for;
    end for;
    match := [];
    for i in [1..#dd[4]] do
      key := [[Integers()!c : c in z] : z in dd[5][i]];
      if not IsDefined(lookup, key) then ok := false; printf "R318_UNMATCHED_POINT p=%o d=%o i=%o\n", p, d, i; continue; end if;
      Append(~match, lookup[key]);
    end for;
    bij := #match eq #dd[4] and #Set(match) eq #xs and #xs eq #dd[4];
    ok and:= bij;
    newcl := [cl(Divisor(xs[match[i]]) - d*Divisor(c0)) : i in [1..#match]];
    eq_entries := #[i : i in [1..#match] | Eltseq(newcl[i]) eq dd[4][i]];
    printf "R318_POINTS p=%o d=%o fresh_places=%o stored_points=%o bijection=%o entrywise_equal_classes=%o time=%o\n",
      p, d, #xs, #dd[4], bij, eq_entries, Cputime(t0);
    oldlist cat:= [oldelt(v) : v in dd[4]]; newlist cat:= newcl;
    oldvec cat:= dd[4]; newvec cat:= [Eltseq(g) : g in newcl];
  end for;

  // (a) exact frame
  printf "R318_EXACT_FRAME p=%o same_relation_matrix_and_cusp_classes=%o all_entries_equal=%o\n",
    p, exact_frame, exact_frame and newvec eq oldvec;

  // (b) frame-independent comparison
  N := #oldlist; Z := FreeAbelianGroup(N);
  hold := hom<Z -> Gold | oldlist>; hnew := hom<Z -> G | newlist>;
  Kold := Kernel(hold); Knew := Kernel(hnew);
  samekernel := Kold eq Knew;
  printf "R318_KERNEL p=%o N=%o image_old=%o image_new=%o same_kernel=%o time=%o\n",
    p, N, Invariants(Image(hold)), Invariants(Image(hnew)), samekernel, Cputime(t0);
  ok and:= samekernel;
  printf "R318_RESULT p=%o %o\n", p, ok select "PASS" else "FAIL";
  allok and:= ok;
end for;

printf "R318_RECOMPUTE_CLASSES %o\n", allok select "PASS" else "FAIL";
quit;
