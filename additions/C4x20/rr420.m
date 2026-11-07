// Z/4 x Z/20 (proof of prop:4x20): Riemann-Roch checks for lem:section.
// X = X_1(4,20) over F_9 (reduction modulo the prime of Q(i) above 3).
// The Hecke sieve leaves 4 non-cuspidal parts B of the reduction of the trace divisor:
// the closed points of degree 4 with orbit indices 4 and 12, and the sums of two closed
// points of degree 2 given by the pairs [1,82] and [2,84] (indices as in
// continuation_20261029/M4N20_quartic_followup_arithmetic.txt and
// continuation_20261029/M4N20_pair_followup_arithmetic.txt, written by
// computations/cont_sieve_4_20_quartic_followup.m and computations/cont_sieve_4_20_pairs_followup.m).
// We check that l(B + c + c') = 1 for all 24 degree 1 places c, c' of X_{F_9} (all of them cusps),
// i.e. 300 unordered pairs (c = c' allowed) for each B, 1200 checks in total.
// Run from the repository root:  magma -b additions/C4x20/rr420.m
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("shared/mdmagma/v2/mdmagma.spec");
Attach("continuation_20261029/deps/fast_hecke.m");
Attach("continuation_20261029/deps/fast_hecke_11.m");
Attach("continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("continuation_20261029/deps/full_aut_4_20.m");
X:=MDX11(4,20,GF(9));
// the degree 1 places are all cusps
assert #NoncuspidalPlaces(X,1) eq 0;
xs:=NoncuspidalPlaces(X,2);rr,labels,keys:=Decic420Orbits(X,xs);
old:=eval Read("continuation_20261029/M4N20_pair_followup_arithmetic.txt");
assert keys eq old[1] and labels eq old[2] and [Index(xs,x):x in rr] eq old[3];
assert old[5] eq [[1,82],[2,84]];
xs4:=NoncuspidalPlaces(X,4);rs4,ll4,kk4:=Decic420Orbits(X,xs4);
old4:=eval Read("continuation_20261029/M4N20_quartic_followup_arithmetic.txt");
assert kk4 eq old4[2] and ll4 eq old4[3] and [Index(xs4,x):x in rs4] eq old4[4];
assert old4[6] eq [4,12];
print "SETUP",#xs,#rr,#xs4,#rs4,Cputime(tt);
CC:=Curve(xs[1]);P1:=Places(CC,1);
print "genus",Genus(CC),"degree 1 places",#P1;
assert Genus(CC) eq 25 and #P1 eq 24;
Bs:=[Divisor(rr[1])+Divisor(xs[82]),Divisor(rr[2])+Divisor(xs[84]),Divisor(rs4[4]),Divisor(rs4[12])];
names:=["[2,2] pair [1,82]","[2,2] pair [2,84]","[4] orbit 4","[4] orbit 12"];
total:=0;
for k in [1..4] do
  B:=Bs[k];assert Degree(B) eq 4;dims:={};n:=0;
  for i in [1..#P1] do for j in [i..#P1] do
    Include(~dims,Dimension(RiemannRochSpace(B+Divisor(P1[i])+Divisor(P1[j]))));n+:=1;
  end for;end for;
  print "RR420",k,names[k],"pairs",n,"dimensions",dims,Cputime(tt);
  assert n eq 300 and dims eq {1};total+:=n;
end for;
print "RR420_DONE all",total,"Riemann-Roch dimensions equal 1",Cputime(tt);
quit;
