/* Fresh exact nine-multiple cusp-diagonal inputs; no completion claim. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);Alarm(1700);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261031b/M3N18/fast_class_group_v4.m";
Attach("work/deg10/continuation_20261031b/M3N18/inputs/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261031b/M3N18/inputs/fast_full_aut_hash_m3_18_root.m");
for p in [13] do
 X:=MDX11(3,18,GF(p));C:=Curve(X);FF:=FunctionField(C);F0,afmap:=AlgorithmicFunctionField(FF);assert Genus(X) eq 10;
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..9] do
  for j in [0..2] do
   A0:=k*Q+j*P;assert A0 ne E!0 and A0[1] ne 0 and A0[2] ne 0;fs cat:=[A0[1],A0[2]];
  end for;
  A0:=k*Q;B0:=A0+P;assert A0 ne B0 and A0 ne -B0;Append(~fs,A0[1]-B0[1]);
 end for;
 cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^18,3^^6*};groups:=AssociativeArray();
 for c in cs do
  sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;blocks:=[&+[Divisor(c):c in groups[s]]:s in sigs];ds:=[Degree(b):b in blocks];assert {*d:d in ds*} eq {*1^^18,3^^6*};
 c0:=Support(blocks[Index(ds,1)])[1];
 print "M3N18_9MULT_CUSP_PROFILES",p,#sigs,ds,Cputime(t0);
 G,unused,phi:=ClassGroup(F0:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
 print "M3N18_9MULT_NATIVE_GROUP_PROVEN",p,Invariants(G),Cputime(t0);
 GG,gm,mg:=MDClassGroupDirectRRFastV4(F0,Sprintf("work/deg10/continuation_20261031b/M3N18/fast_benchmark_p%o_direct_RR_identities.txt",p));assert GG eq G;
 m2:=map<DivisorGroup(C)->G|D:->mg(FunctionFieldDivisor(D))>;
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..#blocks]];H:=sub<G|gs>;Qs:=quo<G|H>;
 print "M3N18_9MULT_CUSP_GROUP",p,Invariants(H),Invariants(Qs),Cputime(t0);
 J:=JH(54,[19]);fp:=FrobeniusPolynomial(J,p);assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!Evaluate(fp,1);
 rat:=[c:c in cs|Degree(c) eq 1];ratcl:=[m2(Divisor(c)-Divisor(c0)):c in rat];
 data:=[*p,sigs,ds,Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],[Eltseq(v):v in gs],Invariants(H),[Eltseq(v):v in ratcl],Coefficients(fp)*];
 Write(Sprintf("work/deg10/continuation_20261031b/M3N18/fast_benchmark_p%o_group_data.m",p),Sprint(data));


 old:=eval Read("work/deg10/continuation_20261031b/M3N18/paired_uniform_s0_partial.m");K:=GF(13^5);R0:=PolynomialRing(GF(13));vals:=[];for i in[1..Minimum(30,#old[1])] do qt:=[Evaluate(R0!v,K.1):v in old[2][i]];pl:=MDPlaceCommonZeros(X`_coordinates,qt);cl:=m2(Divisor(pl)-5*Divisor(c0));assert cl eq G!old[1][i];Append(~vals,Eltseq(cl));print "M318_FAST_MAP_BENCHMARK",i,Cputime(t0);end for;Write("work/deg10/continuation_20261031b/M3N18/fast_map_benchmark.m",Sprint(vals):Overwrite:=true);print "D31B_M318_FAST_MAP_BENCHMARK_COMPLETE",#vals,Cputime(t0);end for;quit;
