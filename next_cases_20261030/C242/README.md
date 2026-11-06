# C2 × C42 / C84: complete decic source census and first quotient sieve

Completed first stage, 2026-10-30: **60971 normalized supports → 4647 survivors**. This is a necessary sieve, not a final exclusion. The subsequent genus5 and X1(21) work is in `geometry/`.

## Source and pure-cusp branch

The source X1(2,42) has genus61, projective index1152, and rational Jacobian rank0. The geometric gonality bound 325·1152/32768>11 gives gonality at least12. Hence a hypothetical rational noncuspidal effective divisor of degree10 cannot reduce entirely to cusps: lift the cusps over an unramified extension at5, apply rank0 and Manin–Drinfeld to the difference, and use injectivity of torsion reduction. This would produce distinct equivalent effective divisors of degree10 below geometric gonality.

The fresh Frobenius computation `next242_frobenius.m` uses JH(84,[43]) at5. Its degree122 polynomial exactly matches the accepted independent polynomial, and gives noncuspidal closed-place counts in degrees1 through10

`0,0,0,126,648,2080,11376,48816,211776,983196`.

The cusp signature at5 is `1^18 2^9 3^6 6^3`. Its geometric degree is72. `frobenius.m` stores the polynomial, complete cusp signature and counts.

## Complete source census

`source4.m` through `source9.m` are exact degree-independent copies of accepted source data. Their original paths and hashes, and those of copied helpers, are pinned in `reuse_manifest.json`. No earlier-degree exclusion is used as a decic exclusion.

Four full degree10 scans partition all976887 Frobenius orbits of j-invariants by the orbit index modulo4. For each curve twist they check full rational2-torsion and rational3- and7-torsion, and construct every cyclic21 subgroup together with an ordered2-torsion basis. A stored point Q has exact order42 and P is independent of21Q. Keys are canonical under all36 rational level automorphisms and every Frobenius power.

The automorphisms are represented by `(P,Q) -> (T_i,T_j+a(22Q))`, where T_i,T_j are distinct nonzero2-torsion points, and a lies in {1,2,4,5,8,10}. Their source/Frobenius orbits are handled with their actual stabilizers; no division by36 is assumed.

`source10_s0.m` through `source10_s3.m` all completed. `next242_source_merge.m` checks disjoint canonical source-key sets, completion of every j scan, identical field polynomials, and exact agreement with every relevant Frobenius count. The merged `source10.m` has27392 records:18 of actual degree5 and27374 of actual degree10. Proper-degree5 records are retained with their correct labels and omitted from degree10 support lists.

Normalized source orbit counts in degrees4 through10 are

`5,18,64,316,1372,5888,27374`.

## Two genus4 quotient maps

The target Y=X1(2,14) has a fresh proven finite Picard presentation in `target14_proven_picard.m`. Its finite Jacobian has invariants [2,2,12,36]; all normalized finite cusp classes generate H=[2,2,6,18], of index4. Rational cusp generation gives the required rational-image containment. The normalization cusp may be replaced by any finite rational cusp because its difference lies in H.

The maps retain `(P,3Q)` and, respectively, quotient by `<14Q>` and retain the images of P,Q. Every projected place is multiplied by its exact residue degree ratio. All finite target cusp classes lie in H, so a source cusp complement does not affect these two tests.

The portable forward class map retains the actual projected prime ideal, its class, and an identity `D=E-R-div(h)` with both finite and infinite maximal orders checked. It never assigns an unproved class-group flag. The finite presentation is freshly proven; the direct Riemann–Roch forward map uses that presentation with exact ideal certificates.

Both classes are computed for every degree4,5,6 support. Normalizing the first summand and retaining every actual second place gives actual second-place counts126,648,2080. For singleton degrees7 through10 a failure of the first map permits skipping the second. Eight degree10 target-filter shards partition global source indices modulo8; the census itself used four shards.

## Complete support coverage

There are no noncuspidal places of degree1,2,3 at5. A noncuspidal part of degree at most10 has one summand of degree4,...,10, or two summands of degrees(4,4),(4,5),(4,6),(5,5), with repetitions. There are no triples. Normalizing the first summand by all36 automorphisms preserves completeness.

|Noncuspidal degrees|Cases|Survivors of both maps|
|---|---:|---:|
|4|5|3|
|5|18|1|
|6|64|14|
|7|316|17|
|8|1372|144|
|9|5888|350|
|10|27374|1832|
|4,4|630|210|
|4,5|3240|180|
|4,6|10400|1128|
|5,5|11664|768|

Case IDs are sequential in this table order. `first_survivors_unified.m` stores `[fields,cases]`, with each case `[caseID,[[degree,key],...]]`. `first_survivors_details.m` preserves summand indices, both class vectors and the pattern table. These exact IDs pass unchanged to later sieves.

## Independent verification and lifecycle

`independent_arithmetic_replay.py` uses plain Python integers. It checks Newton identities, pinned source copies, all four source10 shards and the merged source ordering, the full quotient-lattice index4, and every one of the60971 support decisions. It recovers precisely4647 survivor IDs and source keys; results and input hashes are in `independent_arithmetic_replay.json`.

The independent Sage replay `independent_projection_common.py` reconstructs every source elliptic curve and both explicit quotient formulas. It verifies exact source orders, target coordinate keys, every residue-degree norm, and each assigned class against the recorded target prime. For degrees4,5,6 it also independently reconstructs all36-automorphism source orbits and all actual second places. Results are in `independent_projection_low.json`, `independent_projection_D9.json` and the eight `independent_projection_D10s*.json` files.

Native ideal replays under the `next242_native_*` prefix reconstruct both maximal-order ideals and check every saved identity directly, without repeating Riemann–Roch or computing new class groups. All44819 projected images passed these native checks:44801 explicit Riemann–Roch ideal identities and18 direct factorbase cases. `first_stage_completion.json` records all batches and verifies that every owned source/filter/replay worker exited and its PID is absent. Completion records and logs remain under `../../jobs/`.

Every Magma and Sage worker is single-threaded, has an explicit lifetime bound, and is recorded by `work/deg10/run_job.py` with immutable script hash, PID, owner, start time and actual exit code. The source profile is explicitly incomplete and is not a census certificate. No paused jobs or automatic fallback launchers are retained.

The preparation-only `prepared_manifest.json` is a historical snapshot. The accepted current first-stage hashes are in `certificate_manifest.json`; no files from the independently progressing `geometry/` stage are included there.
