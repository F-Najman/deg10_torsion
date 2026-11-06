C75 over degree-ten fields: complete exclusion certificate.

Proof: resolution.tex; validator: certify75.py.

New finite branch: all6 degree10 diamond orbits on the explicit genus73
quotient are excluded by T7-1-7delta7. Exact TwoGenerators model/divisor
packets are proof_packet_case1.m throughcase6.m, all in one shared model.
Parent rebuilt this model and replayed all6 h0=0 tests, checking36 primes.
Independent Sage verification checks84 actual moduli records, all48 geometric
7kernels and level images, norm factors, and full source Frobenius counts.

New cusp branch: all871630 admissible degree10 quotient images have rank10.
The old common-basis packet covers871520; direct fresh-kernel ranks cover
110 high-multiplicity cases. At the finite base cusp of width25, the ideal
pullback is Q25*Q75^3. A fractional-ideal sum and intersection isolates Q25;
its quotient gives3Q75, avoiding generic cusp decomposition. Exact ranks
are12 at12Q75,10 at10Q25, and15 at12Q75+3delta2^kQ75 for all k=1,...,9.
Surjectivity passes to each required effective subdivisor. All fresh kernels
use their own basis; no basis from another geometric export is combined.

Source cusp fields are not all intact: global octic cusps split into residual
quartics. Quotient cusp fields ARE intact. Tate-label arithmetic checks both
signatures, the residue-degree push factors, and the degree3 NC7 complement.
The native dyadic argument uses the full quotient Jacobian and its full
cotangent space, with the already proved local kernel lemma at2.

All14 tracked jobs (including the independent parent replay) EXITED.
First versions of the three NC10 jobs failed at an absent helper path;
corrected v2 jobs pin the accepted helpers locally. The first cusp verifier
failed only at JSON serialization of a Sage integer; its corrected v2 passes.
The first finite verifier failed at use of len() on a Sage polynomial;
corrected v2 passes. All failures and actual exit statuses remain preserved
and are not treated as certificates. All588 inherited accepted files remain
hash-identical. Accepted C57 and all earlier-degree files were untouched.
