Phi(10): certified partial results, 31 October 2026

STATUS: INCOMPLETE. This directory does not contain a completed classification.

68 exact full torsion groups are proved to occur:
  C_n: n=1,...,40,42,45,50;
  C_2 x C_(2n): n=1,...,15;
  C_3 x C_(3n): n=1,...,5;
  C_4 x C_(4n): n=1,...,3;
  C_6 x C_6, C_6 x C_12.

Of the 109 minimal subgroup obstructions outside this lower list, 99 are
excluded and 10 remain unresolved (6 cyclic, 4 noncyclic). These 10 are
minimal proof tasks, not asserted occurrences or an exhaustive upper list.
Any new occurrence may expose further minimal obstructions.

Read paper/main.pdf (88 pages) for the mathematical arguments and remaining
cases. candidates.txt lists the current exact counts. certified_results.json
is the evidence ledger. candidates.py independently enumerates the frontier.

Accepted inputs:
  - S(10)={2,3,5,7,11,13,17,19,23,29,31,37}, from the completed audited
    prime-order project snapshotted in inputs/. It is not re-proved here.
  - The complete earlier quintic classification and the quadratic theorem.
  - Degree-independent exact rank, winding, newform, and least-CM-degree
    evidence. Every degree-10 inequality is checked afresh.

Highlights:
  - Ten further exclusions: C44, C46, C58, C62, C63, C70, C78, C90, C117, C3 x C39.
    C78 also classifies all noncuspidal rational effective decic divisors
    on X1(39), excluding C3 x C39. C90 uses marked 2-isogeny cycles;
    the possible CM endomorphism fields are excluded by class numbers.
    C63 handles non-global cubic cusp fragments by exact Hecke tests.
    See continuation_20261031a/ for proofs, independent replays and
    immutable certificate manifests.
  - Six further exclusions: C65, C85, C91, C93, C99 and C105.
    C91 uses a relative formal-fibre argument on a genus23 quotient:
    the exact fibre consists of a fixed residual divisor plus pullbacks
    from X0(91), and the fixed part contains a cusp. All120 cusp and12
    mixed CM ranks were replayed independently. C99 needs a genus29
    quotient retaining the9-level to exclude its final Hecke class.
    C65 has7170048 independently checked cusp ranks. C85/C93 use the
    convention-correct full-Jacobian local operator T2mu-1, with genuine
    nonrational cusp arcs. C105 is resolved by six genus25 tests.
    See continuation_20261030c/C{65,85,91,93,99,105}/resolution.tex.
  - C55, C57, C75, C77, C81 and C87 are excluded. Every final test has
    independent exact checks, including the physical isogeny markings and
    residue norm weights. The C57 and C77 proofs use finite Hecke images
    despite positive full-Jacobian rank. C81 uses unique Hecke divisors on
    a genus46 Kummer cover and an impossible finite 2-isogeny-closed j-set.
    The source cusp boundary is covered separately in every case.
    See continuation_20261030b/C{55,57,75,77,81,87}/resolution.tex.
  - C51, C69, C2 x C42 and C84 are now excluded. The source data,
    quotient maps, Hecke norms and final Riemann-Roch tests were replayed
    independently. C51 uses full-Jacobian dyadic formal immersion on all
    9656144 cusp divisors. C242 ends with a unique-section argument on
    its diagonal genus17 quotient (24 patterns, 72 character spaces).
    C69 retains the 23-level generator on a genus67 quotient; its exact
    Hecke divisor is nonprincipal. See next_cases_20261030/{cyclic_further/C51,
    C242/geometry,C69}/resolution.tex.
  - C49 is excluded. The finite Hecke sieve and 2900895 independent
    Cartier restriction checks cover every decic reduction. A complete
    local Q2 formal-group lemma uses the full Jacobian as target.
    See next_cases_20261030/cyclic_primepowers/resolution.tex.
  - C95 and C111 are excluded by the dyadic cuspidal diamond criterion:
    local emptiness and rational cusp lifts give a degree-40 divisor
    comparison below gonality, forcing an impossible order-three orbit.
    See next_cases_20261030/dyadic_95111.tex.
  - C2 x C36 and C72 are excluded. Three finite class sieves reduce
    135282 support patterns to one degree-eight point. An order-four
    quotient and 23 Riemann-Roch tests exclude all quadratic cusp
    complements; the full Frobenius-fixed cusp-unit lattice closes
    the boundary. See C236_resolution_20261030/geometry/resolution.tex.
  - C2 x C34 and C68 are excluded. A complete finite sieve and a full
    rational cusp-unit lattice computation prove gon_Q >= 11 and finish
    the boundary case. See C234_resolution_20261030/resolution.tex.
  - C2 x C32 and C64 are excluded. The involution quotient has genus 13;
    an integral index-two torsion bound and disjoint diamond-3 coset sets
    complete the proof. See C232_resolution_20261030/resolution.tex.
  - A new sporadic C6 x C12 occurrence, with exact full torsion proved by
    a reduction of order 72 at 61 and independently checked in Sage.
    The complete quintic analysis gives one full-level orbit on X1(6,12).
    occurrences/new612_absolute.json gives a compact absolute field and curve.
  - Fifteen modular-unit family certificates and their separate replay.
  - Explicit decic curves, independently checked in Sage.
  - 43 q=3 global exclusions, two q=2 exclusions, and two exclusions by
    descent of a 31-level point to a forbidden quintic field.
  - A new sporadic C3 x C15 example, and a complete proof that it cannot
    grow in degree 10. The Galois cusp calculation proves exact generation
    of J1(3,15)(Q(zeta3)); the integer sieve leaves exactly its 48 level classes.
  - Seven further exclusions: C135, C155, C3 x C21, C3 x C24,
    C4 x C16, C4 x C20, C4 x C24.
    The C155 proof combines quadratic diamond descent, reduction at 2,
    and the Q-gonality of X1(31)/<2>.
  - Six new exclusions: C2 x C38, C2 x C40, C2 x C50, and consequently
    C76, C80, C100. Complete source enumerations, Picard filters and
    Hecke tests are independently checked. See continuation_20261029b/.
  - C2 x C26 is completely resolved: both residual classes have unique
    effective degree-ten representatives, giving two explicit curves with
    full torsion C2 x C26. The independent exact replay identifies their
    classes and excludes cyclic C52 by certified nonsquare root differences.
    See C226_resolution_20261030/.
  - Full 11-torsion is excluded over the forced field Q(zeta_11), using the
    complete 25-point set of X1(11) there and good reduction at 23.

All new work belongs here. Earlier-degree files are protected by the initial
hash manifest. Source scripts, failed diagnostics, logs and actual exit
statuses are retained; a failed or terminated job is not a proof certificate.
Batch coordinators run only when invoked explicitly and have bounded workers.
Keep two host CPU cores free, use one thread per Magma worker, and terminate
obsolete jobs instead of pausing them.

Read-only audit: /usr/bin/python3 work/deg10/verify_project.py
Separate algebraic replays are in computations/replay*.py; their completed
logs are in jobs/. The note states the remaining gaps explicitly.

Final audit: 8431 evidence hashes checked; all 1034 recorded jobs
exited; 25 protected earlier-degree files unchanged. The 88-page paper
compiles without errors or overfull boxes. No degree-ten worker, supervisor
or launcher remains active. See validation.json and lifecycle_final.json.
