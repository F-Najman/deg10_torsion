This packet determines the degree-five points on X_1(3,18) over Q(zeta_3) and gives a new occurrence of C3 × C18 over a degree-ten field. There are exactly two marking orbits under the order-27 modular automorphism group. The represented elliptic curves are related by a 2-isogeny, and both have full torsion C3 × C18. The classification also excludes C3 × C36 and C6 × C18 in degree ten.

The proof is in `resolution.tex`. The displayed relative quintic field and curve are checked exactly in `endomorphism_gluing/global_point13900_absolute.m`; the absolute degree-ten polynomial and all marked coordinates are stored there. `global_RR_char0_reductionfalse.m` contains both exact global Riemann–Roch sections.

The proof uses rank zero over Q(zeta_3), complete finite-field divisor enumerations, a two-prime comparison of three-multiple classes, relative formal fibres, and a final unmultiplied cuspidal-class comparison. It does not assume that rational cuspidal divisors generate the whole Mordell–Weil group.

The accepted proof data, scripts, logs and independent controls are pinned in `certificate_manifest.json`. Failed, superseded and timed-out investigations are retained separately as diagnostics. `lifecycle.json` records actual exits; no paused or obsolete computation remains.
