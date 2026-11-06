# Degree-10 exclusion of C81

`resolution.tex` proves that an elliptic curve over a degree-10 number field has no point of order 81. The final obstruction is a finite Riemann–Roch calculation on the genus-46 quotient `X_1(81)/<28>`. It replaces the insufficient global gonality inequality `gon(X_1(81)) > 30`.

The source curve has no noncuspidal closed points over F2 in degrees 1–8, and has 135 each in degrees 9 and 10. Each degree contributes five diamond orbits. A proven genus-4 `X_0(81)` Jacobian sieve reduces these ten supports to IDs 1, 2, 3, 7, 10. The genus-46 calculation gives `h0(E + 2<2>E) = 1` for all 29 remaining support/cusp combinations. The Hecke relation then forces a finite nonempty set of complex j-invariants closed under all 2-isogenies, a contradiction. Pure-cusp reduction is excluded separately using rank zero, zero rational 2-torsion, and intact global cusps.

All computations here concern degree 10. Earlier-degree evidence was copied into `inputs/` and is pinned in `reuse_manifest.json`; earlier files were not edited.

## Final evidence

- `source81.m`: complete finite source census and ten normalized Tate `(b,c)` keys.
- `parent_independent_source_filter.json`: independent Sage census of all curves, points and cyclic subgroups, manual 81-Vélu j-pair checks, exact norm multiplicities, and an exhaustive global/F2 cusp-orbit comparison.
- `x0_81_proven_picard.m`: proven `J_0(81)(F2) = C3 × C18` presentation.
- `x0_filter.m`: all ten normalized classes and final IDs `[1,2,3,7,10]`.
- `parent_native_x0_filter.m`: fresh native finite/infinite ideal-identity replay.
- `parent_two_torsion_replay.json`: independent binary replay of the exact 380-dimensional integral homology evidence, proving `J_1(81)(Q)[2] = 0`.
- `y46_kummer_model.m`: the genus-13 base and the exact genus-46 Kummer cover, including its identification with the actual dual-isogeny preimage cubic.
- `rr_case{1,2,3,7,10}_summary.m`: all 29 character-space dimensions, matrix sizes and ranks.
- `rr_caseID_cuspCI.m`: complete bases, jet matrices, local prime ideals, uniformizers, finite fields, and all three cubic roots with the physical root selected.
- `rr_caseID_localmeta.m`: explicit absolute field models for independent local-arc replay.
- `parent_native_completion_RR_v2.m`: independent native replay of all 87 base RR spaces and all 29 matrices using 14 genuine completions, with basis membership and independence checks.
- `independent_sage_local_arcs.json`: independent Hensel-series replay to precision 64, explicit source level-point and manual Vélu checks, and entrywise equality of every matrix.
- `certificate_manifest.json`: hashes and actual job outcomes. It distinguishes successful certificates from unsuccessful diagnostics.

For IDs 1, 2, 3, every one of the nine rational cusp complements is tested. The character dimensions are `(15,11,13)` and all 39 × 54 matrices have rank 38. For IDs 7 and 10 there is no complement; character dimensions are `(18,12,15)` and both 45 × 60 matrices have rank 44. All 29 kernels have dimension one.

## Serialization

All `.m` evidence files contain numeric Magma sequences/lists, with no executable statements. Python can parse them by replacing `[*` with `[`, `*]` with `]`, `true` with `True` and `false` with `False`, then using `ast.literal_eval`.

`source81.m = [fields, coverage, records]`. Field rows are `[d, F2_modulus_coefficients]`. Record rows are `[d, b_coefficients + c_coefficients]`, each coefficient block of length d. Source records 1–5 have degree 9; 6–10 have degree 10. They are canonical under Frobenius and all diamond operators. Each actual orbit has 27 closed places.

In the genus-46 snapshot, a coefficient in F4 is `[a,b]`, meaning `a+b*eta` with `eta^2+eta+1=0`. `ER` encodes a rational function in t as `[numerator_coefficients, denominator_coefficients]`. `EF` encodes an element of the degree-8 genus-13 base as a list of ER coefficients in `1,y,...,y^7`.

`y46_kummer_model.m` has ten fields:

1. Coefficients of the genus-13 base equation over F4(t), in ER encoding.
2. Original modular `(x,y)` functions, in EF encoding.
3. The Kummer function f, in EF encoding.
4. Coefficients of the irreducible dual-isogeny x-preimage cubic, in EF encoding.
5. Exponent 1 in the Cardano identity.
6. Function A, in EF encoding.
7. Function B, in EF encoding.
8. Function a in `x_dual = a + A*w + B*w^2`.
9. Divisor-of-f rows `[degree, valuation, finite_flag, ideal_basis_EF]`.
10. Universal Tate-27 elliptic a-invariants, in EF encoding.

For a prime divisor, `PD = [degree, finite_flag, ideal_basis_EF]`.

`rr_caseID_cuspCI.m` has three fields: character bases `[j, EF_functions]` for j=0,1,2; the complete matrix whose entries use the F4 two-bit encoding; and local blocks `[PD, pole_multiplicity, EF_uniformizer, [relative_modulus], roots, actual_root]`. The extra outer brackets on the relative modulus are intentional. Roots use relative-generator coefficient vectors over F4. Matrix columns run through local blocks, the two unselected roots in retained order, jet order 0 then 1 when needed, and the relative residue-field basis. Rows run through character bases in order 0,1,2.

`rr_caseID_localmeta.m = [caseID, d, source_key, delta_permutation, ramified_PDs, blocks]`. Each block is `[PD, pole_multiplicity, absolute_F2_modulus, ABS(eta), ABS(relative_generator), [ABS(t0),ABS(y0)], EF_uniformizer, ABS(all_roots), ABS(actual_root), ABS(actual_dual_x)]`. ABS means the coefficient vector in the explicit absolute generator given by the listed F2 modulus. It removes every finite-field embedding choice from the Sage replay. The native replay reconstructs its own residue-field maps directly from the physical source point.

## Reproduction and lifecycle

All executed scripts are immutable and retained under `work/deg10/computations/b81_*`. Their logs and actual exit records are under `work/deg10/jobs/b81_*`. Use a fresh script basename with the bounded runner rather than overwriting an executed script or status record:

```
python3 work/deg10/run_job.py FRESH_SCRIPT.m --seconds 300
python3 work/deg10/run_job.py FRESH_SCRIPT.py --engine sage --seconds 300
```

Primary successful model construction is `b81_y46_kummer_model_v2.m`. Primary RR jobs are `b81_rr_case{1,2,3,7,10}.m`. Independent successful native replay is `b81_parent_rr_completions_v2.m`; independent successful Sage replay is `b81_sage_local_arcs_v3.py`.

The absolute genus-46 normalization probes `b81_y46_model{,_v2}.m` reached their explicit memory limits and exited 1. Their completed Kummer replacement made further normalization unnecessary. Initial Kummer/native parser diagnostics and Sage API diagnostics are retained with their actual errors and are not certificates, including parser errors for which Magma returned exit 0. The accepted jobs require both exit 0 and their explicit completion marker. No C81 worker is paused or left running, and no automatic launcher is retained to restart obsolete work.
