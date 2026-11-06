Complete resolution of the two decic C2 x C26 divisor classes
30 October 2026

RESULTS
Both residual rational divisor classes have h0=1 over Q. Each unique
effective divisor is a single noncuspidal closed point of degree 10.
The two elliptic curves have full torsion C2 x C26 and distinct quintic
j-invariant orbits. These are the complete two Galois and level orbits
of C2 x C26 curves over decic fields.

The complete classification also excludes cyclic C52: for each of the
three rational 2-isogeny kernels on either representative, the dual
kernel point cannot be halved. All six relevant root differences are
certified nonsquares by exact finite residues.

The Phi(10) ledger now has 68 occurring groups, 64 exclusions and 45
unresolved minimal obstructions (36 cyclic, 9 noncyclic); frontier 109.
The classification of Phi(10) as a whole remains incomplete.

FILES
resolution.tex, coefficient_table.tex: mathematical proof, explicit models,
full torsion reductions, and the C52 argument; included in paper/main.pdf.
models_compact.json: two small field polynomials and polynomial numerators
and denominators for b,c in E:y^2=x^3+c*x^2+b*(1-b-c)*x, P=(0,0),Q=(b,b).
exact_points_30.m: first class, whose F3 reduction has type 5+5.
exact_points_24.m: second class, whose F3 reduction is one degree-10 point.
complete_resolution_inputs.m: complete small input for the native verifier.
complete_native_result.m: final native output for both cases.
independent_class_binding.json: precise Frobenius-class bindings, distinct
quintic j-polynomials and all twelve finite nonsquare certificates.
witness_replay_30.json, witness_replay_24.json: independent Sage full-torsion
checks, including good unramified reductions of order 52 at 41 and 43.
coordinate_fields.json: exact verification that Q(b,c) is the full decic field.
frontier_change.json: exact changes from the preceding 47-case frontier.

SHORT REPRODUCTION
From the project root, the complete native verification is
  magma -b work/deg10/computations/resolve226_complete_native.m
It verifies the field irreducibility, level point orders, full torsion bound,
Frobenius-class bindings, and nonsquare residues. It writes the deterministic
complete_native_result.m; it does not repeat the expensive discovery searches.
Its successful recorded run is jobs/resolve226_complete_native.log.

The independent Sage checks are recorded under jobs/resolve226_replay_witness*,
resolve226_bind_classes and resolve226_coordinate_fields. Their source and
reports are pinned in certified_results.json. Some discovery replays write
reports with elapsed times: use a fresh project copy if rerunning them, to
preserve the hashes of accepted proof evidence.

DISCOVERY AND JOB LIFECYCLE
The large function-field files and exploratory scripts are preserved for
provenance. Several prototype jobs timed out or failed after saving useful
intermediate objects. They are not accepted completed proof certificates.
The final proof instead uses the explicit curves and successful independent
native and Sage verification. The duplicate and final obsolete fallback
were terminated with verified PID/owner/start/command and actual exit -15.
See retired_duplicate.json and retired_obsolete_workers.json. The final
lifecycle audit records no remaining degree-ten workers or launchers.
