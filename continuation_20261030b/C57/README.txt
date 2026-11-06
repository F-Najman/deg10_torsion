C57 over degree-ten fields: complete exclusion certificate.

Proof: resolution.tex. This case does not assume rank zero of J1(57).
The accepted winding decomposition leaves the elliptic factor57a1; diamond
differences and full-Jacobian product operators handle that factor.

Fresh final certificates:
 - all_cusp_FI_sage_certificate.json: all2383311 rank10 cusp configurations;
 - noncuspidal_filter_v2_data.txt: genus9 containing set48, six images retained;
 - source_coverage_verified.json: all12 physical cyclic-subgroup orbits, six
   required orbits, independently matched to the DUAL j-filter via3+19isogenies;
 - good_jets_sage_certificate_v3.json:56 genuine R1/R2 arcs,48 geometric5kernels;
 - ring_RR_all_sage_certificate.json:8 exact nonprincipal point-model divisors,
   six required tests2,3,4,5,7,8 and two harmless extra tests;
 - lifecycle_final.json:23 tracked jobs all EXITED, no live C57 launchers/workers.

The six-case computations are completed intermediate checks, not the final
source coverage. The first independent jet verifier failed because it selected
F1024 as a common field for targets of degrees10 and50. Its script/log/exit1
are preserved and not used as a certificate. Verifier v2 fixed finite-field
alignment; final v3 additionally checks the two missing physical-source orbits.

All170 files of the inherited accepted degree9 packet were hash-verified
unchanged. Their degree-independent geometric, winding, cusp, and degree6/8
point-obstruction inputs are reused; no degree9 total-degree inequality is
transferred without checking degree10. The degree10 cusp enumeration and new
non-global9+1 cusp branch are explicitly recomputed. Physical and dual j-values
are never assumed equal. Necessary point-model source coefficients are[1,4,-5].

Run python3 work/deg10/continuation_20261030b/C57/certify57.py for read-only
hash and structural verification. Native scripts were executed only through
work/deg10/run_job.py with bounded lifetimes and one thread per worker.
