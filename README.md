# deg10_torsion

This repository contains the code, the logs and the data of the computations for the paper
*Torsion of elliptic curves over number fields of degree 10* (F. Najman, 2026). The statement numbers in this
README are those of the version of November 2026 of that paper.

The paper proves (Theorem 1.1) that the torsion group of an elliptic curve over a number field of degree 10 is
one of the following 69 groups, and that each of them occurs:

    Z/n            n = 1, ..., 40, 42, 45, 50
    Z/2 x Z/2n     n = 1, ..., 15
    Z/3 x Z/3n     n = 1, ..., 6
    Z/4 x Z/4n     n = 1, 2, 3
    Z/6 x Z/6n     n = 1, 2

and (Theorem 1.2) determines all the elliptic curves over number fields of degree 10 whose torsion contains
Z/2 x Z/26, Z/3 x Z/15, Z/3 x Z/18 or Z/6 x Z/12. The proof uses that the primes which are orders of points on
elliptic curves over number fields of degree 10 are exactly the primes up to 37 (F. Najman, *Points of prime
order on elliptic curves over number fields of degree 9 to 22*, https://hexagonmath.org/2611.00120v1, code at
https://github.com/F-Najman/Sd), and the methods and some computations of the paper *Torsion of elliptic curves
over number fields of degree 6, 7, 8 and 9* (https://hexagonmath.org/2611.00073v2).

**Which script proves what.** The table in [Which files prove what](#which-files-prove-what) is the map from
each statement of the paper that rests on a computation to the folder, the scripts and the logs. In the digital
version of the paper, the words reporting a computation are links to the script that carries it out, at the
revision of this repository named in the paper.

## Origin and layout

The repository is a copy of the computational record of the proof, the folder `work/deg10` of the workspace in
which the computations were made between 29 and 31 October 2026 on the *Mordell* server (Department of
Mathematics, University of Zagreb; AMD EPYC 9175F, 384 GB RAM), with Magma V2.30-1, SageMath 10.7 and Python 3.11.
The files are the original ones, at their original relative paths; nothing in them was edited for this
repository. Consequently:

* Scripts were run from the root of the record, and paths in them are relative to it (this repository's root).
  Some job records (`jobs/*.status.json`) give the absolute paths of the original workspace,
  `/home/fnajman/Documents/Codex/2026-10-21/on-this-computer-there-are-files/work/deg10/`.
* Folder names record the order in which the cases were done (`continuation_20261031b/C56/` is the last stage
  of the work on Z/56). Where a case reuses computations of the earlier degrees, the files used were copied into
  the case folder (sub-folders `inherited_*`, `inputs/`), with their own manifests.
* Several folders contain more than the proof needs: exploratory runs, diagnostics and runs that failed or were
  stopped. Only the files listed below, and the files they read, are used in the paper. A failed or interrupted
  run is never used as evidence.

Top-level folders:

| Folder | Content |
|---|---|
| `occurrences/` | modular units of Proposition 3.3; data of the sporadic examples |
| `computations/` | scripts of the first stages (explicit examples, global method, modular units, many sieves) |
| `jobs/` | logs (`*.log`) and job records (`*.status.json`: command, script hash, exit code, time) of the scripts in `computations/` and of many others |
| `exclusions/` | first exclusions (global method, levels 115/121/135/155/341/403, full 11-torsion, Z/3 x Z/15, Z/4 x Z/16, Z/4 x Z/20, ...) with their data |
| `continuation_20261029/`, `continuation_20261029b/` | Z/3 x Z/21, Z/3 x Z/24, Z/4 x Z/24, Z/4 x Z/16, Z/4 x Z/20, Z/2 x Z/38, Z/2 x Z/40, Z/2 x Z/50, residual classes on X_1(2,26) |
| `C226_resolution_20261030/` ... `C236_resolution_20261030/` | Z/2 x Z/26 (and Z/52), Z/2 x Z/32, Z/2 x Z/34, Z/2 x Z/36 |
| `next_cases_20261030/` | Z/2 x Z/42, Z/49, Z/51, Z/69, Z/95, Z/111, the dyadic formal immersion criterion |
| `continuation_20261030b/`, `continuation_20261030c/` | Z/55, Z/57, Z/75, Z/77, Z/81, Z/87 and Z/65, Z/85, Z/91, Z/93, Z/99, Z/105 |
| `continuation_20261031a/` | Z/44, Z/46, Z/58, Z/62, Z/63, Z/70, Z/78, Z/90, Z/117, first stages of Z/48, Z/54, Z/56, Z/60, Z/66 |
| `continuation_20261031b/` | Z/48, Z/54, Z/56, Z/60, Z/66, Z/74, Z/3 x Z/18, Z/3 x Z/27, Z/3 x Z/33, Z/4 x Z/28 |
| `inputs/` | inputs from earlier work: van Hoeij's decic points (from the census https://x1mn-sporadic.github.io), LMFDB newform data, least CM degrees of Clark-Genao-Pollack-Saia, archives of the global method and of the octic Z/3 x Z/15 computation, the result S(10) |
| `shared/` | the packages `mdmagma` (Derickx-Sutherland) and `mdsage` (Derickx), at the revisions used, with manifests |
| `notes/` | intermediate status notes of the record |
| `independent_checks/` | checks made for this repository (see below); not part of the original record |

`certified_results.json` is the ledger of the record: for every group, the method and the list of evidence files
with their SHA-256 hashes. `candidates.py` enumerates the minimal groups that have to be ruled out and checks the
ledger, and `verify_project.py` checks the integrity of the whole record; both expect the complete record
(including the files described in the next section and files of the earlier degrees outside this folder), so use
`independent_checks/check_ledger_hashes.py` in this repository.

## Files not in the git repository

* **Large files.** The 243 files of at least 2 MB (1.56 GB in total, mostly exported matrices, presentations of
  Picard groups and checkpoints) are in the archive `deg10_torsion_large_files.tar.zst` attached to the release
  [large-files-2026-11](https://github.com/F-Najman/deg10_torsion/releases/tag/large-files-2026-11) of this
  repository (211 MiB, SHA-256 `932c2ddd33ad3436c1ca73544207f591b0917487362efe8362baa6651da28633`). They are not needed
  to read the proofs; about 30 of them are inputs of scripts linked in the paper. Unpack it at the repository root (`tar --zstd -xf deg10_torsion_large_files.tar.zst`). The list,
  with sizes and SHA-256 hashes, is `independent_checks/large_files_manifest.tsv`.
* **Copies of papers.** The record contained copies of published papers and preprints (LaTeX sources, PDFs and
  text extracts) which the scripts used to look up statements. They are not redistributed here; they are listed in
  `independent_checks/literature_not_distributed.txt`.
* The folders `database_submission*` of the record (submissions of the new points to the census of sporadic points)
  and the working copy of the first draft of the paper.

`independent_checks/check_ledger_hashes.py` checks every evidence file of the ledger that is present against its
recorded hash, and reports the missing ones as belonging to one of the two lists above. On this repository it
reports `verified=10396 mismatched=0 in_release_archive_not_unpacked=216 literature_not_distributed=75
unaccounted=0`.

## Which files prove what

The following table is generated from the LaTeX sources of the paper: for each lemma, proposition or corollary of
the paper whose statement or proof links to this repository, it lists the linked files (scripts, and in a few cases
data or logs). The log of a script `X.m` or `X.py` in `computations/` is `jobs/X.log`, with the job record
`jobs/X.status.json` (command, hash of the script, exit code, time); in the case folders the logs are next to the
scripts. The numbers are those of the version of the paper of November 2026.

| Statement | Label | Files linked in the statement and its proof |
|---|---|---|
| Proposition 3.3 | `prop:families` | `computations/unit_1_23.py`, `occurrences` |
| Proposition 3.4 | `prop:explicit_cyclic` | `computations/witness_C29.m`, `independent_checks/same_pair_check.m`, `independent_checks/vanhoeij_check.m` |
| Proposition 3.5 | `prop:2x28-2x30` | `computations/quintic_full2_28_v2.m`, `computations/quintic_full2_30a_v2.m`, `computations/quintic_full2_30b_v2.m` |
| Lemma 4.1 | `lem:minimal` | `candidates.py`, `independent_checks/minimal_groups.py` |
| Proposition 5.1 | `prop:global` | `computations/replay_new_certificates_v2.py`, `exclusions/global3.json`, `exclusions/global3_additional.json` |
| Lemma 5.2 | `lem:glob-finite` | `inputs/global_newforms_cache.json`, `computations/winding1147.m`, `computations/winding1147_resume.m`, `computations/winding175.m`, `computations/winding185.m`, `computations/winding259.m`, `computations/winding407.m`, `computations/winding481.m`, `computations/winding629.m`, `computations/winding703.m`, `computations/winding851.m`, `computations/winding1073.m`, `computations/winding1369.m` |
| Lemma 5.4 | `lem:glob-odd` | `inputs/global_newforms_cache.json`, `computations/cont_winding135.m`, `computations/image121_no2_v2.m`, `computations/cont_image135_a16_no2.m` |
| Proposition 5.6 | `prop:135` | `computations/cont_image135_a16_no2.m` |
| Proposition 5.8 | `prop:341-403` | `computations/descent31.m` |
| Lemma 5.9 | `lem:glob-Y31` | `computations/cont_C155_local2.m`, `computations/cont_winding31.m` |
| Proposition 5.12 | `prop:155` | `computations/cont_winding155.m` |
| Lemma 5.14 | `lem:glob-95-111` | `computations/next95_winding.m`, `computations/winding111.m` |
| Lemma 6.1 | `lem:mxn-X111` | `inputs/global_newforms_cache.json`, `computations/full11_cyclotomic_v3.m` |
| Proposition 6.2 | `prop:11x11` | `computations/full11_cyclotomic_v3.m` |
| Lemma 6.3 | `lem:mxn-321` | `continuation_20261029/replay_cusp_tate.py`, `computations/cont_local_3_21_q4.m`, `computations/cont_double_cusp_3_21_D5_dt.m` |
| Proposition 6.4 | `prop:3x21` | `continuation_20261029/degree_independent_321/sieve_3_21_D4_p2.m`, `computations/cont_sieve_3_21_D5.m` |
| Lemma 6.5 | `lem:mxn-324-424` | `continuation_20261029/replay_cusp_tate.py`, `computations/cont_local_3_24_q7.m`, `computations/cont_local_4_24_q5.m` |
| Proposition 6.6 | `prop:3x24-4x24` | `computations/cont_sieve_3_24_D5_lower_supports.m`, `computations/cont_sieve_3_24_D5_shard0.m`, `computations/cont_sieve_3_24_cubic31.m`, `computations/cont_sieve_4_24_D5_d4_r1.m`, `computations/cont_sieve_4_24_D5_d5_r1.m` |
| Lemma 6.8 | `lem:mxn-416` | `computations/cont_diagonal_4_16_D5_v2.m`, `computations/cont_local_4_16_q5.m` |
| Proposition 6.9 | `prop:4x16` | `computations/cont_diagonal_4_16_D5_v2.m`, `computations/cont_M4N16_D5_RR13.m` |
| Lemma 6.10 | `lem:mxn-420` | `continuation_20261029/replay_cusp_tate420.py`, `computations/cont_local_4_20_q9.m` |
| Lemma 6.11 | `lem:mxn-420gon` | `computations/cont_gonality_4_20_q9_d6_v2.m` |
| Proposition 6.13 | `prop:4x20` | `computations/cont_sieve_4_20_lower.m`, `computations/cont_sieve_4_20_quartic_followup.m`, `computations/cont_sieve_4_20_quintic_shard0.m`, `computations/cont_sieve_4_20_mixed_shard0.m`, `computations/cont_sieve_4_20_pairs_followup.m` |
| Lemma 6.14 | `lem:mxn-327` | `computations/d31b_m327_torsion_gcd.m`, `computations/d31b_m327_cusp_parity_v2.py`, `computations/d31b_m3n27_frob_preflight.m`, `computations/d31b_m3n27_physical_source.m` |
| Proposition 6.15 | `prop:3x27` | `computations/d31b_m3n27_quotient_A7_v5.m` |
| Lemma 6.16 | `lem:mxn-333image` | `computations/d31b_m333_antih10_winding.m` |
| Lemma 6.19 | `lem:mxn-333` | `independent_checks/cusp_degrees.py`, `computations/d31b_m3n33_frob_preflight.m`, `computations/d31b_m333_kummer_frame.m`, `computations/d31b_m333_analytic_jets_v4.m`, `computations/d31b_m333_prefix_ranks.py`, `computations/d31b_m3n33_physical_source.m` |
| Proposition 6.20 | `prop:3x33` | `computations/d31b_m333_rational_trace_A5.m`, `computations/d31b_m333_rational_trace_3isog_A5.m` |
| Lemma 6.22 | `lem:mxn-428` | `computations/d31b_m428_Y_fourmap.m`, `computations/d31b_Y228_AF_cusp_signature_probe.m`, `independent_checks/cusp_degrees.py` |
| Proposition 6.23 | `prop:4x28` | `computations/d31b_m428_Y_fourmap.m`, `computations/d31b_Y228_A5_v5_s0.m`, `computations/d31b_Y228_A5_common_v5.m`, `computations/d31b_m428_neighbors_A5.m` |
| Lemma 7.3 | `lem:2x2n-32gon` | `computations/next232_gonality_boundary.py` |
| Lemma 7.4 | `lem:2x2n-32Y` | `computations/next232_cusp_signatures.m`, `computations/next232_final_replay_v3.py`, `computations/next232_quotient_cusps.py`, `computations/next232_quotient13_integral_certificate.py`, `computations/next232_quotient13_v2.m`, `computations/next232_quotient13_filter.m` |
| Lemma 7.5 | `lem:2x2n-32sieve` | `computations/next232_source_common.m`, `computations/next232_target_setup_v2.m`, `computations/next232_filter_common.m`, `computations/next232_final_replay_v3.py` |
| Lemma 7.6 | `lem:2x2n-32cosets` | `computations/next232_final_replay_v3.py`, `computations/next232_quotient13_orbit_filter.m` |
| Proposition 7.7 | `prop:2x32` | `computations/new234_cusp_lattice.py` |
| Lemma 7.8 | `lem:2x2n-34sieve` | `computations/new234_source_common.m`, `computations/new234_target_setup_v2.m`, `computations/new234_filter_common.m` |
| Lemma 7.9 | `lem:2x2n-34units` | `computations/new234_cusp_lattice.py`, `computations/new234_short_cusp_units.m` |
| Lemma 7.14 | `lem:2x2n-36sieve` | `computations/new236_source_common.m`, `computations/new236_filter_target18_setup.m`, `computations/new236_filter18_common_v2.m`, `C236_resolution_20261030/geometry/genus9_stage2_common_v2.m`, `computations/new236_geom_cyclic36_setup.m`, `C236_resolution_20261030/geometry/cyclic36_stage3_common.m` |
| Lemma 7.15 | `lem:2x2n-36x` | `computations/new236_geom_cusp2_numeric.m` |
| Lemma 7.17 | `lem:2x2n-36fibres` | `computations/new236_cusp_lattice.py`, `computations/new236_cusp_short.m` |
| Lemma 7.20 | `lem:2x2n-38` | `computations/late238_source.m`, `computations/late238_filter19_v2.m`, `computations/late_cusp_tate.py`, `continuation_20261029b/C238/hecke_template.m` |
| Lemma 7.21 | `lem:2x2n-40` | `continuation_20261029b/C240/enumerate_direct.m`, `continuation_20261029b/C240/filter_quotient20_linear.m`, `continuation_20261029b/C240/filter_quotient20_isogeny.m`, `continuation_20261029b/C240/filter_cyclic40_portable.m` |
| Lemma 7.22 | `lem:2x2n-50` | `computations/late250_source.m`, `computations/late250_filter25.m`, `computations/late_cusp_tate.py`, `continuation_20261029b/C250/hecke_template.m` |
| Lemma 7.26 | `lem:2x2n-42X0` | `computations/next242_geom_x0_42_setup.m` |
| Lemma 7.27 | `lem:2x2n-42sieve` | `computations/next242_source_common.m`, `computations/next242_target14_setup.m`, `computations/next242_filter14_common.m`, `next_cases_20261030/C242/geometry/quotient5_filter_common_v2.m`, `next_cases_20261030/C242/geometry/cyclic21_filter_common.m`, `computations/next242_cheap_stabilizers.py`, `next_cases_20261030/C242/geometry/x0_42_filter_common.m`, `computations/next242_direct_hecke_common_v2.m` |
| Lemma 7.28 | `lem:2x2n-42cubic` | `next_cases_20261030/C242/geometry/x0_42_filter_common.m` |
| Lemma 7.29 | `lem:2x2n-42W` | `computations/next242_case26_phases.m`, `next_cases_20261030/C242/geometry/diagonal17_normalized_model_data.m`, `next_cases_20261030/C242/geometry/diagonal17_parameter_data.m`, `computations/next242_geom_case26_sections.m` |
| Lemma 8.1 | `lem:sv-44` | `continuation_20261031a/C44/new_finite_jets.m`, `continuation_20261031a/C44/new_matrix_sieve.m`, `continuation_20261031a/C44/q7_matrix_sieve.m` |
| Proposition 8.2 | `prop:44` | `continuation_20261031a/C44/cusp_audit.py` |
| Lemma 8.3 | `lem:sv-46` | `continuation_20261031a/C46/new_support_filter.m`, `computations/d31a_46_parent_exceptional_hecke_common.m` |
| Proposition 8.4 | `prop:46` | `continuation_20261031a/C46/cusp_audit.py`, `continuation_20261031b/C48/d31b_C48_low_class_common.m` |
| Lemma 8.5 | `lem:sv-48-cusps` | `continuation_20261031b/C48/d31b_C48_independent_cusp_labels.py`, `continuation_20261031b/C48/d31b_C48_cusp_corrections_v2.m` |
| Lemma 8.6 | `lem:sv-48-sieve` | `continuation_20261031b/C48/inputs/four_characters.m`, `continuation_20261031b/C48/d31b_C48_character_high_common_v3.m`, `continuation_20261031b/C48/d31b_C48_high_class_common.m`, `continuation_20261031b/C48/d31b_C48_low_class_common.m`, `continuation_20261031b/C48/d31b_C48_high_hecke_common.m`, `continuation_20261031b/C48/d31b_C48_low_hecke_common.m` |
| Lemma 8.8 | `lem:sv-54-gon` | `computations/d31a_54_parent_gonality_controls.py` |
| Lemma 8.10 | `lem:sv-54-Z` | `continuation_20261031b/C54/d31b_C54_dyadic_rigidity_v2.m` |
| Lemma 8.13 | `lem:sv-54-dyadic` | `continuation_20261031b/C54/d31b_C54_dyadic_rigidity_v2.m` |
| Lemma 8.14 | `lem:sv-54-five` | `continuation_20261031b/C54/d31b_C54_transport_800_resume.m`, `continuation_20261031b/C54/d31b_C54_enumerate_all_class_members.m`, `continuation_20261031b/C54/d31b_C54_member_two_isogenies_all.py`, `continuation_20261031b/C54/d31b_C54_member_all_classes.m` |
| Proposition 8.15 | `prop:54` | `computations/d31a_C56_target28_setup.m` |
| Lemma 8.16 | `lem:sv-56-sieve` | `continuation_20261031a/C56/fourmap_filter.m`, `computations/d31a_C56_fourmap_filter_v4.m`, `continuation_20261031b/Y228_shared/union_cases.m`, `computations/d31b_Y228_A5_common_v5.m`, `continuation_20261031b/C56/d31b_C56_final_source_subset.py` |
| Lemma 8.17 | `lem:sv-56-known` | `computations/d31b_Y228_known_decic_reduction_v2.m`, `computations/d31b_Y228_known_decic_nonsquare.py` |
| Lemma 8.18 | `lem:sv-56-cover` | `computations/d31b_C56_cover_hecke_common_v4.m` |
| Lemma 8.19 | `lem:sv-56-rr` | `computations/d31b_C56_cover_hecke_common_v4.m`, `computations/d31b_C56_kummer_RR_direct_v4.m` |
| Proposition 8.20 | `prop:56` | `continuation_20261031b/C56/d31b_C56_cusp_A11_audit.py` |
| Lemma 8.21 | `lem:sw-5862-setup` | `continuation_20261031a/C58/inputs/C58_X129/build_presentation.m` |
| Lemma 8.22 | `lem:sw-58-classes` | `continuation_20261031a/C58/inputs/C58_X129/classify58.m`, `continuation_20261031a/C58/degree10_filter.m`, `jobs/d31a_C58_class_filter.log`, `continuation_20261031a/C58/known_lift_binding_job_v3.m` |
| Lemma 8.24 | `lem:sw-62-classes` | `continuation_20261031a/C62/inputs/C62_X131/classify62.m`, `continuation_20261031a/C62/inputs/C62_X131/known31/verify_known31.py`, `continuation_20261031a/C62/known31_global.m`, `continuation_20261031a/C62/degree10_filter_resume_v2.m`, `jobs/d31a_C62_class_resume_v2.log`, `continuation_20261031a/C62/small_support_arithmetic.py` |
| Lemma 8.26 | `lem:sw-60-setup` | `continuation_20261031b/C60/d31b_C60_degree10_series.m`, `continuation_20261031b/C60/d31b_C60_target30_setup.m`, `continuation_20261031b/C60/d31b_C60_merge_three_maps.py` |
| Lemma 8.27 | `lem:sw-60-series` | `continuation_20261031b/C60/d31b_C60_degree10_series.m` |
| Lemma 8.28 | `lem:sw-60-quintic` | `continuation_20261031b/C60/d31b_C60_quintic_2tower.m`, `continuation_20261031b/C60/inputs/quintic_full2_30a_v2.m`, `continuation_20261031b/C60/inputs/quintic_full2_30b_v2.m` |
| Lemma 8.29 | `lem:sw-60-mixed` | `continuation_20261031b/C60/d31b_C60_mixed_fibres_v2.m`, `continuation_20261031b/C60/d31b_C60_three_maps_common_v2.m`, `continuation_20261031b/C60/d31b_C60_merge_three_maps.py`, `continuation_20261031b/C60/d31b_C60_final_halving.m` |
| Lemma 8.32 | `lem:sw-66-F2` | `computations/d31b_C66_F2_frame_v2.m`, `computations/d31b_C66_F2_degree10_classes.py` |
| Lemma 8.33 | `lem:sw-66-F5` | `computations/d31b_C66_F5_extra_RR.m` |
| Lemma 8.34 | `lem:sw-66-pencils` | `computations/d31b_C66_prepare_gonal_cusps.py`, `computations/d31b_C66_two_gonal_orbits_filter.m`, `computations/d31b_C66_gonal_genuine_replay_v3.py` |
| Lemma 8.37 | `lem:sw-70-test` | `independent_checks/cusp_degrees.py`, `computations/d31a_even_p3_census_common.m`, `computations/d31a_C70_A2_common.m` |
| Lemma 8.39 | `lem:sw-74-J` | `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/deck/rational_two_diamond.m`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/rational_two_delta4.log`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/prime_good_branches_64GB.m`, `computations/d31b_C74_formal_export_v4.m`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/good_branches_64GB.log`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/sextic_lifts_control_v2.m`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/finite_nc9_survivors.m` |
| Lemma 8.41 | `lem:sw-74-ranks` | `computations/d31b_C74_formal_export_v4.m`, `computations/d31b_C74_prefix_rank_driver.py`, `computations/d31b_C74_prefix_ranks.cpp` |
| Lemma 9.1 | `lem:fi-49-points` | `next_cases_20261030/cyclic_primepowers/pp49_finite_sieve.m` |
| Lemma 9.2 | `lem:fi-49-cusps` | `next_cases_20261030/cyclic_primepowers/inputs/degree9/cartier_mixed9_export.m`, `next_cases_20261030/cyclic_primepowers/pp49_replay_old_jets_degree10.py`, `independent_checks/rerun_C49/pp49_fast_cartier_export.m`, `next_cases_20261030/cyclic_primepowers/inputs/degree9/cartier_transport_export.m`, `next_cases_20261030/cyclic_primepowers/pp49_native_matrix_replay.py` |
| Lemma 9.4 | `lem:fi-51-points` | `next_cases_20261030/cyclic_further/C51/further51_finite_sieve.m` |
| Lemma 9.5 | `lem:fi-51-cusps` | `next_cases_20261030/cyclic_further/C51/inputs/degree9/generic_cartier_completion_export.m`, `next_cases_20261030/cyclic_further/C51/further51_replay_old_jets.py`, `next_cases_20261030/cyclic_further/C51/further51_cusp_orbits.py`, `next_cases_20261030/cyclic_further/C51/further51_new_jets.m` |
| Lemma 9.7 | `lem:fi-55-inj` | `continuation_20261030b/C55/inputs/N55_quotient_injectivity_D9_v2.m` |
| Lemma 9.8 | `lem:fi-55-points` | `computations/contb55_census.m`, `computations/contb55_source_case1.m`, `computations/contb55_source_case2.m`, `computations/contb55_source_case3.m`, `computations/contb55_source_case4.m`, `computations/contb55_source_case5.m` |
| Proposition 9.9 | `prop:55` | `continuation_20261030b/C55/cusp_audit.py` |
| Lemma 9.10 | `lem:fi-75-points` | `continuation_20261030b/C75/cont75_census10.m`, `continuation_20261030b/C75/cont75_nc10_job1_v2.m`, `continuation_20261030b/C75/cont75_nc10_job2_v2.m`, `continuation_20261030b/C75/cont75_nc10_job3_v2.m` |
| Proposition 9.11 | `prop:75` | `continuation_20261030b/C75/cont75_verify_cusps_v2.py` |
| Lemma 9.12 | `lem:fi-77-Y` | `continuation_20261030b/C77/cusp_audit.py` |
| Lemma 9.13 | `lem:fi-77-cusps` | `computations/contb77_fi10.m` |
| Lemma 9.14 | `lem:fi-77-points` | `computations/contb77_census.m`, `computations/contb77_nc_case1_v3.m`, `computations/contb77_nc_case2_v3.m`, `computations/contb77_nc_case3_v3.m` |
| Lemma 9.16 | `lem:fi-85-MW` | `continuation_20261030c/C85/inherited_degree9/prime_orders/continued_20261026/composite85/rational_cusp_FI_positive.m` |
| Lemma 9.17 | `lem:fi-85-cusps` | `continuation_20261030c/C85/contc85_cusp_blocks.py`, `continuation_20261030c/C85/contc85_cusp_FI.py`, `continuation_20261030c/C85/contc85_arc_q4.m`, `continuation_20261030c/C85/contc85_arc_q8width5.m`, `continuation_20261030c/C85/contc85_arc_q8width1.m` |
| Lemma 9.18 | `lem:fi-85-points` | `computations/contc_parent_source_census.py`, `continuation_20261030c/C85/contc85_jets10.m`, `continuation_20261030c/C85/inherited_degree9/prime_orders/continued_20261026/composite85/verify_stage1_readonly.py`, `continuation_20261030c/C85/contc85_R1_sieve.py` |
| Lemma 9.21 | `lem:fi-91-Y` | `continuation_20261030c/C91/tate7_quotient_common.m`, `computations/contc91_tate7_cusps.m`, `continuation_20261030c/C91/tate7_cusp_audit_v2.py`, `computations/contc91_export_native_v2.m`, `computations/contc_parent_source_census.py`, `continuation_20261030c/C91_CM_control/CM_level7_binding.m` |
| Lemma 9.22 | `lem:fi-91-ranks` | `computations/contc91_tate7_cusps.m`, `computations/contc91_tate7_CM.m` |
| Lemma 9.25 | `lem:fi-93-cusps` | `continuation_20261030c/C93/contc93_cusp_descent.py`, `continuation_20261030c/C93/contc93_cusp_blocks.py`, `continuation_20261030c/C93/contc93_cusp_FI.py`, `continuation_20261030c/C93/contc93_arc_q5.m` |
| Proposition 9.26 | `prop:93` | `computations/contc_parent_source_census.py` |
| Lemma 10.2 | `lem:fa-57-Y` | `continuation_20261030b/C57/cont57_quotient10.m` |
| Lemma 10.3 | `lem:fa-57-filter` | `continuation_20261030b/C57/cont57_census10.m`, `continuation_20261030b/C57/cont57_quotient10.m` |
| Lemma 10.4 | `lem:fa-57-B` | `continuation_20261030b/C57/cont57_dual_match.m`, `continuation_20261030b/C57/cont57_good10_1.m`, `continuation_20261030b/C57/cont57_RR_v2.m` |
| Lemma 10.5 | `lem:fa-57-FI` | `continuation_20261030b/C57/cont57_cusp_FI.py` |
| Lemma 10.8 | `lem:fa-65-B` | `computations/contc65_source_census_v2.m`, `computations/contc65_ordinary10_9.m`, `computations/contc65_R1_sieve.py` |
| Lemma 10.9 | `lem:fa-65-FI` | `computations/contc65_cusp_FI_n1.py`, `computations/contc65_cusp_FI_common.py` |
| Lemma 10.11 | `lem:fa-63-Y` | `continuation_20261031a/C63/inputs/C63/cusp63.py`, `computations/d31a_63_good10.m` |
| Lemma 10.12 | `lem:fa-63-A5` | `computations/d31a_63_good10.m`, `computations/d31a_63_good7.m` |
| Lemma 10.13 | `lem:fa-63-cusp` | `computations/d31a_63_high_1.m`, `computations/d31a_63_existing_boundary.py` |
| Lemma 10.16 | `lem:fa-81-points` | `computations/b81_source_census.m` |
| Lemma 10.17 | `lem:fa-81-X0` | `computations/b81_x0_setup.m`, `computations/b81_x0_filter.m` |
| Lemma 10.18 | `lem:fa-81-model` | `computations/b81_y46_kummer_model_v2.m` |
| Lemma 10.19 | `lem:fa-81-RR` | `computations/b81_rr_case1.m` |
| Lemma 10.21 | `lem:fa-87-points` | `computations/b87_base_test_v3.m` |
| Lemma 10.22 | `lem:fa-87-RR` | `computations/b87_base_test_v3.m` |
| Lemma 10.24 | `lem:fb-W69` | `next_cases_20261030/C69/quotient67/q67_profile_rr_v2.m` |
| Lemma 10.25 | `lem:fb-69red` | `computations/cont69_cover_all_decic_v2.m` |
| Proposition 10.26 | `prop:69` | `computations/cont69_cover_all_decic_v2.m`, `next_cases_20261030/C69/direct_full/direct69_select_input_v2.m`, `next_cases_20261030/C69/quotient67/q67_profile_rr_v2.m` |
| Lemma 10.27 | `lem:fb-x39` | `independent_checks/same_pair_check.m`, `computations/witness_C39.m`, `computations/d31a_78_independent_genuine.py` |
| Lemma 10.28 | `lem:fb-39` | `computations/d31a_78_source.m`, `computations/d31a_78_A5.m`, `computations/d31a_78_independent_genuine.py` |
| Lemma 10.32 | `lem:fb-45` | `computations/d31a_90_source.m`, `computations/d31a_90_quartic_boundary_v2.py` |
| Lemma 10.33 | `lem:fb-45red` | `computations/d31a_90_parent_CM_classnumbers.py` |
| Lemma 10.36 | `lem:fb-99` | `computations/cont_winding99.m`, `computations/contc99_source_census.m` |
| Lemma 10.37 | `lem:fb-V99` | `computations/contc99_independent_Y29_model_v5.m` |
| Proposition 10.38 | `prop:99` | `computations/contc99_quotient33_RR.m`, `computations/contc99_independent_Y29_RR.m` |
| Lemma 10.39 | `lem:fb-105` | `computations/cont_winding105.m`, `computations/contc105_source_census.m` |
| Proposition 10.40 | `prop:105` | `computations/contc105_quotient35_RR.m` |
| Lemma 11.1 | `lem:spor-612-jac` | `computations/cont_local_6_12_q7.m`, `computations/cont_diagonal_6_12_D5_v5.m`, `computations/cont_hecke_torsion_6_12_v2.m` |
| Lemma 11.2 | `lem:spor-612-double` | `jobs/cont_reduce612_quintic_lattice_v2.log`, `computations/cont_diagonal_6_12_D5_v5.m`, `computations/cont_M6N12_D5_all_source_lifts.m`, `computations/cont_RR612_norm13_shard0.m`, `continuation_20261029/degree_independent_612_L8/RR_6_12_L8_char0_six_pencils_root.m`, `computations/cont_replay612_complete_v2.py`, `computations/cont_quintic612_class_orbits_v2.m` |
| Proposition 11.3 | `prop:6x12` | `computations/cont_witness612.m`, `computations/cont_char0_612_norm.m` |
| Lemma 11.4 | `lem:spor-315-jac` | `computations/probe315_galois3.m` |
| Lemma 11.5 | `lem:spor-315-classes` | `computations/diagonal315_D5.m`, `computations/diagonal315_D5_exact.m`, `exclusions/315d5/M3N15_D4_exact_cusp_lifts.txt` |
| Proposition 11.6 | `prop:3x15` | `computations/new315_witness.m`, `computations/simplify315.py` |
| Lemma 11.9 | `lem:spor-318-jac` | `computations/d31b_m318_paired_p7.m`, `computations/d31b_m318_kernel_structure_v2.m`, `computations/d31b_m3n18_frob_preflight.m`, `computations/d31b_m318_paired_p13.m`, `computations/d31b_m318_p13_lean_s0.m` |
| Lemma 11.11 | `lem:spor-318-classes` | `jobs/d31b_m318_short_cusp_divisors.log`, `continuation_20261031b/M3N18/endomorphism_gluing/short_cusp_divisors.m`, `computations/d31b_m318_relative_cubic_rank.m`, `computations/d31b_m318_relative_cubic_rank_p13.m`, `computations/d31b_m318_paired_three_relative_v2.py`, `independent_checks/recompute_318_classes.m`, `computations/d31b_m318_final_Hecke7_v2_s0.m`, `computations/d31b_m318_quartic_phase_binding.m`, `computations/d31b_m318_tau_frame_p7.m`, `computations/d31b_m318_relative_quartic_12.m`, `computations/d31b_m318_relative_quartic_17.m`, `computations/d31b_m318_direct_cusp_membership.py` |
| Proposition 11.12 | `prop:3x18` | `computations/d31b_m318_global_point13900_v2.m`, `computations/d31b_m318_parent_nf_reductions.m`, `independent_checks/sporadic_check.m`, `computations/d31b_m318_global_RR_reductionfalse.m`, `computations/d31b_m318_tau_cusp_class_relation.py` |
| Corollary 11.13 | `cor:3x36-6x18` | `independent_checks/aut_318_check.m`, `computations/late226_picard_setup.m` |
| Lemma 11.15 | `lem:spor-226-j0` | `computations/late226_j0_twists_v2.m` |
| Lemma 11.16 | `lem:spor-226-red` | `computations/late226_source.m`, `computations/late226_filter26.m`, `computations/late226_picard_s0.m`, `continuation_20261029b/C226/picard_template.m` |
| Lemma 11.17 | `lem:spor-226-nine` | `computations/late226_stabilizer_replay.py` |
| Proposition 11.18 | `prop:2x26` | `computations/resolve226_complete_native.m`, `computations/resolve226_bind_classes.py`, `computations/resolve226_coordinate_fields.py`, `computations/late226_two_classes_stabilizers.py` |
| Proposition 11.19 | `prop:52` | `computations/resolve226_complete_native.m` |

### Further checks

Computations that the proofs do not need, mostly repetitions of a computation by an independent program from the
exported data, or by a second method, grouped by the result they check (see Section 13 of the paper):

* Proposition 5.3 (`prop:globaltable`): `computations/global3_characters.m`, `jobs/global3_characters.log`, `exclusions/global3.json`, `computations/cached_kernels.m`, `jobs/cached_kernels.log`, `computations/additional_kernels.m`, `jobs/additional_kernels.log`, `computations/replay_new_certificates_v2.py`, `jobs/replay_new_certificates_v2.log`, `exclusions/global3_additional.json`, `jobs/global_archive_replay.log`, `inputs/global_archive/certify_global_packet.py`
* Lemma 5.4 (`lem:glob-odd`): `inputs/C115_archive/C115_image_even_orders_certified.m`, `inputs/C115_archive/replay_C115_norms.py`, `computations/replay_new_certificates_v2.py`, `computations/cont_replay135.py`, `jobs/cont_replay135.log`, `continuation_20261029/C135_independent_replay.json`, `continuation_20261029/image135_a16_q31_index1_chi9.json`, `continuation_20261029/image135_a16_q31_index2_chi3.json`, `continuation_20261029/image135_a16_q31_index3_chi36.json`, `continuation_20261029/image135_a16_q31_index4_chi12.json`, `continuation_20261029/image135_a16_q31_index5_chi18.json`, `continuation_20261029/image135_a16_q31_index6_chi6.json`, `continuation_20261029/krylov_N135_chi4_order36.json`
* Proposition 5.12 (`prop:155`): `computations/cont_C155_local2.m`, `jobs/cont_C155_local2.log`, `continuation_20261029/C155_local2_arithmetic.txt`, `computations/cont_replay_hecke_cases.py`, `jobs/cont_replay_hecke_cases.log`, `continuation_20261029/new_cases_independent_replay.json`
* Proposition 5.8 (`prop:341-403`): `computations/descent31.m`, `jobs/descent31.log`, `inputs/global_archive/winding341.log`
* Proposition 5.15 (`prop:95-111`): `computations/next95111_finite_controls.m`, `jobs/next95111_finite_controls.log`, `next_cases_20261030/dyadic_local/all_curves.m`, `computations/next95111_independent_finite.py`, `jobs/next95111_independent_finite.log`, `next_cases_20261030/cyclic_hecke/replay.py`
* Proposition 6.2 (`prop:11x11`): `computations/replay_cyclotomic11_v2.py`, `jobs/replay_cyclotomic11_v2.log`, `inputs/global_newforms_cache.json`
* Proposition 6.4 (`prop:3x21`): `computations/cont_replay321.py`, `jobs/cont_replay321.log`, `computations/cont_local_3_21_q4.m`
* Proposition 6.6 (`prop:3x24-4x24`): `computations/cont_local_3_24_q7.m`, `computations/cont_local_4_24_q5.m`, `computations/cont_replay_hecke_cases.py`, `jobs/cont_replay_hecke_cases.log`
* Proposition 6.9 (`prop:4x16`): `computations/cont_replay416_v4.py`, `jobs/cont_replay416_v4.log`
* Proposition 6.13 (`prop:4x20`): `computations/cont_replay420.py`, `jobs/cont_replay420.log`, `computations/cont_local_4_20_q9.m`
* Proposition 6.15 (`prop:3x27`): `computations/d31b_m327_cusp_parity_v2.py`, `computations/d31b_m327_torsion_gcd.m`, `computations/d31b_m327_independent_genuine_v2.py`, `computations/d31b_m3n27_parent_native_v3.m`, `computations/d31b_m3n27_quotient_A7_v4.m`
* Proposition 6.20 (`prop:3x33`): `computations/d31b_m333_winding_independent.py`, `computations/d31b_m333_independent_ranks.py`, `computations/d31b_m333_native_frame_replay.m`, `computations/d31b_m333_independent_genuine35.py`, `computations/d31b_m333_parent_native.m`
* Proposition 6.23 (`prop:4x28`): `computations/d31b_m428_independent_fourmap_v3.py`, `computations/d31b_m428_independent_union.py`, `computations/d31b_m428_parent_class_ideals_v2.m`, `computations/d31b_Y228_independent_genuine_v2_s0.py`, `computations/d31b_Y228_parent_native_v2.m`
* Proposition 7.7 (`prop:2x32`): `computations/next232_replay_sieve_v2.py`, `computations/next232_ideal_replay10_v2.m`, `computations/next232_ideal_replay789_v2.m`, `computations/next232_ideal_replayprobe6_v2.m`, `computations/next232_quotient13_native_replay.m`, `computations/next232_quotient13_cuspidal_saturation.py`, `computations/next232_integral_lattice_replay.py`, `computations/next232_quotient13_torsion_bound.py`, `computations/next232_final_replay_v3.py`
* Proposition 7.12 (`prop:2x34`): `computations/new234_sieve_replay_v2.py`, `computations/new234_ideal_replay0.m`, `computations/new234_lattice_replay.py`
* Proposition 7.19 (`prop:2x36`): `computations/new236_filter18_native_replay_D10s0_v2.m`, `computations/new236_filter18_native_replay_D10s1_v2.m`, `computations/new236_filter18_native_replay_D10s2_v2.m`, `computations/new236_filter18_native_replay_D10s3_v2.m`, `computations/new236_filter18_native_replay_D10s4_v2.m`, `computations/new236_filter18_native_replay_D10s5_v2.m`, `computations/new236_filter18_native_replay_D10s6_v2.m`, `computations/new236_filter18_native_replay_D10s7_v2.m`, `computations/new236_filter18_native_replay_D8_v2.m`, `computations/new236_filter18_native_replay_D9_v2.m`, `computations/new236_filter18_native_replay_common_v2.m`, `computations/new236_filter18_native_replay_ownsmall_v2.m`, `computations/new236_parent_ideal_genus9.m`, `computations/new236_parent_ideal_cyclic36_prefix.m`, `computations/new236_geom_genus9_sage_replay.py`, `computations/new236_geom_cyclic36_sage_replay.py`, `computations/new236_geom_quotient9_setup.m`, `computations/new236_geom_cusp2_replay_v2.py`, `computations/new236_cusp_lattice_replay_v2.py`
* Proposition 7.23 (`prop:2x38-2x40-2x50`): `continuation_20261029b/winding_N76_H39.m`, `computations/late238_control_s0.m`, `computations/late240_replay_v2_s0.m`, `computations/late250_kernel_certificate.m`, `computations/late250_replay_kernels.py`, `computations/late250_control_s0.m`, `computations/late250_control_s1.m`, `computations/late250_control_s2.m`, `continuation_20261029b/six_exclusions_replay.json`
* Proposition 7.30 (`prop:2x42`): `next_cases_20261030/C242/independent_arithmetic_replay.py`, `computations/next242_projection_replay_D10s0.py`, `computations/next242_projection_replay_D10s1.py`, `computations/next242_projection_replay_D10s2.py`, `computations/next242_projection_replay_D10s3.py`, `computations/next242_projection_replay_D10s4.py`, `computations/next242_projection_replay_D10s5.py`, `computations/next242_projection_replay_D10s6.py`, `computations/next242_projection_replay_D10s7.py`, `computations/next242_projection_replay_D9.py`, `computations/next242_projection_replay_low.py`, `computations/next242_native_D10s0.m`, `computations/next242_native_D10s1.m`, `computations/next242_native_D10s2.m`, `computations/next242_native_D10s3.m`, `computations/next242_native_D10s4.m`, `computations/next242_native_D10s5.m`, `computations/next242_native_D10s6.m`, `computations/next242_native_D10s7.m`, `computations/next242_native_D7.m`, `computations/next242_native_D8.m`, `computations/next242_native_D9.m`, `computations/next242_native_common.m`, `computations/next242_native_small.m`, `computations/next242_geom_g5_setup_v2.m`, `computations/next242_geom_g5_sage_replay_v2.py`, `computations/next242_geom_c21_sage_replay_v2.py`, `computations/next242_geom_x0_42_sage_replay_v2.py`, `computations/next242_order3_phases.m`, `computations/next242_independent_norm_A11pure.py`, `computations/next242_geom_case26_sage_replay.py`, `next_cases_20261030/C242/geometry/case26_parent_native_replay.m`
* Proposition 8.2 (`prop:44`): `computations/d31a_C44_independent_matrices.py`, `computations/d31a_C44_independent_jet_geometry.py`, `computations/d31a_C44_independent_isogenies_shard0_v2.py`, `computations/d31a_C44_independent_isogenies_shard1_v2.py`, `computations/d31a_C44_independent_isogenies_shard2_v2.py`, `computations/d31a_C44_independent_isogenies_shard3_v2.py`, `computations/d31a_C44_independent_isogenies_shard4_v2.py`, `computations/d31a_C44_independent_isogenies_shard5_v2.py`, `computations/d31a_C44_independent_isogenies_sharded_common_v2.py`, `computations/d31a_C44_independent_merge.py`, `continuation_20261031a/C44/independent_final_verification.json`, `continuation_20261031a/C44/cusp_audit.py`
* Proposition 8.4 (`prop:46`): `continuation_20261031a/C46/quotient_replay.m`, `computations/d31a_C46_quotient_replay.m`, `jobs/d31a_C46_quotient_replay.log`, `computations/d31a_C46_independent_quotient.py`, `jobs/d31a_C46_independent_quotient.log`, `computations/d31a_46_parent_exceptional_native_common.m`, `computations/d31a_46_parent_exceptional_genuine_common.py`, `continuation_20261031a/C46/parent_exceptional_summary.json`, `computations/d31a_C46_new_supports.m`, `jobs/d31a_C46_new_supports.log`
* Proposition 8.7 (`prop:48`): `continuation_20261031b/C48/cusp_A13_audit.json`, `continuation_20261031b/C48/independent_cusp_A13_labels.json`, `continuation_20261031b/C48/d31b_C48_independent_cusp_principal_v2.m`, `continuation_20261031b/C48/d31b_C48_independent_character_s0.py`, `continuation_20261031b/C48/d31b_C48_independent_character_s1.py`, `continuation_20261031b/C48/d31b_C48_independent_character_s2.py`, `continuation_20261031b/C48/d31b_C48_independent_character_s3.py`, `continuation_20261031b/C48/d31b_C48_independent_character_s3_v2.py`, `computations/d31b_C48_parent_native_q7.m`, `jobs/d31b_C48_parent_native_low_q13.log`, `jobs/d31b_C48_parent_native_low_q7.log`, `jobs/d31b_C48_parent_native_q13.log`, `jobs/d31b_C48_parent_native_q7.log`, `computations/d31b_C48_parent_native_low_q13.m`
* Proposition 8.15 (`prop:54`): `continuation_20261031a/C54/parent_gonality_controls.json`, `continuation_20261031b/C54/inputs/CM54_minimum_check.json`, `continuation_20261031b/C54/d31b_C54_parent_dyadic_replay.m`, `continuation_20261031b/C54/d31b_C54_parent_f5_rr_common.m`
* Proposition 8.20 (`prop:56`): `computations/d31a_C56_independent_fourmap_v3.py`, `continuation_20261031a/C56/independent_fourmap.json`, `computations/d31b_C56_independent_class_ideals.m`, `continuation_20261031b/Y228_shared/union_normalization.json`, `computations/d31b_Y228_known_decic_2descent.m`, `computations/d31b_Y228_known_decic_binding_sage.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case1009_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case271_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case327_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case491_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case578_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case65_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case818_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case837_v3.py`, `continuation_20261031b/C56/d31b_C56_sieve_rr_case978_v3.py`, `continuation_20261031b/C56/cusp_A11_degree2_audit.json`, `computations/d31b_C56_kummer_case491.m`, `computations/d31b_C56_kummer_direct_case491_v3.m`
* Proposition 8.23 (`prop:58`): `continuation_20261031a/C58/independent_decic_replay.py`, `jobs/d31a_C58_independent_replay.log`, `computations/d31a_C58_inherited_replay.py`, `jobs/d31a_C58_inherited_replay.log`, `continuation_20261031a/C58/parent_class_replay.m`, `jobs/d31a_58_parent_class_replay_v2.log`
* Proposition 8.25 (`prop:62`): `continuation_20261031a/C62/independent_decic_replay.py`, `jobs/d31a_C62_independent_replay.log`, `computations/d31a_C62_inherited_replay.py`, `jobs/d31a_C62_inherited_replay.log`, `continuation_20261031a/C62/parent_class_replay.m`, `jobs/d31a_62_parent_class_replay.log`
* Proposition 8.38 (`prop:70`): `computations/d31a_even_p3_counts.m`, `jobs/d31a_even_p3_counts.log`, `computations/d31a_C70_independent_genuine_v2.py`, `jobs/d31a_C70_independent_genuine_v2.log`, `continuation_20261031a/C70/parent_native_common.m`
* Proposition 8.30 (`prop:60`): `continuation_20261031b/C60/d31b_C60_independent_fibres.m`, `jobs/d31b_C60_independent_fibres.log`, `continuation_20261031b/C60/d31b_C60_genuine_source_maps.py`, `jobs/d31b_C60_genuine_source_maps.log`, `continuation_20261031b/C60/d31b_C60_independent_class_ideals.m`, `jobs/d31b_C60_independent_class_ideals.log`, `computations/d31b_C60_parent_halving_v2.py`, `jobs/d31b_C60_parent_halving_v2.log`, `continuation_20261031b/C60/degree10_series.m`
* Proposition 8.35 (`prop:66`): `computations/d31b_C66_F2_independent_v3.py`, `jobs/d31b_C66_F2_independent_v3.log`, `computations/d31b_C66_F5_common_basis58.m`, `computations/d31b_C66_common_basis58_native.m`, `computations/d31b_C66_common_basis58_replay.py`, `jobs/d31b_C66_F5_common_basis58.log`, `jobs/d31b_C66_common_basis58_native.log`, `jobs/d31b_C66_common_basis58_replay.log`, `computations/d31b_C66_gonal_genuine_replay_v3.py`, `jobs/d31b_C66_gonal_genuine_replay_v3.log`, `continuation_20261031b/C66/gonal_fibre_filter.m`
* Proposition 8.42 (`prop:74`): `computations/d31b_C74_parent_census.py`, `jobs/d31b_C74_parent_census.log`, `computations/d31b_C74_parent_prefix_ranks.py`, `jobs/d31b_C74_parent_prefix_ranks.log`, `computations/d31b_C74_parent_native_v2.m`, `jobs/d31b_C74_parent_native_v2.log`, `continuation_20261031b/C74/inputs/work/degree9/prime_orders/continued_20261026/P37/replay_rational_two.py`
* Proposition 9.3 (`prop:49`): `next_cases_20261030/cyclic_primepowers/pp49_finite_controls.m`, `jobs/pp49_finite_controls.log`, `next_cases_20261030/cyclic_primepowers/pp49_completion_control.m`, `next_cases_20261030/cyclic_primepowers/pp49_cartier_transport_replay.py`
* Proposition 9.6 (`prop:51`): `next_cases_20261030/cyclic_further/C51/further51_finite_controls.m`, `next_cases_20261030/cyclic_further/C51/further51_replay_new_jets.py`
* Proposition 9.9 (`prop:55`): `computations/contb55_sage_finite_v2.py`, `computations/contb55_parent_native.m`, `continuation_20261030b/C55/replay_source_injectivity.py`, `computations/contb55_cusp_coupled.m`
* Proposition 9.11 (`prop:75`): `continuation_20261030b/C75/cont75_verify_finite_v2.py`, `computations/cont75_parent_native.m`, `continuation_20261030b/C75/cont75_retained_FI.py`
* Proposition 9.15 (`prop:77`): `computations/contb77_sage_census.py`, `computations/contb77_sage_nc_case1.py`, `computations/contb77_parent_moduli_binding.m`, `computations/contb77_parent_native.m`, `computations/contb77_fi_native_controls.m`, `continuation_20261030b/C77/replay_fi10_v2.py`
* Proposition 9.19 (`prop:85`): `continuation_20261030c/C85/contc85_verify_R1.py`, `continuation_20261030c/C85/contc85_arc_q4_control.m`, `computations/contc_parent_85_93_cusp_replay.py`
* Proposition 9.24 (`prop:91`): `continuation_20261030c/C91/verify_arithmetic.py`, `continuation_20261030c/C91/verify_inherited_R1.py`, `computations/contc91_parent_native_v3.m`, `computations/contc91_independent_CM_level7_binding.m`, `computations/contc91_independent_CM_level7_sage.py`
* Proposition 9.26 (`prop:93`): `continuation_20261030c/C93/contc93_operator_control.py`, `computations/contc_parent_85_93_cusp_replay.py`
* Proposition 10.6 (`prop:57`): `continuation_20261030b/C57/cont57_verify_coverage.py`, `jobs/cont57_verify_coverage.log`, `continuation_20261030b/C57/cont57_verify_jets_v3.py`, `continuation_20261030b/C57/cont57_verify_RR_v2.py`
* Proposition 10.14 (`prop:63`): `computations/d31a_63_native_replay.m`, `computations/d31a_63_independent_genuine.py`, `computations/d31a_63_independent_bind_boundary_v2.py`, `continuation_20261031a/C63/old_kernel_replay.log`, `continuation_20261031a/C63/certify.py`
* Proposition 10.10 (`prop:65`): `computations/contc65_verify_m4ri.py`, `computations/contc65_verify_velu.py`, `computations/contc65_verify_jets_v2.py`, `computations/contc65_source_binding.py`, `computations/contc65_parent_prefix_replay.py`, `computations/contc65_inherited_deep.py`
* Proposition 10.20 (`prop:81`): `computations/b81_parent_two_torsion.py`, `continuation_20261030b/C81/b81_parent_source_filter.py`, `computations/b81_parent_x0_filter.m`, `computations/b81_parent_rr_completions_v2.m`, `computations/b81_sage_local_arcs_v3.py`
* Proposition 10.23 (`prop:87`): `computations/b87_raw_hecke_control.m`, `computations/b87_parent_native_control.m`, `computations/b87_independent_moduli.py`, `continuation_20261030b/C87/independent_source_cusps.py`
* Proposition 10.26 (`prop:69`): `next_cases_20261030/C69/independent_source_cusps.py`, `jobs/independent_source_cusps.log`, `computations/next69_degree10_base.m`, `jobs/next69_degree10_base.log`, `computations/cont69_native_RR_replay_v2.m`, `jobs/cont69_native_RR_replay_v2.log`, `computations/cont69_independent_moduli_v4.py`, `jobs/cont69_independent_moduli_v4.log`, `next_cases_20261030/C69/direct_full/direct69_verify_input_v2.py`, `computations/cont69_q67_parent_native_v2.m`, `jobs/cont69_q67_parent_native_v2.log`, `next_cases_20261030/C69/quotient67/proof_packet_v2.m`, `next_cases_20261030/C69/quotient67/q67_geometry_audit.py`
* Proposition 10.29 (`prop:78-3x39`): `computations/d31a_78_native_replay_v2.m`, `jobs/d31a_78_native_replay_v2.log`, `computations/d31a_78_independent_genuine.py`, `jobs/d31a_78_independent_genuine.log`, `continuation_20261031a/C117/resolution.tex`, `continuation_20261031a/C117/finite_and_inherited_verified.json`, `continuation_20261031a/C78/inputs/work/degree9/cyclic/continuation_20261027/C78_X139/replay_dyadic_injection39.py`, `jobs/witness_C39.log`, `computations/witness_C39.m`
* Proposition 10.34 (`prop:90`): `computations/d31a_90_cycle_and_boundary_controls.py`, `jobs/d31a_90_cycle_and_boundary_controls.log`, `continuation_20261031a/C90/independent_cycle_and_boundary_controls.json`, `computations/d31a_90_independent_genuine.py`, `continuation_20261031a/C90/independent_genuine_and_known_lift.json`, `continuation_20261031a/C90/inputs/old45/N45_quotient_injectivity_control.m`, `computations/d31a_90_A7_part1.m`, `computations/d31a_90_A11_case10.m`, `computations/d31a_90_source.m`
* Proposition 10.38 (`prop:99`): `computations/contc99_parent_native_A5.m`, `jobs/contc99_parent_native_A5.log`, `computations/contc99_independent_genuine.py`, `computations/contc99_independent_native_binding.py`, `computations/contc99_peer_Y29_native.m`, `jobs/contc99_peer_Y29_native.log`, `computations/contc99_independent_Y29_genuine.py`, `computations/contc99_independent_Y29_binding.py`, `computations/contc_parent_source_census.py`
* Proposition 10.40 (`prop:105`): `computations/contc105_parent_native.m`, `jobs/contc105_parent_native.log`, `computations/contc105_sage_velu.py`, `continuation_20261030c/C105/independent_sage_velu.json`, `computations/contc_parent_source_census.py`, `continuation_20261030c/C105/inputs/parent_source_census.json`
* Proposition 11.3 (`prop:6x12`): `computations/cont_replay_witness612_v2.py`, `jobs/cont_replay_witness612_v2.log`, `computations/cont_replay612_absolute.py`, `continuation_20261029/M6N12_absolute_polynomial_coefficients.json`, `computations/cont_replay612_source_v2.py`, `computations/cont_RR612_norm19_shard0.m`, `computations/cont_RR612_norm19_shard1.m`, `computations/cont_RR612_norm19_shard2.m`, `computations/cont_RR612_norm19_shard3.m`, `computations/cont_RR612_norm19_shard4.m`, `computations/cont_RR612_norm19_shard5.m`, `computations/cont_RR612_norm19_shard6.m`, `computations/cont_RR612_norm19_shard7.m`, `computations/cont_replay612_norm.py`, `computations/cont_replay612_complete_v2.py`, `computations/cont_frob612.m`, `computations/cont_validate612_public_basis.m`, `computations/cont_hecke_torsion_6_12_v2.m`
* Proposition 11.6 (`prop:3x15`): `computations/replay315_complete_v2.py`, `jobs/replay315_complete_v2.log`, `exclusions/315d5/independent_replay.json`, `computations/probe315_galois3.m`, `computations/315_frobenius.m`, `inputs/M3N15_octic_archive/noncyclic/M3N15_necessity_frobenius.log`
* Proposition 11.12 (`prop:3x18`): `continuation_20261031b/M3N18/d31b_m318_sieve_genuine_A7.py`, `computations/d31b_m318_parent_relative_p7.m`, `computations/d31b_m318_parent_native_shortclass13_v2.m`, `computations/d31b_m318_parent_global_basis_13900.m`, `continuation_20261031b/M3N18/d31b_m318_pointcheck_exact_v3.py`, `computations/d31b_m318_class_frame_audit.py`, `computations/d31b_m318_sage_geocusp_probe_v3.py`
* Proposition 11.18 (`prop:2x26`): `computations/late226_replay_sieve.py`, `computations/late226_replay_tate.py`, `computations/late226_F5_RR.m`, `computations/late226_second_prime_v3.m`

## Independent checks made for this repository

The folder `independent_checks/` contains checks made after the original computations, for the paper and this
repository, with code written independently of the original scripts. Each script is accompanied by its log.

| Script | What it checks | Paper |
|---|---|---|
| `minimal_groups.py` | enumerates from scratch the 110 minimal groups outside the list of Theorem 1.1 | Lemma 4.1 |
| `occ_torsion_check.m` | `TorsionSubgroup` of the curves of Propositions 3.4 and 3.5 (models of the scripts `computations/witness_C*.m`, `computations/quintic_full2_*_v2.m`) | Props. 3.4, 3.5 |
| `vanhoeij_check.m` | the curves of Proposition 3.4 built from van Hoeij's coordinates: order of (0,0), `TorsionSubgroup`, field | Prop. 3.4 |
| `same_pair_check.m` | the record's models and van Hoeij's curves have isomorphic fields and the same j-invariant | Prop. 3.4, Lemma 10.27 |
| `lmfdb_fields.py`, `lmfdb_sporadic_fields.py` | polredabs and LMFDB labels of the 17 decic fields of Sections 3 and 11 (LMFDB mirror) | Tables, Thm. 1.2 |
| `sporadic_check.m` | the curves of Section 11: `TorsionSubgroup`, minimal polynomials of j (degree 5, not integral, distinct for the two pairs with Z/3 x Z/18 and with Z/2 x Z/26) | Props. 11.3, 11.6, 11.12, 11.18 |
| `aut_318_check.m` | Aut(K_{3,18}) has order 2 and its non-trivial element moves zeta_3 | Remark 11.14 |
| `cusp_degrees.py` | degrees of the closed cusps of X_1(2,28) over Q and modulo 3, of X_1(70) modulo 3 and of X_1(3,33) over F_4 | Lemmas 6.22, 8.37, proof of Prop. 6.20 |
| `recompute_318_classes.m` | recomputes from the model of X_1(3,18) the classes of all non-cuspidal closed points of degree 2 and 3 at the primes above 7 and 13 and compares them with the stored data used for the types [2,2] and [2,3] | Lemma 11.11 |
| `rerun_C49/` | clean reruns of the two exports used for Z/49, compared with the stored data (`COMPARE.txt`) | Prop. 9.3 |
| `check_ledger_hashes.py` | the SHA-256 hashes of the ledger `certified_results.json` against the files present | -- |

All of them agree with the original computations.

## License

MIT, see `LICENSE`. The folders `shared/mdmagma` and `shared/mdsage` (and their copies inside case folders) are
copies of the packages of their authors, under their own licenses.
