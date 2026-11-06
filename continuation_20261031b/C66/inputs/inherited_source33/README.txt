This packet proves C66 is excluded over degree-nine number fields by classifying the rational effective degree-nine divisors on X1(33). The proof is exclusion66.tex.

The authoritative complete F2 computation is rebound/classify_degree9_33.m. It proves the finite Picard group, the full rational cuspidal image, all 155770 class tests, and every known degree7/8/9 witness. Rebound exports the global cusp signatures and all127 finite noncuspidal places through degree9.

The exact global rational Jacobian upper bound is the unconditional theorem and Table2 of DEHMZ, retained in imports/DEHMZ.tex. Equality of this global order with the reduced global cusp image supplies specialization injectivity at2, including the2-primary part. The odd-prime point-count diagnostic outside this manifest is not used to infer that upper bound.

The10 unknown degree9 finite classes are represented by short GLOBAL signed cuspidal degree9 divisors. All have h0=0 modulo5. The shared positive divisor A has degree48; its28 basis functions and15 local restriction tables are checked afresh by F5_rr/verify33.m without any Riemann--Roch or class-group search. F5_rr/replay33.py independently checks basis independence and all ten rank28 restriction matrices. This rules out these global classes by specialization.

full_torsion33.py directly counts all625 affine pairs over F25 for the known nonic Tate curve, obtaining order33 at a good unramified prime. The native global example has a point of exact order33, so its full torsion is33.

Read-only verification:
python3 work/degree9/cyclic/continuation_20261026/C66_quotient17/source33/certify66.py

The model_retained folder preserves the completed F5 model producer and exact metadata. Its source package is byte-identical to the195 retained Magma source/model files. reproduce_model33.m uses the local package copy to reproduce this model independently of the shared package path.

The earlier genus17, Hecke, and finite-character computations are preserved elsewhere under C66_quotient17. They are not dependencies of this exclusion packet. Earlier degree7/8 folders were read only.
