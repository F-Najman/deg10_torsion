C85 over degree-ten fields: COMPLETE EXCLUSION; parent mathematics PASS.

Proof: resolution.tex. Validator: certify85.py. The full source Jacobian
is not assumed to have rank zero. Exact Hecke polynomials k2 and k3 kill
the only possible positive-rank factors. The local full-Jacobian operators
are (T2^mu-1)k2(T2^mu) and (T2^mu-1)k3(T3^mu). On the canonical mu-level
model the dyadic factor is T2-I. Point-level T2-delta2 has a different
normalization, related by Fricke and an invertible diamond.

Every pure cusp branch is covered. The rational-prefix relation minimum
is32. Exact transported quartic and octic arcs give2330064 normalized
mixed patterns, all rank10. Parent independently enumerated the relation
space and checked complete matrices at all mixed leaves in C++.
All native arcs were independently recomputed from different canonical
unit minors and matched coefficientwise. Frobenius and diamond transport
are checked explicitly. Retained precisions are2,1,1 coefficients in the
quartic and two octic types; genuine_cusp_arcs_verified.json is authoritative.
Generic legacy progress labels are not used as precision claims.

Four new NC10 orbits have genuine ordinary85-torsion jets. Every actual
3-isogeny and point image is rebuilt independently in Sage. The four
necessary point-convention divisors have signed degree80 and canonical
multiplication rank225, independently checked with a different quotient
algorithm. The old degree9 obstruction extends with its rational cusp
complement. Order17 occurs quintically; no quintic17 exclusion is used.

All accepted inherited files are locally snapshotted under inherited_degree9.
Their original hashes are unchanged. No Magma workers or launchers remain.
