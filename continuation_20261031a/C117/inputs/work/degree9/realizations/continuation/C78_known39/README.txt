Known nonic39 witnesses for the degree9 C78 investigation.

known39.tex proves the global point order, fulltorsion39, and the four distinct global/nonic diamond trace divisors. Their good unramified prime29 reduction has exactly39 points, so there is no rational point of order2 on any of these known curves.

positive_divisors39.json gives the four reduction codes requested by the C78 worker:
[closed_degree, original_F512_field_polynomial_coefficients, [Tate_b_coefficients,Tate_c_coefficients]].
The original field polynomial is T^9+T^8+1. Every coefficient vector uses its formal root, and each Frobenius-canonical key is the lexicographically smallest concatenation of the two nine-bit coefficient vectors among all nine Frobenius powers. Consequently any explicit field embedding into this fixed F512 yields the same key.

diamond39.jsonseq also records all twelve physical multipliers, their exact global Tate coefficients, primitive-element witnesses, and reduction vectors. The global trace triples are [1,16,17], [2,5,7], [4,10,14], [8,11,19]; each triple describes one global divisor and one finite closed point. Every global trace has degree9 and reduces with multiplicity1.

Read-only verification:
python3 work/degree9/realizations/continuation/C78_known39/certify39.py

This packet supplies positive global witnesses and their full-torsion and reduction bindings. It does not assert completeness of the C78 sieve or exclude C78.
