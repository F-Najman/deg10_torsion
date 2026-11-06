# Final cyclic 42 Hecke tests

Nine of the last ten C242 supports are excluded here. Case26 passes these tests and is resolved separately by the genus17 section argument in `../geometry/`.

The operator A17=T17−1−17⟨17⟩ kills **every** F5-rational cuspidal divisor: on Tate cusp labels σ17=σ5^5 because17≡5^5(mod42), and the boundary formula is Tℓ=σℓ+ℓ⟨ℓ⟩σℓ⁻¹. `boundary17_certificate.json` independently verifies this convention on all 48 geometric X1(42)cusps and every F5 closed-cusp basis divisor. `criterion.tex` gives the proof, including the rational-Jacobian annihilator and specialization. Purely noncuspidal degree 10 supports may use the cheaper A11.

The six source maps retain Q,P+Q,P+22Q, or quotient respectively by P,21Q,P+21Q and retain φQ,φ(P+22Q),φQ. Checks stop after the first nonprincipal norm.

|Case ID|Hecke prime|Selected map|h0|
|---|---:|---:|---:|
|33546|11|1|0|
|39662|11|1|0|
|39669|11|1|0|
|56937|11|1|0|
|59070|11|1|0|
|35579|17|2|0|
|35633|17|1|0|
|35341|17|1|0|
|35601|17|1|0|

`A11pure/`, `A17s0/` and `A17s1/` contain the exact full algebraic function-field model, both maximal-order bases, explicit kernels and level images, both ideal bases for O(Adiv), and native Riemann–Roch output. `final_test_summary.json` and `selected_native_replay_cases.m` are the structured summaries. The original `arithmetic.m` files contain unquoted human-readable diagnostic strings from Magma's Sprint; they are transcripts, not a parser input. All proof objects in the model, kernel and ideal-basis files are numeric and directly readable.

The independent Sage replay checks 16 unique source summands, 78 Frobenius orbits of kernels and 240 geometric cyclic kernels. For each summand the product of the kernel norms equals the entire monic division polynomial, so no kernel is omitted or counted twice. Every Vélu codomain, marked point, cyclic 42 plane equation and residue-degree norm is independently reconstructed. It also reconstructs each final degree-zero divisor; see `independent_norm_*.json` and `divisor_geometry_*.json`.

The parent independently reconstructs the fresh genus 25, degree 15 function field, checks the finite and infinite integral-basis transition matrices in both directions, reconstructs O(Adiv) from the saved ideal bases, and recomputes all nine Riemann–Roch dimensions as zero. The numeric results are `TAG/parent_native_RR_replay.m`.

`nine_exclusions_completion.json` combines these checks and verifies all owned workers have exited. No paused worker or automatic launcher remains. The three original norm workers and all three independent Sage replays completed within their explicit limits.
