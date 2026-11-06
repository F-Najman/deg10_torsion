C236 cuspidal boundary certificate, 30 October 2026

Proved: every degree10 map X1(2,36)/F5 -> P1 has a rational fibre with
noncuspidal support. This is a standalone local statement. It closes the
entirely cuspidal branch once all noncuspidal reductions of rational degree10
divisors with entirely noncuspidal global support have been excluded.
Do not assert that this local certificate alone excludes C2 x C36.

The full Gal(F5)-invariant geometric cusp-unit lattice has rank35 in Z64.
There are36 finite cusps with closed degrees 1^24,2^8,6^4. Weighted norm
on closed coordinates equals ordinary norm on64 geometric coordinates.
Magma ShortVectors<=104 gives24552 nonzero pairs and no positive degree10.
Independent Sage rebuild of the period kernel plus PARI qfminim gives the
same24552 and no positive degree10. Both computations completed normally.

Proof: proof.tex. Data: lattice.json/m,short_vectors.m,independent_replay.json.
Primary scripts: new236_cusp_lattice.py,new236_cusp_short.m,
new236_cusp_lattice_replay_v2.py under computations/. The initial replay
omitted the degree-zero condition and failed; it is not proof evidence.
The successful replay explicitly includes the degree column.
