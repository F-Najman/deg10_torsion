C2 x C32 and C64 over decic fields — resolved, 30 October 2026

RESULT: both subgroup obstructions are excluded. This reduces the global
Phi(10) frontier from 45 unresolved groups to 43. Phi(10) remains incomplete.
The mathematical proof is resolution.tex, with gonality.tex as its first lemma.

FINAL ARGUMENT:
1. A degree-two map to X_Delta(64), Delta=<15,-1>, genus 13, and Orlic's
   gonality bound >=6 prove gon_Q X1(2,32)>=12 by Castelnuovo-Severi.
2. At 3, 9504 noncuspidal closed points through degree ten form 305 full
   level orbits. Four maps to X1(32) eliminate 304. The sole remainder
   has degree six and needs a degree-four cusp complement in every target.
3. Every rational residue cusp is globally rational. Thus the complement
   has no rational cusp. The involution h:(P,Q)->(P,15Q) fixes the remaining
   closed support and all degree-two/four cusps geometrically. Rank-zero
   torsion injectivity and gonality force hD=D. Hence D is pulled back from
   degree five on the genus-13 quotient. Cusp ramification leaves a cubic
   noncuspidal support plus one of eight quadratic quotient cusps.
4. The quotient's finite Jacobian is [2,2,2,2,8,40,40,80], with rational
   cuspidal image [2,2,4,20,40,40], order 512000, index 32.
   Integral modular symbols plus geometric-cuspidal Galois descent prove
   [J(Q):H_Q] divides 2. Two automorphism coset sets, for 1 and diamond 3,
   are disjoint. Their coordinates are in resolution.tex. This excludes
   C2 x C32; the existing isogeny-transfer lemma excludes C64.

COMPACT VERIFICATION (run from the project root):
  /opt/sage/envs/sage/bin/python3 work/deg10/computations/next232_replay_sieve_v2.py
  magma -b work/deg10/computations/next232_ideal_replay10_v2.m
  magma -b work/deg10/computations/next232_ideal_replay789_v2.m
  magma -b work/deg10/computations/next232_ideal_replayprobe6_v2.m
  magma -b work/deg10/computations/next232_quotient13_native_replay.m
  /opt/sage/envs/sage/bin/python3 work/deg10/computations/next232_final_replay_v3.py
  /opt/sage/envs/sage/bin/python3 work/deg10/computations/next232_integral_lattice_replay.py
These replay saved exact evidence. The fresh quotient class-group proof
uses Proof:=true and checks equality of the archived and fresh factor-base
relation kernels, including the free degree coordinate.

QUOTIENT MODEL:
Over X1(2,16), with normalized E:y^2=x^3+c*x^2+b*(1-b-c)*x,
Q=(b,b) and T=8Q, adjoin w with w^2=x(Q)-x(T).
The resulting field has genus 13 and is the quotient by h. The rational
function giving w on a point S of order 32 with Q=2S is
  alpha*((x(S)-x(T))^2-B_T)/(2*y(S)),
where alpha=(x(Q)-x(P))/y(Q) and B_T is the derivative of the cubic at T.
The saved native model uses degree eight over F3(t), rather than the
expensive degree-nineteen source model. Its full Picard proof takes seconds.

IMPORTANT SCOPE:
The exploratory T11 test is principal. The cusp-aware T5 test leaves three
complements; all survive T7. These tests are NOT needed in the final proof:
the absence of rational cusps in the complement already gives h-invariance.
The exponent probe timed out and is not a certificate. The slow genus-37
class-group job and slow moduli-coordinate inversion were terminated when
abandoned; see retired_source_picard.json and retired_slow_inversion.json.
Their logs and actual exits are retained. No uncompleted probe is used as a
proof. The final proof uses only completed, independently checked evidence.
Earlier degree-six through degree-nine work is unchanged. No paused workers
or automatic launchers are needed for this completed case.
