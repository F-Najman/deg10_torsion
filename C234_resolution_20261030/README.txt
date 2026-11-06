C2 x C34 and C68 — completed degree-ten exclusion, 30 October 2026

RESULT: no elliptic curve over a degree-ten number field contains either
C2 x C34 or C68. Main proof: resolution.tex. Do not edit the shared ledger
without checking the parent task's independently coordinated C236 result.

SOURCE SIEVE:
X1(2,34), genus49; noncuspidal counts at3 throughdegree10:
0,0,0,36,48,48,144,738,2160,6504. Full48-element levelgroup yields
2,1,2,3,18,45,138 own-degree orbits. All9678 closedpoints are checked.
The521 supportpatterns are209singles and72+96+96+48 pairs(44,45,46,55).
All521 fail one of the four maps to X1(34). Its rational image is fully
known [8760,595680], index128 in [4,4,4,8760,1191360]. Every finite cusp
class belongs to the rational image. Independent Sage+integer-HNF replay
checks all842 normed summands. Native ideal replays check328 target images.

PURE CUSP BRANCH / GONALITY:
Only gon_F3>=10 is imported from NV. Empty NC sieve forces all four finite
rational fibres of any proposed rational degree10pencil to be cuspidal.
The24 rational cusps distribute over these four fibres. The two largest
rational-support counts satisfy n1>=6 and n2>=ceil((24-n1)/3)>=5. Both
fibres therefore contain only rational cusps (allothercusps are octic).
A degree10 divisor with n distinct rational supports has squarednorm at
most (11-n)^2+n-1; for the two largest fibres the sum is<=60. Their
principal difference lifts to a QQ principal cuspidal divisor by Manin-
Drinfeld and odd unramified torsion injectivity. The full rational-cusp
modular-unit lattice has rank23 inZ24 and has no nonzero vector of norm
<=60, proven independently in Magma and PARI after exact period-kernel
reconstruction. Thus gon_Q>=11. Rankzero then excludes entirelycuspal
reduction of a noncuspidal degree10trace divisor.

VERIFICATION:
- new234_sieve_replay_v2.py: full finite-field source and521-case replay.
- new234_ideal_replay0/1/2.m: exact target class ideal identities.
- new234_lattice_replay.py: exact periodkernel reconstruction andPARI check.
- new234_short_cusp_units.m: native complete ShortVectors bound60 check.
All successful logs have explicit completion markers and exited status.
The initially duplicated Frobenius output from independent source shards
is preserved in frobenius_raw_twice.m; frobenius.m is one identical packet.
See frobenius_deduplication.json. No failed probe is accepted as evidence.

OWNERSHIP:
This task created only C234_resolution_20261030 and computations/new234*.
All case workers have exited; no automatic launchers or paused jobs remain.
Earlier degree-six through degree-nine sources were read but not changed.
