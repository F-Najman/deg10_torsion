# Written for this repository (not part of the original record).
# Degrees of the closed cusps of two modular curves, which the proofs of prop:4x28 (Section 6) and prop:70
# (Section 8) use and for which the original record has no script:
#   * the Galois orbits of cusps of X_1(2,28) over Q (expected 1^12, 2^6, 6^4);
#   * the closed cusps of X_1(70) over F_3, i.e. the orbits of the Frobenius sigma_3 (expected signature below);
#   * the closed cusps of X_1(3,33) over F_4 (expected 1^40, 5^8), used in the proof of prop:3x33 (Section 6).
# X_1(2,28) is modelled, as in the package mdsage, by X_H(56) with H = <-1, 29> (mdsage's Gamma11(2,14)).
# The multisets of orbit sizes do not depend on the convention used for the Galois action on the cusps.
# Run with: sage -python independent_checks/cusp_degrees.py
import sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared/mdsage'))
from sage.all import Gamma1, gcd
from mdsage.maartens_sage_functions import Gamma11, galois_action

def orbit_sizes(G, gens):
    N = G.level()
    cusps = [G.reduce_cusp(c) for c in G.cusps()]
    seen, sizes = set(), []
    for c in cusps:
        if c in seen:
            continue
        orbit, todo = {c}, [c]
        while todo:
            d = todo.pop()
            for t in gens:
                e = G.reduce_cusp(galois_action(d, t, N))
                if e not in orbit:
                    orbit.add(e); todo.append(e)
        seen |= orbit
        sizes.append(len(orbit))
    return sorted(Counter(sizes).items())

G = Gamma11(2, 14)
units = [t for t in range(1, 56) if gcd(t, 56) == 1]
print('CUSP_DEGREES X_1(2,28) over Q:', orbit_sizes(G, units), 'number of cusps', len(G.cusps()))
# How the Galois orbits over Q of X_1(2,28) split modulo 3 (orbits of sigma_3 inside each orbit over Q).
def split(G, units, p):
    N = G.level()
    out = []
    cusps = [G.reduce_cusp(c) for c in G.cusps()]
    done = set()
    for c in cusps:
        if c in done: continue
        orbit, todo = {c}, [c]
        while todo:
            d = todo.pop()
            for t in units:
                e = G.reduce_cusp(galois_action(d, t, N))
                if e not in orbit: orbit.add(e); todo.append(e)
        done |= orbit
        sub, rest = [], set(orbit)
        while rest:
            d = rest.pop(); o = {d}; e = G.reduce_cusp(galois_action(d, p, N))
            while e not in o: o.add(e); e = G.reduce_cusp(galois_action(e, p, N))
            rest -= o; sub.append(len(o))
        out.append((len(orbit), tuple(sorted(sub))))
    return sorted(Counter(out).items())
G = Gamma11(2, 14)
print('CUSP_SPLITTING X_1(2,28) mod 3 (degree over Q, degrees over F_3) -> count:', split(G, units, 3))
G = Gamma1(70)
print('CUSP_DEGREES X_1(70) over Q:', orbit_sizes(G, [t for t in range(1, 70) if gcd(t, 70) == 1]))
print('CUSP_SPLITTING X_1(70) mod 3 (degree over Q, degrees over F_3) -> count:', split(G, [t for t in range(1, 70) if gcd(t, 70) == 1], 3))
print('CUSP_DEGREES X_1(70) over F_3:', orbit_sizes(G, [3]), 'number of cusps', len(G.cusps()))
# X_1(3,33) over F_4 (the prime above 2 of Q(zeta_3)), modelled by X_H(99), H = <-1, 34, 67>: orbits of sigma_4.
G = Gamma11(3, 11)
print('CUSP_DEGREES X_1(3,33) over F_4:', orbit_sizes(G, [4]), 'number of cusps', len(G.cusps()))
print('CUSP_DEGREES_DONE')
