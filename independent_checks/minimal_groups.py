# Written for this repository: independent enumeration of the 110 groups of Lemma 4.1 (minimal groups outside
# the list of Theorem 1.1, with phi(m) | 10 and all primes dividing n at most 37).
# Independent enumeration of the minimal groups Z/m x Z/n (m | n) outside the list L of Theorem 1.1,
# with phi(m) | 10 and all primes dividing n in S(10) = primes up to 37.
from math import gcd
from sympy import primerange, totient, divisors, primefactors
S10 = list(primerange(2, 38))
cyc = set(range(1, 41)) | {42, 45, 50}
L = {(1, n) for n in cyc}
L |= {(2, 2*k) for k in range(1, 16)}
L |= {(3, 3*k) for k in range(1, 7)}
L |= {(4, 4*k) for k in range(1, 4)}
L |= {(6, 6*k) for k in range(1, 3)}
assert len(L) == 69
def subgroups(m, n):
    # isomorphism types of subgroups of Z/m x Z/n: Z/a x Z/b with a | b, a | m, b | n
    return {(a, b) for a in divisors(m) for b in divisors(n) if b % a == 0}
# subgroup closure of L
for (m, n) in L:
    assert subgroups(m, n) <= L, (m, n)
ms = [m for m in range(1, 50) if 10 % totient(m) == 0]
print("possible m:", ms)
maxcyc = max(cyc)
M = []
for m in ms:
    for n in range(m, maxcyc * 37 + 1, m):
        if not all(p in S10 for p in primefactors(n)): continue
        if (m, n) in L: continue
        # minimal: every proper subgroup in L
        if all(s in L for s in subgroups(m, n) if s != (m, n)):
            M.append((m, n))
print(len(M))
for m in ms:
    print(m, [n for (mm, n) in M if mm == m])
