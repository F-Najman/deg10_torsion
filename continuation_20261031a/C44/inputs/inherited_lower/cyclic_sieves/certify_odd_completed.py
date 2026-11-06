#!/usr/bin/env python3
"""Bind completed octic odd-prime sieves after their mathematical review."""
from pathlib import Path
import hashlib
import json
import sys

B = Path(__file__).resolve().parents[1]
n = int(sys.argv[1])
expected = {44:138, 46:384}
assert n in expected
a = json.loads((B/f'cyclic_sieves/N{n}_audit.json').read_text())
assert a['complete'] and a['partition_types'] == 66
assert a['tested_divisors'] == expected[n] and a['survivors'] == 0
assert a['independent_counts_agree'] and a['prime'] == 3
assert hashlib.sha256((B/f'cyclic_sieves/N{n}_p3_smart.log').read_bytes()).hexdigest() == a['log_sha256']
files = [
    f'cyclic_sieves/exclusion{n}.tex',
    f'cyclic_sieves/N{n}_p3_smart.log', f'cyclic_sieves/N{n}_audit.json',
    f'cyclic_sieves/N{n}_audit.log', 'cyclic_sieves/hecke_sieve_smart.m',
    'cyclic_sieves/audit_sieves.py', 'cyclic_sieves/certify_odd_completed.py',
    'cyclic_sieves/odd_prime_criterion.tex', 'cyclic_sieves/independent_counts.m',
    'cyclic_sieves/independent_counts.log',
    'refs/najman_varivoda_2602.03513v4.html',
    'small_cyclic/C34/inputs/DEHMZ.tex',
]
ev = {f:hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files}
ledger = json.loads((B/'certified_results.json').read_text())
reference = next(r for r in ledger['exclusions'] if r['group'] == [1,38])
for r in reference['evidence']:
    if r['path'].startswith('shared/'):
        assert hashlib.sha256((B/r['path']).read_bytes()).hexdigest() == r['sha256']
        ev[r['path']] = r['sha256']
manifest = dict(a, rank_over_Q=0, rational_gonality_lower_bound=10,
    all_cuspidal_branch='Unramified odd-prime torsion injectivity and gonality>8.',
    evidence=[{'path':f,'sha256':h} for f,h in sorted(ev.items())])
mp = B/f'cyclic_sieves/N{n}_manifest.json'
mp.write_text(json.dumps(manifest,indent=2)+'\n')
ev[str(mp.relative_to(B))] = hashlib.sha256(mp.read_bytes()).hexdigest()
assert not any(r['group'] == [1,n] for r in ledger['exclusions'])
ledger['exclusions'].append({
    'degree':8, 'group':[1,n], 'audited':True,
    'proof':f'cyclic_sieves/exclusion{n}.tex',
    'explanation':f'Fresh complete degree8 Hecke sieve at3: all66 partitions, '
    f'{expected[n]} tests, zero survivors. Includes every cusp complement and '
    'repetition; independent Frobenius counts agree. All-cuspidal branch uses '
    'rankQ0 DEHMZ3.1(2), NVTable1 gonality>=10, odd-prime unramified torsion injectivity.',
    'evidence':[{'path':f,'sha256':h} for f,h in sorted(ev.items())],
})
(B/'certified_results.json').write_text(json.dumps(ledger,indent=2)+'\n')
print(f'PROMOTED{n}',len(ev))
