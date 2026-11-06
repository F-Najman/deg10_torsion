#!/usr/bin/env python3
"""Bind the audited fresh octic C45 computation and its independent checks."""
from pathlib import Path
import hashlib
import json
import re

B = Path(__file__).resolve().parents[1]
files = [
    'cyclic_sieves/exclusion45.tex', 'cyclic_sieves/dyadic_criterion.tex',
    'cyclic_sieves/odd_prime_criterion.tex',
    'cyclic_sieves/generic_cartier_completion_export.m',
    'cyclic_sieves/N45_completion.log', 'cyclic_sieves/N45_completion_jets.m',
    'cyclic_sieves/N45_completion_projection_replay.log',
    'cyclic_sieves/N45_completion_rr_crosscheck.log',
    'cyclic_sieves/N45_completion_fi_replay.log',
    'cyclic_sieves/replay_cartier_projection.py',
    'cyclic_sieves/compare_rr_to_completions.py',
    'cyclic_sieves/replay_fi_jets.py',
    'cyclic_sieves/N45_reference_RR_geometry_D7.m',
    'cyclic_sieves/N45_reference_geometry_provenance.json',
    'cyclic_sieves/hecke_sieve_smart.m', 'cyclic_sieves/N45_p2_good_only.log',
    'cyclic_sieves/audit_partial_sieve.py',
    'cyclic_sieves/N45_good_reduction_audit.json',
    'cyclic_sieves/N45_good_reduction_audit.log',
    'cyclic_sieves/independent_counts.m', 'cyclic_sieves/independent_counts.log',
    'small_cyclic/audits/audit45_completions.py',
    'small_cyclic/audits/audit45_completions.log',
    'small_cyclic/audits/audit45_theory.txt',
    'small_cyclic/C34/inputs/DEHMZ.tex', 'cyclic_sieves/certify45.py',
]
for name in ('N45_completion.log', 'N45_p2_good_only.log'):
    log = (B/'cyclic_sieves'/name).read_text()
    assert not re.search(r'Runtime error|Assertion failed|failed memory', log)
good = json.loads((B/'cyclic_sieves/N45_good_reduction_audit.json').read_text())
assert good['noncuspidal_reduction_branch_complete'] and good['tested_divisors'] == 5
fi = (B/'cyclic_sieves/N45_completion_fi_replay.log').read_text()
assert '"enumerated": 148945, "exceptions": 0, "complete": true' in fi
ev = {f: hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files}
ledger = json.loads((B/'certified_results.json').read_text())
reference = next(r for r in ledger['exclusions'] if r['group'] == [1,38])
for r in reference['evidence']:
    if r['path'].startswith('shared/'):
        assert hashlib.sha256((B/r['path']).read_bytes()).hexdigest() == r['sha256']
        ev[r['path']] = r['sha256']
manifest = {
    'degree': 8, 'level': 45, 'status': 'complete exclusion',
    'source_genus': 41, 'rational_jacobian_rank': 0, 'prime': 2,
    'noncuspidal_closed_counts': [0,0,0,0,0,0,0,60],
    'noncuspidal_diamond_representatives': 5, 'cartier_difference_dimension': 39,
    'local_restriction_blocks': 124, 'cuspidal_degree8_divisors': 148945,
    'formal_immersion_exceptions': 0, 'independent_RR_kernel_equalities': 103,
    'independent_agent_audit': 'small_cyclic/audits/audit45_theory.txt',
    'fresh_completion_cpu_seconds': 2280.44, 'fresh_good_branch_cpu_seconds': 662.86,
    'evidence': [{'path':f,'sha256':h} for f,h in sorted(ev.items())],
}
mp = B/'cyclic_sieves/N45_manifest.json'
mp.write_text(json.dumps(manifest,indent=2)+'\n')
ev[str(mp.relative_to(B))] = hashlib.sha256(mp.read_bytes()).hexdigest()
assert not any(r['group'] == [1,45] for r in ledger['exclusions'])
ledger['exclusions'].append({
    'degree':8, 'group':[1,45], 'audited':True,
    'proof':'cyclic_sieves/exclusion45.tex',
    'explanation': 'Fresh degree8 dyadic formal-immersion certificate: all148945 effective '
    'cuspidal divisors rank8 on dimension39 Cartier-difference image. Independent '
    'projection and all103 exact RR kernel cross-checks; agent independently replayed '
    'all124 blocks and148945 ranks. Fresh Hecke tests exclude all5 octic source '
    'orbits, independent count60, no lower noncuspidal places. RankQ0 DEHMZ3.1(2).',
    'evidence':[{'path':f,'sha256':h} for f,h in sorted(ev.items())],
})
(B/'certified_results.json').write_text(json.dumps(ledger,indent=2)+'\n')
print('PROMOTED45', len(ev))
