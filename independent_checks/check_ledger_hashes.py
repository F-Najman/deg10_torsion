#!/usr/bin/env python3
"""Check the SHA-256 hashes recorded in certified_results.json against the files of this repository.

Written for this repository (not part of the original computational record).  Every evidence file of the
ledger is either in the repository, or in the release archive deg10_torsion_large_files.tar.zst (listed in
independent_checks/large_files_manifest.tsv; unpack it at the repository root), or a copy of a published
paper that is not redistributed (listed in independent_checks/literature_not_distributed.txt).
Exit status 0 means that every file that is present has the recorded hash and that every missing file is
accounted for by one of the two lists.
"""
import hashlib, json, sys
from pathlib import Path

H = Path(__file__).resolve().parents[1]
ledger = json.loads((H / 'certified_results.json').read_text())
large = {l.split('\t')[0] for l in (H / 'independent_checks/large_files_manifest.tsv').read_text().splitlines() if l}
lit = set((H / 'independent_checks/literature_not_distributed.txt').read_text().splitlines())

evidence = {}
def walk(x):
    if isinstance(x, dict):
        if 'path' in x and 'sha256' in x:
            evidence.setdefault(x['path'], set()).add(x['sha256'])
        for v in x.values():
            walk(v)
    elif isinstance(x, list):
        for v in x:
            walk(v)
walk(ledger)

ok = bad = in_archive = literature = unaccounted = 0
for path, hashes in sorted(evidence.items()):
    f = H / path
    if f.is_file():
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        if h in hashes and len(hashes) == 1:
            ok += 1
        else:
            bad += 1
            print('HASH_MISMATCH', path)
    elif path in large:
        in_archive += 1
    elif path in lit:
        literature += 1
    else:
        unaccounted += 1
        print('MISSING', path)
print(f'LEDGER_HASHES evidence={len(evidence)} verified={ok} mismatched={bad} '
      f'in_release_archive_not_unpacked={in_archive} literature_not_distributed={literature} unaccounted={unaccounted}')
sys.exit(1 if bad or unaccounted else 0)
