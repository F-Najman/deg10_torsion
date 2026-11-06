from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent;h=p.parents[1];base=h/'parallel_expansion_20261026/N93_source_units';m=json.loads((base/'manifest93.json').read_text());v=json.loads((p/'correction.json').read_text())
assert hashlib.sha256((base/'manifest93.json').read_bytes()).hexdigest()==v['base_manifest_sha256']
for e in m['files']:assert hashlib.sha256((h/e['path']).read_bytes()).hexdigest()==e['sha256']
a=(base/'exclusion93.tex').read_bytes();b=(p/'exclusion93_corrected.tex').read_bytes();assert hashlib.sha256(a).hexdigest()==v['old_proof_sha256']and hashlib.sha256(b).hexdigest()==v['new_proof_sha256']
q=json.loads((base/'canonical_R1_sieve_verified.json').read_text());valid=[r for r in q['results']if r['source']==-1 and r['convention']==2];assert len(valid)==30 and all(r['rank']==261 for r in valid)
assert b'A_5=T_5-1-5' in b and b'B=\\sigma T-T-\\sigma P+P-5\\sigma\\tau P+5\\tau P'in b
print('NONIC93_POINT_PUSH_CORRECTION_PASS all322basehashes=true validB2cases30rank261=true')
