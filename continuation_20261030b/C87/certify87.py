from pathlib import Path
import json,hashlib
B=Path(__file__).resolve().parent;H=B.parents[1]
S=json.loads((B/'final_certificate_summary.json').read_text())
for p,sh in S['file_hashes'].items():assert hashlib.sha256((H/p).read_bytes()).hexdigest()==sh,p
for r in S['completed_jobs']:
 p=H/'jobs'/f"{r['job']}.status.json";s=json.loads(p.read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0
 assert r['marker'] in (H/'jobs'/f"{r['job']}.log").read_text()
 assert hashlib.sha256((H/s['script']).read_bytes()).hexdigest()==s['script_sha256']
assert not [p for p in (H/'jobs').glob('b87*.status.json') if json.loads(p.read_text())['state']!='EXITED']
print('DECIC87_CERTIFICATE_PASS',len(S['file_hashes']),len(S['completed_jobs']))
