"""Read-only finite-fibre binding supplement; no Magma or Sage."""
from pathlib import Path
import json,hashlib,subprocess,sys
p=Path(__file__).resolve().parent;J=lambda q:json.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();m=J(p/'manifest117.json')
for r,v in m['files'].items():
 q=p/r;assert q.stat().st_size==v['bytes']and sha(q)==v['sha256'],r
h=p.parent/'C78_known39';assert sha(h/'manifest39.json')==m['helper_manifest_sha256'];hm=J(h/'manifest39.json')
for r,v in hm['files'].items():assert sha(h/r)==v['sha256']
for cmd in[[sys.executable,str(p/'census_root/census39_F512.py'),'--readonly'],[sys.executable,str(p/'bind117.py'),'--readonly']]:
 r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0,(r.stdout,r.stderr);print(r.stdout.strip())
print('C117_KNOWN39_FINITE_FIBRE_BINDING_CERTIFICATE_PASS')
