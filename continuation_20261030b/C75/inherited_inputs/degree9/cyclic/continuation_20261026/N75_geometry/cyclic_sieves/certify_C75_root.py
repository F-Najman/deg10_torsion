"""Read-only SHA verification and isolated exact binary/elliptic replays."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile,shutil
P=Path(__file__).resolve().parent;B=P.parent
j=json.loads((P/'C75_root_manifest.json').read_text())
for name,digest in j['files'].items():
 p=B/name;assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
print('C75_MANIFEST_SHA_PASS',len(j['files']),flush=True)
with tempfile.TemporaryDirectory(prefix='degree8_C75_readonly_') as d:
 root=Path(d)
 for name in j['files']:
  p=root/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(B/name,p)
 for name,marker in [('replay_C75_good7_root.py','C75_GOOD7_INDEPENDENT_REPLAY_PASS'),('replay_C75_allbad_root.py','C75_ALLBAD_INDEPENDENT_REPLAY_PASS')]:
  r=subprocess.run([sys.executable,str(root/'cyclic_sieves'/name)],cwd=root,text=True,capture_output=True,timeout=120)
  assert r.returncode==0 and marker in r.stdout,(name,r.stdout,r.stderr)
  print('C75_ISOLATED_REPLAY_PASS',name,flush=True)
print('C75_COMPLETE_READONLY_CERTIFICATE_PASS')
