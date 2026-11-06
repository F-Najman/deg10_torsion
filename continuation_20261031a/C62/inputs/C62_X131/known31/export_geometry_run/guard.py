import os,subprocess,time,json,resource,signal
from pathlib import Path
p=Path(__file__).resolve().parent
root=Path.cwd()
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(120,120))
 resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
 os.setsid()
cmd=['nice','-n','10','/usr/local/bin/magma','-b',str(root/'work/degree9/cyclic/continuation_20261026/C62_X131/export_geometry.m')]
with (p/'run.log').open('w') as log:
 t=time.time();pr=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
 stat=Path('/proc')/str(pr.pid)/'stat';a=stat.read_text();fields=a[a.rfind(')')+2:].split()
 info={'pid':pr.pid,'startticks':int(fields[19]),'cmd':cmd,'CPU_budget_seconds':120,'wall_budget_seconds':600,'memory_budget_bytes':3*1024**3,'nice':10}
 (p/'native_status.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
 try: rc=pr.wait(timeout=600)
 except subprocess.TimeoutExpired:
  os.killpg(pr.pid,signal.SIGTERM);rc=pr.wait(timeout=20)
 info.update(actual_exit_code=rc,elapsed_wall_seconds=time.time()-t,completed_marker='NONIC62_GEOMETRY_COUNTS_COMPLETE' in (p/'run.log').read_text())
 (p/'native_status.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
