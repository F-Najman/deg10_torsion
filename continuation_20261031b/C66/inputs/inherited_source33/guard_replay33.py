import os,subprocess,time,json,resource,signal
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path.cwd()
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(1800,1800));resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));os.setsid()
cmd=['nice','-n','10','/opt/sage/envs/sage/bin/sage','-python',str(p/'replay33.py')]
with (p/'replay33.log').open('w') as log:
 t=time.time();pr=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
 a=Path('/proc',str(pr.pid),'stat').read_text();a=a[a.rfind(')')+2:].split()
 info={'pid':pr.pid,'startticks':int(a[19]),'cmd':cmd,'CPU_budget_seconds':1800,'wall_budget_seconds':7200,'memory_budget_bytes':6*1024**3,'native_memory_budget_bytes':4*10**9,'nice':10}
 (p/'replay33.status.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
 try:rc=pr.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  os.killpg(pr.pid,signal.SIGTERM);rc=pr.wait(timeout=20)
 info.update(actual_exit_code=rc,elapsed_wall_seconds=time.time()-t,completed_marker='NONIC33_INDEPENDENT_CLASS_AND_REPRESENTATIVE_REPLAY_PASS' in (p/'replay33.log').read_text())
 (p/'replay33.exit.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
