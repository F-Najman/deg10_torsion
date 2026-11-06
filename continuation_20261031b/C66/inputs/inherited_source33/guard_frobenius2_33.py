import os,subprocess,time,json,resource,signal
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path.cwd()
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(1800,1800));resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));os.setsid()
cmd=['nice','-n','10','/usr/local/bin/magma','-b',str(p/'frobenius2_verify33.m')]
with (p/'frobenius2_verify33.log').open('w') as log:
 t=time.time();pr=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
 a=Path('/proc',str(pr.pid),'stat').read_text();a=a[a.rfind(')')+2:].split()
 info={'pid':pr.pid,'startticks':int(a[19]),'cmd':cmd,'CPU_budget_seconds':1800,'wall_budget_seconds':7200,'memory_budget_bytes':3*1024**3,'native_memory_budget_bytes':2*10**9,'nice':10}
 (p/'frobenius2_verify33.status.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
 rr=subprocess.run(['python3',str(root/'work/degree9/register_resource_job.py'),str(pr.pid),'50',str(3*1024**3),'C66_source33_independent_Frobenius2'],capture_output=True,text=True);print(rr.stdout,rr.stderr,flush=True);assert rr.returncode==0
 try:rc=pr.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  os.killpg(pr.pid,signal.SIGTERM);rc=pr.wait(timeout=20)
 info.update(actual_exit_code=rc,elapsed_wall_seconds=time.time()-t,completed_marker='NONIC33_FRESH_FROBENIUS2_VERIFIED' in (p/'frobenius2_verify33.log').read_text())
 (p/'frobenius2_verify33.exit.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
