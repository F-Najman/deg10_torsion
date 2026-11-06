import os,subprocess,time,json,resource,signal
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path.cwd()
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(1800,1800));resource.setrlimit(resource.RLIMIT_AS,(6000000000,6000000000));os.setsid()
cmd=['nice','-n','10','/usr/local/bin/magma','-b',str(p/'elliptic_hecke5.m')]
with (p/'elliptic_hecke5.log').open('w') as log:
 t=time.time();pr=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
 a=Path('/proc',str(pr.pid),'stat').read_text();a=a[a.rfind(')')+2:].split()
 info={'pid':pr.pid,'startticks':int(a[19]),'cmd':cmd,'CPU_budget_seconds':1800,'wall_budget_seconds':7200,'memory_budget_bytes':6000000000,'native_memory_budget_bytes':4000000000,'nice':10}
 (p/'elliptic_hecke5.status.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
 rr=subprocess.run(['python3',str(root/'work/degree9/register_resource_job.py'),str(pr.pid),'50',str(6000000000),'C93_actual_level_Hecke5_jets'],capture_output=True,text=True);print(rr.stdout,rr.stderr,flush=True);assert rr.returncode==0
 try:rc=pr.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  os.killpg(pr.pid,signal.SIGTERM);rc=pr.wait(timeout=20)
 info.update(actual_exit_code=rc,elapsed_wall_seconds=time.time()-t,completed_marker='NONIC93_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE' in (p/'elliptic_hecke5.log').read_text())
 (p/'elliptic_hecke5.exit.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
