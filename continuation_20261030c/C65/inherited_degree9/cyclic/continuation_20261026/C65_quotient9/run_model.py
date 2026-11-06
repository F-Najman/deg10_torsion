from pathlib import Path
import subprocess,resource,os,time,json,sys,signal
P=Path(__file__).resolve().parent;ROOT=P.parents[4]
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(1770,1800));resource.setrlimit(resource.RLIMIT_AS,(6*10**9,6*10**9));os.nice(10)
t=time.monotonic();timed=False
with (P/'model65.log').open('w') as log:
 proc=subprocess.Popen(['/usr/local/bin/magma','-b',str(P/'model65.m')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits,start_new_session=True)
 (P/'model_launch.json').write_text(json.dumps(dict(supervisor_pid=os.getpid(),launcher_pid=proc.pid,cpu_soft=1770,cpu_hard=1800,memory_budget=4*10**9,address_space_limit=6*10**9,wall_limit=7200,nice=10),indent=2)+'\n')
 time.sleep(.1)
 q=Path('/proc')/str(proc.pid)
 if proc.poll() is None and q.exists():
  st=(q/'stat').read_text().split();j=dict(native_pid=proc.pid,start_ticks=int(st[21]),command=(q/'cmdline').read_bytes().replace(b'\x00',b' ').decode(),memory_budget=4*10**9,cpu_hard=1800)
  (P/'model_native_identity.json').write_text(json.dumps(j,indent=2)+'\n')
  rr=subprocess.run(['python3',str(ROOT/'work/degree9/register_resource_job.py'),str(proc.pid),'60','4000000000','C65_quotient9_model'],cwd=ROOT)
  if rr.returncode and proc.poll() is None:raise RuntimeError('Native registration failed for a live worker')
 try:rc=proc.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  timed=True;os.killpg(proc.pid,signal.SIGTERM)
  try:rc=proc.wait(timeout=15)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
s=(P/'model65.log').read_text();out=dict(actual_exit_code=rc,wall_seconds=time.monotonic()-t,wall_limit=7200,wall_timeout=timed,cpu_hard=1800,memory_budget=4*10**9,complete=rc==0 and not timed and 'NONIC65_MODEL_CERTIFIED' in s)
(P/'model_completion.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);sys.exit(0 if out['complete'] else 1)
