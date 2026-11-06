#!/usr/bin/env python3
"""Run a single-thread Magma task with bounded lifetime; never pause workers."""
from pathlib import Path
import argparse,subprocess,time,json,os,signal,hashlib
p=argparse.ArgumentParser();p.add_argument('script');p.add_argument('--seconds',type=int,default=240);p.add_argument('--engine',choices=['magma','sage'],default='magma');args=p.parse_args()
root=Path(__file__).resolve().parents[2];home=Path(__file__).resolve().parent
script=Path(args.script).resolve();assert home in script.parents
stem=script.stem;log=home/'jobs'/(stem+'.log');status=home/'jobs'/(stem+'.status.json')
assert not status.exists(),status
cmd=['/usr/local/bin/magma','-b',str(script)] if args.engine=='magma' else ['/opt/sage/envs/sage/bin/python3',str(script)]
env=os.environ.copy();env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
t0=time.time()
with log.open('x') as fp:
 proc=subprocess.Popen(cmd,cwd=root,stdout=fp,stderr=subprocess.STDOUT,env=env,start_new_session=True)
 rec={'purpose':stem,'script':str(script.relative_to(home)),'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'command':cmd,'pid':proc.pid,'owner_uid':os.getuid(),'start_time_unix':t0,'state':'RUNNING','time_limit_seconds':args.seconds,'log':str(log.relative_to(home))}
 status.write_text(json.dumps(rec,indent=2)+'\n')
 try:
  code=proc.wait(timeout=args.seconds);reason=None
 except subprocess.TimeoutExpired:
  # The direct child cannot have its PID reused before wait; retire the process group.
  os.killpg(proc.pid,signal.SIGTERM)
  try:code=proc.wait(timeout=5)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
  reason='time limit; computation is not a completed certificate'
 rec.update({'state':'EXITED','actual_exit_code':code,'termination_reason':reason,'elapsed_seconds':time.time()-t0})
 status.write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,sort_keys=True))
raise SystemExit(0 if code==0 else 1)
