#!/usr/bin/env python3
"""Supervise a task-owned Magma process with actual exit recording and finite CPU/wall limits."""
from pathlib import Path
import argparse,subprocess,time,json,os,signal,datetime,hashlib
p=argparse.ArgumentParser();p.add_argument('tag');p.add_argument('script');p.add_argument('args',nargs='*');p.add_argument('--cpu-limit',type=int,default=7200);p.add_argument('--wall-limit',type=int,default=21600);a=p.parse_args()
h=Path(__file__).resolve().parent;log=h/(a.tag+'.log');out=h/(a.tag+'_status.json');cmd=['/usr/local/bin/magma','-b']+a.args+[a.script];hz=os.sysconf('SC_CLK_TCK');started=time.monotonic()
with log.open('w')as f:
 child=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
 def identity():
  raw=Path(f'/proc/{child.pid}/stat').read_text();s=raw[raw.rindex(')')+2:].split();return int(s[19]),(int(s[11])+int(s[12]))/hz,s[0]
 tick,_,_=identity();state={'cmd':cmd,'pid':child.pid,'start_ticks':tick,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cpu_limit':a.cpu_limit,'wall_limit':a.wall_limit,'status':'RUNNING; no proof conclusion yet'};last=-1;cpu=0
 while child.poll()is None:
  try:
   t,cpu,st=identity();assert t==tick
  except FileNotFoundError:break
  wall=time.monotonic()-started
  if cpu>a.cpu_limit or wall>a.wall_limit:
   state['termination_reason']='finite CPU/wall limit reached';child.terminate()
   try:child.wait(timeout=30)
   except subprocess.TimeoutExpired:
    assert identity()[0]==tick;child.kill()
   break
  if int(wall//30)!=last:
   last=int(wall//30);state.update(elapsed_wall_seconds=wall,cpu_seconds=cpu,last_log_lines=log.read_text().splitlines()[-5:]);out.write_text(json.dumps(state,indent=2)+'\n')
  time.sleep(1)
 code=child.wait();state.update(status='EXITED',native_exit_code=code,elapsed_wall_seconds=time.monotonic()-started,cpu_seconds=cpu,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest());out.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state));raise SystemExit(0 if code==0 else 1)
