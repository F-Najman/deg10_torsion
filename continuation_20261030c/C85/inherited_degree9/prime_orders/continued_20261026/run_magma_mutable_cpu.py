#!/usr/bin/env python3
"""Owned actual-exit runner with a finite, root-coordinated CPU budget.
Usage: runner PREFIX WALL_SECONDS INITIAL_CPU_SECONDS [Magma arguments...]
Production: nice10, hard CPU7200, initial effective CPU<=1800, wall<=7200.
Only the root agent changes PREFIX.cpu_budget.json after reviewing progress.
Each applied extension is durably recorded; frozen runners are untouched.
"""
from pathlib import Path
import os,sys,json,time,datetime,subprocess,signal,resource,hashlib
BASE=Path(__file__).resolve().parent

def atomic(path,data):
 t=path.with_name(path.name+'.tmp')
 with t.open('w')as f:json.dump(data,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,path)

def snapshot(pid):
 p=Path('/proc')/str(pid);s=(p/'stat').read_text();a=s[s.rfind(')')+2:].split()
 return {'pid':pid,'startticks':int(a[19]),'ppid':int(a[1]),'state':a[0],
 'cpu':(int(a[11])+int(a[12]))/os.sysconf('SC_CLK_TCK'),
 'uid':int(next(v for v in(p/'status').read_text().splitlines()if v.startswith('Uid:')).split()[1]),
 'cmd':[a.decode()for a in(p/'cmdline').read_bytes().split(b'\0')if a]}

def run_owned(prefix,cmd,wall,cpu,hard=7200,poll=1.0,magma_args=None,hook=None):
 prefix=Path(prefix).resolve();assert prefix.is_relative_to(BASE)
 assert 0<cpu<=1800 and cpu<=hard<=7200 and 0<wall<=7200
 assert os.getpriority(os.PRIO_PROCESS,0)>=10
 paths={k:prefix.with_suffix(s)for k,s in [('log','.log'),('status','.status.json'),('exit','.exit.json'),('budget','.cpu_budget.json'),('history','.cpu_budget_history.jsonl')]}
 assert not any(p.exists()for p in paths.values())
 def history(event,**kw):
  rec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'event':event,**kw}
  with paths['history'].open('a')as f:f.write(json.dumps(rec,sort_keys=True)+'\n');f.flush();os.fsync(f.fileno())
 def hardcap():resource.setrlimit(resource.RLIMIT_CPU,(int(hard),int(hard)))
 def terminate_owned(p,ident):
  if p.poll()is not None:return
  fresh=snapshot(p.pid);assert fresh['startticks']==ident['startticks']and fresh['ppid']==os.getpid()and fresh['uid']==os.getuid()
  if magma_args is None:assert fresh['cmd']==cmd
  else:assert '/magma' in fresh['cmd'][0]and fresh['cmd'][1:]==['-b']+magma_args
  p.terminate()
  if fresh['state']in('T','t'):os.kill(p.pid,signal.SIGCONT)
  try:p.wait(timeout=30)
  except subprocess.TimeoutExpired:
   fresh=snapshot(p.pid);assert fresh['startticks']==ident['startticks']and fresh['ppid']==os.getpid()and fresh['uid']==os.getuid();p.kill();p.wait()
 began=datetime.datetime.now(datetime.timezone.utc).isoformat();start=time.monotonic();term=None;limit=cpu;last_bad=None
 owner=snapshot(os.getpid());data={}
 with paths['log'].open('wb')as log:
  p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=hardcap)
  ident=None
  try:
   ident=snapshot(p.pid);assert ident['ppid']==os.getpid()and ident['uid']==os.getuid()
   budget={'native_pid':p.pid,'native_startticks':ident['startticks'],'owner_pid':owner['pid'],'owner_startticks':owner['startticks'],'cpu_limit_seconds':cpu,'hard_cpu_backstop_seconds':hard,'wall_limit_seconds':wall,'requested_by':'root','approval_note':'Initial allocated finite budget'}
   atomic(paths['budget'],budget);history('INITIAL_BUDGET',**budget)
   while p.poll()is None:
    z=snapshot(p.pid);assert z['startticks']==ident['startticks']and z['ppid']==os.getpid()and z['uid']==os.getuid()
    elapsed=time.monotonic()-start
    if hook is not None:hook(paths['budget'],z,elapsed)
    try:
     raw=paths['budget'].read_bytes();b=json.loads(raw);new=b['cpu_limit_seconds']
     assert b['native_pid']==p.pid and b['native_startticks']==ident['startticks']and b['owner_pid']==owner['pid']and b['owner_startticks']==owner['startticks']
     assert isinstance(new,(int,float))and 0<new<=hard and b['hard_cpu_backstop_seconds']==hard and b['wall_limit_seconds']==wall and b['requested_by']=='root'and isinstance(b['approval_note'],str)and b['approval_note'].strip()
     if new!=limit:history('APPLIED_ROOT_BUDGET_CHANGE',old_cpu_limit_seconds=limit,new_cpu_limit_seconds=new,observed_cpu_seconds=z['cpu'],budget_sha256=hashlib.sha256(raw).hexdigest(),approval_note=b['approval_note']);limit=new
    except (OSError,ValueError,KeyError,TypeError,AssertionError)as e:
     bad=repr(e)
     if bad!=last_bad:history('INVALID_BUDGET_IGNORED',error=bad,retained_cpu_limit_seconds=limit);last_bad=bad
    data={'started_utc':began,'worker_pid':p.pid,'worker_start_ticks':ident['startticks'],'owner_pid':owner['pid'],'owner_startticks':owner['startticks'],'command':cmd,'actual_cmdline':z['cmd'],'wall_seconds':elapsed,'cpu_seconds':z['cpu'],'state':z['state'],'wall_cap_seconds':wall,'cpu_budget_seconds':limit,'hard_cpu_backstop_seconds':hard,'status':'RUNNING; not accepted evidence'}
    atomic(paths['status'],data)
    reason='finite CPU budget expired'if z['cpu']>=limit else('finite wall guard expired'if elapsed>=wall else None)
    if reason and p.poll()is None:
     term=reason;history('GUARD_TERMINATING_OWN_CHILD',reason=reason,native_pid=p.pid,native_startticks=ident['startticks']);terminate_owned(p,ident)
     break
    time.sleep(poll)
   rc=p.wait()
  except BaseException as err:
   term='unexpected runner error: '+repr(err)
   if p.poll()is None:
    if ident is None:ident=snapshot(p.pid)
    terminate_owned(p,ident)
   rc=p.wait()
   try:history('UNEXPECTED_RUNNER_ERROR',error=repr(err),native_exit_code=rc)
   except OSError:pass

 if ident is None:ident={'startticks':None}
 if not data:data={'started_utc':began,'worker_pid':p.pid,'worker_start_ticks':ident['startticks'],'owner_pid':owner['pid'],'owner_startticks':owner['startticks'],'command':cmd,'wall_cap_seconds':wall,'cpu_budget_seconds':limit,'hard_cpu_backstop_seconds':hard}
 data.update(status='COMPLETED'if rc==0 and term is None else'FAILED OR GUARDED; NO EXCLUSION INFERRED',native_exit_code=rc,wall_seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),termination_reason=term)
 atomic(paths['exit'],data);atomic(paths['status'],data);history('ACTUAL_POPEN_WAIT',native_exit_code=rc,termination_reason=term)
 return data

def selftest():
 root=BASE/'mutable_guard_selftest_20261026_v2';root.mkdir(exist_ok=False)
 code="import time\ns=time.process_time()\nwhile time.process_time()-s<0.35: pass\n"
 changed=[False]
 def extend(path,z,elapsed):
  if z['cpu']>.08 and not changed[0]:
   b=json.loads(path.read_text());b['cpu_limit_seconds']=.6;b['approval_note']='Self-test root-authorized finite extension';atomic(path,b);changed[0]=True
 a=run_owned(root/'extension',[sys.executable,'-c',code],2,.2,hard=2,poll=.01,hook=extend);assert a['native_exit_code']==0 and changed[0]and a['cpu_budget_seconds']==.6
 a=run_owned(root/'cpu',[sys.executable,'-c','while True: pass'],2,.15,hard=2,poll=.01);assert a['native_exit_code']<0 and a['termination_reason']=='finite CPU budget expired'
 a=run_owned(root/'wall',[sys.executable,'-c','import time;time.sleep(3)'],.2,.5,hard=2,poll=.01);assert a['native_exit_code']<0 and a['termination_reason']=='finite wall guard expired'
 def fail_stopped(path,z,elapsed):
  if z['cpu']>.04:
   os.kill(z['pid'],signal.SIGSTOP)
   for _ in range(100):
    if snapshot(z['pid'])['state']=='T':break
    time.sleep(.001)
   assert snapshot(z['pid'])['state']=='T';raise RuntimeError('intentional cleanup self-test')
 a=run_owned(root/'error_stopped',[sys.executable,'-c','while True: pass'],2,.5,hard=2,poll=.01,hook=fail_stopped);assert a['native_exit_code']<0 and 'unexpected runner error'in a['termination_reason']and not Path('/proc',str(a['worker_pid'])).exists()
 print('MUTABLE_OWNED_GUARD_SELFTEST_PASS actual_exit0_after_extension=true finite_cpu_guard=true finite_wall_guard=true own_children_only=true unexpected_error_stopped_child_reaped=true')

if __name__=='__main__':
 if sys.argv[1:]==['--self-test']:selftest()
 else:
  prefix=Path(sys.argv[1]);wall=int(sys.argv[2]);cpu=int(sys.argv[3]);args=sys.argv[4:];assert args
  z=run_owned(prefix,['/usr/local/bin/magma','-b']+args,wall,cpu,magma_args=args)
  print('OWNED_MAGMA_ACTUAL_EXIT',json.dumps(z),flush=True);sys.exit(0 if z['native_exit_code']==0 and z['termination_reason']is None else 1)
