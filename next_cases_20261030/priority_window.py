from pathlib import Path
import psutil,os,time,json,signal
H=Path(__file__).resolve().parent; S=H/'priority_window.json'; release=H/'priority_window.release'
p=psutil.Process(2008682); expected=['python3','work/degree20/run_characters_hybrid24.py']; born=1793365555.01
assert p.username()=='fnajman' and p.cmdline()==expected and abs(p.create_time()-born)<0.01
r={'purpose':'Bounded temporary launcher drain for user-prioritized Phi(10); useful existing workers are preserved','pid':p.pid,'owner':p.username(),'create_time':born,'command':expected,'max_seconds':900,'state':'VERIFIED','started':time.time()}
def save():S.write_text(json.dumps(r,indent=2)+'\n')
def interrupted(sig,frame):raise KeyboardInterrupt
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
save();stopped=False
try:
 p.send_signal(signal.SIGSTOP);stopped=True;r['state']='LAUNCHER_TEMPORARILY_STOPPED';save();print('PHI10_PRIORITY_WINDOW_STARTED; useful children preserved; automatic restore within900seconds',flush=True)
 deadline=time.monotonic()+900
 while not release.exists() and time.monotonic()<deadline:time.sleep(2)
 r['release_reason']='parent finished priority window' if release.exists() else 'bounded lease expired'
finally:
 if stopped:
  try:
   if p.is_running() and abs(p.create_time()-born)<0.01 and p.cmdline()==expected:p.send_signal(signal.SIGCONT);r['restored']=True
   else:r['restored']=False;r['reason']='original launcher no longer exists'
  except psutil.NoSuchProcess:r['restored']=False;r['reason']='launcher exited'
 r['state']='EXITED';r['actual_exit_code']=0;r['finished']=time.time();save();print('PHI10_PRIORITY_WINDOW_RESTORED',r,flush=True)
