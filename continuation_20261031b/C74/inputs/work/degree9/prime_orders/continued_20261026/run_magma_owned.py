#!/usr/bin/env python3
"""Own a Magma worker, record actual exit, enforce a finitewall/idle guard.
Usage: run_magma_owned.py output_prefix cap_seconds [magma arguments...]
All output files are under the new continuation directory. Only the exact
Popen child owned by this runner can be signalled. A stopped/silent but
CPU-active RR computation is not an idle failure.
"""
from pathlib import Path
import subprocess,sys,time,os,json,datetime,signal
prefix=Path(sys.argv[1]).resolve();cap=int(sys.argv[2]);args=sys.argv[3:]
assert prefix.is_relative_to(Path(__file__).resolve().parent)
cmd=['/usr/local/bin/magma','-b']+args
log=prefix.with_suffix('.log');status=prefix.with_suffix('.status.json');ex=prefix.with_suffix('.exit.json')
assert not log.exists() and not status.exists() and not ex.exists()
start=time.monotonic();beg=datetime.datetime.now(datetime.timezone.utc).isoformat();term=None
with log.open('wb')as f:
 p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
 while p.poll()is None:
  elapsed=time.monotonic()-start
  proc=Path(f'/proc/{p.pid}')
  try:
   raw=(proc/'stat').read_text();s=raw[raw.rindex(')')+2:].split();ticks=int(s[11])+int(s[12]);st=int(s[19]);state=s[0]
   pcmd=[a.decode()for a in(proc/'cmdline').read_bytes().split(b'\0')if a]
  except FileNotFoundError:continue
  data={'started_utc':beg,'worker_pid':p.pid,'worker_start_ticks':st,'command':cmd,'actual_cmdline':pcmd,'wall_seconds':elapsed,'cpu_seconds':ticks/os.sysconf('SC_CLK_TCK'),'state':state,'wall_cap_seconds':cap,'status':'RUNNING; not accepted evidence'}
  status.write_text(json.dumps(data,indent=2)+'\n')
  if elapsed>=cap:
   assert p.poll()is None and '/magma' in pcmd[0] and pcmd[1:]==['-b']+args,pcmd
   term='finite wall guard expired';p.terminate()
   try:p.wait(timeout=30)
   except subprocess.TimeoutExpired:p.kill();p.wait()
   break
  time.sleep(10)
 rc=p.wait()
data.update(status='COMPLETED'if rc==0 else'FAILED OR GUARDED; NO EXCLUSION INFERRED',native_exit_code=rc,wall_seconds=time.monotonic()-start,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),termination_reason=term)
ex.write_text(json.dumps(data,indent=2)+'\n');status.write_text(json.dumps(data,indent=2)+'\n')
print('OWNED_MAGMA_ACTUAL_EXIT',json.dumps(data),flush=True)
sys.exit(0 if rc==0 else 1)
