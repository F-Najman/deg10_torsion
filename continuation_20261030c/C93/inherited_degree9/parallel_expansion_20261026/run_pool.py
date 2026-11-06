#!/usr/bin/env python3
"""Ten independent single-core census workers, guarded by total-host headroom."""
from pathlib import Path
import concurrent.futures,subprocess,time,json,os,signal,hashlib,datetime
H=Path(__file__).resolve().parent;ROOT=H.parents[2];SCRIPT=H/'odd_census.m';DEG=H.parent

def run(pair):
 n,d=pair;p=H/f'N{n}_f{d}';p.mkdir(exist_ok=True);log=p/'run.log';out=p/'records.m';args=['/usr/local/bin/magma','-b',f'N:={n}',f'FieldDegree:={d}',f'OutFile:={out}',str(SCRIPT)]
 shell='ulimit -t 125\nulimit -v 2097152\nexec nice -n 10 '+ ' '.join(args)
 start=time.monotonic();ident=None
 with log.open('w')as f:
  child=subprocess.Popen(['bash','-c',shell],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  for _ in range(50):
   if child.poll()is not None:break
   try:
    raw=Path(f'/proc/{child.pid}/stat').read_text();a=raw[raw.rfind(')')+2:].split();cmd=Path(f'/proc/{child.pid}/cmdline').read_bytes()
    if b'magma' in cmd:
     ident={'pid':child.pid,'startticks':int(a[19])};reg=subprocess.run(['python3',str(DEG/'register_resource_job.py'),str(child.pid),'60','1000000000',f'parallel_odd_N{n}_f{d}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
     if reg.returncode==0 or child.poll()is not None:break
   except FileNotFoundError:break
   time.sleep(.03)
  state={'N':n,'field_degree':d,'command':args,'identity':ident,'status':'RUNNING','CPU_limit_seconds':120,'kernel_CPU_backstop_seconds':125,'wall_limit_seconds':600,'native_memory_bytes':10**9,'AS_bytes':2*1024**3,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(p/'status.json').write_text(json.dumps(state,indent=2)+'\n')
  try:code=child.wait(timeout=600)
  except subprocess.TimeoutExpired:
   if ident:
    raw=Path(f'/proc/{child.pid}/stat').read_text();a=raw[raw.rfind(')')+2:].split();assert int(a[19])==ident['startticks']
   child.kill();code=child.wait()
 s=log.read_text();complete=code==0 and f'NONIC_ODD_CENSUS_COMPLETE N={n} field_degree={d}' in s and all(v not in s for v in ['Runtime error','User error','Assertion failed','Segmentation fault'])
 state.update(status='COMPLETE'if complete else'FAILED_OR_LIMIT',actual_exit_code=code,completed_marker=complete,elapsed_wall_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest());(p/'status.json').write_text(json.dumps(state,indent=2)+'\n');print(n,d,state['status'],flush=True);return state
allpairs=[(n,d)for n in[57,65,85,91,93,117]for d in range(6,10)]
pairs=[(n,d)for n,d in allpairs if not (H/f'N{n}_f{d}'/'status.json').exists() or not json.loads((H/f'N{n}_f{d}'/'status.json').read_text()).get('completed_marker')]
for n,d in pairs:
 p=H/f'N{n}_f{d}'
 for name in ['run.log','records.m']:
  q=p/name
  if q.exists():q.rename(p/('initial_unmonitored_'+name))
with concurrent.futures.ThreadPoolExecutor(max_workers=10)as pool:results=list(pool.map(run,pairs))
results=[json.loads((H/f'N{n}_f{d}'/'status.json').read_text())for n,d in allpairs]
(H/'pool_complete.json').write_text(json.dumps({'jobs':24,'max_workers':10,'results':results,'all_complete':all(s['completed_marker']for s in results)},indent=2)+'\n')
