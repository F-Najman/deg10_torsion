#!/usr/bin/env python3
from pathlib import Path
import concurrent.futures,subprocess,time,json,os,resource,datetime,hashlib,signal
H=Path(__file__).resolve().parent;ROOT=H.parents[2];DEG=H.parent

def prep():
 os.nice(10);resource.setrlimit(resource.RLIMIT_CPU,(1800,1810));resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
def run(n):
 p=H/f'N{n}_source_units';script=p/'products.m';log=p/'products.log';status=p/'products_status.json';assert not status.exists(),str(status)
 args=['/usr/local/bin/magma','-b',f'N:={n}',f'OutDir:={p.relative_to(ROOT)}',str(script.relative_to(ROOT))]
 started=time.monotonic();ident=None
 with log.open('w')as out:
  child=subprocess.Popen(args,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=prep)
  for _ in range(100):
   if child.poll()is not None:break
   try:
    raw=Path(f'/proc/{child.pid}/stat').read_text();a=raw[raw.rfind(')')+2:].split();cmd=Path(f'/proc/{child.pid}/cmdline').read_bytes()
    if b'magma' in cmd:
     ident={'pid':child.pid,'start_ticks':int(a[19])}
     r=subprocess.run(['python3',str(DEG/'register_resource_job.py'),str(child.pid),'60',str(4*10**9),f'source_units_N{n}'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
     if r.returncode==0 or child.poll()is not None:break
   except FileNotFoundError:break
   time.sleep(.03)
  state={'N':n,'status':'RUNNING','identity':ident,'args':args,'nice':10,'native_memory_bytes':4*10**9,'AS_bytes':6*1024**3,'CPU_limit_seconds':1800,'hard_CPU_seconds':1810,'wall_limit_seconds':7200,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
  status.write_text(json.dumps(state,indent=2)+'\n')
  try:code=child.wait(timeout=7200)
  except subprocess.TimeoutExpired:
   raw=Path(f'/proc/{child.pid}/stat').read_text();a=raw[raw.rfind(')')+2:].split();assert ident and int(a[19])==ident['start_ticks'];child.terminate()
   try:code=child.wait(timeout=10)
   except subprocess.TimeoutExpired:child.kill();code=child.wait()
 text=log.read_text();ok=code==0 and f'SOURCE_PRODUCT_COMPLETE level={n}'in text and all(t not in text for t in ['Runtime error','User error','Assertion failed'])
 state.update(status='COMPLETE'if ok else'FAILED_OR_GUARDED',actual_exit_code=code,completion_marker=ok,elapsed_wall_seconds=time.monotonic()-started,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest());status.write_text(json.dumps(state,indent=2)+'\n');print(n,state['status'],flush=True);return state
with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool: results=list(pool.map(run,[65,93,91,117]))
(H/'source_products_pool_complete.json').write_text(json.dumps({'max_concurrent':4,'memory_budget_bytes':16*10**9,'jobs':results,'all_complete':all(r['completion_marker']for r in results)},indent=2)+'\n')
