from pathlib import Path
import concurrent.futures,subprocess,json,os,time,threading
H=Path(__file__).resolve().parents[1];ROOT=H.parents[1];P=Path(__file__).resolve().parent
rows=[];launch_lock=threading.Lock();record_lock=threading.Lock()
def one(k):
 stem='late240_replay_s'+str(k);t=time.monotonic()
 while True:
  with launch_lock:
   ls=subprocess.check_output(['ps','-eo','stat=,comm='],text=True).splitlines()
   n=sum(len(s.split())==2 and s.split()[1]=='magma.exe' and not s.split()[0].startswith(('T','Z')) for s in ls)
   if n<min(28,(os.cpu_count() or 2)-2):
    proc=subprocess.Popen([os.sys.executable,str(H/'run_job.py'),str(H/'computations'/(stem+'.m')),'--seconds','600'],cwd=ROOT,stdout=subprocess.PIPE,text=True)
    time.sleep(0.25)
    break
  if time.monotonic()-t>300:raise RuntimeError('No safe core capacity; no worker launched')
  time.sleep(3)
 proc.communicate();rec=json.loads((H/'jobs'/(stem+'.status.json')).read_text());rec['host_magma_before_admission']=n
 with record_lock:
  rows.append(rec);(P/'C240/native_replay_completion.json').write_text(json.dumps(rows,indent=2)+'\n')
 print('C240_REPLAY_SHARD_EXIT',k,proc.returncode,rec['elapsed_seconds'],flush=True)
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:result=list(ex.map(one,range(8)))
(P/'C240/native_replay_completion.json').write_text(json.dumps(result,indent=2)+'\n')
