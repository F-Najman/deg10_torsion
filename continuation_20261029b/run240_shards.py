from pathlib import Path
import concurrent.futures,subprocess,json,time,os
H=Path(__file__).resolve().parents[1];ROOT=H.parents[1];P=Path(__file__).resolve().parent
records=[]
def one(k):
 stem='late240_picard_s'+str(k)
 t=time.monotonic()
 while True:
  ls=subprocess.check_output(['ps','-eo','stat=,comm='],text=True).splitlines()
  n=sum(len(s.split())==2 and s.split()[1]=='magma.exe' and not s.split()[0].startswith(('T','Z')) for s in ls)
  if n<29:break
  if time.monotonic()-t>180:raise RuntimeError('No safe core capacity')
  time.sleep(3)
 r=subprocess.run([os.sys.executable,str(H/'run_job.py'),str(H/'computations'/(stem+'.m')),'--seconds','1200'],cwd=ROOT,stdout=subprocess.PIPE,text=True)
 rec=json.loads((H/'jobs'/(stem+'.status.json')).read_text());records.append(rec);(P/'C240/shard_completion.json').write_text(json.dumps(records,indent=2)+'\n')
 print('C240_PICARD_SHARD_EXIT',k,r.returncode,rec['elapsed_seconds'],flush=True)
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:result=list(ex.map(one,range(8)))
(P/'C240/shard_completion.json').write_text(json.dumps(result,indent=2)+'\n')
print('ALL_EIGHT_SHARDS_EXITED',flush=True)
