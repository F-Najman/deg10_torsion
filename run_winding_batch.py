from pathlib import Path
import concurrent.futures,subprocess,json
p=Path(__file__).resolve().parent
ss=sorted((p/'computations').glob('winding*.m'),key=lambda s:int(s.stem[7:]))
def run(s):
 r=subprocess.run(['/usr/bin/python3',str(p/'run_job.py'),str(s),'--seconds','600'],capture_output=True,text=True)
 print(s.stem,'supervisor_exit',r.returncode,flush=True)
 return {'script':str(s.relative_to(p)),'supervisor_exit':r.returncode}
print('WINDING_BATCH_START',len(ss),'workers',8,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool: rr=list(pool.map(run,ss))
(p/'jobs/winding_batch_result.json').write_text(json.dumps(rr,indent=2)+'\n')
print('WINDING_BATCH_COMPLETE',flush=True)
