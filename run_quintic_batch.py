from pathlib import Path
import concurrent.futures,subprocess,json
p=Path(__file__).resolve().parent
ss=sorted((p/'computations').glob('quintic_*.m'))
def run(s):
 r=subprocess.run(['/usr/bin/python3',str(p/'run_job.py'),str(s),'--seconds','180'],capture_output=True,text=True)
 print(s.stem,'supervisor_exit',r.returncode,flush=True)
 return {'script':str(s.relative_to(p)),'supervisor_exit':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: rr=list(pool.map(run,ss))
(p/'jobs/quintic_batch_result.json').write_text(json.dumps(rr,indent=2)+'\n')
print('QUINTIC_BATCH_COMPLETE',flush=True)
