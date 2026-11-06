from pathlib import Path
import concurrent.futures,subprocess,json
p=Path(__file__).resolve().parent
ss=[p/'computations'/s for s in ['winding111.m','winding117.m','winding851.m','winding1147_resume.m']]
def run(s):
 r=subprocess.run(['/usr/bin/python3',str(p/'run_job.py'),str(s),'--seconds','600'],capture_output=True,text=True)
 print(s.stem,'supervisor_exit',r.returncode,flush=True)
 return {'script':str(s.relative_to(p)),'supervisor_exit':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: rr=list(pool.map(run,ss))
(p/'jobs/third_batch_result.json').write_text(json.dumps(rr,indent=2)+'\n')
print('THIRD_BATCH_COMPLETE',flush=True)
