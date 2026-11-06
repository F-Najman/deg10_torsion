from pathlib import Path
import concurrent.futures,subprocess,json
p=Path(__file__).resolve().parent
scripts=sorted((p/'computations').glob('unit_*.py'))+sorted((p/'computations').glob('function_*.m'))
print('SECOND_BATCH_START tasks',len(scripts),'workers',8,flush=True)
def run(s):
 r=subprocess.run(['/usr/bin/python3',str(p/'run_job.py'),str(s),'--seconds','180','--engine','sage' if s.suffix=='.py' else 'magma'],capture_output=True,text=True)
 print(s.stem,'supervisor_exit',r.returncode,flush=True)
 return {'script':str(s.relative_to(p)),'supervisor_exit':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool: rr=list(pool.map(run,scripts))
(p/'jobs/second_batch_result.json').write_text(json.dumps(rr,indent=2)+'\n')
print('SECOND_BATCH_COMPLETE',flush=True)
