from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import subprocess,json,time
H=Path(__file__).resolve().parents[1];root=H.parents[1]
local=[f'cont_local_{m}_{n}_q{q}' for m,n,q in [(3,18,7),(4,16,5),(6,12,7),(3,21,4),(3,24,7),(4,24,5),(3,27,4),(4,20,9),(4,28,9)]]
wind=[f'cont_winding{n}' for n in [66,68,69,70,72,74,76,78,80,81,84,87,90,93,99,100,105,135]]
# 8 one-thread workers, well below the checked capacity (32 host CPUs, four existing busy workers).
def run(stem):
 sec=600 if stem in local[:3] else 240
 r=subprocess.run(['python3',str(H/'run_job.py'),str(H/'computations'/(stem+'.m')),'--seconds',str(sec)],capture_output=True,text=True)
 status=H/'jobs'/(stem+'.status.json')
 if status.exists():rec=json.loads(status.read_text());rec={k:rec[k] for k in ['purpose','state','actual_exit_code','elapsed_seconds','termination_reason']}
 else:rec={'purpose':stem,'runner_error':r.stderr[-1500:]}
 print(json.dumps(rec),flush=True);return rec
with ThreadPoolExecutor(max_workers=8) as pool:
 rows=[f.result() for f in as_completed([pool.submit(run,s) for s in local+wind])]
(H/'continuation_20261029/initial_batch_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('INITIAL_BATCH_EXITED',flush=True)
