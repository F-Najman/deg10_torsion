from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import subprocess,json
H=Path(__file__).resolve().parents[1];jobs=json.loads((H/'continuation_20261029/shard_jobs.json').read_text())
def run(t):
 n,sec=t;r=subprocess.run(['python3',str(H/'run_job.py'),str(H/'computations'/(n+'.m')),'--seconds',str(sec)],capture_output=True,text=True)
 p=H/'jobs'/(n+'.status.json');d=json.loads(p.read_text()) if p.exists() else {'runner_error':r.stderr}
 out={k:d[k] for k in ['purpose','state','actual_exit_code','elapsed_seconds','termination_reason'] if k in d};print(json.dumps(out),flush=True);return out
# Current host count: 9 Magma workers, including 4 ours. 18 additional one-thread workers leave >=5/32 CPUs free.
with ThreadPoolExecutor(max_workers=18) as pool:rows=[f.result() for f in as_completed([pool.submit(run,t) for t in jobs])]
(H/'continuation_20261029/shard_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('ALL_SHARDS_EXITED',flush=True)
