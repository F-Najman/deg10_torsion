from pathlib import Path
import concurrent.futures,subprocess,json,time
H=Path(__file__).resolve().parents[1];root=H.parents[1]
scripts=sorted((H/'computations').glob('cont_sieve_4_20_*_shard*.m'))
assert len(scripts)==20

def one(p):
 start=time.time();r=subprocess.run(['/usr/bin/python3',str(H/'run_job.py'),str(p),'--seconds','900'],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 s=json.loads((H/'jobs'/(p.stem+'.status.json')).read_text());print(json.dumps({k:s.get(k) for k in ['purpose','state','actual_exit_code','elapsed_seconds']}),flush=True);return s
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool:rows=list(pool.map(one,scripts))
(H/'continuation_20261029/420_shard_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('DECIC420_ALL_SHARDS_FINISHED',len(rows),flush=True)
