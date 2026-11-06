from pathlib import Path
import subprocess,concurrent.futures,json
H=Path(__file__).resolve().parents[1];root=H.parents[1]
scripts=sorted((H/'computations').glob('cont_source612_recert_v2_shard*.m'));assert len(scripts)==8

def one(p):
 r=subprocess.run(['/usr/bin/python3',str(H/'run_job.py'),str(p),'--seconds','300'],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 s=json.loads((H/'jobs'/(p.stem+'.status.json')).read_text());print(json.dumps({k:s.get(k) for k in ['purpose','state','actual_exit_code','elapsed_seconds']}),flush=True);return s
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(one,scripts))
(H/'continuation_20261029/612_source_recert_v2_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('DECIC612_SOURCE_RECERT_V2_FINISHED',len(rows),flush=True)
