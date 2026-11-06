from pathlib import Path
import subprocess,concurrent.futures
H=Path(__file__).resolve().parents[1]
def run(p):
 r=subprocess.run(['/usr/bin/python3',str(H/'run_job.py'),str(p),'--seconds','120','--engine','sage'],capture_output=True,text=True)
 print(p.stem,'exit',r.returncode,flush=True)
ss=sorted((H/'computations').glob('sporadic_check_unit_*.py'))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run,ss))
