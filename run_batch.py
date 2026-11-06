#!/usr/bin/env python3
from pathlib import Path
import concurrent.futures,subprocess,json,time,os
p=Path(__file__).resolve().parent
scripts=sorted((p/'computations').glob('witness_*.m'))+sorted((p/'computations').glob('family_*_v2.m'))
scripts=[s for s in scripts if not(p/'jobs'/(s.stem+'.status.json')).exists()]
# At launch the measured host use was 11 cores on 32 cores. Twelve single-thread
# workers leave substantial headroom beyond the user's two-core reserve.
print('DECIC_BATCH_START tasks='+str(len(scripts))+' max_workers=12',flush=True)
def run(s):
 r=subprocess.run(['/usr/bin/python3',str(p/'run_job.py'),str(s),'--seconds','120'],capture_output=True,text=True)
 print(s.stem+' supervisor_exit='+str(r.returncode)+' '+r.stdout.strip(),flush=True)
 return {'script':str(s.relative_to(p)),'supervisor_exit':r.returncode}
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:results=list(pool.map(run,scripts))
(p/'jobs/initial_batch_result.json').write_text(json.dumps({'purpose':'Exact decic occurrences and family probes','workers':12,'results':results},indent=2)+'\n')
print('DECIC_BATCH_COMPLETE',flush=True)
