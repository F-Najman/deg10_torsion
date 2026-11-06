#!/usr/bin/env python3
import concurrent.futures,json,os,pathlib,subprocess,threading,time
BASE=pathlib.Path.cwd();P=BASE/'work/deg10';DEST=P/'C236_resolution_20261030/enumeration';lock=threading.Lock();results=[]
tags=['ownsmall','D8','D9']+[f'D10s{i}' for i in range(1,8)]
def one(tag):
    end=time.monotonic()+90
    while True:
        with lock:
            ps=subprocess.check_output(['ps','-eo','stat=,comm='],text=True).splitlines()
            active=sum(len(s.split())==2 and s.split()[1]=='magma.exe' and not s.split()[0].startswith(('T','Z')) for s in ps)
            if active<min(30,(os.cpu_count() or 2)-2):
                stem=f'new236_filter18_native_replay_{tag}_v2'
                proc=subprocess.Popen(['python3',str(P/'run_job.py'),str(P/'computations'/f'{stem}.m'),'--seconds','240'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                time.sleep(.1)
                break
        if time.monotonic()>end:raise RuntimeError('No safe core admission')
        time.sleep(1)
    out,_=proc.communicate();rec=json.loads((P/'jobs'/f'{stem}.status.json').read_text());rec['host_magma_before_admission']=active
    assert rec['actual_exit_code']==0 and 'C236_IDEAL_REPLAY_COMPLETE' in (P/'jobs'/f'{stem}.log').read_text(),out
    with lock:
        results.append(rec);(DEST/'native_replay_completion.json').write_text(json.dumps(results,indent=2)+'\n')
    print('C236_REPLAY_EXIT',tag,rec['elapsed_seconds'],flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(one,tags))
# The independently launched benchmark is the remaining eleventh shard.
p=P/'jobs/new236_filter18_native_replay_D10s0_v2.status.json'
rec=json.loads(p.read_text());assert rec['state']=='EXITED' and rec['actual_exit_code']==0
assert 'C236_IDEAL_REPLAY_COMPLETE' in (P/'jobs/new236_filter18_native_replay_D10s0_v2.log').read_text()
results.append(rec);(DEST/'native_replay_completion.json').write_text(json.dumps(results,indent=2)+'\n')
print('C236_ALL_NATIVE_REPLAYS_COMPLETE',len(results),flush=True)
