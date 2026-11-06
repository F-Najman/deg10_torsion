from pathlib import Path
import subprocess,concurrent.futures,threading,time,json,os
H=Path(__file__).resolve().parents[1];ROOT=H.parents[1];lock=threading.Lock()
def run(p):
 stem=f'resolve226_finite_p{p}'
 with lock:
  lines=subprocess.check_output(['ps','-eo','stat=,comm='],text=True).splitlines();n=sum(len(x.split())==2 and x.split()[1]=='magma.exe' and not x.split()[0].startswith(('T','Z')) for x in lines);assert n<min(28,os.cpu_count()-2)
  proc=subprocess.Popen([os.sys.executable,str(H/'run_job.py'),str(H/'computations'/f'{stem}.m'),'--seconds','180'],cwd=ROOT,stdout=subprocess.PIPE,text=True);time.sleep(.25)
 proc.communicate();s=json.loads((H/'jobs'/f'{stem}.status.json').read_text());print(p,s['actual_exit_code'],s['elapsed_seconds'],flush=True);return s
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(run,[17,23,29,31]))
Path(__file__).with_name('medium_prime_completion.json').write_text(json.dumps(rows,indent=2)+'\n')
