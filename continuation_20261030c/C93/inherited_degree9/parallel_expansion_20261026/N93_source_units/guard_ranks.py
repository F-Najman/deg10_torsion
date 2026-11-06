from pathlib import Path
import subprocess,resource,os,time,json,hashlib
p=Path(__file__).resolve().parent
def limits():
 os.nice(10);resource.setrlimit(resource.RLIMIT_CPU,(300,310));resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
t=time.monotonic()
with (p/'ranks_m4ri.log').open('w')as f:
 child=subprocess.Popen(['python3',str(p/'verify_R1_ranks_m4ri.py')],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits)
 try:code=child.wait(timeout=1200)
 except subprocess.TimeoutExpired:child.kill();code=child.wait()
state={'actual_exit_code':code,'CPU_limit':300,'wall_limit':1200,'AS_bytes':8*1024**3,'nice':10,'elapsed_seconds':time.monotonic()-t,'log_sha256':hashlib.sha256((p/'ranks_m4ri.log').read_bytes()).hexdigest()};(p/'ranks_m4ri.exit.json').write_text(json.dumps(state,indent=2)+'\n');print(state)
