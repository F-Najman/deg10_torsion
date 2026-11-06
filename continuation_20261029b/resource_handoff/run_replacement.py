from pathlib import Path
import sys,subprocess,time,json,os,hashlib
script=Path(sys.argv[1]);out=Path(sys.argv[2]);status=Path(sys.argv[3]);t=time.time()
with out.open('x') as fp:
 p=subprocess.Popen([sys.executable,str(script)],cwd=script.parents[2],stdout=fp,stderr=subprocess.STDOUT)
 rec={'purpose':'Continue existing degree-18 computations at user-authorized lower concurrency','state':'RUNNING','pid':p.pid,'owner_uid':os.getuid(),'start_time_unix':t,'script':str(script),'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()};status.write_text(json.dumps(rec,indent=2)+'\n')
 code=p.wait();rec.update(state='EXITED',actual_exit_code=code,elapsed_seconds=time.time()-t);status.write_text(json.dumps(rec,indent=2)+'\n')
