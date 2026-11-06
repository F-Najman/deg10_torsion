from pathlib import Path
import subprocess,time,json,datetime,sys
H=Path(__file__).resolve().parents[1]
ROOT=H.parents[1]
report=Path(__file__).resolve().parent/'admission.json'
steps=[('late240_isogeny',420),('late240_cyclic',900)]
rows=[]
for stem,cap in steps:
 t=time.monotonic()
 while True:
  lines=subprocess.check_output(['ps','-eo','stat=,comm='],text=True).splitlines()
  active=sum(len(z.split())==2 and z.split()[1]=='magma.exe' and not z.split()[0].startswith(('T','Z')) for z in lines)
  if active<min(30,(__import__('os').cpu_count() or 2)-2):break
  record={'pending':stem,'host_active_magma':active,'wait_seconds':round(time.monotonic()-t,1),'steps':rows,'state':'WAITING_FOR_CAPACITY'}
  report.write_text(json.dumps(record,indent=2)+'\n')
  if time.monotonic()-t>600:
   record['state']='EXITED_WITHOUT_CAPACITY';report.write_text(json.dumps(record,indent=2)+'\n');print(record,flush=True);sys.exit(2)
  time.sleep(10)
 cmd=[sys.executable,str(H/'run_job.py'),str(H/'computations'/(stem+'.m')),'--seconds',str(cap)]
 r=subprocess.run(cmd,cwd=ROOT)
 status=json.loads((H/'jobs'/(stem+'.status.json')).read_text())
 rows.append({'step':stem,'host_active_magma_before_start':active,'wrapper_exit':r.returncode,'job_status':status})
 report.write_text(json.dumps({'state':'ACTIVE','steps':rows},indent=2)+'\n')
 if r.returncode:sys.exit(r.returncode)
 if stem=='late240_cyclic':
  txt=(H/'jobs'/(stem+'.log')).read_text()
  if 'CYCLIC40_DONE' not in txt:raise RuntimeError('No completed certificate marker')
report.write_text(json.dumps({'state':'EXITED','steps':rows},indent=2)+'\n')
