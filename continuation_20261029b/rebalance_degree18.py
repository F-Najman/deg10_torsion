from pathlib import Path
import os,signal,subprocess,time,json,hashlib,sys
P=Path(__file__).resolve().parent;A=P/'resource_handoff';A.mkdir(exist_ok=True)
W=Path('/home/fnajman/Documents/Codex/2026-10-28/on-this-computer-there-are-a/work/degree18')
targets=[(840893,'large_order_batch.py',16,6),(850950,'run_new_ordinary.py',12,4),(868598,'run_ordinary131.py',6,2)]
def stat(pid):
 p=Path('/proc')/str(pid);s=(p/'stat').read_text().rsplit(')',1)[1].split()
 return {'pid':pid,'uid':p.stat().st_uid,'state':s[0],'ppid':int(s[1]),'start_ticks':int(s[19]),'wait_status':int(s[49]),'command':(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()}
records=[]
for pid,name,old,new in targets:
 s=stat(pid);assert s['uid']==os.getuid()==1001 and s['command'].strip()=='python work/degree18/'+name
 assert s['state'] not in ['T','Z'];data=(W/name).read_bytes();assert f'max_workers={old}'.encode() in data
 snap=A/(name+'.before');assert not snap.exists();snap.write_bytes(data)
 records.append({'launcher':s,'script':str(W/name),'before_sha256':hashlib.sha256(data).hexdigest(),'old_concurrency':old,'new_concurrency':new,'children':[]})
for r in records:
 s=stat(r['launcher']['pid']);assert s['start_ticks']==r['launcher']['start_ticks'] and s['command']==r['launcher']['command']
 os.kill(s['pid'],signal.SIGSTOP)
for r in records:
 pid=r['launcher']['pid'];assert stat(pid)['state']=='T'
 kids=subprocess.check_output(['ps','--ppid',str(pid),'-o','pid='],text=True).split()
 for kid in kids:
  s=stat(int(kid));assert s['uid']==1001 and s['ppid']==pid
  assert 'magma.exe' in s['command'] or s['state']=='Z'
  r['children'].append(s)
 path=Path(r['script']);data=path.read_text();data=data.replace(f'max_workers={r["old_concurrency"]}',f'max_workers={r["new_concurrency"]}')
 compile(data,str(path),'exec');path.write_text(data);r['after_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
report={'authorization':'User explicitly prioritized Phi(10) and authorized reducing degree-18 launcher concurrency while preserving useful computations.','state':'DRAINING_EXISTING_WORKERS','records':records}
def save(): (A/'handoff.json').write_text(json.dumps(report,indent=2)+'\n')
save();print('Paused only three Python launchers; all their Magma workers continue to completion.',flush=True)
t0=time.monotonic()
while True:
 active=0
 for r in records:
  for c in r['children']:
   if 'actual_exit_code' in c:continue
   try:s=stat(c['pid'])
   except FileNotFoundError: c['state']='GONE';c['actual_exit_code']=None;continue
   assert s['start_ticks']==c['start_ticks'] and s['uid']==c['uid']
   if s['state']=='Z':c['state']='EXITED';c['actual_exit_code']=os.waitstatus_to_exitcode(s['wait_status'])
   else:active+=1
 save()
 if not active:break
 if time.monotonic()-t0>900:raise RuntimeError('Useful old workers still running; preserve them and inspect handoff before any further signals.')
 time.sleep(2)
for r in records:
 pid=r['launcher']['pid'];s=stat(pid);assert s['start_ticks']==r['launcher']['start_ticks'] and s['state']=='T'
 os.kill(pid,signal.SIGTERM);os.kill(pid,signal.SIGCONT);r['retired_with']='SIGTERM after every captured child completed; SIGCONT delivered pending termination'
report['state']='OLD_WORKERS_COMPLETED_LAUNCHERS_RETIRED';save()
wrapper=A/'run_replacement.py'
wrapper.write_text('''from pathlib import Path
import sys,subprocess,time,json,os,hashlib
script=Path(sys.argv[1]);out=Path(sys.argv[2]);status=Path(sys.argv[3]);t=time.time()
with out.open('x') as fp:
 p=subprocess.Popen([sys.executable,str(script)],cwd=script.parents[2],stdout=fp,stderr=subprocess.STDOUT)
 rec={'purpose':'Continue existing degree-18 computations at user-authorized lower concurrency','state':'RUNNING','pid':p.pid,'owner_uid':os.getuid(),'start_time_unix':t,'script':str(script),'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()};status.write_text(json.dumps(rec,indent=2)+'\\n')
 code=p.wait();rec.update(state='EXITED',actual_exit_code=code,elapsed_seconds=time.time()-t);status.write_text(json.dumps(rec,indent=2)+'\\n')
''')
compile(wrapper.read_text(),str(wrapper),'exec')
for r in records:
 name=Path(r['script']).stem
 log=A/(name+'.replacement.log');status=A/(name+'.replacement.status.json')
 proc=subprocess.Popen([sys.executable,str(wrapper),r['script'],str(log),str(status)],stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
 r['replacement_supervisor_pid']=proc.pid;r['replacement_status']=str(status)
report['state']='REBALANCED';report['degree18_concurrency_after']=14;report['unchanged_special_launcher_concurrency']=2;save()
(W/'PHI10_RESOURCE_PRIORITY_NOTICE.json').write_text(json.dumps({'request':'User prioritized Phi(10); keep two host cores idle.','degree18_limits':{Path(r['script']).name:r['new_concurrency'] for r in records},'special_limit':2,'audit_report':str(A/'handoff.json'),'note':'Old workers completed; completed logs are reused. Do not restore previous oversubscribed concurrency while Phi(10) work is active.'},indent=2)+'\n')
print('DEGREE18_REBALANCED_WITH_ALL_CAPTURED_WORKERS_FINISHED',[(Path(r['script']).name,len(r['children']),[c['actual_exit_code'] for c in r['children']]) for r in records],flush=True)
