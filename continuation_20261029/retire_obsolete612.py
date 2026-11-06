from pathlib import Path
import os,json,time,signal
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
probe=json.loads((H/'jobs/cont_RR612_norm_probe.status.json').read_text());assert probe['state']=='EXITED' and probe['actual_exit_code']==0
uid=os.getuid();boot=int(next(s.split()[1] for s in Path('/proc/stat').read_text().splitlines() if s.startswith('btime ')));ticks=os.sysconf('SC_CLK_TCK');records=[]
def inspect(pid):
 p=Path('/proc')/str(pid);st=p.stat();cmd=p.joinpath('cmdline').read_bytes().split(b'\0');cmd=[s.decode() for s in cmd if s];stat=p.joinpath('stat').read_text().rsplit(')',1)[1].split();start=boot+int(stat[19])/ticks
 return st.st_uid,cmd,start
# Retire the finite coordinator first, so queued obsolete workers cannot start.
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:u,cmd,start=inspect(int(p.name))
 except (FileNotFoundError,ProcessLookupError,PermissionError):continue
 if u==uid and any(s in [str(D/'run612_RR13_shards.py'),'work/deg10/continuation_20261029/run612_RR13_shards.py'] for s in cmd):
  os.kill(int(p.name),signal.SIGTERM);records.append({'kind':'retired coordinator','pid':int(p.name),'owner_uid':u,'command':cmd,'start_time_unix':start,'signal':'SIGTERM'})
for path in sorted((H/'jobs').glob('*.status.json')):
 s=json.loads(path.read_text());name=s['purpose']
 if s['state']!='RUNNING' or not (name.startswith('cont_RR612_source_p13_shard') or name=='cont_diagonal_6_12_D5_v5'):continue
 pid=s['pid']
 try:u,cmd,start=inspect(pid)
 except (FileNotFoundError,ProcessLookupError):continue
 assert u==uid==s['owner_uid'] and abs(start-s['start_time_unix'])<3 and cmd[-1]==str((H/s['script']).resolve()) and os.getpgid(pid)==pid,(s,cmd,start)
 os.killpg(pid,signal.SIGTERM);records.append({'kind':'retired obsolete worker','purpose':name,'pid':pid,'owner_uid':u,'command':cmd,'start_time_unix':start,'signal':'SIGTERM'})
time.sleep(2)
for r in records:
 if 'purpose' in r:
  p=H/'jobs'/(r['purpose']+'.status.json');s=json.loads(p.read_text());r['recorded_state']=s['state'];r['actual_exit_code']=s.get('actual_exit_code')
  if s['state']=='EXITED' and s.get('actual_exit_code')==-15:s['termination_reason']='obsolete method retired after verified norm-polynomial parity replacement; not a completed certificate';p.write_text(json.dumps(s,indent=2)+'\n')
(D/'retired_obsolete612.json').write_text(json.dumps({'reason':'Verified norm-polynomial parity test replaces full divisor factorization and target orbit enumeration. Preserve all partial logs; do not promote terminated runs as complete proofs.','actions':records},indent=2)+'\n')
print('Retired obsolete coordinator/workers:',[(r.get('purpose',r['kind']),r['pid'],r.get('actual_exit_code')) for r in records])
