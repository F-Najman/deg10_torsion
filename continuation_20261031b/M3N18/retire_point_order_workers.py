from pathlib import Path
import json,os,signal,time,shutil,hashlib
W=Path('work/deg10').resolve();assert 'D31B_M318_DIRECT_TORSION_ASSERTIONS_COMPLETE' in (W/'jobs/d31b_m318_torsion_assertions.log').read_text();H=W/'continuation_20261031b/M3N18';out=[]
for h in range(4):
 p=H/f'paired_fast_s{h}_partial.m'
 if p.exists():shutil.copy2(p,H/f'preserved_fast_s{h}_partial.m')

for stem in [f'd31b_m318_p13_fast_s{i}' for i in range(4)]:
 status=W/'jobs'/(stem+'.status.json');s=json.loads(status.read_text());rec={'stem':stem,'before':s,'reason':'Retire repeated SEA/group-order computations; equivalent direct torsion and Weil-pairing assertions passed on every cached source point'}
 if s['state']=='RUNNING':
  pid=s['pid'];q=Path(f'/proc/{pid}');assert q.stat().st_uid==s['owner_uid']==os.getuid();cmd=q.joinpath('cmdline').read_bytes().split(b'\0');assert str((W/s['script']).resolve()).encode() in cmd;exe=q.joinpath('exe').resolve();assert str(exe).startswith('/opt/magma-') and exe.name in ['magma.exe','magma.avx2.exe'];raw=q.joinpath('stat').read_text().rsplit(') ',1)[1].split();bt=int(next(x.split()[1] for x in Path('/proc/stat').read_text().splitlines() if x.startswith('btime ')));start=bt+int(raw[19])/os.sysconf('SC_CLK_TCK');assert abs(start-s['start_time_unix'])<3 and os.getpgid(pid)==pid
  rec['verified_identity']={'pid':pid,'owner_uid':s['owner_uid'],'start_time_unix':start,'command':[a.decode() for a in cmd if a],'executable':str(exe)};os.killpg(pid,signal.SIGTERM)
  for _ in range(100):
   time.sleep(.1);now=json.loads(status.read_text())
   if now['state']=='EXITED':break
  assert now['state']=='EXITED';rec['after']=now
 else:rec['after']=s
 out.append(rec)
(H/'point_order_workers_retirement.json').write_text(json.dumps({'records':out,'verified_replacement_marker':'D31B_M318_FAST_MAP_BENCHMARK_COMPLETE 30'},indent=2)+'\n');print('M318_UNBALANCED_WORKERS_RETIRED',[(a['stem'],a['after']['actual_exit_code']) for a in out])
