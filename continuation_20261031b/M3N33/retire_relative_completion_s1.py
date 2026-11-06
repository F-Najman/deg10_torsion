from pathlib import Path
import json,os,signal,time,hashlib
W=Path('work/deg10').resolve();stem='d31b_m333_cusp_jets_s1';status=W/'jobs'/(stem+'.status.json');s=json.loads(status.read_text());p=s['pid'];report={'purpose':'Abandon stalled relative-field Completion; verified base-field series plus explicit tame Kummer Hensel lifting replace the algorithm without changing the boundary criterion','status_before':s}
if s['state']=='RUNNING':
 q=Path(f'/proc/{p}');assert q.stat().st_uid==s['owner_uid']==os.getuid();cmd=q.joinpath('cmdline').read_bytes().split(b'\0');assert str((W/s['script']).resolve()).encode() in cmd;exe=q.joinpath('exe').resolve();assert str(exe).startswith('/opt/magma-') and exe.name in ['magma.exe','magma.avx2.exe']
 raw=q.joinpath('stat').read_text().rsplit(') ',1)[1].split();bt=int(next(x.split()[1] for x in Path('/proc/stat').read_text().splitlines() if x.startswith('btime ')));start=bt+int(raw[19])/os.sysconf('SC_CLK_TCK');assert abs(start-s['start_time_unix'])<3 and os.getpgid(p)==p
 replacement=(W/'jobs/d31b_m333_base_completion_probe.log').read_text();assert 'D31B_M333_BASE_COMPLETION_PROBE_COMPLETE' in replacement
 report['verified_identity']={'pid':p,'uid':s['owner_uid'],'start_time_unix':start,'command':[a.decode() for a in cmd if a],'exe':str(exe)};os.killpg(p,signal.SIGTERM)
 for _ in range(80):
  time.sleep(.1);now=json.loads(status.read_text())
  if now['state']=='EXITED':break
 assert now['state']=='EXITED';report['actual_exit_code']=now['actual_exit_code'];report['status_after']=now
else:report['already_exited']=True;report['actual_exit_code']=s['actual_exit_code']
Path('work/deg10/continuation_20261031b/M3N33/obsolete_relative_completion_s1_retirement.json').write_text(json.dumps(report,indent=2)+'\n');print('M333_OBSOLETE_FACTOR_RETIRED',report['actual_exit_code'])
