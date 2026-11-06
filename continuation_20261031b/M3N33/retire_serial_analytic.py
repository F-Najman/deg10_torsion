from pathlib import Path
import json,os,signal,time
W=Path('work/deg10').resolve();stem='d31b_m333_independent_analytic_v2';st=W/'jobs'/(stem+'.status.json');s=json.loads(st.read_text());p=s['pid'];report={'purpose':'Retire serial independent series fallback after disjoint replacement shards verified base cusps1 and4, including the original low-precision failure','status_before':s}
assert 'M333_INDEPENDENT_ANALYTIC_BASE 1 'in(W/'jobs/d31b_m333_independent_analytic_v3_s0.log').read_text();assert 'M333_INDEPENDENT_ANALYTIC_BASE 4 'in(W/'jobs/d31b_m333_independent_analytic_v3_s3.log').read_text()
if s['state']=='RUNNING':
 q=Path(f'/proc/{p}');assert q.stat().st_uid==s['owner_uid']==os.getuid();cmd=q.joinpath('cmdline').read_bytes().split(b'\0');assert str((W/s['script']).resolve()).encode()in cmd;exe=q.joinpath('exe').resolve();assert str(exe).startswith('/opt/sage/') and 'python' in exe.name
 raw=q.joinpath('stat').read_text().rsplit(') ',1)[1].split();bt=int(next(x.split()[1]for x in Path('/proc/stat').read_text().splitlines()if x.startswith('btime ')));start=bt+int(raw[19])/os.sysconf('SC_CLK_TCK');assert abs(start-s['start_time_unix'])<3 and os.getpgid(p)==p;report['verified_identity']={'pid':p,'uid':s['owner_uid'],'start_time_unix':start,'command':[a.decode()for a in cmd if a],'exe':str(exe)};os.killpg(p,signal.SIGTERM)
 for _ in range(100):
  time.sleep(.1);now=json.loads(st.read_text())
  if now['state']=='EXITED':break
 assert now['state']=='EXITED';report['status_after']=now;report['actual_exit_code']=now['actual_exit_code']
else:report['already_exited']=True;report['actual_exit_code']=s['actual_exit_code']
(W/'continuation_20261031b/M3N33/serial_analytic_retirement.json').write_text(json.dumps(report,indent=2)+'\n');print('M333_SERIAL_ANALYTIC_RETIRED',report['actual_exit_code'])
