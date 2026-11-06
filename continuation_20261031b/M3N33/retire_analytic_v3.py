from pathlib import Path
import json,os,signal,time
W=Path('work/deg10').resolve();checks=[]
assert 'TypeError: object of type' in(W/'jobs/d31b_m333_independent_analytic_v3_s3.log').read_text()
for j in range(4):
 stem=f'd31b_m333_independent_analytic_v3_s{j}';st=W/'jobs'/(stem+'.status.json');s=json.loads(st.read_text());p=s['pid'];r={'status_before':s,'reason':'All identical v3 scripts contain the demonstrated GF(2^10) polynomial len() TypeError, so none can finish; rational-cusp partial checks remain diagnostic'}
 assert '2*d-len((th**i*om**j).polynomial())'in(W/s['script']).read_text()
 if s['state']=='RUNNING':
  q=Path(f'/proc/{p}');assert q.stat().st_uid==s['owner_uid']==os.getuid();cmd=q.joinpath('cmdline').read_bytes().split(b'\0');assert str((W/s['script']).resolve()).encode()in cmd;exe=q.joinpath('exe').resolve();assert str(exe).startswith('/opt/sage/') and 'python'in exe.name
  raw=q.joinpath('stat').read_text().rsplit(') ',1)[1].split();bt=int(next(x.split()[1]for x in Path('/proc/stat').read_text().splitlines()if x.startswith('btime ')));start=bt+int(raw[19])/os.sysconf('SC_CLK_TCK');assert abs(start-s['start_time_unix'])<3 and os.getpgid(p)==p;r['verified_identity']={'pid':p,'uid':s['owner_uid'],'start_time_unix':start,'command':[a.decode()for a in cmd if a],'exe':str(exe)};os.killpg(p,signal.SIGTERM)
  for _ in range(100):
   time.sleep(.1);now=json.loads(st.read_text())
   if now['state']=='EXITED':break
  assert now['state']=='EXITED';r['status_after']=now;r['actual_exit_code']=now['actual_exit_code']
 else:r['already_exited']=True;r['actual_exit_code']=s['actual_exit_code']
 checks.append(r)
(W/'continuation_20261031b/M3N33/analytic_v3_retirement.json').write_text(json.dumps(checks,indent=2)+'\n');print('M333_ANALYTIC_V3_BUGGY_JOBS_RETIRED',[r['actual_exit_code']for r in checks])
