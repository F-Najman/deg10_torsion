from pathlib import Path
import json,os,signal,time
H=Path(__file__).resolve().parents[1];C=H/'continuation_20261031a/C44';final=json.loads((C/'independent_final_verification.json').read_text());assert final['status']=='PASS'and final['unique_genuine_geometric_kernels']==3916
sp=H/'jobs/d31a_C44_independent_isogenies.status.json';st=json.loads(sp.read_text());report={'replacement':'Six complete disjoint independently checked genuine-isogeny shards','replacement_summary':'independent_final_verification.json','old_status_path':str(sp.relative_to(H))}
if st['state']=='EXITED':report.update(action='already exited; no signal sent',actual_exit_code=st['actual_exit_code'])
else:
 assert st['purpose']=='d31a_C44_independent_isogenies'and st['owner_uid']==os.getuid();pid=st['pid'];fd=os.pidfd_open(pid)
 try:
  pp=Path('/proc')/str(pid);assert pp.stat().st_uid==st['owner_uid'];cmd=[v.decode()for v in(pp/'cmdline').read_bytes().split(b'\0')if v];assert cmd==st['command']
  stat=(pp/'stat').read_text();v=stat[stat.rfind(')')+2:].split();ticks=int(v[19]);boot=int(next(s.split()[1]for s in Path('/proc/stat').read_text().splitlines()if s.startswith('btime ')));started=boot+ticks/os.sysconf('SC_CLK_TCK');assert abs(started-st['start_time_unix'])<5
  children=(pp/'task'/str(pid)/'children').read_text().strip();assert not children
  assert os.getpgid(pid)==pid;report.update(verified_pid=pid,verified_owner_uid=st['owner_uid'],verified_start_time_unix=started,verified_start_ticks=ticks,verified_command=cmd,verified_no_child_workers=True)
  signal.pidfd_send_signal(fd,signal.SIGTERM);report['action']='SIGTERM via verified pidfd after complete replacement verification'
 finally:os.close(fd)
 for i in range(100):
  now=json.loads(sp.read_text())
  if now['state']=='EXITED':report['actual_exit_code']=now['actual_exit_code'];break
  time.sleep(.1)
 else:raise AssertionError('bounded runner has not yet recorded retired worker exit')
(C/'serial_backup_retirement.json').write_text(json.dumps(report,indent=2)+'\n');print('D31A_C44_SERIAL_BACKUP_RETIRED',report,flush=True)
