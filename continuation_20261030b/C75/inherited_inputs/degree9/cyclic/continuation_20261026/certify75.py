from pathlib import Path
import hashlib,json,re,subprocess,sys
H=Path(__file__).resolve().parent;B=H.parents[1];m=json.loads((H/'packet75.json').read_text());assert m['degree']==9 and m['complete']
for e in m['files']:assert hashlib.sha256((B/e['path']).read_bytes()).hexdigest()==e['sha256'],e['path']
for tag in ['N75_pure9_kernel','independent_counts75']:
 s=(H/(tag+'.log')).read_text();t=json.loads((H/(tag+'_status.json')).read_text());assert t['native_exit_code']==0 and t['log_sha256']==hashlib.sha256(s.encode()).hexdigest()
 assert not re.search('Runtime error|User error|Assertion failed|syntax error',s,re.I)
s=(H/'N75_pure9_kernel.log').read_text();assert re.search(r'NONIC75_PURE9_RR_COMPLETE genus73=true image71=true width75=true rank=9 covered10=true time=[0-9.]+\s*$',s)
assert 'N75_GON5_CANONICAL_GENUS73_MODEL_VERIFIED'in s and 'N75_D2_AUTOMORPHISM_VERIFIED'in s
assert 'closed_place_counts=[ 0, 0, 0, 0, 0, 0, 20, 0, 0 ]'in(H/'independent_counts75.log').read_text()
for script,args,marker in [('N75_replay_degree9.py',['--readonly'],'NONIC75_RETAINED_FI9_REPLAY'),('N75_cusp_complements_D9.py',[],'NONIC75_CUSP_COMPLEMENT_DEGREE2_PASS'),('N75_replay_good7_readonly.py',[],'C75_GOOD7_INDEPENDENT_REPLAY_PASS')]:
 r=subprocess.run([sys.executable,str(H/script)]+args,capture_output=True,text=True,check=True);assert marker in r.stdout;print(marker)
v=json.loads((H/'N75_degree9_FI_replay.json').read_text());assert v['tested_divisors']==347778 and v['exceptions_count']==0 and v['missing_pure9_width75_rational']==10 and v['all_admissible_degree9_divisors']==347788
print('NONIC75_READONLY_CERTIFICATE_PASS genus145 quotient73 good7_complete all347788_admissible_FI9_rank9 complete=true')
