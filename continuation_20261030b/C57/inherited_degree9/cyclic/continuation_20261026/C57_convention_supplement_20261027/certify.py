from pathlib import Path
import hashlib,json,subprocess,sys,os,time
S=Path(__file__).resolve().parent;D=S.parents[2];P=S.parent/'C57_source85';m=json.loads((S/'manifest.json').read_text())
assert m['degree']==9 and m['complete_exclusion']and m['excluded_groups']==[[1,57]]and m['necessary_point_coefficients']==[1,4,-5]and m['also_tested_coefficients']==[5,-4,-1]
for r in m['files']:
 p=D/r['path'];assert p.is_file()and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
assert hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest()=='5e604ee91fd1cc8633e1b7c5c8a5fb682e5bc21dde3b34bd0607e70ee2f8f6d4'
a=json.loads((P/'rr24_v3_completion.json').read_text());assert a['complete']and a['actual_exit_code']==0 and'NONIC57_SOURCE85_RR24_ALL_EXCLUSIONS_COMPLETE'in(P/'rr24_v3.log').read_text()
a=json.loads((P/'ring_RR_sage_certificate.json').read_text());assert a['complete']and[r['positive_degree']for r in a['nonprincipal_divisors']]==[42,42,88]and all(r['final_criterion_rank']==252 for r in a['nonprincipal_divisors'])
assert'enumerate([1,4,-5])'in(P/'assemble_RR_divisors_sage.py').read_text()
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args=[sys.executable,str(P/'certify.py')]+(['--deep']if'--deep'in sys.argv else[]);r=subprocess.run(args,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);assert r.returncode==0 and'C57_D9_COMPLETE_READ_ONLY_PACKET_PASS'in r.stdout,r.stdout;print(r.stdout,end='',flush=True)
if'--deep'in sys.argv:
 for script,marker in[('assemble_RR_divisors_sage.py','C57_NECESSARY_RR_DIVISORS_ASSEMBLED'),('replay_ring_RR_sage.py','C57_CANONICAL_RING_AND_ALL_RR_OBSTRUCTIONS_INDEPENDENT_SAGE_PASS')]:
  t=time.monotonic();r=subprocess.run([sys.executable,str(P/script),'--readonly'],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);assert r.returncode==0 and marker in r.stdout,(script,r.stdout);print('POINT_NECESSITY_READ_ONLY_REPLAY_PASS',script,round(time.monotonic()-t,2),flush=True)
print('C57_CONVENTION_SUPPLEMENT_COMPLETE_PASS both_arrangements_checked=true point_coefficients=1,4,-5 files='+str(len(m['files'])),flush=True)
