"""Freeze complete C55 evidence after the independent native replay passes."""
from pathlib import Path
import ast,datetime,hashlib,json,os
H=Path(__file__).resolve().parent;ROOT=H.parents[3];BASE=ROOT/'work/deg10'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,data):(H/name).write_text(json.dumps(data,indent=2)+'\n')
def read(p):return ast.literal_eval(Path(p).read_text().replace('[*','[').replace('*]',']'))
native=BASE/'jobs/contb55_parent_native.status.json';d=json.loads(native.read_text());assert d['state']=='EXITED'and d['actual_exit_code']==0
assert 'C55_PARENT_NATIVE_COMPLETE 5 'in(BASE/d['log']).read_text()
a=read(H/'parent_native_replay.m');assert len(a)==5 and all(r[1:]==[10,10,40,0]for r in a)
for i in range(1,6):
 r=read(H/f'source_case{i}_result.m');assert r[0]==i and r[1]==0
for name in ['cusp_audit.json','source_torsion_injectivity.json','independent_finite_geometry.json']:
 assert json.loads((H/name).read_text())['status'].startswith('PASS')
reuse=json.loads((H/'reuse_manifest.json').read_text())
for row in reuse['files']:
 assert sha(ROOT/row['source'])==sha(ROOT/row['snapshot'])==row['sha256'],row
stems=['contb55_census','contb55_cusp_coupled']+[f'contb55_source_case{i}'for i in range(1,6)]+['contb55_sage_finite','contb55_sage_finite_v2']
records=[]
for stem in stems:
 p=BASE/'jobs'/f'{stem}.status.json';r=json.loads(p.read_text());assert r['state']=='EXITED'and sha(BASE/r['script'])==r['script_sha256']
 role=('superseded Sage finite-field API trial; not a certificate'if stem=='contb55_sage_finite'else
       'completed corroborating cusp calculation; not needed by final proof'if stem=='contb55_cusp_coupled'else'accepted')
 if role=='accepted':assert r['actual_exit_code']==0
 records.append(dict(purpose=stem,owner_uid=r['owner_uid'],pid=r['pid'],command=r['command'],start_time_unix=r['start_time_unix'],state=r['state'],actual_exit_code=r['actual_exit_code'],termination_reason=r['termination_reason'],role=role))
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  if any(stem in cmd for stem in stems):
   fields=(p/'stat').read_text().split(') ')[1].split();active.append(dict(pid=int(p.name),owner_uid=p.stat().st_uid,command=cmd,state=fields[0],start_ticks=int(fields[19])))
 except(OSError,IndexError):pass
assert not active,active
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
write('lifecycle.json',dict(status='PASS',finished_utc=now,owned_jobs=records,
 all_owned_jobs_exited=True,active_owned_processes=active,paused_owned_workers=0,automatic_launchers=0,
 parent_native_job='contb55_parent_native',parent_native_state='EXITED0'))
write('completion.json',dict(status='PASS: C55 excluded over degree10 number fields',finished_utc=now,
 source_genus=81,source_projective_index=1440,source_rank_Q=0,
 source_closed_NC_counts=[0,0,0,0,0,0,20,0,40,100],source_diamond_orbits_by_degree={'7':1,'9':2,'10':5},
 accepted_degree_independent_A3_tests_reused=3,fresh_degree10_A3_tests=5,all_fresh_h0=0,
 independent_new_marked_cyclic3_kernels=20,independent_norm_orbits=12,
 source_model_absolute_degree_over_F2x=37,source_prime_serialization='TwoGenerators, with generated-ideal equality and exact function roundtrips',
 fresh_native_replayed_tests=5,all_native_h0=0,
 source_Q_2primary_invariants=[2,16],quotient_global_cusp_2primary_image_order=32,
 quotient_degree=5,full_source_torsion_reduction_at2_injective=True,
 source_global_and_F2_cusp_signature={'1':20,'4':5,'10':2,'20':1},
 every_source_cusp_globally_intact_at2=True,
 pure_cusp_proof='Rational torsion reduction injectivity gives D~C; gonality>=14625/1024>10 forces D=C.',
 proof_scope='Arbitrary rational effective degree10 divisors with entirely noncuspidal support; reducibility is allowed.',
 corroborating_cusp_calculation=dict(patterns=382,Cartier_exceptions=32,all_exception_h0=1,remaining=0,used_in_final_proof=False)))
files={}
def add(p):
 p=Path(p);assert p.is_file();files[str(p.relative_to(ROOT))]=dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p))
for p in H.rglob('*'):
 if p.is_file()and p.name!='certificate_manifest.json'and '__pycache__'not in p.parts:add(p)
for stem in stems+['contb55_parent_native']:
 p=BASE/'jobs'/f'{stem}.status.json';r=json.loads(p.read_text());assert sha(BASE/r['script'])==r['script_sha256'];add(p);add(BASE/r['script']);add(BASE/r['log'])
for p in(BASE/'shared/mdmagma/v2').rglob('*'):
 if p.is_file():add(p)
for row in reuse['files']:add(ROOT/row['source'])
write('certificate_manifest.json',dict(status='PASS',scope='Complete C55 degree10 exclusion with source torsion-injectivity transfer, immutable degree-independent controls, fresh census and five source A3 tests, independent finite geometry and native RR, plus preserved corroboration and job provenance',files=sorted(files.values(),key=lambda r:r['path'])))
print('C55_COMPLETE_PACKET_FROZEN',len(files),'files',len(stems),'owned jobs all EXITED')
