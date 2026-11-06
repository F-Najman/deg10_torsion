"""Freeze the complete C77 exclusion only after all independent controls pass."""
from pathlib import Path
import ast,datetime,hashlib,json,os,re
H=Path(__file__).resolve().parent;ROOT=H.parents[3];BASE=ROOT/'work/deg10'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return ast.literal_eval(Path(p).read_text().replace('[*','[').replace('*]',']'))
def write(name,data):(H/name).write_text(json.dumps(data,indent=2)+'\n')
parents=['contb77_parent_native','contb77_parent_moduli_binding']
for stem in parents:
 d=json.loads((BASE/'jobs'/f'{stem}.status.json').read_text());assert d['state']=='EXITED'and d['actual_exit_code']==0
nr=read(H/'parent_native_RR.m');assert len(nr[0])==24 and all(r[2]==0 for r in nr[0])
binding=read(H/'parent_T3_moduli_binding.m');assert len(binding)==24
for ci,n in [(1,4),(2,16),(3,4)]:
 p=read(H/f'nc_case{ci}_v3_packet.m');assert len(p[9])==n and all(r[4]==0 for r in p[9])
 assert [r[1]for r in nr[0]if r[0]==ci]==list(range(1,n+1))
 assert len([r for r in binding if r[0]==ci])==n
for name in ['independent_census.json','independent_fi10_v2.json','independent_nc_geometry.json','fi_completion.json']:
 assert json.loads((H/name).read_text())['status']=='PASS'
a=read(H/'fi_native_controls.m');assert a==[[10,27,10],[10,27,10]]
reuse=json.loads((H/'reuse_manifest.json').read_text())
for row in reuse['files']:assert sha(ROOT/row['source'])==sha(ROOT/row['snapshot'])==row['sha256']
winding=(H/'inputs/winding77.log').read_text();rows=[tuple(map(int,m.groups()))for m in re.finditer(r'order=(\d+) conductor=(\d+) relative_dim=(\d+) winding_dim=(\d+) coefficient_degree=(\d+)',winding)]
assert sum(r[2]*r[4]for r in rows)==181 and all(r[2]==r[3]for r in rows if r[0]!=1)
assert next(r for r in rows if r[0]==1)[2:]==(7,6,1)and 'POSITIVE_UPPER_T2_EXACTLY_ZERO [ 0 ]'in winding
accepted=['contb77_census','contb77_fi10','contb77_projection_data','contb77_sage_census','contb77_fi_native_controls']+[f'contb77_nc_case{i}_v3'for i in range(1,4)]+[f'contb77_sage_nc_case{i}'for i in range(1,4)]
alljobs=[];stems=[]
for p in sorted((BASE/'jobs').glob('contb77*.status.json')):
 d=json.loads(p.read_text());stem=d['purpose']
 if stem in parents:continue
 assert d['state']=='EXITED'and sha(BASE/d['script'])==d['script_sha256']
 if stem in accepted:assert d['actual_exit_code']==0;role='accepted'
 elif 'absolute_model'in stem:role='retired obsolete preflight or fallback; not a certificate'
 elif stem.endswith('_v2'):role='memory-capped default RR trial; not a certificate'
 else:role='superseded serialization preflight; not a certificate'
 alljobs.append(dict(purpose=stem,role=role,state=d['state'],actual_exit_code=d['actual_exit_code'],pid=d['pid'],owner_uid=d['owner_uid'],start_time_unix=d['start_time_unix'],command=d['command'],termination_reason=d['termination_reason']));stems.append(stem)
assert len(alljobs)==19 and len(accepted)==11
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  if any(s in cmd for s in stems):
   fields=(p/'stat').read_text().split(') ')[1].split();active.append(dict(pid=int(p.name),owner_uid=p.stat().st_uid,command=cmd,state=fields[0],start_ticks=int(fields[19])))
 except(OSError,IndexError):pass
assert not active,active
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
write('lifecycle.json',dict(status='PASS',finished_utc=now,owned_job_records=alljobs,owned_status_count=19,accepted_status_count=11,all_owned_jobs_exited=True,active_owned_processes=active,paused_owned_workers=0,automatic_launchers=0,unregistered_failed_launcher='absolute_model_v2_launcher_failure.json',absolute_fallback_retirements=['absolute_model_trial_retirement.json','absolute_fallback_retired.json'],parent_native_jobs_all_exited0=parents))
write('completion.json',dict(status='PASS: C77 excluded over degree10 number fields',finished_utc=now,
 source_level=77,source_genus=181,source_index=2880,quotient_genus=37,source_quotient_degree=5,quotient_to_X0_degree=6,
 finite_image_operator='T2',source_winding_dimension_accounted=181,nonwinding_dimension=1,T2_zero_on_nonwinding_quotient=True,
 source_closed_NC_counts=[0,0,0,0,0,0,0,0,0,90],new_source_diamond_orbits=3,source_curve_group='C1078',
 pure_cusp=dict(patterns=4992,ordinary_patterns=4989,special_5Q_patterns=3,Cartier_image='C(C-1)',image_dimension=32,restriction_rank=10,RR_blocks=69,retained_blocks_reidentified=60,new_native_geometric_controls=2),
 noncuspidal=dict(operator='A3*T2',delta2_physical_action='t->(t+1)/t, s->s, z->z+1',delta3_is_delta2_inverse=True,divisor_formula='T+2delta(T)-7P-2delta(P)-3delta^(-1)(P)',positive_degree_before_cancellation=120,AS_lift_allocations=[4,16,4],total_allocations=24,all_primary_h0=0,all_fresh_native_h0=0,positive_T3_moduli_bindings=24,independent_marked_Velu3_kernels=12,independent_marked_Velu11_quotients=15,all_norm_weights=1),
 prime_encoding='TwoGenerators with exact generated-ideal and function roundtrip checks; all native RR uses Reduction=false,Simplification=None',
 proof='continuation_20261030b/C77/resolution.tex',proof_scope='Arbitrary rational effective degree10 divisors with entirely noncuspidal support; local cubic cusp traces are allowed.'))
files={}
def add(p):
 p=Path(p);assert p.is_file();files[str(p.relative_to(ROOT))]=dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p))
for p in H.rglob('*'):
 if p.is_file()and p.name!='certificate_manifest.json'and '__pycache__'not in p.parts:add(p)
for stem in stems+parents:
 p=BASE/'jobs'/f'{stem}.status.json';d=json.loads(p.read_text());assert sha(BASE/d['script'])==d['script_sha256'];add(p);add(BASE/d['script']);add(BASE/d['log'])
add(BASE/'jobs/contb77_absolute_model_v2.log')
for p in(BASE/'shared/mdmagma/v2').rglob('*'):
 if p.is_file():add(p)
for row in reuse['files']:add(ROOT/row['source'])
write('certificate_manifest.json',dict(status='PASS',scope='Complete C77 exclusion with retained degree-independent geometry and finite-image winding, all4992 cusp patterns, all24 NC RR allocations, independent finite geometry/native RR/moduli binding, plus preserved lifecycle provenance',files=sorted(files.values(),key=lambda r:r['path'])))
print('C77_COMPLETE_PACKET_FROZEN',len(files),'files',len(alljobs),'owned jobs all EXITED')
