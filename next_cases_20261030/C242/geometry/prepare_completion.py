"""Package accepted C242 evidence; do not run mathematical computations."""
import datetime, hashlib, json, os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'work/deg10'
PACKET=BASE/'next_cases_20261030/C242'
HERE=PACKET/'geometry'
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,data): (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
accepted_stems=[
 'next242_geom_g5_setup_v2','next242_geom_g5_s0','next242_geom_g5_s1',
 'next242_geom_g5_s2','next242_geom_g5_merge','next242_geom_g5_sage_replay_v2',
 'next242_geom_c21_setup','next242_geom_c21_full','next242_geom_c21_sage_replay_v2',
 'next242_geom_x0_42_setup','next242_geom_x0_42_cusp_audit',
 'next242_geom_x0_42_coordinates','next242_geom_x0_42_full',
 'next242_geom_x0_42_sage_replay_v2',
 'next242_geom_case26_sections','next242_geom_case26_sage_replay']
all_jobs=[]
for p in sorted((BASE/'jobs').glob('next242_geom*.status.json')):
 d=json.loads(p.read_text());assert d['state']=='EXITED', (p,d)
 assert digest(BASE/d['script'])==d['script_sha256'],('modified historical script',p)
 label=('accepted' if d['purpose'] in accepted_stems else
        'abandoned bounded experiment; not a completed certificate' if d['actual_exit_code']==-15 else
        'superseded or benchmark; not used as a certificate')
 if d['purpose'] in accepted_stems: assert d['actual_exit_code']==0,(p,d)
 all_jobs.append(dict(purpose=d['purpose'],pid=d['pid'],owner_uid=d['owner_uid'],
   start_time_unix=d['start_time_unix'],command=d['command'],state=d['state'],
   actual_exit_code=d['actual_exit_code'],termination_reason=d['termination_reason'],
   disposition=label,script=d['script'],status=str(p.relative_to(ROOT)),
   log='work/deg10/'+d['log']))
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit(): continue
 try:
  if p.stat().st_uid!=os.getuid(): continue
  cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace').strip()
  if 'next242_geom' not in cmd: continue
  stat=(p/'stat').read_text();tail=stat[stat.rindex(')')+2:].split()
  active.append(dict(pid=int(p.name),owner_uid=p.stat().st_uid,command=cmd,
                     state=tail[0],start_ticks=int(tail[19])))
 except (OSError,ValueError): pass
assert not active, active
write('lifecycle_and_completion.json',dict(status='PASS',finished_utc=now,
 scope='Owned C242 geometry workers, supervisors and launchers',jobs=all_jobs,
 owned_worker_count=len(all_jobs),accepted_job_count=len(accepted_stems),
 all_owned_jobs_exited=True,active_owned_processes=active,paused_owned_workers=0,
 automatic_launchers=0,notes=[
 'The two V9 model-normalization trials exited -15 at their 300-second limits. They are not proof certificates.',
 'The successful X0(42) route superseded the V9 trials; no fallback is left paused or pending.',
 'Earlier failed parsing/metadata trials and successful benchmarks are preserved, but are not used as accepted certificates.']))
summary=dict(status='PASS: C2 x C42 and C84 excluded in degree 10',finished_utc=now,
 degree=10,groups=['C2 x C42','C84'],source_genus=61,source_geometric_gonality_lower_bound=12,
 proof='next_cases_20261030/C242/geometry/resolution.tex',
 support_pipeline=[dict(stage=s,patterns=n) for s,n in [
 ('complete normalized source census',60971),('two genus4 quotients',4647),
 ('genus5 C3xC3 quotient',752),('four X1(21) images',111),
 ('three order3 stabilizer exclusions',108),('six X0(42) images',10),
 ('nine exact full X1(42) Hecke exclusions',1),('final genus17 orbit-section obstruction',0)]],
 small_quotient_supported_prime_images=5346,small_quotient_exact_native_RR_identities=5057,
 small_quotient_smooth_factorbase_controls=289,unsupported_genus5_j0_rows_retained=129,
 final_Hecke_source_summands=16,final_Hecke_kernel_orbits=78,final_Hecke_geometric_kernels=240,
 final_Hecke_RR_zero_tests=9,case26=dict(chosen_h='(P+21Q,P+4Q)',h_Frobenius_phase=2,
 old_genus17_generator='h^2=(21Q,P+37Q)',old_generator_Frobenius_phase=4,
 source_C4_options=108,quotient_cusp_types=12,quadratic_quotient_candidates=2,
 RR_patterns=24,character_spaces=72,all_character_dimensions=[1,0,0]),
 proof_scope='Arbitrary rational effective degree-10 divisors with entirely noncuspidal support; reducibility is allowed.',
 global_cuspidal_descent='DEHMZ Theorems 3.1 and 4.13: global Galois-invariant cuspidal divisors generate J1(42)(Q). Finite cubic source cusps are never silently substituted for global rational divisors.',
 native_replays=['parent_native_replays.json','parent_final_checks.json'],
 lifecycle='lifecycle_and_completion.json')
write('completion.json',summary)
files={}
def add(p,expected=None):
 p=Path(p)
 if not p.is_absolute(): p=ROOT/p
 assert p.is_file(),p
 h=digest(p)
 if expected is not None: assert h==expected,(p,h,expected)
 files[str(p.relative_to(ROOT))]=dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=h)
def nested_manifest(p):
 add(p);d=json.loads(p.read_text())
 for row in d['files']: add(row['path'],row['sha256'])
for name in [
 'resolution.tex','quotient5_geometry.tex','x0_42_image.tex','case26_section_obstruction.tex',
 'quotient5_proven_picard_v2.m','quotient5_A11_numeric.m','quotient5_A11_arithmetic_v2.m',
 'quotient5_helpers.m','quotient5_filter_common_v2.m','genus5_merged_arithmetic.m',
 'genus5_independent_replay.json','replay_genus5_v2.py',
 'quotient5_setup_v2_RR.log','quotient5_setup_v2_RR.log.witnesses.py',
 'cyclic21_proven_picard.m','cyclic21_filter_common.m','cyclic21_full_arithmetic.m',
 'cyclic21_full_RR.log','cyclic21_full_RR.log.witnesses.py',
 'cyclic21_cusp_RR.log','cyclic21_cusp_RR.log.witnesses.py',
 'cyclic21_independent_replay.json','replay_cyclic21_v2.py',
 'x0_42_proven_picard.m','x0_42_proven_picard_v2.m','x0_42_filter_common.m',
 'x0_42_full_arithmetic.m','x0_42_full_RR.log','x0_42_full_RR.log.witnesses.py',
 'x0_42_setup_RR.log','x0_42_setup_RR.log.witnesses.py',
 'x0_42_direct_cusp_pushforward.m','x0_42_independent_replay.json','replay_x0_42_v2.py',
 'DEHMZ_global_cuspidal_definition_and_theorem.tex','diagonal17_normalized_model_data.m',
 'diagonal17_parameter_data.m','case26_section_model.m','case26_section_RR.m',
 'case26_section_independent_replay.json','replay_case26_sections.py',
 'reuse_manifest.json','prepare_completion.py','lifecycle_and_completion.json','completion.json']:
 add(HERE/name)
for i in range(3):
 for suffix in ['arithmetic.m','RR.log','RR.log.witnesses.py']: add(HERE/f'g5_s{i}_{suffix}')
for folder in ['deps','accepted_quotient5']:
 for p in (HERE/folder).rglob('*'):
  if p.is_file(): add(p)
for stem in accepted_stems:
 p=BASE/'jobs'/f'{stem}.status.json';d=json.loads(p.read_text())
 add(p);add(BASE/d['log']);add(BASE/d['script'],d['script_sha256'])
for p in (BASE/'shared/mdmagma/v2').rglob('*'):
 if p.is_file(): add(p)
for manifest in [PACKET/'certificate_manifest.json',PACKET/'direct_hecke17/certificate_manifest.json']:
 nested_manifest(manifest)
for name in [
 'first_stage_completion.json','first_survivors_unified.m','first_survivors_details.m',
 'CM_j0_degree10.tex','order3_stabilizer_degree10.tex','after_cheap_stabilizers_unified.m',
 'cheap_stabilizers.json','degree9_order3_native_phases.m','case26_diagonal_phases.m',
 'case26_cubic_complement.tex','case26_cubic_complement.json','case26_cusp_omega_types.json',
 'direct_hecke17/nine_exclusions_completion.json']:
 add(PACKET/name)
for name in ['parent_native_replays.json','parent_final_checks.json']:
 p=HERE/name;add(p)
 for row in json.loads(p.read_text())['evidence']: add(BASE/row['path'],row['sha256'])
# Lifecycle records preserve failed/obsolete scripts and logs without asserting certificates.
for d in all_jobs:
 for key in ['script','log','status']:
  p=BASE/d[key] if key=='script' else ROOT/d[key]
  add(p)
write('certificate_manifest.json',dict(status='PASS',scope='Complete C242/C84 exclusion, including source census, all quotient and full Hecke controls, final genus17 section obstruction, and preserved job provenance',
 files=sorted(files.values(),key=lambda d:d['path'])))
print('C242_GEOMETRY_PACKET_COMPLETE',len(files),'files',len(all_jobs),'EXITED jobs',len(accepted_stems),'accepted geometry jobs')
