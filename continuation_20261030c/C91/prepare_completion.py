"""Freeze and verify the complete C91 decic proof packet."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[4]; H=Path(__file__).resolve().parent
# All paths are relative to the shared repository root.
def rel(p):return str(p.relative_to(R))
def ent(p):return dict(path=rel(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def add(p):
 p=p if p.is_absolute()else R/p;assert p.is_file();accepted.add(p)
accepted=set();diag=set()
local=['resolution.tex','relative_formal_fibre.tex','tate7_quotient_common.m','tate7_cusp_filter.m','tate7_CM_filter.m',
 'tate7_cusp_blocks.m','tate7_CM_blocks.m','tate7_differential_matrices.m','export_native_v2.m','native_packet.m',
 'parent_native_replay.m','tate7_cusp_audit_v2.py','tate7_cusp_audit.json','verify_arithmetic.py','independent_arithmetic_replay.json',
 'verify_inherited_R1.py','inherited_R1_rank_replay.json','lifecycle_and_completion.json','completion_summary.json','prepare_completion.py']
for s in local:add(H/s)
for p in (H/'inputs').rglob('*'):
 if p.is_file():add(p)
prov=json.loads((H/'inputs/inherited_provenance.json').read_text())
for a in prov['files']:
 assert hashlib.sha256((R/a['snapshot']).read_bytes()).hexdigest()==a['sha256']
 assert hashlib.sha256((R/a['source']).read_bytes()).hexdigest()==a['sha256']
cmfile=R/'work/deg10/continuation_20261030c/C91_CM_control/level7_binding_manifest.json';cm=json.loads(cmfile.read_text());assert cm['status']=='PASS';add(cmfile)
for a in cm['files']:
 q=R/a['path'];assert hashlib.sha256(q.read_bytes()).hexdigest()==a['sha256'];add(q)
for stem,marker in [('contc91_tate7_cusps','C91_TATE7_CUSP_FILTER_COMPLETE patterns=120 exceptions=90'),
 ('contc91_tate7_CM','C91_TATE7_CM_FILTER_COMPLETE tested=12 exceptions=[]'),
 ('contc91_export_native_v2','C91_NATIVE_PACKET_EXPORT_COMPLETE primes=6 kernels=45'),
 ('contc91_parent_native_v3','C91_PARENT_NATIVE_COMPLETE 6 45 120 12 4.450')]:
 for q in [R/'work/deg10/computations'/f'{stem}.m',R/'work/deg10/jobs'/f'{stem}.log',R/'work/deg10/jobs'/f'{stem}.status.json']:add(q)
 x=json.loads((R/'work/deg10/jobs'/f'{stem}.status.json').read_text());assert x['state']=='EXITED'and x['actual_exit_code']==0
 assert marker in (R/'work/deg10/jobs'/f'{stem}.log').read_text()
for q in ['work/deg10/continuation_20261030c/parent_source_census.json','work/deg10/computations/contc_parent_source_census.py',
 'work/deg10/jobs/contc_parent_source_census.log','work/deg10/jobs/contc_parent_source_census.status.json',
 'work/deg10/next_cases_20261030/dyadic_local/all_curves.m','work/deg10/next_cases_20261030/dyadic_local/independent_point_counts.json',
 'work/deg10/next_cases_20261030/dyadic_formal_kernel.tex','work/deg10/run_job.py']:add(Path(q))
for q in H.rglob('*'):
 if q.is_file()and q.name!='certificate_manifest.json'and q not in accepted:diag.add(q)
# Preserve all failed/retired/backup job diagnostics without using them as proof evidence.
for directory,pattern in [('work/deg10/computations','contc91*'),('work/deg10/jobs','contc91*')]:
 for q in (R/directory).glob(pattern):
  if q.is_file()and q not in accepted:diag.add(q)
assert not (accepted&diag)
life=json.loads((H/'lifecycle_and_completion.json').read_text());assert life['status']=='PASS'and not life['active_owned_magma']
manifest=dict(status='PASS',scope='Complete degree10 C91 exclusion: genus23 relative full-Jacobian dyadic formal fibres (120 cuspidal,12 CM cases), exact source CM binding, independent native45-kernel replay, and inherited ordinary-support canonical nonprincipality with all288 ranks independently replayed',
 files=[ent(q)for q in sorted(accepted)],diagnostic_scope='Preserved exploratory C4 computations, retired fallback, failed serializers/parsers, and full-source backup; none is required mathematical proof evidence',
 diagnostic_files=[ent(q)for q in sorted(diag)])
(H/'certificate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for a in manifest['files']+manifest['diagnostic_files']:
 assert hashlib.sha256((R/a['path']).read_bytes()).hexdigest()==a['sha256']
print('C91_CERTIFICATE_FROZEN',len(manifest['files']),len(manifest['diagnostic_files']),hashlib.sha256((H/'certificate_manifest.json').read_bytes()).hexdigest())
