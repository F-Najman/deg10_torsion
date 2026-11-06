from pathlib import Path
from collections import Counter
import json,hashlib,os,time
BASE=Path(__file__).resolve().parent;HOME=BASE.parent;REPO=HOME.parent.parent
GROUPS=[48,54,56,60,66]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return json.loads(p.read_text().replace('[*','[').replace('*]',']'))
markers={'d31a_even_p5_counts':'D31A_P5_SOURCE_FROBENIUS_COMPLETE','d31a_even_p5_torsion_low':'D31A_P5_TORSION_RUN_FINISHED','d31a_C54_target27_setup_v3':'D31A_SMALLTARGET_PROVEN_COMPLETE','d31a_C54_low_actual':'D31A_LOW_ACTUAL_SOURCE_COMPLETE','d31a_C66_low_actual':'D31A_LOW_ACTUAL_SOURCE_COMPLETE','d31a_C56_target28_setup':'D31A_C56_TARGET28_PROVEN_COMPLETE','d31a_C56_fourmap_filter_v4':'D31A_C56_FOURMAP_FILTER_COMPLETE','d31a_C56_independent_fourmap_v3':'D31A_C56_INDEPENDENT_FOURMAP_COMPLETE'}
for i in range(6):markers[f'd31a_even_p5_torsion_high_s{i}']='D31A_P5_TORSION_RUN_FINISHED'
for name in['d31a_even_p3_low','d31a_even_p3_d8','d31a_even_p3_d9','d31a_even_p3_d10']:markers[name]='D31A_P3_ALL_SOURCE_CENSUS_COMPLETE'
markers['d31a_even_p3_counts']='D31A_SOURCE_FROBENIUS_COMPLETE'
def ours(name):return name.startswith(('d31a_even_p3_','d31a_even_p5_','d31a_C48_','d31a_C54_','d31a_C56_','d31a_C60_','d31a_C66_'))
jobs=[]
for p in sorted((HOME/'jobs').glob('*.status.json')):
 name=p.name.removesuffix('.status.json')
 if not ours(name):continue
 a=json.loads(p.read_text());assert a['state']=='EXITED',(name,a['state']);assert sha(HOME/a['script'])==a['script_sha256'];log=(HOME/a['log']).read_text();success=name in markers
 if success:assert a['actual_exit_code']==0 and markers[name]in log and not any(z in log for z in['User error:','Runtime error','Traceback (most recent call last):']),name
 if name=='d31a_C66_target33_setup_v3':success=a['actual_exit_code']==0 and 'D31A_SMALLTARGET_PROVEN_COMPLETE'in log
 if name=='d31a_C56_A5_pure':success=a['actual_exit_code']==0 and 'D31A_C56_A5_PURE_COMPLETE'in log
 jobs.append({'purpose':name,'state':a['state'],'actual_exit_code':a['actual_exit_code'],'complete_certificate':success,'partial_or_diagnostic':not success,'script':str((HOME/a['script']).relative_to(REPO)),'log':str((HOME/a['log']).relative_to(REPO)),'status':str(p.relative_to(REPO))})
active=[]
for pp in Path('/proc').iterdir():
 if not pp.name.isdigit():continue
 try:
  cmd=[v.decode(errors='replace')for v in(pp/'cmdline').read_bytes().split(b'\0')if v]
  if cmd and(Path(cmd[0]).name.startswith('magma')or Path(cmd[0]).name.startswith('python'))and any(ours(Path(v).stem)for v in cmd[1:]):active.append({'pid':int(pp.name),'argv':cmd})
 except OSError:pass
assert not active,active
for n in GROUPS:
 for rec in json.loads((BASE/f'C{n}/reuse_manifest.json').read_text()):assert sha(REPO/rec['source'])==sha(REPO/rec['copy'])==rec['sha256']
for n in[48,54,66]:
 summ=json.loads((BASE/f'C{n}/source_census_verified.json').read_text());src=rd(BASE/f'C{n}/source{n}_p5.m');assert summ['status']=='PASS'and src[1]==summ['noncuspidal_closed_counts']and len(src[2])==summ['source_diamond_orbits_total']
 for rec in summ['inputs']:assert sha(REPO/rec['path'])==rec['sha256']
a=rd(BASE/'C56/fourmap_filter.m');ss=json.loads((BASE/'C56/independent_fourmap.json').read_text());assert ss['status']=='PASS'and len(a[2])==1229 and len(a[3])==135==ss['remaining_support_patterns']
summary={'status':'PARTIAL BACKUP WORK; NO ADDITIONAL GROUP EXCLUSION CLAIMED','groups':GROUPS,'complete_source_census_groups':[48,54,56,66],'C54_proven_rational_gonality_lower_bound':12,'C54_pure_cusp_and_bad_dyadic_boundary_resolved':True,'C56_fourmap_tested_patterns':1229,'C56_fourmap_remaining_patterns':135,'C56_independent_physical_maps':ss['genuine_marked_projection_maps'],'C56_remaining_pattern_types':ss['remaining_pattern_types'],'complete_C70_exclusion_separate_frozen_packet':True,'lifecycle':{'active_owned_workers':active,'paused_workers':0,'automatic_launchers':0,'all_recorded_jobs_exited':True},'jobs':jobs}
(BASE/'even_backup_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={BASE/'freeze_even_backups.py',BASE/'merge_F5_torsion_census.py',BASE/'even_backup_summary.json',HOME/'run_job.py'}
for n in GROUPS:paths.update(p for p in(BASE/f'C{n}').rglob('*')if p.is_file()and'__pycache__'not in p.parts)
paths.update(p for p in(HOME/'computations').glob('d31a_*')if ours(p.stem)or p.name in['d31a_low_actual_common.m','d31a_smalltarget_setup_common.m','d31a_smalltarget_setup_common_v2.m','d31a_smalltarget_setup_common_v3.m','d31a_54_parent_gonality_controls.py'])
for rec in jobs:
 for k in['script','log','status']:paths.add(REPO/rec[k])
paths.update(p for p in(HOME/'shared/mdmagma/v2').rglob('*')if p.is_file())
files=[{'path':str(p.relative_to(REPO)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
out={'status':'PASS; PARTIAL MATHEMATICAL WORK','summary_sha256':sha(BASE/'even_backup_summary.json'),'files':files,'jobs':jobs}
(BASE/'even_backup_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print('D31A_EVEN_BACKUP_MANIFEST_COMPLETE',len(files),len(jobs),sha(BASE/'even_backup_manifest.json'))
