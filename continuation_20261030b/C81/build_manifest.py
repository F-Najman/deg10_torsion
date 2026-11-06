"""Validate and hash the final C81 packet; no Magma jobs are launched."""
from pathlib import Path
import ast, hashlib, json, os

REPO=Path.cwd()
HOME=REPO/'work/deg10'
ROOT=HOME/'continuation_20261030b/C81'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))

markers={
 'b81_source_census':'C81_DYADIC_SOURCE_CENSUS_COMPLETE',
 'b81_x0_setup':'C81_X0_PROVEN_PICARD_COMPLETE',
 'b81_x0_filter':'C81_X0_CLASS_FILTER_COMPLETE',
 'b81_y46_kummer_model_v2':'C81_Y46_KUMMER_MODEL_COMPLETE',
 'b81_parent_source_filter':'C81_PARENT_INDEPENDENT_SOURCE_FILTER_PASS',
 'b81_parent_x0_filter':'C81_PARENT_NATIVE_X0_FILTER_COMPLETE',
 'b81_parent_two_torsion':'C81_PARENT_RATIONAL_TWO_TORSION_ZERO_PASS',
 'b81_parent_rr_completions_v2':'C81_PARENT_COMPLETION_RR_V2_COMPLETE',
 'b81_sage_local_arcs_v3':'C81_SAGE_LOCAL_ARCS_ALL_29_PASS',
}
for i in [1,2,3,7,10]:
 markers['b81_rr_case%d'%i]='C81_Y46_RR_JET_COMPLETE'
 markers['b81_rr_localmeta%d'%i]='C81_RR_LOCAL_ABSOLUTE_META_COMPLETE'
diagnostics={
 'b81_parent_rr_completions':'Magma parser diagnostic; actual exit 0 is not a certificate',
 'b81_y46_kummer_model':'Magma serialization/parser diagnostic; actual exit 0 is not a certificate',
 'b81_y46_model':'Absolute genus-46 normalization reached explicit memory bound; exited 1; replaced by compact Kummer/character method',
 'b81_y46_model_v2':'Absolute genus-46 normalization reached explicit memory bound; exited 1; replaced by compact Kummer/character method',
 'b81_sage_local_arcs':'Finite-field API diagnostic; exited 1; superseded by successful v3',
 'b81_sage_local_arcs_v2':'NTL Laurent-series coercion diagnostic; exited 1; superseded by successful v3',
}
jobs=[]
for p in sorted((HOME/'jobs').glob('b81*.status.json')):
 a=json.loads(p.read_text());name=a['purpose'];script=HOME/a['script'];log=HOME/a['log']
 assert a['state']=='EXITED',a
 assert digest(script)==a['script_sha256'],name
 text=log.read_text()
 if name in markers:
  assert a['actual_exit_code']==0 and markers[name] in text,name
  assert not any(t in text for t in ['User error:','Runtime error','Traceback (most recent call last):','System Error:']),name
  category='successful certificate'
 elif name in diagnostics:category='unsuccessful diagnostic, not certificate'
 else:category='completed auxiliary assessment, not required for final proof'
 jobs.append(dict(purpose=name,category=category,state=a['state'],actual_exit_code=a['actual_exit_code'],
                  expected_completion_marker=markers.get(name),diagnostic=diagnostics.get(name),
                  elapsed_seconds=a['elapsed_seconds'],script=str(script.relative_to(REPO)),
                  status=str(p.relative_to(REPO)),log=str(log.relative_to(REPO))))
assert set(markers)<=set(a['purpose'] for a in jobs)

# Inspect actual command lines, avoiding both stale PID assumptions and textual
# matches against a shell that happens to mention a script name.
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  cmd=(p/'cmdline').read_bytes().split(b'\0');argv=[v.decode(errors='replace') for v in cmd if v]
  if not argv:continue
  executable=Path(argv[0]).name
  relevant=(executable.startswith('magma') and any('b81_' in v for v in argv[1:])) or \
           (executable.startswith('python') and any(Path(v).name.startswith('b81_') for v in argv[1:])) or \
           (executable.startswith('python') and any(Path(v).name=='run_job.py' for v in argv[1:]) and any('b81_' in v for v in argv[1:]))
  if relevant and int(p.name)!=os.getpid():active.append({'pid':int(p.name),'argv':argv})
 except (OSError,ProcessLookupError):pass
assert not active,active

src=read(ROOT/'source81.m');assert len(src[2])==10
assert [r[5] for r in src[1]]==[0]*8+[135,135]
x0=read(ROOT/'x0_filter.m');assert x0[2]==[1,2,3,7,10]
tests=[]
for i in [1,2,3,7,10]:tests+=read(ROOT/('rr_case%d_summary.m'%i))
assert len(tests)==29 and all(r[6]==1 for r in tests)
native=read(ROOT/'parent_native_completion_RR_v2.m');assert native[0]==tests and native[1:]==[87,14]
sage=json.loads((ROOT/'independent_sage_local_arcs.json').read_text());assert sage['status']=='PASS' and len(sage['tests'])==29
assert all(a['matrix_entrywise_equal'] and a['h0']==1 for a in sage['tests'])
sourcecontrol=json.loads((ROOT/'parent_independent_source_filter.json').read_text());assert sourcecontrol['status']=='PASS' and sourcecontrol['all_cusps_globally_intact']
parentchecks=json.loads((ROOT/'parent_checks.json').read_text())
for item in parentchecks['evidence']:
 assert digest(HOME/item['path'])==item['sha256'],item['path']
for item in json.loads((ROOT/'reuse_manifest.json').read_text()):
 assert digest(REPO/item['source'])==item['sha256'] and digest(REPO/item['copy'])==item['sha256'],item['copy']
summary={
 'status':'COMPLETE: C81 excluded over degree-10 number fields',
 'degree':10,'group':'C81','source_genus':190,'source_rank':0,'source_rational_2_torsion_dimension':0,
 'finite_NC_closed_counts_through10':[0]*8+[135,135],
 'source_diamond_orbits_degree9':5,'source_diamond_orbits_degree10':5,
 'source_candidate_supports':10,'after_proven_X0_81_odd_coset_sieve':[1,2,3,7,10],
 'RR_target':'X1(81)/<28>','RR_target_genus':46,'base_curve':'X1(27)','base_genus':13,
 'Kummer_extension_degree':3,'rational_branched_cusp_count':9,
 'RR_tests':29,'native_base_RR_spaces':87,'native_genuine_completions':14,
 'all_RR_dimensions':1,'native_completion_replay':'PASS','independent_Sage_Hensel_arcs':'PASS',
 'parent_nested_input_hashes_verified':len(parentchecks['evidence']),
 'Sage_series_precision':64,'remaining_supports':[],
 'lifecycle':{'all_records_exited':True,'active_C81_workers_or_supervisors':active,'obsolete_paused_workers':0,'automatic_launchers':0,'diagnostics_preserved':list(diagnostics)},
 'primary_proof':'resolution.tex',
}
(ROOT/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

paths={p for p in ROOT.rglob('*') if p.is_file() and p.name!='certificate_manifest.json' and '__pycache__' not in p.parts}
paths.update((HOME/'computations').glob('b81_*'))
paths.update((HOME/'jobs').glob('b81*.log'))
paths.update((HOME/'jobs').glob('b81*.status.json'))
paths.update(p for p in (HOME/'shared/mdmagma/v2').rglob('*') if p.is_file())
paths.add(HOME/'run_job.py')
files=[{'path':str(p.relative_to(REPO)),'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(paths)]
manifest={'status':'PASS','format':'C81 degree-10 computational certificate packet v1',
          'summary_sha256':digest(ROOT/'certificate_summary.json'),'file_count':len(files),'files':files,'jobs':jobs}
(ROOT/'certificate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('C81_CERTIFICATE_PACKET_HASHED',len(files),'files;',len(markers),'successful jobs;',len(jobs),'total exited records')
print('manifest_sha256',digest(ROOT/'certificate_manifest.json'))
