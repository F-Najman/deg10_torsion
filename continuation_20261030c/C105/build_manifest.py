"""Hash and verify this completed degree-ten exclusion packet."""
from pathlib import Path
import ast,json,hashlib,os
ROOT=Path(__file__).resolve().parent
N=int(ROOT.name[1:]);HOME=ROOT.parents[1];REPO=HOME.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
common={'contc_parent_source_census':'CONT_C_PARENT_SOURCE_CENSUS_COMPLETE'}
if N==65:
 markers={'contc65_source_census_v2':'C65_DYADIC_SOURCE_CENSUS_COMPLETE','contc65_R1_sieve':'DECIC65_R1_SIEVE_RESULT 240 excluded 0 unresolved','contc65_verify_jets_v2':'DECIC65_ORDINARY_STAGE2_PASS','contc65_verify_velu':'DECIC65_SAGE_KERNEL_ISOGENY_TRANSPORT_PASS','contc65_verify_m4ri':'DECIC65_R1_M4RI_RANK_PASS cases240','contc65_source_binding':'C65_DECIC_PHYSICAL_SOURCE_BINDING_PASS','contc65_inherited_deep':'C65_INHERITED_DEGREE_INDEPENDENT_EVIDENCE_ALL_PASS','contc65_parent_prefix_replay':'C65_PARENT_PREFIX_REPLAY_COMPLETE'}
 for i in range(9,14):markers['contc65_ordinary10_'+str(i)]='C65_DECIC_ELLIPTIC_JETS_COMPLETE'
 for tag in ['n1','n2','n3','n4plus']:markers['contc65_cusp_FI_'+tag]='C65_DECIC_CUSP_FI_SHARD_PASS'
 diagnostics={'contc65_source_census':'Initial wrong new-directory path; exit1, replaced by source_census_v2','contc65_verify_jets':'Inherited degree whitelist excluded10/20; exit1, replaced by verify_jets_v2; all arithmetic checks unchanged'}
 expected=[0,0,0,0,0,2,48,48,72,108]
else:
 assert N==105
 markers={'contc105_source_census':'C105_DYADIC_SOURCE_CENSUS_COMPLETE','contc105_quotient35_RR':'C105_TARGET35_A3_RR_COMPLETE','contc105_kernel_export':'C105_EXPLICIT_KERNEL_EXPORT_COMPLETE','contc105_sage_velu':'C105_SAGE_VELU_ALL_SIX_PASS','contc105_parent_native':'C105_PARENT_NATIVE_COMPLETE'}
 diagnostics={};expected=[0]*9+[144]
markers.update(common)
source=rd(ROOT/('source'+str(N)+'.m'));assert [r[5]for r in source[1]]==expected
independent=json.loads((ROOT/'inputs/parent_source_census.json').read_text());assert independent['status']=='PASS' and independent['cases'][str(N)]['nc_closed_counts']==expected
for nm in ['reuse_manifest.json']+(['inherited_snapshots.json']if N==65 else[]):
 for v in json.loads((ROOT/nm).read_text()):
  assert sha(REPO/v['copy'])==v['sha256'] and sha(REPO/v['source'])==v['sha256'],v
if N==65:
 r=json.loads((ROOT/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==240 and not r['unresolved'] and all(v['rank']==121 for v in r['results'])
 assert json.loads((ROOT/'independent_m4ri.json').read_text())['all_ranks']==121
 assert json.loads((ROOT/'independent_source_binding.json').read_text())['diamond_orbit_sizes']==[24,24,24,12,24]
 assert json.loads((ROOT/'stage2_independent_verified.json').read_text())['ordinary65_torsion_families']==15
 assert json.loads((ROOT/'hecke_transport_independent_verified.json').read_text())['independent_actual3_isogenies']==20
 f=[json.loads((ROOT/('cusp_FI_'+tag+'.json')).read_text())for tag in ['n1','n2','n3','n4plus']];assert sum(v['all_normalized_configurations']for v in f)==7170048 and all(v['all_FI_ranks_equal10']for v in f)
 pp=json.loads((ROOT/'parent_prefix_replay.json').read_text());assert pp['status']=='PASS' and pp['total_configurations']==7170048 and pp['all_prefix_matrices_have_rank10']
 assert sha(HOME/'computations/contc65_parent_prefix_native.cpp')==pp['cpp_sha256'] and sha(ROOT/'parent_prefix_vectors.txt')==pp['vector_data_sha256']
 summary={'status':'COMPLETE EXCLUSION','degree':10,'group':'C65','source_genus':121,'full_source_rank_zero_assumed':False,'finite_NC_closed_counts':expected,'new_degree10_orbit_sizes':[24,24,24,12,24],'new_RR_tests':240,'RR_domain_dimension':121,'RR_product_dimension':360,'all_final_ranks':121,'all_new_positive_negative_degrees':70,'new_genuine_jet_families':15,'independent_actual3_isogenies':20,'new_global_cusp_prefix_tests':7170048,'all_global_cusp_FI_ranks':10,'parent_independent_full_prefix_replay':'PASS','old417_pinned_files_and_full_deep_replay':'PASS','non_global_sextic_cusp_obstruction':'same norm6 class obstruction, globally intact degree4 complement','remaining_cases':0}
else:
 rr=rd(ROOT/'quotient35_RR.m');assert len(rr[2])==6 and all(v[6]==0 for v in rr[2])
 native=rd(ROOT/'parent_native_RR.m');assert len(native)==6 and all(v[4]==0 for v in native)
 vel=json.loads((ROOT/'independent_sage_velu.json').read_text());assert vel['status']=='PASS'and vel['actual_kernel_orbits']==12 and vel['exact_native_RR_prime_bindings']
 summary={'status':'COMPLETE EXCLUSION','degree':10,'group':'C105','source_genus':289,'full_source_rank_zero_assumed':False,'finite_image_operator':'diamond22-1','pure_cusp_comparison_degree':40,'gonality_bound':[1008,25],'free_diamond_order':4,'finite_NC_closed_counts':expected,'source_diamond_orbits':6,'RR_curve':'X1(35)','RR_curve_genus':25,'RR_curve_rank':0,'RR_tests':6,'all_RR_h0':0,'all_native_primes_rebuilt':24,'independent_actual_kernel_orbits':12,'all_source_and_marked_transport_bindings':'PASS','remaining_cases':0}
records=[]
statuspaths=set((HOME/'jobs').glob('contc'+str(N)+'_*.status.json'))|{HOME/'jobs/contc_parent_source_census.status.json'}
for p in sorted(statuspaths):
 a=json.loads(p.read_text());name=a['purpose'];assert a['state']=='EXITED',name
 assert sha(HOME/a['script'])==a['script_sha256'],name
 text=(HOME/a['log']).read_text()
 if name in markers:
  assert a['actual_exit_code']==0 and markers[name]in text,name
  assert not any(v in text for v in ['User error:','Runtime error','System Error:','Traceback (most recent call last):']),name
  category='successful certificate'
 else:assert name in diagnostics,name;category='unsuccessful diagnostic, not certificate'
 records.append({'purpose':name,'category':category,'state':'EXITED','actual_exit_code':a['actual_exit_code'],'elapsed_seconds':a['elapsed_seconds'],'required_marker':markers.get(name),'diagnostic':diagnostics.get(name),'status':str(p.relative_to(REPO)),'script':a['script'],'log':a['log']})
assert set(markers)=={r['purpose']for r in records if r['category']=='successful certificate'}
active=[]
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:
  argv=[v.decode(errors='replace')for v in(proc/'cmdline').read_bytes().split(b'\0')if v]
  if not argv:continue
  exe=Path(argv[0]).name
  relevant=(exe.startswith('magma')or exe.startswith('python'))and any(Path(v).name.startswith('contc'+str(N)+'_')for v in argv[1:])
  if relevant and int(proc.name)!=os.getpid():active.append({'pid':int(proc.name),'argv':argv})
 except(OSError,ProcessLookupError):pass
assert not active,active
summary['lifecycle']={'all_new_job_records_exited':True,'successful_jobs':len(markers),'unsuccessful_diagnostics':diagnostics,'active_workers_supervisors_or_launchers':active,'paused_workers':0,'automatic_launchers':0}
(ROOT/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={p for p in ROOT.rglob('*')if p.is_file()and p.name!='certificate_manifest.json'and'__pycache__'not in p.parts and not p.read_bytes().startswith(b'\x7fELF')}
paths.update((HOME/'computations').glob('contc'+str(N)+'_*'))
paths.add(HOME/'computations/contc_parent_source_census.py');paths.add(HOME/'run_job.py')
for p in statuspaths:
 paths.add(p);a=json.loads(p.read_text());paths.add(HOME/a['log']);paths.add(HOME/a['script'])
paths.update(p for p in(HOME/'shared/mdmagma/v2').rglob('*')if p.is_file())
files=[{'path':str(p.relative_to(REPO)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
out={'status':'PASS','group':'C'+str(N),'degree':10,'file_count':len(files),'summary_sha256':sha(ROOT/'certificate_summary.json'),'files':files,'new_jobs':records,'historical_evidence_lifecycle':'Copied historical records preserve their original outcomes; actual new computations and bounded supervisors are separately verified exited. No historical launcher is invoked.'}
(ROOT/'certificate_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
for item in files:assert sha(REPO/item['path'])==item['sha256']
print('C'+str(N)+'_FINAL_PACKET_PASS',len(files),'files',len(markers),'successful jobs',len(diagnostics),'diagnostics')
print('manifest_sha256',sha(ROOT/'certificate_manifest.json'))
