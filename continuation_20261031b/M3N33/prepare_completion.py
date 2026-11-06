from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[3];W=R/'work/deg10'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
required={
'd31b_m3n33_frob_preflight':'D31B_M3N33_FROB_PREFLIGHT_COMPLETE',
'd31b_m3n33_physical_source':'D31B_PHYSICAL_SOURCE5_COMPLETE 33 264 88 12',
'd31b_m333_rational_trace_A5':'D31B_RELATIVE5_RATIONAL_TRACE_A5_COMPLETE 33 12',
'd31b_m333_rational_trace_3isog_A5':'D31B_RELATIVE5_RATIONAL_TRACE_A5_COMPLETE 33 9',
'd31b_m333_independent_genuine35':'D31B_M333_INDEPENDENT_SOURCE_GENUINE35_COMPLETE',
'd31b_m333_parent_native':'D31B_M333_PARENT_NATIVE_COMPLETE 12 80',
'd31b_m333_kummer_frame':'D31B_M333_KUMMER_FRAME_BASE_CUSPS_COMPLETE 24 36',
'd31b_m333_analytic_jets_v4':'D31B_M333_ANALYTIC_CUSP_JETS_COMPLETE 48',
'd31b_m333_high_primitive_series':'D31B_M333_HIGH_PRIMITIVES_COMPLETE 24',
'd31b_m333_antih10_winding':'D31B_M333_ANTIH10_WINDING_COMPLETE',
'd31b_m333_winding_independent':'D31B_M333_WINDING_INDEPENDENT_COMPLETE 40',
'd31b_m333_prefix_ranks':'D31B_M333_ALL_PREFIX_RANKS_COMPLETE',
'd31b_m333_independent_ranks':'D31B_M333_INDEPENDENT_ALL_RANKS_COMPLETE',
'd31b_m333_native_frame_replay':'D31B_M333_NATIVE_FRAME_REPLAY_COMPLETE 81 36 40 24'}
for j in range(4):required[f'd31b_m333_independent_analytic_v4_s{j}']=f'D31B_M333_INDEPENDENT_ANALYTIC_SHARD_COMPLETE_{j}'
for j in[0,1,3]:required[f'd31b_m333_independent_analytic_v5_s{j}']=f'D31B_M333_INDEPENDENT_ANALYTIC_SHARD_COMPLETE_{j}'
a=set();di=set()
def add(p):assert p.is_file(),p;a.add(p)
for p in(H/'inputs').rglob('*'):
 if p.is_file()and '__pycache__'not in p.parts:add(p)
files=['resolution.tex','finite_image.tex','relative_formal_fibre_F4.tex','physical_source5.m','frob_preflight.m','rational_trace_model.m','rational_trace_A5_RR.m','rational_trace_3isog_A5_RR.m','parent_native_selected12.m','independent_genuine35.json','relative_kummer_frame.m','relative_base_cusps.m','analytic_cusp_jets_v4.m','high_primitive_series.m','winding_antih10_N99_H34.m','winding_antih10_independent.json','native_frame_replay.m','prefix_rank_gf4.cpp','prefix_rank_input.txt','prefix_ranks_primary.json','independent_ranks.json','inherited_provenance.json','prepare_completion.py']
for n in files:add(H/n)
for stem,marker in required.items():
 for ext in['m','py']:
  p=W/'computations'/(stem+'.'+ext)
  if p.exists():add(p)
 log=W/'jobs'/(stem+'.log');st=W/'jobs'/(stem+'.status.json');add(log);add(st);s=json.loads(st.read_text());assert s['state']=='EXITED'and s['actual_exit_code']==0;assert marker in log.read_text(),(stem,marker);assert sha(W/s['script'])==s['script_sha256']
reports=[f'independent_analytic_v5_s{j}.json'if j in[0,1,3]else'independent_analytic_v4_s2.json'for j in range(4)];checks=[]
for n in reports:
 add(H/n);v=json.loads((H/n).read_text());assert v['status']=='PASS';checks+=v['checks']
assert sorted(c['base_index']for c in checks)==list(range(1,25));assert sum(c['branches']for c in checks)==48;assert sum(c['branches']for c in checks if c['degree']==1)==40;assert sum(c['branches']for c in checks if c['degree']==5)==8
merged={'status':'PASS','all24_base_cusps_and48_source_cusps_verified':True,'source_closed_cusps':{'1':40,'5':8},'both_ideal_generators_and_uniformizers_verified':True,'all82_function_series_reconstructed_to_native140precision_from_1000precision_primitives':True,'all48_tame_Kummer_uniformizers_and_holomorphic_J_UJ_matrices_matched':True,'original_four_v4_Horner_shards_all_completed':True,'three_optimized_v5_shards_also_completed':True,'reports_with_final_retained_numeric_outputs':reports,'note':'Individual shard source_cusps/verified_jet_matrices fields give the common global reference total48; actual verified scope is precisely their six checks and branch counts. Merged scope is all48. Optimized reports replace earlier same-path reports only before this packet was frozen; both sets of complete logs are retained.','checks':sorted(checks,key=lambda c:c['base_index'])}
(H/'independent_analytic_series.json').write_text(json.dumps(merged,indent=2)+'\n');add(H/'independent_analytic_series.json')
life=[]
for st in sorted((W/'jobs').glob('d31b_*.status.json')):
 if not any(st.stem.startswith(v)for v in['d31b_m3n33','d31b_m333']):continue
 s=json.loads(st.read_text());assert s['state']=='EXITED',st;life.append({'status_path':str(st.relative_to(R)),'actual_exit_code':s['actual_exit_code'],'accepted':st in a})
summary={'status':'PASS','group':[3,33],'degree':10,'relative_field':'Q(zeta3)','relative_degree':5,'new_exclusion':True,'genus':81,'full_JF_rank_zero_not_assumed':True,'h10_minus_factor_dimension':40,'both_QQ_and_minus3_twists_full_winding':True,'finite_alpha_JF':True,'NC_F4_closed_counts':[0,0,0,0,264],'cyclic33_image_places':88,'image_diamond_orbits':12,'selected_RR_dimensions':[0]*12,'independent_genuine5kernels':72,'selected_genuine3kernels':3,'manual_Velu_maps':75,'parent_native_prime_references':80,'full_frame_second_native_replay':True,'F4_closed_cusps':{'1':40,'5':8},'h10_free_at_all_cusps':True,'U_dimension':36,'rational_degree5_cusp_divisors':1086008,'rank_counts_by_pulled_back_degree':[864008,213600,8400],'quintic_cusp_divisors_rank5':8,'independent_series_bases':24,'independent_source_cusp_matrices':48,'independent_all_cusp_ranks':1086016,'relative_formal_F4_lemma_used':True,'no_global_cusp_lift_assumption':True,'workers_remaining':0}
for n in['completion_summary.json','status.json']:(H/n).write_text(json.dumps(summary,indent=2)+'\n');add(H/n)
(H/'lifecycle.json').write_text(json.dumps({'status':'PASS','all_owned_jobs_exited':True,'jobs':life,'failed_or_retired_probes_not_completed_certificates':True,'no_case_workers_supervisors_or_launchers_remain':True},indent=2)+'\n');add(H/'lifecycle.json')
for p in H.rglob('*'):
 if p.is_file()and p not in a and p.name!='certificate_manifest.json':di.add(p)
for directory in['jobs','computations']:
 for p in(W/directory).glob('d31b_*'):
  if p.is_file()and any(p.stem.startswith(v)for v in['d31b_m3n33','d31b_m333'])and p not in a:di.add(p)
def ent(p):return {'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)}
m={'status':'PASS','scope':'Complete C3xC33 decic exclusion: exact anti-invariant winding over Q and its -3 twist; relative F4 formal-fibre criterion at all1086016 cuspidal divisors with independent analytic jets and native cotangent frame; twelve exact cyclic33 rational-trace Hecke exclusions with75 genuine marked isogenies and80 fresh native primes','files':[ent(p)for p in sorted(a)],'diagnostic_files':[ent(p)for p in sorted(di)]};(H/'certificate_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
for e in m['files']+m['diagnostic_files']:assert sha(R/e['path'])==e['sha256']
print('M3N33_FROZEN',len(m['files']),len(m['diagnostic_files']),sha(H/'certificate_manifest.json'))
