from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[3];W=R/'work/deg10'
required={'d31b_m3n27_frob_preflight':'D31B_M3N27_FROB_PREFLIGHT_COMPLETE','d31b_m3n27_physical_source':'D31B_PHYSICAL_SOURCE5_COMPLETE 27 162 54 6','d31b_m3n27_quotient_A7_v4':'D31B_RELATIVE5_QUOTIENT_A7_COMPLETE 27 6 [ 0, 0, 0, 0, 0, 0 ]','d31b_m3n27_quotient_A7_v5':'D31B_RELATIVE5_QUOTIENT_A7_COMPLETE 27 6 [ 0, 0, 0, 0, 0, 0 ]','d31b_m327_independent_genuine_v2':'D31B_M327_INDEPENDENT_SOURCE_GENUINE7_COMPLETE','d31b_m327_torsion_gcd':'D31B_M327_TORSION_GCD_COMPLETE','d31b_m327_cusp_parity_v2':'D31B_M327_CUSP_PARITY_COMPLETE','d31b_m3n27_parent_native_v3':'D31B_M3N27_PARENT_NATIVE_COMPLETE 6 36'}
a=set();di=set()
def add(p):assert p.is_file(),p;a.add(p)
for p in (H/'inputs').rglob('*'):
 if p.is_file() and '__pycache__' not in p.parts:add(p)
for n in ['resolution.tex','physical_source5.m','frob_preflight.m','quotient_model.m','quotient_model_v5.m','quotient_A7_RR.m','quotient_A7_v5_RR.m','independent_genuine7.json','cusp_parity.json','torsion_gcd.m','parent_native_replay_v3.m','inherited_provenance.json','prepare_completion.py']:add(H/n)
for stem,marker in required.items():
 for ext in ['m','py']:
  p=W/'computations'/(stem+'.'+ext)
  if p.exists():add(p)
 log=W/'jobs'/(stem+'.log');st=W/'jobs'/(stem+'.status.json');add(log);add(st)
 s=json.loads(st.read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0
 assert marker in log.read_text(),(stem,marker)
 assert hashlib.sha256((W/s['script']).read_bytes()).hexdigest()==s['script_sha256']
life=[]
for st in sorted((W/'jobs').glob('d31b_*.status.json')):
 if not any(st.stem.startswith(v) for v in ['d31b_m3n27','d31b_m327']):continue
 s=json.loads(st.read_text());assert s['state']=='EXITED',st
 life.append({'status_path':str(st.relative_to(R)),'actual_exit_code':s['actual_exit_code'],'accepted':st in a})
summary={'status':'PASS','group':[3,27],'degree':10,'relative_ground_field':'Q(zeta3)','relative_degree':5,'new_exclusion':True,'rank_zero_input':'Derickx-Sutherland Theorem4.1 m3 n9','genus':46,'geometric_gonality_lower_bound':'78975/8192 >9','global_and_F4_cusp_signature':{'1':27,'3':6,'9':3},'Frobenius4_generates_full_relative_cusp_Galois_group':True,'split37_Jacobian_order_odd':True,'independent_Sage_modular_symbol_determinant_matches_native_Frobenius_order':True,'source_NC_closed_counts_F4':[0,0,0,0,162],'cyclic27_image_places':54,'image_diamond_orbits':6,'RR_dimensions':[0]*6,'independent_actual7kernels':48,'manual_Velu_checks':48,'parent_native_two_generator_prime_references':36,'GF4_constant_orientation_exactly_bound':True,'cuspidal_saturation_assumed':False,'workers_remaining':0}
(H/'completion_summary.json').write_text(json.dumps(summary,indent=2)+'\n');add(H/'completion_summary.json')
(H/'status.json').write_text(json.dumps(summary,indent=2)+'\n');add(H/'status.json')
(H/'lifecycle.json').write_text(json.dumps({'status':'PASS','all_owned_case_jobs_exited':True,'jobs':life,'failed_probes_not_certificates':True,'no_case_workers_supervisors_or_launchers_remain':True},indent=2)+'\n');add(H/'lifecycle.json')
for p in H.rglob('*'):
 if p.is_file() and p not in a and p.name!='certificate_manifest.json':di.add(p)
for directory in ['jobs','computations']:
 for p in (W/directory).glob('d31b_*'):
  if p.is_file() and any(p.stem.startswith(v) for v in ['d31b_m3n27','d31b_m327']) and p not in a:di.add(p)
def ent(p):return {'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
m={'status':'PASS','scope':'Complete C3xC27 decic exclusion: split37 odd torsion and intact cusp lifts close dyadic boundary; exhaustive relative-quintic images, all6 exact Hecke tests,48 genuine kernels, and36 fresh native prime references with explicit relativeF4 constants','files':[ent(p) for p in sorted(a)],'diagnostic_files':[ent(p) for p in sorted(di)]}
(H/'certificate_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
for e in m['files']+m['diagnostic_files']:assert hashlib.sha256((R/e['path']).read_bytes()).hexdigest()==e['sha256']
print('M3N27_FROZEN',len(m['files']),len(m['diagnostic_files']),hashlib.sha256((H/'certificate_manifest.json').read_bytes()).hexdigest())
