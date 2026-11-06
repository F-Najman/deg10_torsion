#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,ast,subprocess
R=Path.cwd();D=R/'work/deg10';H=D/'next_cases_20261030/C69';J=D/'jobs'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_evidence(d):
 for e in d['evidence']:assert sha(D/e['path'])==e['sha256'],e['path']
small=json.loads((H/'small_cover_verified_partial.json').read_text());parent=json.loads((H/'parent_final_checks.json').read_text());check_evidence(small);check_evidence(parent)
assert parent['status'].startswith('PASS') and parent['native_RR_dimension']==0
required={'direct69_select_input_v2':'DIRECT69_SELECTED_INPUT_COMPLETE','direct69_native_hecke_control':'DIRECT69_NATIVE_HECKE_CONTROL_COMPLETE','direct69_verify_input_v2':'DIRECT69_SELECTED_INPUT_REPLAY_COMPLETE','q67_local_preflight':'Q67_LOCAL_PREFLIGHT_COMPLETE','q67_profile_rr_v2':'Q67_V2_A5_RR_COMPLETE base_case=2 principal=false','q67_geometry_audit':'Q67_GEOMETRY_PACKET_AUDIT_COMPLETE','cont69_q67_parent_native_v2':'C69_Q67_PARENT_NATIVE_RR_V2_COMPLETE 4 0'}
checks=[];extras=[]
for stem,marker in required.items():
 s=json.loads((J/(stem+'.status.json')).read_text());log=(J/(stem+'.log')).read_text();assert s['state']=='EXITED' and s['actual_exit_code']==0 and s['termination_reason'] is None and marker in log
 assert 'Runtime error'not in log and 'Traceback'not in log and sha(D/s['script'])==s['script_sha256'];checks.append({'job':stem,'marker':marker,'actual_exit_code':0});extras.extend([J/(stem+'.status.json'),J/(stem+'.log'),D/s['script']])
failed=[]
for stem,code in [('direct69_select_input',1),('direct69_verify_input',1),('direct69_profile_and_rr',-14),('q67_profile_rr',-14)]:
 s=json.loads((J/(stem+'.status.json')).read_text());assert s['state']=='EXITED' and s['actual_exit_code']==code;assert sha(D/s['script'])==s['script_sha256'];failed.append({'job':stem,'actual_exit_code':code,'used_as_proof':False});extras.extend([J/(stem+'.status.json'),J/(stem+'.log'),D/s['script']])
rr=json.loads((H/'quotient67/rr_result_v2.json').read_text());assert rr['completed'] and rr['principal']is False and rr['operator']=='T5-1-5<5>'
native=ast.literal_eval((H/'quotient67/parent_native_replay.m').read_text().replace('[*','[').replace('*]',']'));assert native[1:4]==[60,0,0]
geo=json.loads((H/'quotient67/geometry_verified.json').read_text());assert geo['status']=='PASS' and geo['quotient_genus']==67 and geo['packet_prime_generators_per_place']==2
inp=json.loads((H/'direct_full/selected_input_verified.json').read_text());assert inp['status']=='PASS' and inp['point_order']==69 and inp['curve_order']==966
snap=json.loads((H/'proof_dependencies/snapshot_manifest.json').read_text())
for p in snap['snapshots']:assert sha(R/p['source'])==p['sha256']==sha(R/p['snapshot'])
prose=(H/'resolution.tex').read_text();assert 'cite[Theorem~3.1]{DEHMZ}'in prose and 'gon}_{\\overline{\\mathbb Q}}'in prose and 'exact winding computation' not in prose
active=[]
for line in subprocess.check_output(['ps','-eo','uid,pid,args'],text=True).splitlines()[1:]:
 a=line.split(None,2)
 if len(a)<3:continue
 uid,pid,cmd=a;v=cmd.split();name=Path(v[0]).name if v else '';scope='/C69/'in cmd or any(x in cmd for x in ['cont69_q67','direct69_','q67_profile'])
 worker=name.startswith('magma')or(name.startswith('python')and len(v)>1 and (v[1].endswith('run_job.py')or Path(v[1]).name.startswith(('direct69_','q67_','cont69_q67'))))
 if scope and worker:active.append({'owner_uid':uid,'pid':pid,'command':cmd})
assert not active,active
life={'all_case_workers_and_supervisors_exited':True,'active_workers_supervisors_launchers':active,'paused_jobs':[],'automatic_launchers_created':False,'successful_new_checks':checks,'unsuccessful_attempts_preserved_not_used_as_proof':failed}
(H/'lifecycle_final.json').write_text(json.dumps(life,indent=2)+'\n')
files=[H/'resolution.tex',H/'certify69.py',H/'lifecycle_final.json',H/'parent_final_checks.json',H/'small_cover_verified_partial.json',H/'degree10_base.m',H/'base_j_functions.m']
for d in ['direct_full','quotient67','proof_dependencies']:files.extend(p for p in(H/d).rglob('*')if p.is_file()and'__pycache__'not in p.parts)
for d in[small,parent]:files.extend(D/e['path']for e in d['evidence'])
files.extend(extras);files=sorted(set(files));hashes={str(p.relative_to(R)):sha(p)for p in files}
summary={'status':'COMPLETE; parent mathematical and independent geometry audits PASS','degree':10,'excluded_group':'C69','claim':'No elliptic curve over a degree10 number field has a point of order69','proof':str((H/'resolution.tex').relative_to(R)),'rank_zero_input':'DEHMZ Theorem3.1','source_genus':133,'geometric_gonality_bound':'10725/512 > 20','noncuspidal_closed_point_counts':small['noncuspidal_closed_points_through_degree10'],'small_cover_completed_intermediate_stage':{'status_at_that_stage':small['status'],'base_tests':9,'cover_tests':18,'residual_then':'one degree10 diamond orbit','resolved_by_current_genus67_test':True,'source_sha256':sha(H/'small_cover_verified_partial.json')},'final_quotient_genus':67,'quartic_over_X1_23_genus12':'v^4+v^3+1/j','exact_constants_proof':'simple root v=1 at a rational base cusp; Hensel','actual_source_and_Hecke_lifts':'unique same-degree root, factors[1,3], all exact prime restrictions independently replayed','correct_operator':'T5-1-5<5>','Hecke_degrees':[10,50],'Hecke_observations':[1,5],'Hecke_norm_weights':[1,1],'positive_Hecke_degree':60,'primary_principal':False,'independent_native_RR_dimension':0,'source_and_norm_inputs_pinned_locally':True,'terminated_attempts_not_used_as_proof':True,'lifecycle':life,'file_hashes':hashes}
(H/'final_certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('DECIC69_FINAL_CERTIFICATE_COMPLETE',json.dumps({'group':'C69','parent_audits':'PASS','fresh_native_RR_dimension':0,'hashed_files':len(hashes),'active_jobs':0},sort_keys=True))
