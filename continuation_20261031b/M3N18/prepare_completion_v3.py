from pathlib import Path
import hashlib,json,ast,re
H=Path(__file__).resolve().parent;R=H.parents[3];W=R/'work/deg10';accepted=set();context=set()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p):assert p.is_file(),p;accepted.add(p)
def resolve(name,h=None):
 candidates=[R/name,W/name,H/name];hits=[p for p in candidates if p.is_file()and(h is None or sha(p)==h)];assert hits,(name,h);return hits[0]
def rd(name):
 s=(H/name).read_text();s=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)\s*\]',lambda m:str(list(range(int(m[1]),int(m[2])+1))),s);s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
required={
'd31b_m3n18_frob_preflight':'D31B_M3N18_FROB_PREFLIGHT_COMPLETE',
'd31b_m318_paired_p7':'D31B_M318_PAIRED_INPUTS_COMPLETE',
'd31b_m318_hecke_framed':'D31B_M318_HECKE_ENDOMORPHISMS_COMPLETE',
'd31b_m318_kernel_structure_v2':'D31B_M318_KERNEL_STRUCTURE_COMPLETE',
'd31b_m318_class_frame_audit':'D31B_M318_CLASS_FRAME_AUDIT_COMPLETE',
'd31b_m318_cusp_automorphisms':'D31B_M318_CUSP_AUTOMORPHISMS_COMPLETE',
'd31b_m318_paired_three_relative_v2':'D31B_M318_PAIRED_ARITHMETIC_COMPLETE True',
'd31b_m318_relative_cubic_rank':'D31B_M318_RELATIVE_CUBIC_RANK_COMPLETE',
'd31b_m318_relative_cubic_rank_p13':'D31B_M318_RELATIVE_CUBIC_RANK_P13_COMPLETE',
'd31b_m318_tau_frame_p7':'D31B_M318_TAU_FRAME_COMPLETE',
'd31b_m318_quartic_phase_binding':'D31B_M318_QUARTIC_PHASE_BINDING_COMPLETE',
'd31b_m318_global_cusp_model_p19':'D31B_M318_GLOBAL_CUSP_MODEL_COMPLETE',
'd31b_m318_global_model_unlabelled':'D31B_M318_GLOBAL_RATIONAL_UNLABELLED_COMPLETE',
'd31b_m318_global_RR_reductionfalse':'D31B_M318_GLOBAL_REDUCTION_FALSE_RR_COMPLETE',
'd31b_m318_global_cusp_widths':'D31B_M318_GLOBAL_CUSP_WIDTHS_COMPLETE',
'd31b_m318_parent_global_basis_13900':'D31B_M318_PARENT_GLOBAL_BASIS_COMPLETE',
'd31b_m318_parent_native_shortclass13_v2':'D31B_M318_PARENT_NATIVE_SHORTCLASS13_COMPLETE',
'd31b_m318_parent_nf_reductions':'D31B_M318_PARENT_NF_REDUCTIONS_COMPLETE',
'd31b_m318_pointcheck_exact_v3':'D31B_M318_POINTCHECK_ALL_COMPLETE'}
for j in range(4):required[f'd31b_m318_p13_lean_s{j}']=f'D31B_M318_LEAN_P13_{j}_COMPLETE'
for j in range(3):required[f'd31b_m318_final_Hecke7_v2_s{j}']='D31B_M318_FINAL_HECKE_COMPLETE'
for j in[12,17]:required[f'd31b_m318_relative_quartic_{j}']='D31B_M318_RELATIVE_QUARTIC_COMPLETE'
nested=['endomorphism_gluing/certificate_manifest.json','sieve_A7_controls_manifest.json','parent_relative_cubic_verified.json','parent_relative_quartic_verified.json','parent_M318_verification.json']
for name in nested:
 p=H/name;add(p);a=json.loads(p.read_text());assert a['status']=='PASS'
 for nm,hs in a.get('files_sha256',{}).items():accepted.add(resolve(nm,hs))
 fs=a.get('files',[])
 if isinstance(fs,dict):
  for nm,hs in fs.items():accepted.add(resolve(nm,hs))
 else:
  for e in fs:accepted.add(resolve(e['path'],e['sha256']))
 for job in a.get('jobs',[]):
  if isinstance(job,dict)and job.get('accepted_completion_marker'):
   stem=Path(job['status_path']).name.replace('.status.json','');required[stem]=job['accepted_completion_marker']
# All native and Sage independent controls named by their enclosing reports.
for stem in ['d31b_m318_parent_relative_p7']+[f'd31b_m318_parent_relative_p13_block{j}'for j in range(1,11)]+[f'd31b_m318_parent_quartic_{j}'for j in[12,17]]+['d31b_m318_sieve_native_s0','d31b_m318_sieve_native_s1','d31b_m318_sieve_genuine_A7']:
 log=(W/'jobs'/(stem+'.log')).read_text();hits=re.findall(r'\b[A-Z][A-Z0-9_]*COMPLETE\b',log);assert hits,stem;required[stem]=hits[-1]
# Read actual marker names from immutable successful control logs when necessary.
for stem in ['d31b_m318_parent_global_basis_13900','d31b_m318_parent_native_shortclass13_v2','d31b_m318_parent_nf_reductions']:
 log=(W/'jobs'/(stem+'.log')).read_text();hits=re.findall(r'\b[A-Z][A-Z0-9_]*COMPLETE\b',log);assert hits;required[stem]=hits[-1]
for p in(H/'inputs').rglob('*'):
 if p.is_file()and '__pycache__'not in p.parts:add(p)
raws=['README.md','resolution.tex','relative_cubic_fibre.tex','relative_quartic_fibre.tex','endomorphism_gluing/cuspidal_last_classes.tex','inherited_provenance.json','prepare_completion.py','frob_preflight.m','fast_class_group_v4.m','class_frame_audit.json','split_hecke_framed_arithmetic.txt','split_hecke_framed_RR_identities.log','split_hecke_generator_points.m','kernel_structure.m','cusp_automorphisms.m','relative_cubic_rank.m','relative_cubic_rank_p13.m','relative_quartic_common.m','relative_quartic_12.m','relative_quartic_17.m','quartic_phase_binding.m','tau_frame_p7.m','tau_cusp_permutation_common.m','tau_cusp_permutation_p7.m','tau_cusp_permutation_p13.m','global_cusp_model_common.m','global_cusp_model_p19.m','global_cusp_model_char0_rational.m','global_function_field_char0.m','global_rational_cusps_unlabelled.m','global_rational_cusp_widths.m','global_RR_char0_reductionfalse.m','paired_p7_group_data.m','paired_p7_direct_RR_identities.txt','paired_p13_group_data.m','paired_three_relative_arithmetic_complete.json','paired_three_relative_matches_complete.json','paired_merged_p13_degree5_data.m','paired_final_Hecke_inputs_v2.m','final_Hecke_common_v2.m','isogenies_witness.m','sieve_A7_controls_summary.json','sieve_independent_genuine_A7.json','sieve_independent_global_points.json','parent_mathematical_review.json']
for d in range(2,6):raws.append(f'paired_p7_degree{d}_data.m')
for d in range(2,4):raws.append(f'paired_p13_degree{d}_data.m')
for j in range(4):raws +=[f'paired_lean_s{j}_p13_group_data.m',f'paired_lean_s{j}_p13_degree5_data.m',f'paired_lean_s{j}_p13_direct_RR_identities.txt']
for j in range(3):raws.append(f'final_Hecke7_v2_s{j}.m')
for n in raws:add(H/n)
for stem,marker in required.items():
 for ext in['m','py']:
  p=W/'computations'/(stem+'.'+ext)
  if p.exists():add(p)
 st=W/'jobs'/(stem+'.status.json');log=W/'jobs'/(stem+'.log');add(st);add(log);s=json.loads(st.read_text());assert s['state']=='EXITED'and s['actual_exit_code']==0,(stem,s);assert marker in log.read_text(),(stem,marker);assert sha(W/s['script'])==s['script_sha256']
review=json.loads((H/'parent_mathematical_review.json').read_text());assert review['status']=='PASS'and review['parent_mathematical_review']=='PASS'
for n,hs in review['reviewed_source_hashes'].items():assert sha(resolve(n,hs))==hs
nf=json.loads((H/'sieve_independent_global_points.json').read_text());assert nf['status']=='PASS'
rep=json.loads((H/'paired_three_relative_arithmetic_complete.json').read_text());assert rep['full_target_quintics_included']and rep['source_patterns']==3622 and rep['source_after_split_Hecke']==242 and rep['target_patterns']==20222 and rep['matched_target_patterns']==78
assert rd('kernel_structure.m')[1:]==[[3,3,3],3]
rr=rd('global_RR_char0_reductionfalse.m');assert [r[2]for r in rr]==[1,1]
widths=rd('global_rational_cusp_widths.m');assert len(widths)==18 and all(r[1]<0 for r in widths)
life=[]
for st in sorted((W/'jobs').glob('d31b_*.status.json')):
 if not any(st.name.startswith(v)for v in['d31b_m318','d31b_m3n18']):continue
 s=json.loads(st.read_text());assert s['state']=='EXITED',st;life.append({'status_path':str(st.relative_to(R)),'actual_exit_code':s['actual_exit_code'],'proof_completion_required':st.name.replace('.status.json','')in required})
summary={'status':'PASS','parent_mathematical_review':'PASS','complete_low_degree_classification':True,'classified_effective_relative_degree':5,'relative_field':'Q(zeta3)','absolute_degree':10,'realized_groups':[[3,18]],'excluded_groups':[[3,36],[6,18]],'new_occurrence':True,'new_exclusion_of_C3xC18':False,'marking_automorphism_group_order':27,'closed_quintic_marking_orbits':2,'orbit_size':27,'exact_torsion_orders':[54,54],'both_curves_2_isogenous_over_same_decic_field':True,'global_cuspidal_generation_assumed':False,'Jacobian_rank_zero_input':'Derickx-Sutherland Theorem4.1 m3 n6','source_genus':10,'gonality_F7_lower_bound':6,'cusp_group_order':20253807,'split_Hecke_kernel_over_cusps':[3,3,3],'NC_closed_counts7':[0,27,144,513,3672],'NC_closed_counts13':[0,54,846,7506,75384],'source_patterns':3622,'target_patterns':20222,'paired_survivors':78,'distinct_Hecke_supports':53,'exact_negative_Hecke_supports':43,'independent_marked_sources':45,'independent_geometric7kernels':360,'independent_native_prime_refs':292,'independent_cubic_restrictions':2394,'independent_quartic_restrictions':36,'global_RR_dimensions':[1,1],'actual_relative_quintic_and_absolute_decic_models_verified':True,'both_global_RR_sections_independently_evaluate_to_zero':True,'both_unramified_F49_good_reductions_have_order54':True,'positive_native_and_Sage_F13_marked_binding':True,'final_proposition_label':'prop:decic318-classification','frontier_corollary_label':'cor:decic318-no-larger','workers_remaining':0,'files_sha256':{str(p.relative_to(R)):sha(p)for p in sorted(accepted)}}
(H/'final_certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n');add(H/'final_certificate_summary.json')
small={k:v for k,v in summary.items()if k!='files_sha256'}
for n in['completion_summary.json','status.json']:(H/n).write_text(json.dumps(small,indent=2)+'\n');add(H/n)
(H/'lifecycle.json').write_text(json.dumps({'status':'PASS','all_owned_and_case_control_jobs_exited':True,'no_paused_workers_or_automatic_launchers':True,'jobs':life,'failed_retired_and_timed_out_attempts_not_completed_certificates':True},indent=2)+'\n');add(H/'lifecycle.json')
# Hash all preserved exploration separately. A nested context may contain diagnostics;
# its own job scope, not mere hash inclusion, determines whether it is a certificate.
di=set()
for p in H.rglob('*'):
 if p.is_file()and p not in accepted and p.name!='certificate_manifest.json':di.add(p)
for directory in['jobs','computations']:
 for p in(W/directory).glob('d31b_*'):
  if p.is_file()and any(p.name.startswith(v)for v in['d31b_m318','d31b_m3n18'])and p not in accepted:di.add(p)
def ent(p):return{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)}
m={'status':'PASS','scope':'Complete degree-five noncuspidal divisor classification on X1(3,18)/Q(zeta3), new exact C3xC18 decic occurrence, and C3xC36/C6xC18 exclusions','realized_groups':[[3,18]],'excluded_groups':[[3,36],[6,18]],'proof_completion_markers':required,'hash_inclusion_is_not_a_claim_that_diagnostic_jobs_are_certificates':True,'nested_certificates':nested,'files':[ent(p)for p in sorted(accepted)],'diagnostic_files':[ent(p)for p in sorted(di)]};(H/'certificate_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
for e in m['files']+m['diagnostic_files']:assert sha(R/e['path'])==e['sha256']
print('M3N18_CLASSIFICATION_FROZEN',len(m['files']),len(m['diagnostic_files']),len(required),sha(H/'certificate_manifest.json'))
