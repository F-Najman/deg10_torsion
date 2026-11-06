from pathlib import Path
import json,hashlib
N=Path(__file__).resolve().parent;P=N.parent;H=P.parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();files=set()
def add(p):assert p.is_file(),p;files.add(p)
accepted_names=['framed_prefix.snapshot.py','direct_cusp_membership_report.json','short_cusp_divisors.m','short_cusp_divisors.json','tau_cusp_class_relation.json','actual_target13_points.m','actual_target13_points.json','global_RR_section13900.snapshot.m','global_section13900_provenance.json','global_point13900_absolute.m','global_point13900_extraction.m','compact_realization.json','cuspidal_last_classes.tex','realization_and_larger_groups.tex','certify_subpacket.py','freeze_subpacket.py']
for n in accepted_names:add(N/n)
for n in ['kernel_structure.m','paired_p7_group_data.m','paired_p13_group_data.m','split_hecke_framed_arithmetic.txt','split_hecke_generator_points.m','relative_cubic_rank.m','relative_cubic_rank_p13.m','cusp_automorphisms.m','paired_final_survivors.json','paired_merged_p13_degree5_data.m','tau_cusp_permutation_p7.m','tau_cusp_permutation_p13.m','global_function_field_char0.m','global_cusp_model_p19.m','global_cusp_model_char0_rational.m','global_RR_char0_reductionfalse.m']:add(P/n)
for d in range(2,6):add(P/f'paired_p7_degree{d}_data.m')
for d in range(2,5):add(P/f'paired_p13_degree{d}_data.m')
for p in(P/'inputs').rglob('*'):
 if p.is_file():add(p)
accepted={
'd31b_m318_direct_cusp_membership':'D31B_M318_DIRECT_CUSP_MEMBERSHIP_COMPLETE',
'd31b_m318_short_cusp_divisors':'D31B_M318_SHORT_CUSP_DIVISORS_COMPLETE',
'd31b_m318_tau_cusp_class_relation':'D31B_M318_TAU_CUSP_CLASS_RELATION_COMPLETE',
'd31b_m318_global_point13900_v2':'D31B_M318_GLOBAL_POINT13900_COMPLETE',
'd31b_m318_cusp_automorphisms':'D31B_M318_CUSP_AUTOMORPHISMS_COMPLETE',
'd31b_m318_tau_cusps_p7':'D31B_M318_TAU_CUSP_PERMUTATION_COMPLETE',
'd31b_m318_tau_cusps_p13':'D31B_M318_TAU_CUSP_PERMUTATION_COMPLETE',
'd31b_m318_global_RR_reductionfalse':'D31B_M318_GLOBAL_REDUCTION_FALSE_RR_COMPLETE'}
diagnostic=['d31b_m318_orbit_gluing','d31b_m318_orbit_gluing_tau','d31b_m318_actual_cusp_gluing','d31b_m318_global_point13900','d31b_m318_effective_cusp_p19','d31b_m318_effective_cusp_p19_v2','d31b_m318_effective_cusp_p19_v3','d31b_m318_effective_cusp_p19_v4','d31b_m318_effective_cusp_p19_v5','d31b_m318_effective_cusp_p31','d31b_m318_cusp_model_p37','d31b_m318_effective_cusp_p37']
jobs=[]
for stem in list(accepted)+diagnostic:
 sp=H/'jobs'/(stem+'.status.json');st=json.loads(sp.read_text());assert st['state']=='EXITED',sp
 marker=accepted.get(stem)
 if marker:assert st['actual_exit_code']==0 and marker in(H/st['log']).read_text(),stem
 for p in[sp,H/st['log'],H/st['script']]:add(p)
 assert sha(H/st['script'])==st['script_sha256'],stem
 jobs.append({'status_path':str(sp.relative_to(H)),'actual_exit_code':st['actual_exit_code'],'accepted_completion_marker':marker,'scope':'completed proof component'if marker else'discovery, diagnostic or independently superseded experiment; not needed for the final proof'})
# Keep discovery dependencies identifiable, without treating them as exclusions.
for n in ['full_aut_actions.json','orbit_relations_report.json','actual_cusp_gluing_report.json','effective_cusp_divisors_p19.m','effective_cusp_divisors_p31.m','global_cusp_model_p37.m','effective_cusp_divisors_p37.m']:
 if(N/n).exists():add(N/n)
add(H/'run_job.py')
m={'status':'PASS','scope':'Owned arithmetic and exact number-field subcertificate; enclosing classification separately supplies exhaustive finite/relative filters and independent number-field verification','not_a_standalone_classification_certificate':True,'degree':10,'new_full_torsion_group':[3,18],'larger_groups_excluded_by_enclosing_classification':[[3,36],[6,18]],'conditional_final_class_reduction':{'four_normalized_pairs':4,'excluded_by_unmultiplied_cuspidal_membership':2,'remaining_global_degree5_classes':2,'tau_global_relation':'tau(A1) ~ g(A2), g=[2,0,2,5]','exact_realization':'global_point13900_absolute.m'},'files_sha256':{str(p.relative_to(H)):sha(p)for p in sorted(files)},'jobs':jobs,'upstream_partial_reports_are_historical':True,'independent_NF_and_full_torsion_report':'Integrated separately by the enclosing parent proof packet','diagnostics_are_not_certificates':True}
p=N/'certificate_manifest.json';assert not p.exists();p.write_text(json.dumps(m,indent=2)+'\n');print('M318_OWNED_SUBPACKET_FROZEN',len(files),sha(p))
