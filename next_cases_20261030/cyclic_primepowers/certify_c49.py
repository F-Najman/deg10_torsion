#!/usr/bin/env python3
"""Close the C49 degree-ten proof with explicit successful and partial provenance."""
from pathlib import Path
import json,hashlib,subprocess,os,re
H=Path(__file__).resolve().parent; root=Path.cwd();jobs=root/'work/deg10/jobs'
required={
 'pp49_finite_sieve':'DECIC_C49_GOOD_BRANCHES_COMPLETED',
 'pp49_finite_controls':'DECIC_FINITE_AND_HECKE_CONTROLS_COMPLETE',
 'pp49_replay_old_jets_degree10':'DECIC_C49_OLD_JETS_DEGREE10_REPLAY_COMPLETE',
 'pp49_native_matrix_replay':'DECIC_C49_CARTIER_TRANSPORT_REPLAY_COMPLETE'}
completed=[]
for stem,marker in required.items():
 s=json.loads((jobs/(stem+'.status.json')).read_text());log=(jobs/(stem+'.log')).read_text()
 assert s['state']=='EXITED' and s['actual_exit_code']==0 and s['termination_reason'] is None
 assert marker in log and 'Runtime error' not in log and 'Traceback' not in log
 script=root/'work/deg10'/s['script'];assert hashlib.sha256(script.read_bytes()).hexdigest()==s['script_sha256']
 completed.append({'stem':stem,'marker':marker,'status':s})
r1=json.loads((H/'c49_old_jets_degree10_verified.json').read_text());r2=json.loads((H/'c49_cartier_transport_verified.json').read_text())
assert r1['status']==r2['status']=='PASS' and r1['rank_all']==r2['rank_all']==10
assert r1['tested_divisors']==2203209 and r2['seven_distinct_rational_support_cases']==697680 and r2['sevenfold_fixed_rational_cusp_cases']==6
assert r2['Cartier_minimal_polynomial_degree']==69 and r2['known_pair_Cartier_span_rank']==69
assert r2['basis_alignment_single_evaluator_Cartier_cyclic_rank']==69 and r2['diamond_order']==21
assert r2['diamond_checks_all_old_rational_evaluators']==18 and r2['cubic_covariance_orders_checked']==[1,2,3]
partial_specs={
 'pp49_fast_cartier_export':('DECIC_C49_FAST_CARTIER_JET_DONE i=1 ', 'only the completed order-seven first-rational-cusp block and basis alignment data; the full export is not a completed certificate'),
 'pp49_completion_control':('DECIC_C49_COMPLETION_RATIONAL_ALL_SEVEN_KERNELS_AGREE', 'only the completed seven rational restriction-kernel comparisons with native Completion; the cubic and combined checks did not complete'),
 'pp49_cusp_permutation':(None,'no output used')}
partials=[]
for stem,(marker,scope)in partial_specs.items():
 s=json.loads((jobs/(stem+'.status.json')).read_text());log=(jobs/(stem+'.log')).read_text()
 assert s['state']=='EXITED' and s['actual_exit_code']==-15 and s['termination_reason']
 if marker:assert marker in log
 script=root/'work/deg10'/s['script'];assert hashlib.sha256(script.read_bytes()).hexdigest()==s['script_sha256']
 partials.append({'stem':stem,'whole_job_is_complete_certificate':False,'accepted_completed_subcalculation':scope,'marker':marker,'status':s})
old=root/'work/degree9/prime_orders/continued_20261026/P49'
oldsrc=(old/'cartier_transport_export.m').read_text();oldlog=(old/'cartier_transport_export.log').read_text();oldexit=json.loads((old/'cartier_transport_export.exit.json').read_text())
assert oldexit['native_exit_code']==1 and 'Runtime error' in oldlog and 'Incompatible sequences' in oldlog
formula='MD:=Matrix(GF(2),[Eltseq((dg(w/eta)*deta)@@om):w in B])'
assert formula in oldsrc
assert oldsrc.index(formula)<oldsrc.index('MD^21 eq I')<oldsrc.index('Eltseq(MD)')<oldsrc.index('blocks:=')
assert 'NONIC49_CARTIER_TRANSPORT_SETUP genus=69 Udim=68 diamond2_exact_order=21' in oldlog
native_origin={'source_job_is_complete_certificate':False,'source_native_exit_code':oldexit['native_exit_code'],'later_failure':'heterogeneous matrix-sequence concatenation after matrix export','accepted_completed_subcalculation':'all69 actual differential pullbacks via dg(w/du)*d(dg(u)), exact coordinates via the inverse native differential map, followed by order21/Cartier/U checks and export','basis_invariance_verified':'the exported diamond equals an explicit polynomial in the common cyclic Cartier matrix','actual_source_formula':formula}
active=[]
for line in subprocess.check_output(['ps','-eo','uid,pid,args'],text=True).splitlines()[1:]:
 bits=line.split(None,2)
 if len(bits)<3:continue
 uid,pid,args=bits;argv=args.split();first=Path(argv[0]).name if argv else ''
 own=str(H)in args or str(H.relative_to(root))in args
 if own and (first.startswith('magma')or(first.startswith('python')and len(argv)>1 and (argv[1].endswith('run_job.py')or'/pp49_'in argv[1]))):
  active.append({'uid':int(uid),'pid':int(pid),'command':args})
assert not active,active
lifecycle={'all_new_C49_jobs_exited':True,'new_job_count':7,'successful_jobs':list(required),'terminated_jobs':list(partial_specs),'active_workers_supervisors_launchers':active,'paused_workers':[],'automatic_launchers_created':False,'no_new_C49_job_is_to_be_restarted_without_a_new_purpose':True}
(H/'lifecycle_final.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
extras=[old/x for x in ['cartier_mixed9_jets.m','cartier_mixed9_export.m','cartier_mixed9_export.log','cartier_mixed9_export.exit.json','cartier_mixed9_export.status.json','cartier_transport_jets.m','cartier_transport_export.m','cartier_transport_export.log','cartier_transport_export.exit.json','mixed9_RR_control.m','mixed9_RR_control.log','refs/DEHMZ_source.tex']]
extras.append(root/'work/deg10/next_cases_20261030/dyadic_formal_kernel.tex')
for stem in list(required)+list(partial_specs):extras +=[jobs/(stem+'.status.json'),jobs/(stem+'.log')]
owned=[p for p in H.rglob('*')if p.is_file() and '__pycache__' not in p.parts and p.name!='final_certificate_summary.json']
files=sorted(set(owned+extras));hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files}
summary={'degree':10,'groups_excluded':['C49'],'claim':'No elliptic curve over a degree10 number field has a point of order49','status':'COMPLETE; parent mathematical audit PASS','proof':str((H/'resolution.tex').relative_to(root)),'proof_dependencies':['work/deg10/next_cases_20261030/dyadic_formal_kernel.tex'],'noncuspidal_reduction_counts':{'degree9':{'points':42,'diamond_orbits':2},'degree10':{'points':126,'diamond_orbits':6}},'Cartier_ranks':{'old_selected18_checks':2203209,'seven_distinct_support_checks':697680,'sevenfold_base_cusp_checks':6,'total_checks':2900895,'all_rank':10},'successful_jobs':completed,'terminated_jobs_with_explicit_scopes':partials,'completed_old_native_matrix_subcalculation':native_origin,'lifecycle':lifecycle,'file_hashes':hashes,'C81':'NOT EXCLUDED; rank0 and rational2torsion0 known, but the previous Hecke argument needs gonality>30 and the verified analytic bound is only236925/8192'}
(H/'final_certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('DECIC_C49_FINAL_CERTIFICATE_COMPLETE',json.dumps({'excluded':'C49','total_rank_checks':2900895,'hashed_files':len(hashes),'new_jobs_exited':7,'active_jobs':0},sort_keys=True))
