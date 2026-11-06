#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess
H=Path(__file__).resolve().parent;root=Path.cwd();jobs=root/'work/deg10/jobs'
subprocess.run(['python3',str(H/'further51_cusp_orbits.py')],check=True,stdout=subprocess.PIPE,text=True)
required={
'further51_finite_sieve':'DECIC51_GOOD_BRANCHES_COMPLETED',
'further51_finite_controls':'DECIC51_FINITE_AND_HECKE_CONTROLS_COMPLETE',
'further51_replay_old_jets':'DECIC51_DEGREE10_CUSP_REPLAY_DONE tested=9655880 exceptions=0 complete=true',
'further51_new_jets':'DECIC51_NEW_JETS_COMPLETE',
'further51_replay_new_jets':'DECIC51_NEW_JETS_REPLAY_COMPLETE'}
records=[]
for stem,marker in required.items():
 status=json.loads((jobs/(stem+'.status.json')).read_text());log=(jobs/(stem+'.log')).read_text()
 assert status['state']=='EXITED' and status['actual_exit_code']==0 and status['termination_reason'] is None
 assert marker in log and 'Runtime error' not in log and 'Traceback' not in log
 script=root/'work/deg10'/status['script'];assert hashlib.sha256(script.read_bytes()).hexdigest()==status['script_sha256']
 records.append({'stem':stem,'marker':marker,'status':status})
flog=(jobs/'further51_finite_sieve.log').read_text();assert flog.count('principal=false')==9 and 'summary=[ <9, 80, 5, 0>, <10, 64, 4, 0> ]'in flog
clog=(jobs/'further51_finite_controls.log').read_text();assert 'degree=9 all_closed_points=80 source_reps=5'in clog and 'degree=10 all_closed_points=64 source_reps=4'in clog and clog.count('native_maps=6')==9
old=json.loads((H/'old_jets_degree10_replay.json').read_text());new=json.loads((H/'new_jets_verified.json').read_text());orbits=json.loads((H/'cusp_orbits_verified.json').read_text())
assert old['exceptions_count']==0 and old['tested_divisors']==9655880 and old['complete_enumeration']
assert new['status']=='PASS'and new['rank_all']==10 and new['basis_transition_dimension']==1 and new['basis_transition_invertible']
assert new['all_cuspidal_divisors']==9656144 and new['new_high_multiplicity_divisors']==264 and new['nine_plus_one_ranks']==[10]*15
assert orbits['status']=='PASS'and orbits['global_rational_cusps']==16 and orbits['rational_cusps_transitive']and orbits['quadratic_cusps_transitive']
active=[]
for line in subprocess.check_output(['ps','-eo','uid,pid,args'],text=True).splitlines()[1:]:
 pieces=line.split(None,2)
 if len(pieces)<3:continue
 uid,pid,args=pieces;argv=args.split();first=Path(argv[0]).name if argv else ''
 own=str(H)in args or str(H.relative_to(root))in args
 if own and (first.startswith('magma')or(first.startswith('python')and len(argv)>1 and (argv[1].endswith('run_job.py')or'/further51_'in argv[1]))):active.append({'uid':int(uid),'pid':int(pid),'command':args})
assert not active,active
life={'new_jobs':5,'all_jobs_exited':True,'all_actual_exit_codes':0,'active_workers_supervisors_launchers':active,'paused_jobs':[],'terminated_jobs':[],'automatic_launchers_created':False}
(H/'lifecycle_final.json').write_text(json.dumps(life,indent=2)+'\n')
base=root/'work/degree9/cyclic/continuation_20261026/N51_geometry/cyclic_sieves'
extra=[base/p for p in ['N51_frozen_jets.m','generic_cartier_completion_export.m','N51_frozen.log','N51_projection_replay.log','N51_full_RR_crosscheck_root.log','N51_prefix27_RR_control_root.m','N51_prefix27_RR_control_root_data.m','N51_prefix27_RR_control_root.log']]
extra+=[root/'work/deg10/next_cases_20261030/dyadic_formal_kernel.tex',root/'work/degree9/cyclic/continuation_20261026/exclusion51.tex',root/'work/degree9/prime_orders/continued_20261026/P49/refs/DEHMZ_source.tex',root/'work/degree9/shared/fast_hecke.m']
for stem in required:extra+=[jobs/(stem+'.status.json'),jobs/(stem+'.log')]
owned=[p for p in H.rglob('*')if p.is_file() and '__pycache__' not in p.parts and p.name!='final_certificate_summary.json']
files=sorted(set(owned+extra));hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files}
summary={'degree':10,'group_excluded':'C51','claim':'No elliptic curve over a degree10 number field has a point of order51','status':'COMPLETE; parent mathematical review pending','proof':str((H/'resolution.tex').relative_to(root)),'proof_dependencies':['work/deg10/next_cases_20261030/dyadic_formal_kernel.tex'],'rank_zero_source':'DEHMZ Theorem3.1(2)','noncuspidal_closed_points':{'degrees1_to8':0,'degree9':80,'degree10':64},'noncuspidal_diamond_orbits':{'degree9':5,'degree10':4},'A5_nonprincipal_tests':9,'cusp_signature_at2':{'1':16,'2':8,'8':4},'Cartier_image_dimension':62,'old_individual_caps_new_total_degree_replayed':9655880,'new_high_multiplicity_divisors':264,'new_native_representative_cases':17,'all_cuspidal_degree10_divisors':9656144,'all_restriction_ranks':10,'basis_alignment':{'Cartier_intertwiner_and_degree8_anchor_dimension':1,'invertible':True,'actual_geometric_transition_identified_up_to_diamond':True},'earlier_degree_exclusions_not_transferred':True,'earlier_degree_independent_jet_data_used_at_fresh_degree10':True,'successful_jobs':records,'lifecycle':life,'file_hashes':hashes}
(H/'final_certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('DECIC51_FINAL_CERTIFICATE_COMPLETE',json.dumps({'excluded':'C51','all_cuspidal_divisors':9656144,'new_native_cases':17,'A5_nonprincipal_tests':9,'hashed_files':len(hashes),'jobs_exited':5,'active_jobs':0},sort_keys=True))
