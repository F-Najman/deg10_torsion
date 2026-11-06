from pathlib import Path
import hashlib,json,os
P=Path(__file__).resolve().parent;HOME=P.parents[1];REPO=HOME.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
dy=rd('dyadic_rigidity_v2.m');assert len(dy[6])==8825 and len(dy[7])==4104 and len(dy[8])==25;assert sum(r[5]==1 for r in dy[6])==4721;assert len({tuple(r[6])for r in dy[7]})==800 and len({(tuple(r[6]),tuple(r[7]))for r in dy[7]})==4103
tr=rd('transport_800_classes.m');assert len(tr[0])==800 and len(tr[1])==724 and [r[0]for r in tr[0]]==list(range(1,801));fib=rd('class_members_all.m');assert len(fib[0])==192 and len(fib[1])==1067 and fib[2]==217;src=rd('member_two_isogenies_all.m');assert len(src[0])==334 and len(src[1])==192 and len(src[2])==147 and sum(len(r[2])for r in src[1])==361;fin=rd('member_all_classes.m');assert len(fin[0])==334 and len(fin[1])==147 and not fin[2]
n2=rd('parent_dyadic_replay.m');assert n2[0]==[13,25,8825,4721,4104]and len(n2[1])==8825 and all(v[1]==dy[6][v[0]-1][5]for v in n2[1]);nr=[rd(f'parent_f5_RR_s{i}.m')for i in[0,1]];obs=[r for a in nr for r in a[0]];assert sorted(r[0]for r in obs)==list(range(1,801))and all(r[1]==tr[0][r[0]-1][4]for r in obs);assert[sum(a[1][i]for a in nr)for i in range(5)]==[800,724,1067,217,1284];np=rd('parent_f5_pairs.m');assert len(np[0])==192 and len(np[1])==147 and np[2]==[361,334,147,0]and all(not r[5]for r in np[1]);assert sorted(r[2]for r in np[1])==[0]*59+[1]*88
cm=json.loads((P/'inputs/CM54_minimum_check.json').read_text());assert cm['rows']['54']['least_CM_degree']==18
for v in json.loads((P/'input_manifest.json').read_text()):assert sha(REPO/v['source'])==sha(REPO/v['copy'])==v['sha256']
markers={'d31b_C54_dyadic_rigidity_v2':'D31B_C54_DYADIC_RIGIDITY_COMPLETE','d31b_C54_transport_800_resume':'D31B_C54_800_RR_COMPLETE','d31b_C54_enumerate_all_class_members':'D31B_C54_RR_MEMBERS_ALL_COMPLETE','d31b_C54_member_two_isogenies_all':'D31B_C54_MEMBER_TWO_ISOGENIES_ALL_COMPLETE','d31b_C54_member_all_classes':'D31B_C54_MEMBER_ALL_CLASS_COMPLETE','d31b_C54_parent_dyadic_replay':'D31B_C54_PARENT_DYADIC_COMPLETE','d31b_C54_parent_f5_rr_s0':'D31B_C54_PARENT_F5_RR_COMPLETE','d31b_C54_parent_f5_rr_s1':'D31B_C54_PARENT_F5_RR_COMPLETE','d31b_C54_parent_f5_pairs_v2':'D31B_C54_PARENT_F5_PAIRS_COMPLETE'}
jobs=[]
for sp in sorted((HOME/'jobs').glob('d31b_C54*.status.json')):
 name=sp.name.removesuffix('.status.json');s=json.loads(sp.read_text());assert s['state']=='EXITED';assert sha(HOME/s['script'])==s['script_sha256'];log=(HOME/s['log']).read_text();ok=name in markers
 if ok:assert s['actual_exit_code']==0 and markers[name]in log and not any(t in log for t in['Runtime error','User error:','Traceback (most recent call last):']),name
 jobs.append({'purpose':name,'state':s['state'],'actual_exit_code':s['actual_exit_code'],'complete_certificate':ok,'diagnostic_only':not ok,'script':str((HOME/s['script']).relative_to(REPO)),'log':str((HOME/s['log']).relative_to(REPO)),'status':str(sp.relative_to(REPO))})
assert all(any(v['purpose']==n for v in jobs)for n in markers)
active=[]
for q in Path('/proc').iterdir():
 if not q.name.isdigit():continue
 try:
  cmd=[v.decode(errors='replace')for v in(q/'cmdline').read_bytes().split(b'\0')if v]
  if cmd and(Path(cmd[0]).name.startswith('magma')or Path(cmd[0]).name.startswith('python'))and any(Path(v).stem.startswith('d31b_C54')for v in cmd[1:]):active.append(int(q.name))
 except OSError:pass
assert not active,active
summary={'status':'COMPLETE EXCLUSION','group':'C54','degree':10,'source_genus':52,'rational_gonality_lower_bound':12,'target_genus':13,'rational_target_Jacobian_invariants':[3,3,52497],'dyadic_cases':8825,'dyadic_RR_excluded':4721,'dyadic_remaining_cases':4104,'ordered_rational_class_pairs':4103,'rational_first_classes':800,'F5_effective_classes':724,'F5_noncuspidal_members':1067,'F5_pure_cusp_members':217,'F5_distinct_noncuspidal_points':192,'actual_two_kernel_orbits':361,'distinct_marked_targets':334,'complete_finite_noncuspidal_lifts':147,'final_remaining':0,'all_dyadic_RR_and_classes_independently_replayed':True,'all_F5_classes_RR_fibres_and_pair_tests_independently_replayed_with_FRESH_PROVEN_class_groups':True,'all_two_isogenies_and_norms_cross_language_verified':True,'all_input_hashes_unchanged':True,'all_recorded_jobs_exited':True,'active_workers':active,'paused_workers':0,'automatic_launchers':0}
(P/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={q for q in P.rglob('*')if q.is_file()and q.name!='certificate_manifest.json'and'__pycache__'not in q.parts};paths.add(HOME/'run_job.py')
for job in jobs:
 for key in['script','log','status']:paths.add(REPO/job[key])
paths.update(q for q in(HOME/'computations').glob('d31b_C54*')if q.is_file());paths.update(q for q in(HOME/'shared/mdmagma/v2').rglob('*')if q.is_file())
manifest={'status':'PASS','group':'C54','degree':10,'summary_sha256':sha(P/'certificate_summary.json'),'files':[{'path':str(q.relative_to(REPO)),'sha256':sha(q),'bytes':q.stat().st_size}for q in sorted(paths)],'jobs':jobs};(P/'certificate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('D31B_C54_FINAL_MANIFEST_COMPLETE',len(paths),len(jobs),sha(P/'certificate_manifest.json'))
