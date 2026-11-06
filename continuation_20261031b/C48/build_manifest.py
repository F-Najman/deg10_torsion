from pathlib import Path
import hashlib,json,os
P=Path(__file__).resolve().parent;HOME=P.parents[1];REPO=HOME.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
hc=rd('high_character_merged.json');assert len(hc[1])==149508 and sum(r[5]for r in hc[1])==12604 and sum(not r[4]for r in hc[1])==1 and len(hc[2])==5619;hl=rd('high_class_merged.json');assert len(hl[0])==5619 and len(hl[1])==107
lo=rd('low_geometry.m');assert len(lo[1])==20404 and len(lo[2])==3550;ls=json.loads((P/'low_patterns_summary.json').read_text());assert ls['status']=='PASS'and ls['support_cases']==119986 and ls['remaining_count']==243;ded=json.loads((P/'low_source_orbits.json').read_text());assert ded['status']=='PASS'and ded['original_cases']==243 and ded['source_diamond_orbits']==201
nclass=0;nwit=0;ndirect=0
for kind,n in[('high',4),('low',3)]:
 for sh in range(n):
  a=rd(f'{kind}_class_native_s{sh}.m');assert len(a[0])==a[1][0]and a[1][0]==a[1][1]+a[1][2];nclass+=a[1][0];nwit+=a[1][1];ndirect+=a[1][2]
assert(nclass,nwit,ndirect)==(9169,9078,91)
nhigh=0;nlower=0;nrefs=0;ncontexts=0
for q,wanted in[(7,57),(13,50)]:
 native=rd(f'parent_native_high_q{q}.m');primary=[r for sh in[0,1]for r in rd(f'hecke{q}_s{sh}.m')[0]];assert len(native)==len(primary)==wanted and sorted(r[0]for r in native)==sorted(r[0]for r in primary)and all(r[8]==0 for r in native+primary);nhigh+=wanted;nrefs+=sum(2+len(r[7])for r in native);ncontexts+=wanted
for q,wanted,pts in[(7,50,73),(13,151,156)]:
 native=rd(f'parent_native_low_q{q}.m');primary=rd(f'low_hecke{q}_merged.m');assert len(native[0])==pts and len(native[1])==wanted and native[1]==primary[3]and all(r[4]==0 for r in native[1]);assert sorted(r[0]for r in native[0])==sorted(r[0]for r in primary[2]);nlower+=wanted;nrefs+=sum(2+len(r[7])for r in native[0]);ncontexts+=pts
assert(nhigh,nlower,nrefs,ncontexts)==(107,201,1681,336)
genuine_count=0;kernel_count=0
for kind,q,n in[('high',7,1),('high',13,2),('low',7,1),('low',13,2)]:
 seen=set()
 for sh in range(n):
  a=json.loads((P/f'independent_genuine_{kind}_q{q}_s{sh}.json').read_text());assert a['status']=='PASS'and a['hecke']==q
  for r in a['cases']:assert r['source']not in seen and r['all_geometric_kernels']==q+1;seen.add(r['source']);kernel_count+=q+1
  genuine_count+=len(a['cases'])
 native=rd(f'parent_native_{"high"if kind=="high"else"low"}_q{q}.m');native=native if kind=='high'else native[0];assert seen=={r[0]for r in native}
assert genuine_count==336 and kernel_count==3924
cc=[json.loads((P/f'independent_character_s{i}.json').read_text())for i in range(4)];assert all(v['status']=='PASS'for v in cc)and sum(v['source_orbits']for v in cc)==149508 and sum(v['admitted_orbits']for v in cc)==12604 and sum(v.get('conservatively_retained_poles',0)for v in cc)==1
cache=json.loads((P/'independent_moduli_cache.json').read_text());assert cache['status']=='PASS'and cache['low_actual_source_points']==20404 and cache['high_retained_source_points']==12604 and cache['all_mixed_normalizing_automorphism_labels_bound_to_actual_projected_points']
cs=rd('independent_cusp_principal.m');assert cs==[[1,2,12,1],[2,4,12,1]];lab=json.loads((P/'independent_cusp_A13_labels.json').read_text());assert lab['status']=='PASS'and lab['geometric_cusps']==56 and len(lab['nonzero_A13_columns'])==4
for v in json.loads((P/'input_manifest.json').read_text()):assert sha(REPO/v['source'])==sha(REPO/v['copy'])==v['sha256']
markers={}
for i in[0,1]:markers[f'd31b_C48_character_high_s{i}_v3']='D31B_C48_HIGH_CHARACTER_COMPLETE'
for kind,n,marker in[('high',4,'D31B_C48_HIGH_CLASS_COMPLETE'),('low',3,'D31B_C48_LOW_CLASS_COMPLETE')]:
 for i in range(n):markers[f'd31b_C48_{kind}_class_s{i}']=marker;markers[f'd31b_C48_class_native_{kind}_s{i}_v4']='D31B_C48_CLASS_IDEAL_REPLAY_COMPLETE'
markers.update({'d31b_C48_low_geometry':'D31B_C48_LOW_GEOMETRY_COMPLETE','d31b_C48_label_order':'D31B_C48_ACTION_LABELS_COMPLETE','d31b_C48_merge_low_patterns':'D31B_C48_LOW_PATTERN_SIEVE_COMPLETE','d31b_C48_dedup_low':'D31B_C48_LOW_ORBITS_COMPLETE','d31b_C48_cusp_corrections_v2':'D31B_C48_CUSP_CORRECTIONS_COMPLETE','d31b_C48_independent_cusp_principal_v2':'D31B_C48_NATIVE_CUSP_PRINCIPAL_COMPLETE','d31b_C48_independent_moduli_cache':'D31B_C48_INDEPENDENT_ALL_MODULI_CACHE_COMPLETE'})
for q in[7,13]:
 for i in[0,1]:markers[f'd31b_C48_hecke{q}_s{i}']='D31B_C48_HEK_BATCH_COMPLETE'
 markers[f'd31b_C48_parent_native_q{q}']='D31B_C48_PARENT_NATIVE_COMPLETE';markers[f'd31b_C48_parent_native_low_q{q}']='D31B_C48_PARENT_NATIVE_LOW_COMPLETE'
for q,shards in[(7,[0,1]),(13,[0,1,3])]:
 for i in shards:markers[f'd31b_C48_low_hecke{q}_s{i}']='D31B_C48_LOW_HEK_BATCH_COMPLETE'
markers['d31b_C48_low_hecke13_tail41265']='D31B_C48_LOW_HEK_BATCH_COMPLETE'
for kind,q,n in[('high',7,1),('high',13,2),('low',7,1),('low',13,2)]:
 for i in range(n):markers[f'd31b_C48_genuine_{kind}_q{q}_s{i}']=f'D31B_C48_GENUINE_{kind.upper()}_COMPLETE'
for i in range(4):markers[f'd31b_C48_independent_character_s{i}'+('_v2'if i==3 else'')]='D31B_C48_INDEPENDENT_CHARACTER_COMPLETE'
jobs=[]
for sp in sorted((HOME/'jobs').glob('d31b_C48*.status.json')):
 name=sp.name.removesuffix('.status.json');s=json.loads(sp.read_text());assert s['state']=='EXITED';assert sha(HOME/s['script'])==s['script_sha256'];log=(HOME/s['log']).read_text();ok=name in markers
 if ok:assert s['actual_exit_code']==0 and markers[name]in log and not any(t in log for t in['Runtime error','User error:','Traceback (most recent call last):']),name
 jobs.append({'purpose':name,'state':s['state'],'actual_exit_code':s['actual_exit_code'],'complete_certificate':ok,'diagnostic_only':not ok,'script':str((HOME/s['script']).relative_to(REPO)),'log':str((HOME/s['log']).relative_to(REPO)),'status':str(sp.relative_to(REPO))})
assert all(any(v['purpose']==n for v in jobs)for n in markers)
active=[]
for q in Path('/proc').iterdir():
 if not q.name.isdigit():continue
 try:
  cmd=[v.decode(errors='replace')for v in(q/'cmdline').read_bytes().split(b'\0')if v]
  if cmd and(Path(cmd[0]).name.startswith('magma')or Path(cmd[0]).name.startswith('python'))and any(Path(v).stem.startswith('d31b_C48')for v in cmd[1:]):active.append(int(q.name))
 except OSError:pass
assert not active,active
summary={'status':'COMPLETE EXCLUSION','group':'C48','degree':10,'source_genus':37,'rational_gonality':16,'target_genus':13,'rational_target_Jacobian_invariants':[2,4,4,120,240],'finite_quotient_invariants':[2,2,24,48],'high_character_source_orbits':149508,'high_character_retained':12604,'high_conservative_pole_fallbacks':1,'target_class_queries':9169,'independent_both_order_smoothing_identities':9078,'independent_direct_factorbase_hits':91,'lower_support_patterns':119986,'lower_class_survivors':243,'lower_simultaneous_diamond_orbits':201,'final_high_Hecke_cases':107,'final_lower_Hecke_cases':201,'final_total_Hecke_cases':308,'independent_full_source_prime_references':1681,'genuine_marked_isogeny_source_contexts':336,'actual_geometric_cyclic_kernels':3924,'final_remaining':0,'exceptional_A13_cusp_pair_differences_independently_principal':2,'all_mixed_aut_labels_and_high_cache_bindings_independently_verified':True,'all_input_hashes_unchanged':True,'all_recorded_jobs_exited':True,'active_workers':active,'paused_workers':0,'automatic_launchers':0}
(P/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={q for q in P.rglob('*')if q.is_file()and q.name!='certificate_manifest.json'and'__pycache__'not in q.parts};paths.add(HOME/'run_job.py')
for job in jobs:
 for key in['script','log','status']:paths.add(REPO/job[key])
paths.update(q for q in(HOME/'computations').glob('d31b_C48*')if q.is_file());paths.update(q for q in(HOME/'shared/mdmagma/v2').rglob('*')if q.is_file())
manifest={'status':'PASS','group':'C48','degree':10,'summary_sha256':sha(P/'certificate_summary.json'),'files':[{'path':str(q.relative_to(REPO)),'sha256':sha(q),'bytes':q.stat().st_size}for q in sorted(paths)],'jobs':jobs};(P/'certificate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('D31B_C48_FINAL_MANIFEST_COMPLETE',len(paths),len(jobs),sha(P/'certificate_manifest.json'))
