from pathlib import Path
import json,re,hashlib,os
P=Path(__file__).resolve().parent;HOME=P.parents[1];REPO=HOME.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return json.loads(p.read_text().replace('[*','[').replace('*]',']'))
def assigned(p,n):return json.loads(re.search(n+r'\s*:=\s*(.*?);',p.read_text(),re.S)[1].replace('[*','[').replace('*]',']'))
def fibres(p):
 ps=assigned(p,'C230Points');cs=assigned(p,'C230Candidates');fs=assigned(p,'C230Fibres');ks=[json.dumps(v,separators=(',',':'))for v in ps];return ps,cs,fs,set(ks),{tuple(sorted((ks[i-1],m)for i,m in c))for c in cs}
x=fibres(P/'exceptional_fibres_data.m');y=fibres(P/'independent/exceptional_fibres_data.m');assert x[3:]==y[3:]and(len(x[0]),len(x[1]),len(x[2]))==(1674,3414,1848)==(len(y[0]),len(y[1]),len(y[2]))
a=rd(P/'degree10_series.m');assert len(a[5])==256 and len(a[6])==84 and len(a[7])==4;fixed=[v for v in a[5]if v[6]==0];assert len(fixed)==168 and all(v[5]==0 and v[9]in[2,4,6]for v in fixed)
sm=json.loads((P/'three_maps_summary.json').read_text());assert sm['status']=='PASS'and sm['source_patterns']==3414 and sm['remaining_count']==24
hh=rd(P/'final_halving.m');assert len(hh[0])==72 and len(hh[1])==24 and not hh[2]
sg=json.loads((P/'independent_genuine_source_maps.json').read_text());assert sg['status']=='PASS'and sg['source_points']==1674 and sg['marked_target_maps']==5022
nr=rd(P/'independent_native_class_ideals.m');assert nr[1]==[[0,837,741,96],[1,837,741,96],[2,837,741,96]]
pp=json.loads((P/'parent_independent_halving.json').read_text());assert pp['status']=='PASS'
for v in json.loads((P/'input_manifest.json').read_text()):assert sha(REPO/v['source'])==sha(REPO/v['copy'])==v['sha256']
markers={'d31b_C60_degree10_series':'D31B_C60_D10_SERIES_COMPLETE','d31b_C60_mixed_fibres_v2':'D31B_C60_MIXED_FIBRES_COMPLETE','d31b_C60_quintic_2tower':'D31B_C60_QUINTIC_2TOWER_COMPLETE','d31b_C60_target30_setup':'D31B_C60_TARGET30_PROVEN_COMPLETE','d31b_C60_merge_three_maps':'D31B_C60_THREE_MAP_MERGE_COMPLETE','d31b_C60_final_halving':'D31B_C60_FINAL_HALF_COMPLETE','d31b_C60_genuine_source_maps':'D31B_C60_GENUINE_SOURCE_MAPS_COMPLETE','d31b_C60_independent_fibres':'D31B_C60_INDEPENDENT_FIBRES_COMPLETE','d31b_C60_independent_class_ideals':'D31B_C60_NATIVE_IDEAL_REPLAY_COMPLETE','d31b_C60_parent_halving_v2':'D31B_C60_PARENT_HALVING_COMPLETE'}
for i in range(3):markers[f'd31b_C60_three_maps_s{i}_v2']='D31B_C60_THREE_MAPS_COMPLETE'
jobs=[]
for sp in sorted((HOME/'jobs').glob('d31b_C60*.status.json')):
 name=sp.name.removesuffix('.status.json');s=json.loads(sp.read_text());assert s['state']=='EXITED';assert sha(HOME/s['script'])==s['script_sha256'];log=(HOME/s['log']).read_text();ok=name in markers
 if ok:assert s['actual_exit_code']==0 and markers[name]in log and not any(z in log for z in['User error:','Runtime error','Traceback (most recent call last):']),name
 jobs.append({'purpose':name,'state':s['state'],'actual_exit_code':s['actual_exit_code'],'complete_certificate':ok,'diagnostic_only':not ok,'script':str((HOME/s['script']).relative_to(REPO)),'log':str((HOME/s['log']).relative_to(REPO)),'status':str(sp.relative_to(REPO))})
assert all(any(v['purpose']==n for v in jobs)for n in markers)
active=[]
for z in Path('/proc').iterdir():
 if not z.name.isdigit():continue
 try:
  cmd=[v.decode(errors='replace')for v in(z/'cmdline').read_bytes().split(b'\0')if v]
  if cmd and(Path(cmd[0]).name.startswith('magma')or Path(cmd[0]).name.startswith('python'))and any(Path(v).stem.startswith('d31b_C60')for v in cmd[1:]):active.append(int(z.name))
 except OSError:pass
assert not active,active
summary={'status':'COMPLETE EXCLUSION','group':'C60','degree':10,'genus9_cusp_image_classes':256,'global_fixed_cusp_classes':168,'quintic_pullback_families':4,'mixed_character_series':84,'finite_projective_members':1848,'physical_full2_source_places':1674,'full2_candidate_divisors':3414,'target_level_maps':5022,'class_sieve_remaining':24,'final_marked_half_supports':72,'final_remaining':0,'independent_native_target_primes':2511,'independent_both_order_smoothing_identities':2223,'independent_direct_factorbase_classes':288,'all_input_hashes_unchanged':True,'all_recorded_jobs_exited':True,'active_workers':active,'paused_workers':0,'automatic_launchers':0}
(P/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={z for z in P.rglob('*')if z.is_file()and z.name!='certificate_manifest.json'and'__pycache__'not in z.parts};paths.add(HOME/'run_job.py')
for job in jobs:
 for k in['script','log','status']:paths.add(REPO/job[k])
paths.update(z for z in(HOME/'shared/mdmagma/v2').rglob('*')if z.is_file())
manifest={'status':'PASS','group':'C60','degree':10,'summary_sha256':sha(P/'certificate_summary.json'),'files':[{'path':str(z.relative_to(REPO)),'sha256':sha(z),'bytes':z.stat().st_size}for z in sorted(paths)],'jobs':jobs}
(P/'certificate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('D31B_C60_FINAL_MANIFEST_COMPLETE',len(paths),len(jobs),sha(P/'certificate_manifest.json'))
