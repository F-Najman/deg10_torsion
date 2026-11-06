from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def job(stem):
 status=json.loads((H/'jobs'/(stem+'.status.json')).read_text());assert status['state']=='EXITED' and status['actual_exit_code']==0
 return [status['script'],'jobs/'+stem+'.log','jobs/'+stem+'.status.json']
def evidence(paths):
 return [{'path':s,'sha256':sha(H/s)} for s in sorted(set(paths))]
prefix='continuation_20261029/'
proof=['exclusions/continuation_20261029.tex']
common=proof+[prefix+'cusp_tate_replay.json',prefix+'replay_cusp_tate.py',prefix+'new_cases_independent_replay.json']+job('cont_replay_hecke_cases')
helpers=[prefix+'deps/'+n for n in ['fast_hecke.m','fast_hecke_11.m','fast_diamond_hash_m2.m']]
ledger=json.loads((H/'certified_results.json').read_text());present={tuple(r['group']) for r in ledger['exclusions']}
rows=[]
files=proof+job('cont_winding135')+job('cont_image135_a16_no2')+job('cont_replay135')+[prefix+'C135_independent_replay.json']
files += [str(p.relative_to(H)) for p in D.glob('image135_a16_q31_*.json')]+[prefix+'krylov_N135_chi4_order36.json']
rows.append(([1,135],'Dyadic global Hecke criterion with diamond 16 of order 9; its finite image has odd reduction order at 31',files))
files=common+helpers+[prefix+'deps/'+n for n in ['fast_full_aut_hash_m3_15.m','fast_full_aut_m3_21.m','dt_fast_differential_expansions.m']]
for stem in ['cont_local_3_21_q4','cont_sieve_3_21_D5','cont_double_cusp_3_21_D5_dt','cont_replay321']:files+=job(stem)
files += [prefix+s for s in ['cont_local_3_21_q4_census.m','sieve_3_21_D5_arithmetic.txt','double_cusp_3_21_D5_dt_arithmetic.txt','M3N21_independent_replay.json','degree_independent_321_manifest.json']]
files += [str(p.relative_to(H)) for p in (D/'degree_independent_321').iterdir()]
rows.append(([3,21],'Relative-quintic rank-zero Hecke sieve at F4 with complete doubled-cusp Riemann-Roch exclusion',files))
for m in [3,4]:
 files=common+helpers
 stem=f'cont_local_{m}_24_q'+('7' if m==3 else '5');files+=job(stem)+[prefix+stem+'_census.m']
 if m==3:
  files += [prefix+'deps/fast_full_aut_hash_m3_15.m',prefix+'deps/fast_full_aut_hash_m3_24.m']
  stems=[f'cont_sieve_3_24_D5_shard{i}' for i in range(6)]+['cont_sieve_3_24_D5_lower_supports','cont_sieve_3_24_cubic31']
 else:
  files += [prefix+'deps/fast_full_aut_hash_m4_24_root.m']
  stems=[f'cont_sieve_4_24_D5_d{d}_r{i}' for d,c in [(4,4),(5,11)] for i in range(1,c+1)]
 for stem in stems:files+=job(stem)+[prefix+stem+'_arithmetic.txt']
 rows.append(([m,24],'Relative-quintic rank-zero Hecke sieve; complete full-level automorphism orbits, and odd-prime cuspidal lifting. For m=3, T31 handles every quadratic cusp complement.',files))
files=common+job('cont_winding31')+job('cont_winding155')+job('cont_C155_local2')
files += [str(p.relative_to(H)) for p in D.glob('krylov_N155_*.json')]
files += [prefix+s for s in ['Orlic_2407.14512.pdf','Orlic_2407.14512.txt','Orlic_source.json','C155_local2_arithmetic.txt']]
rows.append(([1,155],'Diamond descent to a quadratic curve with 5-torsion and a 31-isogeny; good reduction at 2 is impossible, and Q-gonality at least five on X1(31)/<2> excludes cuspidal quadratic reduction.',files))
assert not any(tuple(g) in present for g,method,files in rows)
for g,method,files in rows:ledger['exclusions'].append({'group':g,'verified':True,'method':method,'evidence':evidence(files)})
(H/'certified_results.json').write_text(json.dumps(ledger,indent=2)+'\n')
print('Registered five verified exclusions:',[g for g,m,f in rows])
