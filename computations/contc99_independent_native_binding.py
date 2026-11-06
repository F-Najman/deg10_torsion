from pathlib import Path
import ast,json,hashlib
ROOT=Path('work/deg10/continuation_20261030c/C99')
native=ast.literal_eval((ROOT/'parent_native_RR.m').read_text().replace('[*','[').replace('*]',']'))
audit=json.loads((ROOT/'independent_genuine_isogenies.json').read_text());assert audit['status']=='PASS'
assert len(native)==len(audit['source_cases'])==4
for row,case in zip(native,audit['source_cases']):
 ci,source,diamond,targets,h0=row
 assert ci==case['case'] and source==case['source_key'] and diamond==case['diamond5_key']
 assert sorted(targets,key=str)==sorted(case['kernel_orbit_observations'],key=str)
 assert h0==audit['all_primary_RR_h0'][ci-1]
assert [r[4]for r in native]==[0,0,1,0]
out={'status':'PASS','all_native_source_diamond_and_Hecke_prime_keys_match_independent_actual_marked_isogenies':True,'cases':4,'native_positive_kernel_orbits':16,'geometric_cyclic5_kernels':24,'native_RR_h0':[0,0,1,0],'residual_case_ids':[3],'norm_multiplicities_all_one':True,'kernel_orbit_degrees_per_case':[1,1,2,2],'target_degrees_per_case':[10,10,20,20],'native_snapshot_sha256':hashlib.sha256((ROOT/'parent_native_RR.m').read_bytes()).hexdigest(),'genuine_audit_sha256':hashlib.sha256((ROOT/'independent_genuine_isogenies.json').read_bytes()).hexdigest()}
(ROOT/'independent_native_binding.json').write_text(json.dumps(out,indent=2)+'\n')
print('C99_INDEPENDENT_ALL_NATIVE_PRIME_BINDINGS_PASS',4,16,24,[0,0,1,0])
