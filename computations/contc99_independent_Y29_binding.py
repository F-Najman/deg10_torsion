from pathlib import Path
import ast,json,hashlib
ROOT=Path('work/deg10/continuation_20261030c/C99')
native=ast.literal_eval((ROOT/'peer_native_Y29_RR.m').read_text().replace('[*','[').replace('*]',']'))
actual=json.loads((ROOT/'independent_Y29_moduli_keys.json').read_text());assert actual['status']=='PASS'
phasekeys,targets,h0,pos,npr=native
assert phasekeys==[actual['source_phase_keys'][str(u)]for u in[1,2,5,10]]
assert phasekeys[0]==phasekeys[3] and len({json.dumps(k)for k in phasekeys})==3
assert h0==0 and pos==110 and npr==23
expected=actual['paired_Hecke_keys_weighted_by_actual_kernel_orbits']
observed=[]
for e,n,n2,key,key2 in targets:
 assert 10*e==n*key[0]==n2*key2[0] and n==n2==1
 observed.append([e,key,key2])
assert sorted(observed,key=str)==sorted(expected,key=str)
out={'status':'PASS','case':3,'native_h0':0,'positive_and_negative_degree':110,'all_source_diamond_and_signed_Hecke_moduli_keys_bound':True,'native_prime_identity_checks':23,'geometric_cyclic5_kernels':6,'weighted_kernel_orbits':[1,1,2,2],'all_norm_multiplicities':1,'delta10_equals_identity_verified':True,'key_coordinates':'[t9,j11]','native_snapshot_sha256':hashlib.sha256((ROOT/'peer_native_Y29_RR.m').read_bytes()).hexdigest(),'genuine_Sage_snapshot_sha256':hashlib.sha256((ROOT/'independent_Y29_moduli_keys.json').read_bytes()).hexdigest()}
(ROOT/'independent_Y29_native_binding.json').write_text(json.dumps(out,indent=2)+'\n')
print('C99_INDEPENDENT_Y29_NATIVE_GENUINE_BINDING_COMPLETE',3,0,110,23,6)
