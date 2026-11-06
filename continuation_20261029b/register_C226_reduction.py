from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];P=Path(__file__).resolve().parent;C=P/'C226';pre='continuation_20261029b/C226/'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def job(stem,marker=None):
 p=H/'jobs'/(stem+'.status.json');s=json.loads(p.read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0
 assert sha(H/s['script'])==s['script_sha256']
 if marker:assert marker in (H/s['log']).read_text()
 return [s['script'],s['log'],str(p.relative_to(H))]
files=[pre+s for s in ['residual_classes.tex','frobenius.m','filter26.m','target26_RR.log','target26_RR.log.witnesses.py','proven_picard.m','picard_template.m','picard_completion.json','native_replay_template.m','all_picard_witnesses.m','native_replay_completion.json','cusp_preimages.m','small_cusp_divisors.m','small_cusp_divisors_v2.m','cusp_lattice_replay.json','independent_sieve_replay.json','stabilizer_j0_replay.json','j0_twists.m','second_prime_RR.m','F5_RR.m','F5_effective_24.txt','F5_effective_30.txt','two_classes_stabilizers.json']]
files+=['continuation_20261029b/literature_inputs.json','continuation_20261029b/cusp_tate.json','continuation_20261029b/independent_sources_all.json','continuation_20261029b/C240/portable/portable_forward.m']
files += [str(p.relative_to(H)) for p in C.glob('source*.m')]+[str(p.relative_to(H)) for p in (C/'tate_alignment_inputs').iterdir()]+[str(p.relative_to(H)) for p in (P/'deps').glob('*.m')]
for stem in ['late226_source','late226_filter26','late226_picard_setup','late226_cusp_preimages','late226_reduce_cusps','late226_replay_sieve','late226_j0_twists_v2','late226_stabilizer_replay','late226_second_prime_v3','late226_F5_RR','late226_two_classes_stabilizers','late226_replay_tate']:files+=job(stem)
for k in range(8):
 files+=job(f'late226_picard_s{k}','C226_PICARD_SHARD_COMPLETE')+job(f'late226_native_replay_s{k}','C226_ALL_PICARD_IDEAL_REPLAYS_COMPLETE')
 files+=[pre+s for s in [f'picard{k}.m',f'picard{k}_RR.log',f'picard{k}_RR.log.witnesses.py',f'native_ideal_replay_shard{k}.m']]
assert json.loads((C/'independent_sieve_replay.json').read_text())['status']=='PASS'
assert json.loads((C/'stabilizer_j0_replay.json').read_text())['status']=='PASS'
L=json.loads((H/'certified_results.json').read_text());name='decic_C226_reduction_to_two_rational_divisor_classes';assert not any(r['name']==name for r in L['auxiliary_results'])
L['auxiliary_results'].append({'name':name,'verified':True,'statement':'A hypothetical decic C2 x C26 trace divisor, modulo level automorphisms, belongs to one of two explicit rational divisor classes; the degree-nine-plus-cusp branch is excluded. Effectivity of the two classes is unresolved.','remaining_classes':2,'does_not_exclude_or_prove_occurrence':[2,26],'evidence':[{'path':s,'sha256':sha(H/s)} for s in sorted(set(files))]})
(H/'certified_results.json').write_text(json.dumps(L,indent=2)+'\n');print('REGISTERED_C226_TWO_CLASS_REDUCTION_NO_GROUP_EXCLUSION')
