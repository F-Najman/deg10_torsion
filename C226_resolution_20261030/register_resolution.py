from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
L=json.loads((H/'certified_results.json').read_text());before=json.loads((H/'candidates.json').read_text());(P/'frontier_before.json').write_text(json.dumps(before,indent=2)+'\n')
changed={'continuation_20261029b/C226/residual_classes.tex','exclusions/isogeny_transfer.tex'}
for row in L['realizations']+L['exclusions']+L['auxiliary_results']:
 for e in row['evidence']:
  if e['path'] in changed:e['sha256']=sha(H/e['path'])
old=next(r for r in L['auxiliary_results'] if r['name']=='decic_C226_reduction_to_two_rational_divisor_classes')
old['statement']='Every decic C2 x C26 trace divisor, modulo level automorphisms, belongs to one of two explicit rational classes; the degree-nine-plus-cusp branch is excluded. Effectivity and the two representatives are determined in complete_decic_C226_case.'
old['residual_classes_before_resolution']=old.pop('remaining_classes')
def job(stem,marker=None):
 p=H/'jobs'/f'{stem}.status.json';s=json.loads(p.read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0;assert sha(H/s['script'])==s['script_sha256']
 if marker:assert marker in (H/s['log']).read_text(),stem
 return [s['script'],s['log'],str(p.relative_to(H))]
pre='C226_resolution_20261030/'
files=[pre+s for s in ['resolution.tex','coefficient_table.tex','models_compact.json','exact_points_30.m','exact_points_24.m','complete_resolution_inputs.m','complete_native_result.m','independent_class_binding.json','coordinate_fields.json','witness_replay_30.json','witness_replay_24.json','witness_replay.py']]
files+=job('resolve226_complete_native','Q226_COMPLETE_RESOLUTION_CERTIFICATE_DONE')
for i in [30,24]:files+=job(f'resolve226_replay_witness{i}','C226_EXACT_WITNESS_REPLAY_COMPLETE')
files+=job('resolve226_bind_classes')+job('resolve226_coordinate_fields')
files+=['exclusions/isogeny_transfer.tex']
for e in old['evidence']:assert sha(H/e['path'])==e['sha256']
allfiles=sorted(set(files+[e['path'] for e in old['evidence']]))
def evidence(ps):return [{'path':s,'sha256':sha(H/s)} for s in sorted(set(ps))]
assert not any(r['group']==[2,26] for r in L['realizations']);assert not any(r['group']==[1,52] for r in L['exclusions'])
assert json.loads((P/'independent_class_binding.json').read_text())['status']=='PASS'
for i in [24,30]:assert json.loads((P/f'witness_replay_{i}.json').read_text())['status']=='PASS'
L['realizations'].append({'group':[2,26],'verified':True,'method':'Two explicit irreducible degree-ten curves with independent exact torsion checks and good reductions of order 52 at 41 and 43','number_of_Galois_and_level_orbits':2,'evidence':evidence(files)})
L['exclusions'].append({'group':[1,52],'verified':True,'method':'Complete decic X1(2,26) classification; six certified nonsquare root differences on each representative prevent every dual 2-isogeny kernel from having a rational half-point','evidence':evidence(allfiles)})
L['auxiliary_results'].append({'name':'complete_decic_C226_case','verified':True,'statement':'Both residual rational classes have h0=1 and a unique irreducible noncuspidal degree-ten representative. Exactly two Galois and level orbits occur, each with full torsion C2 x C26 and quintic j-invariant. No dual 2-isogeny kernel can be halved over either field.','effective_classes':2,'full_torsion_group':[2,26],'minimal_level_field_degree':10,'j_invariant_degree':5,'cyclic52_excluded':True,'evidence':evidence(allfiles)})
(H/'certified_results.json').write_text(json.dumps(L,indent=2)+'\n');print('REGISTERED_C226_OCCURRENCE_TWO_CLASS_CLASSIFICATION_AND_C52_EXCLUSION')
