from pathlib import Path
import hashlib,json
C=Path(__file__).resolve().parent;H=C.parents[1];R=C.parent;S=R/'Y228_shared';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
jobs={
'd31b_m428_Y_fourmap':'D31B_M428_Y_FOURMAP_COMPLETE',
'd31b_Y228_prepare_union':'D31B_Y228_UNION_NORMALIZED',
'd31b_m428_parent_class_ideals_v2':'D31B_M428_NATIVE_IDEAL_COMPLETE',
'd31b_m428_independent_fourmap_v3':'D31B_M428_INDEPENDENT_FOURMAP_COMPLETE',
'd31b_m428_independent_union':'D31B_M428_INDEPENDENT_UNION_COMPLETE',
'd31b_Y228_parent_native_v2':'D31B_Y228_PARENT_NATIVE_COMPLETE',
'd31b_m428_neighbors_A5':'D31B_M428_NEIGHBOR_A5_COMPLETE',
'd31b_Y228_independent_genuine_v2_s0':'D31B_Y228_INDEPENDENT_GENUINE_COMPLETE',
'd31b_Y228_independent_genuine_v2_s1':'D31B_Y228_INDEPENDENT_GENUINE_COMPLETE',
**{f'd31b_Y228_A5_v5_s{i}':'D31B_Y228_A5_COMPLETE'for i in range(6)}}
files=set()
def add(p):assert p.is_file(),p;files.add(p)
for n in ['resolution.tex','certify.py','build_manifest.py','fourmap_Y28.m','fourmap_native_class.log','fourmap_native_class.log.witnesses.py','fourmap_class_witnesses_v2.m','parent_native_class_ideals_v2.m','parent_independent_fourmap_v3.json','parent_independent_union.json','mathematical_review.json']:add(C/n)
for p in(C/'inputs').rglob('*'):
 if p.is_file():add(p)
for n in ['union_cases.m','union_normalization.json','A5_v5_summary.json','fast_primitive_inverse.m','fast_hecke_11.m','fast_hecke_11.sig','parent_native_all45_and_neighbors_v2.m','m428_neighbor_A5.m','model_m428_neighbors.m','independent_genuine_v2_s0.json','independent_genuine_v2_s1.json']:add(S/n)
for i in range(6):
 for n in[f'model_v5_s{i}.m',f'A5_v5_s{i}.m']:add(S/n)
for n in ['d31b_Y228_A5_common_v5.m','d31b_Y228_independent_genuine_common_v2.py']:add(H/'computations'/n)
add(R/'C56/projected_Y28.m');add(H/'run_job.py')
for p in(H/'shared/mdmagma/v2').rglob('*'):
 if p.is_file()and p.suffix in['.m','.spec','.txt']:add(p)
recs=[]
for stem,marker in jobs.items():
 p=H/'jobs'/(stem+'.status.json');s=json.loads(p.read_text());assert s['state']=='EXITED'and s['actual_exit_code']==0;assert marker in(H/s['log']).read_text();assert sha(H/s['script'])==s['script_sha256'];add(p);add(H/s['log']);add(H/s['script']);recs.append({'status_path':str(p.relative_to(H)),'actual_exit_code':0,'accepted_completion_marker':marker})
m={'status':'PASS','degree':10,'excluded_groups':[[4,28]],'parent_mathematical_review':'PASS','method':'Full-four-torsion residue parity, four X1(28) Picard tests, exact X1(2,28) Hecke tests and rational 2-isogeny quotients; independent marked Velu and native ideal/Riemann--Roch verification','files_sha256':{str(p.relative_to(H)):sha(p)for p in sorted(files)},'jobs':recs,'diagnostics_are_not_certificates':True}
p=C/'final_certificate_summary.json';assert not p.exists();p.write_text(json.dumps(m,indent=2)+'\n');print('M428_MANIFEST_FROZEN',len(files),sha(p))
