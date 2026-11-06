"""Freeze the complete order46 exclusion after mathematical and independent reviews."""
from pathlib import Path
import hashlib,json,re
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;C=H/'continuation_20261031a/C46';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();mf=C/'certificate_manifest.json';assert not mf.exists()
source=json.loads((C/'source_cases.json').read_text());rows=source['cases'];assert len(rows)==1023 and sum(r['quotient_reject']for r in rows)==986
q=json.loads((C/'independent_quotient_arithmetic.json').read_text());assert q['status']=='PASS'and q['support_cases']==1023 and q['quotient_rejections']==986 and q['source_components']==676 and q['source_degree10_points']==6292 and q['all_marked_source_and_target_points_verified']and q['all_residue_norms_and_class_sums_verified']
r=json.loads((C/'parent_exceptional_summary.json').read_text());assert r['status']=='PASS'and r['source_supports']==37 and r['source_components']==54 and r['native_two_generator_prime_references']==292 and r['genuine_geometric5_kernels']==324 and r['all_h0_zero']and r['original_source_projective_coordinate_bindings_verified']
provenance=0
for p in(C/'inputs').glob('*provenance.json'):
 d=json.loads(p.read_text());assert d['status']=='PASS'
 for e in d['files']:
  a,b=ROOT/e['source'],ROOT/e['snapshot'];assert b.is_relative_to(C)and sha(a)==sha(b)==e['sha256'];provenance+=1
assert provenance>=292
pdf=C/'inputs/Varivoda_2610.00219v1.pdf';assert sha(pdf)=='246a8196fe7d15acee9322297c2d5d4f618ad6bcfed6e6b2e0a30c5f3f9ba481'
tex=(C/'resolution.tex').read_text();assert 'Varivoda26'in tex and 'Corollary~4.14'not in tex and 'Table~2'in tex
paths=set(p for p in C.rglob('*')if p.is_file());paths.add(Path(__file__).resolve());jobs=[]
for p in sorted((H/'jobs').glob('*.status.json')):
 d=json.loads(p.read_text());stem=d['purpose']
 if not(stem.startswith('d31a_C46')or stem.startswith('d31a_46')):continue
 assert d['state']=='EXITED';script=H/d['script'];assert sha(script)==d['script_sha256'];log=H/d['log'];paths.update([p,script,log]);code=d['actual_exit_code'];diag=stem.startswith('d31a_46_boundary')
 assert (diag and code==-14)or(not diag and code==0)
 jobs.append({'purpose':stem,'actual_exit_code':code,'certificate':not diag,'diagnostic_reason':'Superseded cusp computation reached native alarm; not used in the gonality proof'if diag else None,'status_path':str(p.relative_to(H))})
assert len(jobs)==16
markers={'d31a_C46_new_supports':'D31A_C46_NEW_SUPPORTS_COMPLETE all_excluded=true','d31a_C46_quotient_replay':'D31A_C46_QUOTIENT_REPLAY_COMPLETE cases=1023 points=676 quotient_rejects=986 exceptions=37','d31a_C46_independent_source_counts':'D31A_C46_INDEPENDENT_COUNTS_DONE complete=true','d31a_46_parent_merge_exceptional_v2':'C46_PARENT_MERGED_EXCEPTIONAL_COMPLETE 37 54 292 324'}
for stem,marker in markers.items():assert marker in(H/f'jobs/{stem}.log').read_text()
for p in(H/'computations').glob('d31a_46_parent_*'):paths.add(p)
summary={'status':'PASS','parent_mathematical_review':'PASS','group':[1,46],'degree':10,'complete_exclusion':True,'method':'Rank-zero odd-prime cusp lifting and rational gonality19; complete finite source support enumeration; exact X1(23) rational class quotientC23 and independently reconstructed marked A5 nonprincipality tests on the37 residual supports','source_genus':45,'rational_J_rank':0,'rational_gonality':19,'prime':3,'noncuspidal_closed_counts':[0,0,0,55,44,55,264,781,1958,6292],'cusp_signature':'1^22 11^2 globally and overF3','inherited_negative_support_tests':782,'new_support_tests':1023,'quotient_rejections':986,'A5_nonprincipal_exceptions':37,'native_prime_references':292,'genuine5_kernels':324,'source_components_in_Hecke_tests':54,'inherited_source_snapshot_bindings_checked':provenance,'old_cusp_alarms_are_not_premises':True,'initial_unused_prepared_scripts_and_failed_merge_are_not_certificates':True,'all_case_jobs_exited':True,'classification_complete':False}
sp=C/'certificate_summary.json';sp.write_text(json.dumps(summary,indent=2)+'\n');paths.add(sp)
entries=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)];mf.write_text(json.dumps({'status':'PASS','parent_mathematical_review':'PASS','degree':10,'group':[1,46],'scope':summary['method'],'summary_sha256':sha(sp),'files':entries,'jobs':jobs,'diagnostic_scope':'Boundary partial completions and native-alarm failures are retained but are not premises; unused prepared unsplit scripts and the first merge parser trial are not accepted computations'},indent=2)+'\n')
for e in entries:assert sha(ROOT/e['path'])==e['sha256']
print('D31A_C46_FINAL_CERTIFICATE_COMPLETE',len(entries),len(jobs),sha(mf))
