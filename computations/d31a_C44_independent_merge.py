from pathlib import Path
from collections import Counter
import json,hashlib,time
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C44';H=C.parents[1];start=time.monotonic();all_audits={};casekeys={};q7keys={};coverage10=set();coverage5=set();reps10=set();reps5=set();pairs=[];inputs={};shards=[]
def tag(x):return json.dumps(x,separators=(',',':'))
for i in range(6):
 p=C/('independent_genuine_shard%d_v2.json'%i);r=json.loads(p.read_text());assert r['status']=='PASS'and r['shard']==i and r['shards']==6
 stp=H/'jobs'/('d31a_C44_independent_isogenies_shard%d_v2.status.json'%i);st=json.loads(stp.read_text());assert st['state']=='EXITED'and st['actual_exit_code']==0;assert'D31A_C44_INDEPENDENT_GENUINE_SHARD_COMPLETE'in(H/st['log']).read_text()
 cp=C/('independent_genuine_shard%d_v2.checkpoint.jsonl'%i);rows=[json.loads(s)for s in cp.read_text().splitlines()];assert rows[0]['type']=='header'and rows[-1]['type']=='complete'and rows[-1]['shard']==i
 ar=[{k:v for k,v in row.items()if k!='type'}for row in rows if row['type']=='unique_source_check'];assert ar==r['all_bound_moduli_keys'];src_rows=[row for row in rows if row['type']=='case_check']
 for a in ar:
  ell,d=a['ell'],a['source_degree'];key=(ell,tag(a['source_key']));assert d==a['source_key'][0]and a['diamond_key'][0]==d
  counts={tag(z):n for z,n in a['geometric_kernel_target_counts']};ns={tag(z):n for z,n in a['target_norms']};assert sum(counts.values())==ell+1
  assert set(counts)==set(ns)
  for z,n in ns.items():assert json.loads(z)[0]*n==d*counts[z]
  assert sum(json.loads(z)[0]*n for z,n in ns.items())==d*(ell+1)
  if key in all_audits:assert all_audits[key]==a
  else:all_audits[key]=a
 for row in src_rows:
  dic=casekeys if row['stage']=='newA5'else q7keys;ci=row['case'];assert ci not in dic and row['all_target_norms_and_negative_terms_bound'];dic[ci]=row['source_keys']
 for v in r['stages']:
  if v['stage']=='newA5':
   assert v['case_ids']==[ci for ci in range(1,957)if(ci-1)%6==i];coverage10.update(v['diamond_coverage10_keys']);coverage5.update(v['diamond_coverage5_keys']);reps10.update(v['representative_source10_keys']);reps5.update(v['representative_source5_keys']);pairs.extend(map(tuple,v['normalized_pair_keys']))
  else:assert i==0 and v['stage']=='residualA7'and v['case_ids']==[1,2]
 for q in[p,cp,stp,H/st['log'],H/st['script']]:inputs[str(q.relative_to(H))]=hashlib.sha256(q.read_bytes()).hexdigest()
 shards.append({'shard':i,'case_count':len([x for x in src_rows if x['stage']=='newA5']),'unique_operator_source_pairs':len(ar),'seconds':r['elapsed_seconds']})
assert set(casekeys)==set(range(1,957))and set(q7keys)=={1,2};assert q7keys[1]==casekeys[510]and q7keys[2]==casekeys[511]
assert len(coverage10)==5854 and len(coverage5)==60 and len(reps10)==590 and len(reps5)==6
assert len(pairs)==360 and set(pairs)=={(a,b)for a in reps5 for b in coverage5};assert all((a,a)in pairs for a in reps5)
assert len(all_audits)==652 and sum(ell+1 for ell,ky in all_audits)==3916 and Counter(ell for ell,ky in all_audits)=={5:650,7:2}
assert len({ky for ell,ky in all_audits})==650
for ci,keys in casekeys.items():
 for ky in keys:assert(5,ky)in all_audits
for ci,keys in q7keys.items():
 for ky in keys:assert(7,ky)in all_audits
frob=json.loads((C/'degree10_frobenius_audit.json').read_text());assert frob['complete']and frob['independent_rational_matrix_replay']and frob['noncuspidal_closed_counts']==[0,0,0,0,60,160,240,920,2260,5854]
mat=json.loads((C/'independent_matrices_verified.json').read_text());geo=json.loads((C/'independent_jet_geometry_verified.json').read_text());assert mat['status']==geo['status']=='PASS'
assert[(v['stage'],v['cases'],v['rank105'],v['survivors'])for v in mat['stages']]==[('newA5',956,954,[510,511]),('residualA7',2,2,[]),('oldA5',226,226,[])]
assert pow(3,8,44)==5
res={'status':'PASS','new_A5_cases':956,'new_A5_rejections':954,'A5_surviving_cases_refined_by_A7':[510,511],'new_A7_rejections':2,'old_A5_rejections_independently_replayed':226,'all_entrywise_canonical_matrices_and_quotient_product_identities_verified':True,'all_saved_higher_jet_geometry_verified':True,'distinct_marked_source_points':650,'operator_source_pairs_checked':652,'unique_genuine_geometric_kernels':3916,'manual_Velu_coordinate_checks_unique':3916,'all_positive_and_negative_divisor_norm_multiplicities_bound':True,'source_degree10_closed_places_covered':5854,'source_degree5_closed_places_covered':60,'degree5_pair_patterns':360,'all_six_diagonal_pairs_included':True,'quintic_cusp_complements_killed_by_A5':'5=3^8mod44, so Galois5 fixes every closedF3cusp orbit','shards':shards,'audit_files_sha256':inputs,'elapsed_seconds':time.monotonic()-start}
(C/'independent_final_verification.json').write_text(json.dumps(res,indent=2)+'\n');print('D31A_C44_INDEPENDENT_FINAL_VERIFICATION_COMPLETE',956,2,226,652,3916,time.monotonic()-start,flush=True)
