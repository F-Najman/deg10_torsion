#!/usr/bin/env python3
"""Read-only frozen certificate validation; no arithmetic workers are launched."""
from pathlib import Path
import ast,hashlib,json
P=Path(__file__).resolve().parent;ROOT=P.parents[3]
a=json.loads((P/'final_certificate_summary.json').read_text());N=a['level'];assert N in[85,93] and a['complete_exclusion'] and a['parent_mathematical_review']=='PASS'
for ent in a['files']:
 p=ROOT/ent['path'];assert p.is_file()and hashlib.sha256(p.read_bytes()).hexdigest()==ent['sha256'],p
snap=json.loads((P/'inherited_snapshot.json').read_text());assert snap['status']=='PASS'
for ent in snap['files']:
 for key in['path','source']:
  p=ROOT/ent[key];assert hashlib.sha256(p.read_bytes()).hexdigest()==ent['sha256'],p
primary=json.loads((P/'cusp_FI_primary.json').read_text());parent=json.loads((P/'parent_cusp_rank_replay.json').read_text());ops=json.loads((P/'operator_and_frobenius_control.json').read_text());arcs=json.loads((P/'genuine_cusp_arcs_verified.json').read_text())
assert primary['complete']and primary['all_tested_ranks10']and parent['status']=='PASS'
assert primary['minimum_rational_prefix_cost']==parent['rational']['minimum_normalized_prefix_cost']==(32 if N==85 else 30)
assert primary['rational_relation_dimension']==parent['rational']['relation_dimension']==(13 if N==85 else 5)
assert not primary['mixed_failures']and sum(row['configurations']for row in primary['mixed_families'])==parent['total_mixed_patterns']==(2330064 if N==85 else 278262)
assert all(row['all_full_rank']for row in parent['mixed_families'])
assert hashlib.sha256((P/'cusp_blocks.json').read_bytes()).hexdigest()==parent['source_blocks_sha256']
assert ops['mu_local_killer']=='T2-I'and ops['canonical_cusp_frobenius_verified']and ops['diamond_stability_verified']and ops['all_blocks_recomputed_equal']
assert len(ops['W'])==(198 if N==85 else 248)and arcs['status']=='PASS'
def lit(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
for row in arcs['arcs']:
 p=ROOT/row['primary_path'];q=ROOT/row['control_path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['primary_sha256']and hashlib.sha256(q.read_bytes()).hexdigest()==row['control_sha256'];u,v=lit(p),lit(q)
 assert u[6]==v[6]and all(len(z)==row['retained_coefficients']for z in u[6])and row['all_retained_coefficients_identical']
if N==85:
 p=json.loads((P/'positive_polynomials_exact.json').read_text());assert p['k2']==[-1,0,4,4,1]and p['k3']==[4,12,12,6,1]
 def product(f,g):return[sum(f[i]*g[k-i]for i in range(len(f))if 0<=k-i<len(g))for k in range(len(f)+len(g)-1)]
 assert product(*[v['T2_characteristic']for v in p['only_possible_positive_factors']])==p['k2']
 assert product(*[v['T3_characteristic']for v in p['only_possible_positive_factors']])==p['k3']
 r=json.loads((P/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==4 and not r['unresolved']and r['source_closed_points']==128
 assert all(v['positive_degree']==80 and v['rank']==225 and v['convention']==2 and v['target_diamond_translate']==0 for v in r['results'])
 s=json.loads((P/'R1_independent_sage.json').read_text());assert s['status']=='PASS'and s['cases']==4 and all(v['independent_joint_rank']==225 for v in s['reports'])
 s=json.loads((P/'jets10_independent.json').read_text());assert s['status']=='PASS'and s['ordinary_85_torsion_families']==12 and s['all_finite_unit_ranks']==224
 s=json.loads((P/'isogenies10_independent.json').read_text());assert s['status']=='PASS'and s['geometric_kernels']==16 and s['physical_isogenies_independently_rebuilt']
else:
 s=json.loads((P/'cusp_descent_verified.json').read_text());assert s['status']=='PASS'and s['h_scalar']==32 and s['h_fixes_every_quadratic_and_decic_closed_cusp']and s['h_stability_descends_31_level_to_quintic']and s['prime31_absent_from_accepted_Phi5']and s['prime17_is_present_in_Phi5']
life=json.loads((P/'lifecycle_final.json').read_text());assert life['status']=='PASS'and life['all_jobs_exited']and not life['active_workers']
for row in life['jobs']:
 st=json.loads((ROOT/row['status_path']).read_text());assert st['state']=='EXITED'and st['actual_exit_code']==0
 assert hashlib.sha256((ROOT/'work/deg10'/st['script']).read_bytes()).hexdigest()==st['script_sha256']
print('C%d_DEG10_FROZEN_CERTIFICATE_READONLY_PASS'%N,len(a['files']),'evidence hashes;',len(snap['files']),'inherited originals unchanged;',len(life['jobs']),'jobs EXITED')
