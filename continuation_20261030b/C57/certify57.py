#!/usr/bin/env python3
"""Read-only structural/hash validation of the complete C57 decic certificate."""
from pathlib import Path
import ast,hashlib,json
P=Path(__file__).resolve().parent
ROOT=P.parents[3]
assert (ROOT/'work/deg10').is_dir(),ROOT
summary=json.loads((P/'final_certificate_summary.json').read_text())
assert summary['parent_mathematical_review']=='PASS' and summary['complete_exclusion']
for ent in summary['files']:
 path=ROOT/ent['path']
 assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==ent['sha256'],path
old=json.loads((P/'inherited_evidence_verified.json').read_text())
assert old['complete'] and old['accepted_packet_files']==170
for ent in old['files']:
 path=ROOT/ent['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==ent['sha256'],path
cusp=json.loads((P/'all_cusp_FI_sage_certificate.json').read_text())
assert cusp['complete'] and cusp['all_FI_ranks_equal10'] and cusp['all_normalized_configurations']==2383311
assert cusp['all_quadratic_configurations']==1287
assert sum(r['configurations']for r in cusp['by_maximum_rational_multiplicity'])+1287==2383311
coverage=json.loads((P/'source_coverage_verified.json').read_text())
assert coverage['complete'] and coverage['all_points_enumerated'] and coverage['independent_dual_j_by_3_then_19_isogenies']
assert coverage['all_physical_subgroup_orbits']==12 and coverage['required_test_indices']==[2,3,4,5,7,8]
jets=json.loads((P/'good_jets_sage_certificate_v3.json').read_text())
assert jets['complete'] and jets['independent_genuine_R1_R2_arcs']==56 and jets['geometric_cyclic5_kernels']==48
assert jets['all_division_polynomial_norms_exact'] and jets['all_actual_P_5P_25P_images_verified']
assert jets['all_unit_recovery_ranks84'] and jets['all_canonical_tangents_nonzero']
rr=json.loads((P/'ring_RR_all_sage_certificate.json').read_text())
assert rr['complete'] and rr['necessary_divisor_source_coefficients']==[1,4,-5]
assert rr['geometric_push_annihilator']=='T5-1-5sigma' and rr['all_q_product_identities_exact'] and rr['global_R1_and_R2_unit_identities']
assert len(rr['nonprincipal_divisors'])==8
for row in rr['nonprincipal_divisors']:
 assert row['positive_degree']==row['negative_degree']==row['positive_vanishing_rank']==row['negative_vanishing_rank']==110
 assert row['multiplication_by_F_rank']==row['final_criterion_rank']==252
oldrr=json.loads((P/'proof_dependencies/old_point_RR_certificate.json').read_text())
assert oldrr['complete'] and [r['positive_degree']for r in oldrr['nonprincipal_divisors']]==[42,42,88]
assert all(r['final_criterion_rank']==252 for r in oldrr['nonprincipal_divisors'])
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']').replace('true','True').replace('false','False'))
filtered=lit((P/'noncuspidal_filter_v2_data.txt').read_text())
assert len(filtered[0])==48 and len(filtered[1])==6 and len({tuple(r[3])for r in filtered[1]})==3
life=json.loads((P/'lifecycle_final.json').read_text());assert life['complete'] and life['all_tracked_jobs_exited'] and not life['active_case_workers']
assert len(life['jobs'])==23
failed=[]
for row in life['jobs']:
 status=json.loads((ROOT/row['status_path']).read_text());assert status['state']=='EXITED'
 assert hashlib.sha256((ROOT/'work/deg10'/status['script']).read_bytes()).hexdigest()==status['script_sha256']
 if status['actual_exit_code']!=0:failed.append(status['purpose'])
assert failed==['cont57_verify_jets'] # preserved failed verifier; superseded by the successful v3
assert 'C57_DECIC_SOURCE_COVERAGE_INDEPENDENT_SAGE_PASS' in (ROOT/'work/deg10/jobs/cont57_verify_coverage.log').read_text()
assert 'DECIC57_NONGLOBAL_CUSP_9_PLUS_1_EXCLUDED' in (ROOT/'work/deg10/jobs/cont57_quotient10.log').read_text()
assert 'DECIC57_QUOTIENT9_NC_EXPORT_COMPLETE' in (ROOT/'work/deg10/jobs/cont57_quotient10.log').read_text()
print('C57_DECIC_COMPLETE_CERTIFICATE_READONLY_PASS',len(summary['files']),'new evidence hashes;',len(old['files']),'inherited hashes;',len(life['jobs']),'jobs EXITED')
