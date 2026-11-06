#!/usr/bin/env python3
"""Read-only validation of the complete degree-ten C75 proof packet."""
from pathlib import Path
import ast,hashlib,json,re
P=Path(__file__).resolve().parent;ROOT=P.parents[3]
a=json.loads((P/'final_certificate_summary.json').read_text());assert a['complete_exclusion']
assert a['parent_mathematical_review']=='PASS'
for ent in a['files']:
 p=ROOT/ent['path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==ent['sha256'],p
old=json.loads((P/'inherited_evidence_verified.json').read_text());assert old['complete'] and old['accepted_packet_files']==588
for ent in old['files']:
 p=ROOT/ent['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==ent['sha256'],p
c=json.loads((P/'cusp_degree10_verified.json').read_text());assert c['complete'] and c['quotient_all_cusps_globally_intact']
assert c['admissible_degree10_cusp_images']==871630 and c['retained_common_basis_tests']==871520 and c['new_high_cases']==110
assert c['rank_at_12_width75_cusp']==12 and c['rank_at_10_width25_cusp']==10 and c['all_nine_rank_at_12_plus_3_width75_cusps']==[15]*9
assert c['all_degree10_cusp_images_rank10'] and c['no_cross_export_basis_alignment_used']
f=json.loads((P/'finite_degree10_verified.json').read_text());assert f['complete'] and f['actual_moduli_records_checked']==84 and f['all_geometric7_kernels_independently_verified']==48
assert f['source_NC_counts_degrees1to10']==[0,0,0,0,0,0,20,0,0,120]
assert f['all_actual_level_images_independently_verified'] and f['all_norm_multiplicities_checked']
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
packets=[lit((P/f'proof_packet_case{i}.m').read_text())for i in range(1,7)]
assert all(p[1:4]==packets[0][1:4]for p in packets)
for i,p in enumerate(packets,1):
 assert p[0]==i
 assert sum(pr[1]*n for pr,n in p[4])==10 and sum(pr[1]*n for pr,n in p[5])==10
 assert sum(pr[1]*n for pr,n in p[6])==80 and sum(pr[1]*n for pr,n in p[7])==0
 assert sorted(p[8])==[[1,10,1],[1,10,1],[3,30,1],[3,30,1]]
life=json.loads((P/'lifecycle_final.json').read_text());assert life['complete'] and life['all_tracked_jobs_exited'] and not life['active_case_workers']
assert len(life['jobs'])==14
failures=[]
for row in life['jobs']:
 st=json.loads((ROOT/row['status_path']).read_text());assert st['state']=='EXITED'
 assert hashlib.sha256((ROOT/'work/deg10'/st['script']).read_bytes()).hexdigest()==st['script_sha256']
 if st['actual_exit_code']!=0:failures.append(st['purpose'])
assert sorted(failures)==['cont75_nc10_job1','cont75_nc10_job2','cont75_nc10_job3','cont75_verify_cusps','cont75_verify_finite']
for job in range(1,4):
 text=(ROOT/f'work/deg10/jobs/cont75_nc10_job{job}_v2.log').read_text();assert 'DECIC75_NC10_JOB_COMPLETE' in text and text.count('DECIC75_A7_RR_COMPLETE')==2 and text.count('DECIC75_NATIVE_GEOMETRIC_HECKE_CONTROL_PASS')==2
assert 'C75_PARENT_NATIVE_COMPLETE 6 36' in (ROOT/'work/deg10/jobs/cont75_parent_native.log').read_text()
assert 'DECIC75_HIGH_CUSP_COMPLETE' in (ROOT/'work/deg10/jobs/cont75_high_cusps.log').read_text()
assert 'C75_RANK_ZERO_WINDING_COMPLETE 145' in (P/'proof_dependencies/accepted_rank0.log').read_text()
assert 'C75_GOOD7_CERTIFICATE_DONE complete=true' in (P/'proof_dependencies/accepted_good7_native.log').read_text()
print('C75_DECIC_COMPLETE_CERTIFICATE_READONLY_PASS',len(a['files']),'new hashes;',len(old['files']),'inherited hashes;',len(life['jobs']),'jobs EXITED')
