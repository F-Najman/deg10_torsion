#!/usr/bin/env python3
"""One-shot freeze after mathematical review; never overwrites a manifest."""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent;H=P.parents[2];out=P/'manifest44.json';assert not out.exists()
B=P.parent;old=B/'inherited_lower/cyclic_sieves/N44_manifest.json';O=json.loads(old.read_text());paths={old}
for e in O['evidence']:
 q=old.parent.parent/e['path'];assert hashlib.sha256(q.read_bytes()).hexdigest()==e['sha256'];paths.add(q)
files=['build_canonical_v2.m','canonical_data.m','build_tensor22.m','multiplication_22_data.m','export_finite_jets_v2.m','finite_all_v2_driver.m','finite_jets_all_v2.raw.jsonseq','finite_functions_v2.m','native_control.m','native_control_inputs.m','prepare_sieve_inputs.py','sieve_all_inputs.m','sieve_all_inputs.normalized.jsonl','matrix22_sieve_template.m','matrix22_all.m','matrix22_all.raw.jsonseq','replay_matrices22_sage.py','matrix22_all_independent.json','matrix22_all_independent.log','audit_cusps44.py','cusp_audit44.json','audit_frobenius44.py','frobenius_audit44.json','exclusion44.tex','certify44.py','freeze44.py']
paths.update(P/s for s in files)
paths.update(B/s for s in ['N44_source_counts.m','N44_frobenius3_source.json','N44_lower_D9_transfer.json','N44_partial_reduction_D9.tex','C44_pair_af_place_linear.m','run_native_job.py'])
paths.update([H/'shared/fast_hecke.m',H/'shared/fast_hecke.sig',H/'shared/snapshot_manifest.json',H/'noncyclic/C230/refs/DEHMZ_2007.13929.tex',H/'prime_orders/continued_20261026/P89/quotient27/canonical_ring_criterion.tex'])
jobs=[]
for tag,marker in [('N44_source_counts','NONIC44_SOURCE_FROBENIUS_COUNTS_PASS'),('C44_canonical_build_v2','NONIC44_CANONICAL_MODEL_CERTIFIED genus36=true degree70=true quadrics561=true'),('C44_finite_jets_all_v2','NONIC44_FINITE_JETS_COMPLETE cases=226 source_distinct_keys=2260 all226=true'),('C44_native_control','NONIC44_NATIVE_RR_CONTROL_DONE dimension=0 no_class_cache=true'),('C44_tensor22','NONIC44_RING22_COMPLETE R2=105 R4=245 products11025=true'),('C44_matrix22_all','NONIC44_CANONICAL22_COMPLETE cases=226 survivors=[] all_rank105=true')]:
 st=B/(tag+'_status.json');log=B/(tag+'.log');s=json.loads(st.read_text());assert s['native_exit_code']==0 and marker in log.read_text();paths.update([st,log]);jobs.append({'status':st.relative_to(H).as_posix(),'log':log.relative_to(H).as_posix(),'completion_marker':marker})
meta={'degree':9,'excluded_group':[1,44],'complete':True,'method':'rank-zero, unramified cuspidal lift, lower138+degree4complement audit, R2xR2 canonical multiplication for226 pure9cases','source_genus':36,'prime':3,'source_noncuspidal_counts':[0,0,0,0,60,160,240,920,2260],'new_cases':226,'matrix_size':[105,194],'all_rank':105,'historical_evidence_used_only_as_finite_arithmetic':True,'no_unfinished_native_worker_required':True,'new_native_jobs':jobs,'files':[{'path':q.relative_to(H).as_posix(),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in sorted(paths)]}
out.write_text(json.dumps(meta,indent=2)+'\n');print('C44_FROZEN files='+str(len(meta['files']))+' sha256='+hashlib.sha256(out.read_bytes()).hexdigest())
