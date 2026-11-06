#!/usr/bin/env python3
"""Freeze only completed C85 proof dependencies; preserve previous diagnostics."""
from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;R=C.parents[1];H=C/'composite85';assert not(H/'manifest.json').exists()
fset={R/a['path']for a in json.loads((R/'shared/snapshot_manifest.json').read_text())['files']};fset.add(R/'shared/snapshot_manifest.json')
for stem in['canonical_unit_stage1_v2','elliptic_hecke_jets_stage2','rational_cusp_FI_positive']:
 for ext in['.m','.log','.exit.json']:fset.add(H/(stem+ext))
 assert json.loads((H/(stem+'.exit.json')).read_text())['native_exit_code']==0
for name in['canonical_unit_stage1_data.m','elliptic_hecke_jets_stage2_data.m','FI_positive_rows.txt','FI_positive_data.m','stage1_independent_verified.json','stage2_independent_verified.json','canonical_ring_R1_tensor.json','canonical_R1_sieve_verified.json','canonical_R1_matrices.json.gz','FI_positive_independent_verified.json','census_independent_verified.json','cusp_action_verified.json','cusp_pole_degrees_verified.json','hecke_transport_independent_verified.json','exclusion.tex','reviewed_exclusion.tex','provenance.tex','unit_arithmetic.py','verify_stage1_readonly.py','verify_stage2_readonly.py','build_R1_product_tensor.py','canonical_R1_sieve.py','verify_R1_ranks_m4ri.py','verify_FI_positive_readonly.py','verify_census_readonly.py','verify_cusp_action_readonly.py','verify_cusp_pole_degrees.py','verify_hecke_transport_readonly.py']:fset.add(H/name)
for f in range(6,10):
 for name in['records.m','status.json','run.log']:fset.add(R/'parallel_expansion_20261026'/('N85_f'+str(f))/name)
for name in['odd_census.m','character_probes/probe_winding.m','character_probes/N85_winding_data.m']:fset.add(R/'parallel_expansion_20261026'/name)
for name in['winding_N85.log','winding_N85_status.json']:fset.add(R/'cyclic/continuation_20261026'/name)
for name in['certify_85_readonly.py','finish85_packet.py','run_magma_mutable_cpu.py']:fset.add(C/name)
for name in['DvH_1307.5719v3.pdf','DvH_1307.5719v3.txt']:fset.add(C/'P83/pole_bound_supplement'/name)
files=[]
for p in sorted(fset):assert p.is_file(),p;files.append({'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
m={'degree':9,'group':[1,85],'status':'COMPLETE EXCLUSION','claim':'An elliptic curve over a degree9 number field has no point of order85','path_base':'work/degree9','proof':'prime_orders/continued_20261026/composite85/reviewed_exclusion.tex','independent_checker':'prime_orders/continued_20261026/certify_85_readonly.py','native_exits_selected_zero':8,'binary_cases':320,'binary_shape':[225,519],'binary_rank':225,'rational_cusp_rows':133,'rational_cusp_rank':122,'rational_cusp_kernel_dimension':11,'minimum_cusp_prefix_cost':32,'normalized_cusp_configurations':19070333,'files':files};p=H/'manifest.json';p.write_text(json.dumps(m,indent=2)+'\n');print('C85_COMPLETE_PACKET_FROZEN files',len(files),'SHA',hashlib.sha256(p.read_bytes()).hexdigest())
