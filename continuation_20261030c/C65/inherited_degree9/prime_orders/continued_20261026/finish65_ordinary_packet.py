#!/usr/bin/env python3
"""Freeze only the completed ordinary7/8/9 C65 branches, with explicit scope."""
from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;R=C.parents[1];H=C/'composite65';U=R/'parallel_expansion_20261026/N65_source_units';assert not(H/'ordinary_manifest.json').exists();fset={R/a['path']for a in json.loads((R/'shared/snapshot_manifest.json').read_text())['files']};fset.add(R/'shared/snapshot_manifest.json')
for p in H.iterdir():
 if p.is_file()and p.suffix in['.py','.json','.m','.tex','.gz']:fset.add(p)
for f,k in[(7,2),(8,2),(9,3)]:
 for i in range(1,k+1):
  for name in['run.m','run.log','run.exit.json','elliptic_hecke3_data.m']:fset.add(H/f'hecke3_f{f}_source{i}'/name)
for name in['geometry_v2.m','status_v2.json','run_v2.log','canonical_unit_data.m']:fset.add(U/name)
for f in range(6,10):
 for name in['records.m','status.json','run.log']:fset.add(R/'parallel_expansion_20261026'/f'N65_f{f}'/name)
for name in['odd_census.m','character_probes/probe_winding.m','character_probes/N65_winding_data.m']:fset.add(R/'parallel_expansion_20261026'/name)
for name in['winding_N65.log','winding_N65_status.json']:fset.add(R/'cyclic/continuation_20261026'/name)
for name in['certify_65_ordinary_readonly.py','finish65_ordinary_packet.py','run_magma_mutable_cpu.py']:fset.add(C/name)
for name in['DvH_1307.5719v3.pdf','DvH_1307.5719v3.txt']:fset.add(C/'P83/pole_bound_supplement'/name)
files=[]
for p in sorted(fset):assert p.is_file(),p;files.append({'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
m={'degree':9,'group':[1,65],'status':'COMPLETE EXCLUSION OF ORDINARY DEGREE7_8_9 REDUCTION BRANCHES','claim':'A noncuspidal effective rational degree9 divisor on X1(65) cannot reduce to an ordinary degree7/8/9 point plus its cusp complement','not_claimed':'Full C65 exclusion; supersingular degree6 and entirely cuspidal reduction are separate','path_base':'work/degree9','proof':'prime_orders/continued_20261026/composite65/ordinary_branches.tex','independent_checker':'prime_orders/continued_20261026/certify_65_ordinary_readonly.py','native_exits_selected_zero':13,'binary_cases':336,'binary_shapes_counts':[[121,288,96],[121,295,96],[121,302,144]],'binary_rank':121,'files':files};p=H/'ordinary_manifest.json';p.write_text(json.dumps(m,indent=2)+'\n');print('C65_ORDINARY_PACKET_FROZEN files',len(files),'SHA',hashlib.sha256(p.read_bytes()).hexdigest())
