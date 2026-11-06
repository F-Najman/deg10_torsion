#!/usr/bin/env python3
"""Freeze the C91 exclusion after every exact read-only component passes."""
from pathlib import Path
import json,hashlib
C=Path(__file__).resolve().parent;R=C.parents[1];H=C/'composite91';U=R/'parallel_expansion_20261026/N91_source_units';assert not(H/'manifest.json').exists()
r=json.loads((H/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==288 and r['source_diamond_orbits']==4 and not r['unresolved'];assert{(v['positive_degree'],v['columns'],v['rank'])for v in r['results']}=={(63,590,265)}
fset={R/a['path']for a in json.loads((R/'shared/snapshot_manifest.json').read_text())['files']};fset.add(R/'shared/snapshot_manifest.json')
for p in H.iterdir():
 if p.is_file()and p.suffix in['.py','.json','.m','.tex','.gz']:fset.add(p)
for name in['geometry_v2.m','status_v2.json','run_v2.log','canonical_unit_data.m']:fset.add(U/name)
for i in range(1,5):
 for name in['run.m','run.log','run.exit.json','guard.py','elliptic_hecke3_data.m']:fset.add(U/f'hecke3_source{i}'/name)
for f in range(6,10):
 for name in['records.m','status.json','run.log']:fset.add(R/'parallel_expansion_20261026'/f'N91_f{f}'/name)
for name in['odd_census.m','character_probes/probe_winding.m','character_probes/N91_winding_data.m']:fset.add(R/'parallel_expansion_20261026'/name)
for name in['winding_N91.log','winding_N91_status.json']:fset.add(R/'cyclic/continuation_20261026'/name)
for name in['certify_91_readonly.py','finish91_packet.py']:fset.add(C/name)
for name in['DvH_1307.5719v3.pdf','DvH_1307.5719v3.txt']:fset.add(C/'P83/pole_bound_supplement'/name)
files=[]
for p in sorted(fset):assert p.is_file(),p;files.append({'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
m={'degree':9,'group':[1,91],'status':'COMPLETE EXCLUSION','claim':'An elliptic curve over a degree9 number field has no point of order91','path_base':'work/degree9','proof':'prime_orders/continued_20261026/composite91/reviewed_exclusion.tex','independent_checker':'prime_orders/continued_20261026/certify_91_readonly.py','native_exits_selected_zero':10,'binary_cases':288,'binary_shape':[265,590],'binary_rank':265,'special_degree8_closed_points':3,'special_involution_scalar':27,'files':files};p=H/'manifest.json';p.write_text(json.dumps(m,indent=2)+'\n');print('C91_COMPLETE_PACKET_FROZEN files',len(files),'SHA',hashlib.sha256(p.read_bytes()).hexdigest())
