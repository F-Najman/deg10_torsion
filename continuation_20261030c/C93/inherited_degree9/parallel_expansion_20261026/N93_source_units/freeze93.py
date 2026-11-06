from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent;D=H.parents[1];assert not(H/'manifest93.json').exists()
files={p for p in H.rglob('*')if p.is_file()and '__pycache__'not in p.parts and p.suffix not in['.pyc','.tmp']and not p.name.startswith('core')}
base=json.loads((D/'prime_orders/continued_20261026/composite85/manifest.json').read_text())
files.update(D/e['path']for e in base['files']if e['path'].startswith('shared/'))
for n in[6,7,8,9]:files.update(p for p in(D/f'parallel_expansion_20261026/N93_f{n}').glob('*')if p.is_file())
extras=['parallel_expansion_20261026/odd_census.m','parallel_expansion_20261026/run_pool.py','parallel_expansion_20261026/run_source_units_v2.py','parallel_expansion_20261026/run_source_products.py','parallel_expansion_20261026/pole_audit_squarefree_v3.py','parallel_expansion_20261026/character_probes/probe_winding.m','parallel_expansion_20261026/character_probes/N93_winding_data.m','cyclic/continuation_20261026/winding_N93.log','cyclic/continuation_20261026/winding_N93_status.json','cyclic/continuation_20261026/run_native_job.py','noncyclic/C230/refs/DEHMZ_2007.13929.tex','noncyclic/C230/refs/NV_2602.03513v4.html']
files.update(D/e for e in extras)
jobs=[(H/'status_v2.json','actual_exit_code'),(H/'products_status.json','actual_exit_code'),(H/'elliptic_hecke5_v2.exit.json','actual_exit_code'),(D/'cyclic/continuation_20261026/winding_N93_status.json','native_exit_code')]
jobs +=[(D/f'parallel_expansion_20261026/N93_f{n}/status.json','actual_exit_code')for n in[6,7,8,9]]
jobs +=[(H/e,'actual_exit_code')for e in['tensor.exit.json','sieve93.exit.json','stage1_independent_v3.exit.json','stage2_independent_v2.exit.json','ranks_m4ri.exit.json']]
for p,k in jobs:assert json.loads(p.read_text())[k]==0
out={'degree':9,'group':[1,93],'complete_exclusion':True,'proof':str((H/'exclusion93.tex').relative_to(D)),'independent_checker':str((H/'certify93.py').relative_to(D)),'successful_jobs':[{'path':str(p.relative_to(D)),'exit_key':k}for p,k in jobs],'diagnostics_not_used':'Earlier failed setup/precision/virtual-address-limit scripts and their logs are retained; only selected successful records support the proof.','files':[{'path':str(p.relative_to(D)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in sorted(files)]}
(H/'manifest93.json').write_text(json.dumps(out,indent=2)+'\n');print('FROZEN_C93',len(files),hashlib.sha256((H/'manifest93.json').read_bytes()).hexdigest())
