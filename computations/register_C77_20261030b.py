#!/usr/bin/env python3
"""Register the frozen C77 proof and local immutable inherited snapshots."""
from pathlib import Path
import hashlib,json,shutil
H=Path(__file__).resolve().parents[1]; ROOT=H.parent.parent; C=H/'continuation_20261030b/C77'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=C/'certificate_manifest.json'; M=json.loads(manifest.read_text())
assert M['status']=='PASS' and len(M['files'])==232
for stem,marker in [('contb77_parent_native','C77_PARENT_NATIVE_COMPLETE 24 36'),('contb77_parent_moduli_binding','C77_PARENT_T3_MODULI_BINDING_COMPLETE 24 16')]:
 s=json.loads((H/f'jobs/{stem}.status.json').read_text());assert s['state']=='EXITED' and s['actual_exit_code']==0
 assert marker in (H/f'jobs/{stem}.log').read_text()
evidence={}; inherited=[]
def add(p,expected=None):
 digest=sha(p)
 if expected is not None:assert digest==expected,p
 key=str(p.relative_to(H));assert key not in evidence or evidence[key]==digest
 evidence[key]=digest
for e in M['files']:
 p=ROOT/e['path'];assert p.stat().st_size==e['bytes'] and sha(p)==e['sha256'],p
 if not p.is_relative_to(H):
  rel=p.relative_to(ROOT/'work/degree9');dest=C/'inherited_degree9'/rel
  dest.parent.mkdir(parents=True,exist_ok=True)
  if not dest.exists():shutil.copy2(p,dest)
  assert sha(dest)==e['sha256'];inherited.append({'source':e['path'],'snapshot':str(dest.relative_to(ROOT)),'sha256':e['sha256']});p=dest
 add(p,e['sha256'])
snap=C/'inherited_snapshots.json'
snap.write_text(json.dumps({'scope':'Immutable local copies of all external files named in the frozen C77 manifest; historical paths and originals are unchanged. Only degree-independent geometry and the explicitly accepted proof inputs are transferred; no degree9 exclusion is assumed.','files':inherited},indent=2)+'\n')
for p in [manifest,snap,Path(__file__).resolve()]:add(p)
ev=[{'path':p,'sha256':s} for p,s in sorted(evidence.items())]
p=H/'certified_results.json'; L=json.loads(p.read_text());assert not any(r['group']==[1,77] for r in L['exclusions'])
method='Finite rational T2-image on the genus37 quotient; full-Jacobian dyadic formal immersion for all4992 cusp patterns; all24 possible quadratic T3-lift allocations have nonprincipal A3T2 divisors, with independent native RR and marked-isogeny/norm binding'
L['exclusions'].append({'group':[1,77],'verified':True,'method':method,'proof':'continuation_20261030b/C77/resolution.tex','evidence':ev})
L['auxiliary_results'].append({'name':'decic_C77_finite_T2_image_genus37','verified':True,'statement':method+'. The three source diamond orbits exhaust all90 noncuspidal decic places. Positive rank of the full source Jacobian is retained; no rank-zero assumption is used. All failed or retired trials are provenance, not certificates.','evidence':ev})
p.write_text(json.dumps(L,indent=2)+'\n')
print(json.dumps({'registered':'C77 exclusion','evidence':len(ev),'external_local_snapshots':len(inherited),'frozen_manifest_sha256':sha(manifest)},indent=2))
