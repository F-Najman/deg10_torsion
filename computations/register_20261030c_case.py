#!/usr/bin/env python3
"""Register a frozen, mathematically reviewed decic exclusion packet."""
from pathlib import Path
import json,hashlib,sys
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent
methods={91:'Relative full-Jacobian dyadic formal fibres on genus23 X1(7) x X0(13); all120 cusp ranks and12 CM ranks equal the expected codimension of the embedded pullback family, forcing a fixed cusp; corrected inherited ordinary Hecke obstructions independently rechecked',99:'Finite diamond image and gonality exclude pure cusps; three A5 obstructions on rank-zero X1(33), and the final alpha2 A5 nonprincipality test on genus29 X1(9) x X0(11), with fresh native RR and exact marked-isogeny norm bindings',65:'Finite-image Hecke products on the full Jacobian; all7170048 global cusp restrictions and the split-cusp finite class filter; complete old degree6-9 tests and240 new marked canonical-multiplication obstructions with independent controls',105:'Finite-image free order4 diamond and gonality exclude pure cusps; six exact A3 nonprincipality tests on rank-zero genus25 X1(35), independently replayed and bound to genuine marked Velu isogenies',85:'Finite-image integral Hecke polynomials with the full-Jacobian local multiplicative-level operator T2-1; exact rational and nonrational cusp restrictions; all new and inherited necessary canonical-multiplication tests independently verified',93:'Finite diamond images; full-Jacobian dyadic formal immersion for rational and quintic cusp supports, with order2 descent to quintic31 excluding the other boundary; independently verified marked noncuspidal Hecke obstructions'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for n in map(int,sys.argv[1:]):
 assert n in methods
 C=H/f'continuation_20261030c/C{n}';mf=C/'certificate_manifest.json'
 if not mf.exists():mf=C/'final_certificate_summary.json'
 d=json.loads(mf.read_text());assert d.get('status')=='PASS' or (d.get('complete_exclusion') is True and d.get('parent_mathematical_review')=='PASS' and d.get('independent_arithmetic')=='PASS')
 ev={}
 for e in d['files']:
  p=ROOT/e['path'];assert p.is_relative_to(H),p;assert sha(p)==e['sha256'],p
  if 'bytes'in e:assert p.stat().st_size==e['bytes']
  ev[str(p.relative_to(H))]=e['sha256']
 review=C/'registration_review.json';review.write_text(json.dumps({'status':'PASS','group':[1,n],'degree':10,'parent_mathematical_review':'PASS','frozen_manifest_sha256':sha(mf),'method':methods[n],'independent_arithmetic_checks':'PASS','classification_complete':False},indent=2)+'\n')
 for p in[mf,review]:ev[str(p.relative_to(H))]=sha(p)
 evidence=[{'path':p,'sha256':s}for p,s in sorted(ev.items())]
 p=H/'certified_results.json';l=json.loads(p.read_text());assert not any(r['group']==[1,n]for r in l['exclusions'])
 l['exclusions'].append({'group':[1,n],'verified':True,'method':methods[n],'proof':str((C/'resolution.tex').relative_to(H)),'evidence':evidence})
 l['auxiliary_results'].append({'name':f'decic_C{n}_continuation_20261030c','verified':True,'statement':methods[n]+'. Parent mathematical review and independent exact controls passed. Failed or terminated diagnostics are not exclusion certificates.','evidence':evidence})
 p.write_text(json.dumps(l,indent=2)+'\n');print('REGISTERED',n,len(evidence))
