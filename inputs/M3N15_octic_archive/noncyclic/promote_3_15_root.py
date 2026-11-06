#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
B=Path(__file__).resolve().parents[1];P=B/'noncyclic';m=json.loads((P/'M3N15_D4_necessity_manifest.json').read_text());ev={}
for f,h in m['files'].items():
 p=(P/f).resolve();assert hashlib.sha256(p.read_bytes()).hexdigest()==h,f;ev[str(p.relative_to(B))]=h
assert len(ev)==274
s=(P/'realization_3_15_direct_v4_root.log').read_text();assert 'M3N15_DIRECT_FULL_TORSION_3_15_CERTIFIED' in s and 'M3N15_FINE_MODULI_FIELD_RELATIVE_DEGREE 4' in s and 'error' not in s.lower()
a=json.loads((P/'replay_3_15_direct_root.json').read_text());assert a['full_global_torsion_exact']==[3,15] and a['independent_reduction_subgroup_3_15_order']==45
assert 'M3N15_FULL_SUPPORT_AND_60_PAIR_COVERAGE_REPLAY_DONE' in (P/'M3N15_necessity_root_replay.log').read_text()
for f in ['M3N15_D4_necessity_manifest.json','octic_3_15_classification.tex','M3N15_root_final_review.txt','realization_3_15_direct_v4_root.m','realization_3_15_direct_v4_root.log','realization_3_15_direct_v4_root_arithmetic.txt','replay_3_15_direct_root.py','replay_3_15_direct_root.json','replay_3_15_direct_root.log','M3N15_necessity_root_replay.log','root_315_necessity_replay/replay_M3N15_support_coverage.json','char0_reconstruct_3_15_cubic_v9_root.m','char0_reconstruct_3_15_cubic_v9_root.log','M3N15_char0_RR_encoded_witness_root.m','promote_3_15_root.py']:
 ev['noncyclic/'+f]=hashlib.sha256((P/f).read_bytes()).hexdigest()
f='small_cyclic/audits/audit_3_15_occurrence_logic.txt';ev[f]=hashlib.sha256((B/f).read_bytes()).hexdigest()
manifest={'degree':8,'group':[3,15],'status':'COMPLETE OCCURRENCE AND CLASSIFICATION','new_exclusion':[3,30],'relative_field_polynomial':[3,-6,8,-5,1],'field':'Q(zeta3,w)','j_invariant':'1331/8','full_torsion':[3,15],'closed_moduli_degree_over_Q':8,'method':'Standalone explicit exact example plus full12-class necessity and one level-only Aut48 orbit. Gonality>4 gives uniqueness.','RR_role':'Auxiliary discovery provenance only; not needed for the accepted occurrence or completeness proof.','files':ev}
p=P/'M3N15_complete_root_manifest.json';assert not p.exists();p.write_text(json.dumps(manifest,indent=2)+'\n');ev['noncyclic/'+p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
r=json.loads((B/'certified_results.json').read_text())
for kind,g,why in [('realizations',[3,15],'Exact explicit octic example over Q(zeta3,w), w4-5w3+8w2-6w+3=0; primitive3-Weil pairing, exact point orders3and15, and good unramified prime43 reduction of order45 give full torsion C3xC15. Independent polynomial irreducibility and finite group-law replay passed.'),('exclusions',[3,30],'Every octic curve containing C3xC15 has exactly that torsion: complete two-prime support/class enumeration leaves12 global classes in one level-only Aut48 orbit, and gonality>4 gives unique effective representatives. The explicit example has closed octic moduli field, so no extra field extension can grow its torsion. This excludes C3xC30.')]:
 assert not any(x['group']==g for x in r[kind])
 r[kind].append({'degree':8,'group':g,'audited':True,'proof':'noncyclic/octic_3_15_classification.tex','explanation':why,'evidence':[{'path':f,'sha256':h} for f,h in sorted(ev.items())]})
(B/'certified_results.json').write_text(json.dumps(r,indent=2)+'\n');print('C315_OCCURRENCE_AND_C330_EXCLUSION_ACCEPTED',len(ev))
