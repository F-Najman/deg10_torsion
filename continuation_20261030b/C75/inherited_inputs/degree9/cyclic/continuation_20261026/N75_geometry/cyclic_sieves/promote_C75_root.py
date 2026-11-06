"""Freeze and accept the complete octic C75 exclusion."""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent;B=P.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files={}
for row in json.loads((B/'shared/degree8_snapshot_manifest.json').read_text())['files']:
 p=B/'shared'/row['path'];assert sha(p)==row['sha256'];files[str(p.relative_to(B))]=sha(p)
names=['C75_exclusion_root.tex','C75_root_review.txt','promote_C75_root.py','certify_C75_root.py',
'C75_rank_zero_root.m','C75_rank_zero_root.log','C75_generic_canonical3_root.m','C75_generic_canonical3_root.log',
'C75_source_frob_root.m','C75_source_frob_root.log','C75_source_frob_root_arithmetic.txt',
'C75_good7_certificate_root_v2.m','C75_good7_certificate_root_v2.log','C75_good7_moduli_records_root_v2.txt',
'C75_binary_fields_root.py','replay_C75_good7_root.py','C75_good7_independent_root.json','C75_good7_independent_root.log',
'N75_gon5_allbad_orbits_root.m','N75_gon5_allbad_orbits_root.log',
'N75_allbad_orbit_matrices_root.m','N75_allbad_orbit_kernels_root.m','N75_allbad_orbit_result_root.txt',
'replay_C75_allbad_root.py','C75_allbad_independent_root.json','C75_allbad_independent_root.log']
for n in names:
 p=P/n;assert p.is_file(),p;files[str(p.relative_to(B))]=sha(p)
for n in ['noncyclic/deps/pair_af_place_linear_C48.m',
'large_levels/audits/audit75_tate_parities.py','large_levels/audits/audit75_tate_parities.json',
'large_levels/audits/audit75_tate_parities.log','large_levels/audits/audit75_theory.txt']:
 p=B/n;assert p.is_file();files[n]=sha(p)
assert 'C75_GOOD7_CERTIFICATE_DONE complete=true' in (P/'C75_good7_certificate_root_v2.log').read_text()
assert 'N75_ALLBAD_FI_DONE 131790 0 complete true' in (P/'N75_gon5_allbad_orbits_root.log').read_text()
for n,marker in [('C75_good7_independent_root.log','C75_GOOD7_INDEPENDENT_REPLAY_PASS'),('C75_allbad_independent_root.log','C75_ALLBAD_INDEPENDENT_REPLAY_PASS')]:
 assert marker in (P/n).read_text()
mp=P/'C75_root_manifest.json';assert not mp.exists()
mp.write_text(json.dumps({'degree':8,'group':[1,75],'complete_exclusion':True,'files':files},indent=2)+'\n')
files[str(mp.relative_to(B))]=sha(mp)
lp=B/'certified_results.json';r=json.loads(lp.read_text());assert not any(x['group']==[1,75] for x in r['exclusions'])
r['exclusions'].append({'degree':8,'group':[1,75],'audited':True,'proof':'cyclic_sieves/C75_exclusion_root.tex',
'explanation':'Rankzero J1(75), a good genus73 point25/cyclic3 quotient and the complete Frobenius/cusp counts reduce the good branch to one degree7 diamond orbit, excluded by actual A7 with full-isogeny reference. The allbad branch has131790 admissible cusp divisors; all have full restriction rank8 on a71-dimensional Cartier image. Independent elliptic enumeration,128-block binary replay and separate theory audit passed.',
'evidence':[{'path':n,'sha256':h} for n,h in sorted(files.items())]})
lp.write_text(json.dumps(r,indent=2)+'\n');print('C75_COMPLETE_EXCLUSION_ACCEPTED',len(files))
