#!/usr/bin/env python3
"""Strict coverage audit and fixed-scope freeze for the completed C39 proof."""
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
SMALL=HERE.parent
DEGREE=SMALL.parent
log=(HERE/'hecke_good_39.log').read_text()
assert not re.search(r'Runtime error|User error|Assertion failed',log)
partitions=re.findall(r'C39_D8_GOOD_PARTITION degree=(\d+) places=(\d+) representatives=(\d+) exact_cover=true',log)
assert partitions==[('6','12','1'),('8','49','5')]
witnesses=re.findall(r'C39_D8_GOOD_WITNESS degree=(\d+) index=(\d+)/(\d+) q=(\d+) principal=(true|false)',log)
assert len(witnesses)==5 and all(q=='5' and p=='false' for _,_,_,q,p in witnesses)
known=re.findall(r'C39_D8_GOOD_CASE degree=8 index=(\d+)/(\d+) known_exact_reduction=true',log)
assert len(known)==1 and known[0][1]=='5'
assert {(d,i) for d,i,_,_,_ in witnesses}|{('8',known[0][0])}=={('6','1')}|{('8',str(i)) for i in range(1,6)}
assert 'full_diamond_orbit=1' in log
assert 'representatives=[ 1, 5 ] nonknown_tested=5 known_representatives=1 unknown=0 complete=true' in log
dyadic=json.loads((HERE/'backup/C4_independent_replay.json').read_text())
assert dyadic['source_torsion_specialization_injective']
assert dyadic['target_2primary_reduction_image_order']==32
positive=json.loads((SMALL/'realizations_manifest.json').read_text())
assert positive['complete'] and positive['examples']==10
for name,digest in positive['files'].items():
    assert sha256((SMALL/name).read_bytes()).hexdigest()==digest
report={'degree':8,'level':39,'complete':True,'source_genus':33,
        'rational_gonality':14,'source_rational_2primary_invariants':[2,16],
        'torsion_specialization_at2_injective':True,
        'noncuspidal_counts':[0,0,0,0,0,12,0,49],
        'sextic_orbits':1,'octic_orbits':5,'nonprincipality_witnesses':5,
        'known_octic_reduction_orbits':1,'known_closed_octic_points':1,
        'all_known_full_torsion':[39],'unknown':0}
(HERE/'quotient_coverage_replay.json').write_text(json.dumps(report,indent=2)+'\n')
names=['hecke_good_39.m','hecke_good_39.log','reproduce_good_39.m','reproduce_dyadic_39.m',
       'octic_completeness_39.tex','CONTROL_STATUS.txt','freeze_quotient_certificate.py',
       'quotient_coverage_replay.json','inputs/DEHMZ.tex','inputs/gonality',
       'deps/fast_diamond_hash_m2.m','deps/fast_diamond_hash_cyclic_bc.m',
       'backup/C4_dyadic_certificate_v2.m','backup/C4_dyadic_certificate_v2.log',
       'backup/C4_cusp_arithmetic.py','backup/replay_C4_arithmetic.py',
       'backup/C4_independent_replay.log','backup/C4_independent_replay.json']
paths=[HERE/name for name in names]+[SMALL/'realizations_manifest.json']
paths += [SMALL/name for name in positive['files']]
paths += sorted(p for p in (DEGREE/'shared/mdmagma').rglob('*') if p.is_file() and p.suffix != '.sig')
paths += [DEGREE/'shared'/name for name in ['fast_hecke.m','fast_hecke_11.m','fast_hecke_cyclic_commonzero.m']]
assert len(paths)==len(set(paths))
files={str(p.relative_to(DEGREE)):sha256(p.read_bytes()).hexdigest() for p in paths}
manifest={'degree':8,'level':39,'status':'complete','complete':True,
          'method':'odd-degree quotient proves dyadic torsion injectivity; fresh source Hecke and gonality',
          'known_full_torsion':[39],'excluded_minimal_enlargements':[[78],[117],[3,39]],
          'source_rational_Jacobian_invariants':[7,31122,3236688],
          'source_rational_2primary_invariants':[2,16],
          'quotient_genus':9,'quotient_degree':3,'target_F2_Jacobian_invariants':[6,336],
          'target_rational_cusp_2primary_image':32,'good_nonprincipality_witnesses':5,
          'allbad_cusps_lift_Q':True,'rational_gonality':14,
          'optional_direct_class_group_control_used':False,
          'Magma_version':'2.30-1',
          'reproduction_commands':[
              'magma -b work/degree8/small_cyclic/C39/reproduce_dyadic_39.m',
              'magma -b work/degree8/small_cyclic/C39/reproduce_good_39.m'],
          'files':files}
(HERE/'quotient_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('C39_QUOTIENT_CERTIFICATE_FROZEN complete=true good_witnesses=5 known_octic=1 files='+str(len(files))+' sha256='+sha256((HERE/'quotient_manifest.json').read_bytes()).hexdigest())
