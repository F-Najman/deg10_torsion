#!/usr/bin/env python3
"""Bind the completed octic55 proof; no excluded experiment is accepted."""
from pathlib import Path
import ast
import hashlib
import json
import re

P = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


rank_input = json.loads((P / 'C55_inputs/newforms_55.json').read_text())
weighted = 0
forms = []
for m in (1, 5, 11, 55):
    multiplicity = sum((55 // m) % d == 0 for d in range(1, 55 // m + 1))
    for f in rank_input[str(m)]:
        assert f['analytic_rank_proved'] is True and f['analytic_rank'] == 0
        weighted += multiplicity * f['dim']
        forms.append({'label': f['label'], 'dimension': f['dim'],
                      'multiplicity': multiplicity, 'proved_analytic_rank': 0})
assert weighted == 81
rank = {'level': 55, 'genus': 81, 'weighted_newform_dimension': weighted,
        'all_constituents_have_proved_analytic_rank_zero': True,
        'Jacobian_rank_over_Q': 0, 'quotient_rank_over_Q': 0, 'forms': forms}
(P / 'C55_rank_certificate.json').write_text(json.dumps(rank, indent=2) + '\n')

counts = (P / 'C55_inputs/independent_counts.log').read_text()
line = next(s for s in counts.splitlines() if s.startswith('DEGREE8_FAST_COUNT m=1 n=55 '))
nc = ast.literal_eval(re.search(r'closed_place_counts=(\[.*?\])', line)[1])
assert nc == [0, 0, 0, 0, 0, 0, 20, 0]
cusps = next(a for a in json.loads((P / 'cusps_order5_quotients.json').read_text())
             if a['level'] == 55)
assert cusps['H_order'] == 5 and cusps['H_free_on_cusps']
assert cusps['all_rational_F2_cusps_lift_toQ']
assert cusps['source_F2_signature'] == {'1': 20, '4': 5, '10': 2, '20': 1}
assert cusps['quotient_F2_signature'] == {'1': 4, '2': 2, '4': 2}
assert cusps['admissible_F2_signature'] == {'1': 4, '4': 1}
assert cusps['admissible_degree8_divisors'] == 201
assert cusps['Frobenius_exponents_relative_to_diamond2'] == {'55': 0, '11': 1, '5': 2, '1': 3}

good = (P / 'C55_inputs/N55_D8_NC7_q3.log').read_text()
assert 'N55_D8_NC7_PLUS_CUSP_SETUP 20 1 ' in good
assert 'N55_D8_NC7_PLUS_CUSP_CASE 1 principal false ' in good
assert 'N55_D8_NC7_PLUS_CUSP_DONE total=1 remaining=0 complete=true ' in good

for name in ('C55_cover_model_v4.log', 'C55_cartier_reproduce.log'):
    s = (P / name).read_text()
    for bad in ('Runtime error', 'System error', 'SystemError', 'Assertion failed'):
        assert bad not in s, (name, bad)
    assert 'OCTIC55_SIMPLE_WITT_MODEL exact_inverse_maps=true cyclic4=true ' in s
    assert 'OCTIC55_WITT_CHARACTER_VERIFIED samples=2 orientation_units={ 3 } all_Picard_mod4=true ' in s
    assert 'OCTIC55_COVER_GEOMETRY genus17=true split_rational_cusp=true unramified_by_RH=true ' in s
    assert 'RANK0_CARTIER_DONE N=55 D=8 tested=201 total=201 exceptions=0 complete=true ' in s
    assert len(re.findall(r'^C55_RR_KERNEL_DONE ', s, re.M)) == 34
assert digest(P / 'C55_cartier_kernels.m') == digest(P / 'C55_cartier_reproduced_kernels.m')
assert "'exact_kernel_equalities': 34" in (P / 'C55_rr_replay.log').read_text()
assert "'projected_jet_matrices': 34" in (P / 'C55_projection_replay.log').read_text()
replay = (P / 'C55_rank_replay.log').read_text()
r = json.loads(next(s.removeprefix('FI_REPLAY_DONE ') for s in replay.splitlines()
                    if s.startswith('FI_REPLAY_DONE ')))
assert r['complete'] and r['representatives'] == r['covered'] == r['enumerated'] == 201
assert r['exceptions'] == 0 and r['degree'] == 8

files = [
    'C55_exclusion.tex', 'certify55.py', 'C55_rank_certificate.json',
    'C55_builtin_model.m', 'C55_cover_model_v4.m', 'C55_cover_model_v4.log',
    'C55_cartier_continue.m', 'C55_cartier_reproduce.m', 'C55_cartier_reproduce.log',
    'C55_cartier_kernels.m', 'C55_cartier_reproduced_kernels.m',
    'C55_rr_replay.log', 'C55_projection_replay.log', 'C55_rank_replay.log',
    'cusps_order5_quotients.py', 'cusps_order5_quotients.json', 'cusps_order5_quotients.log',
    *['C55_inputs/' + n for n in (
        'newforms_55.json', 'dyadic_criterion.tex', 'independent_counts.m',
        'independent_counts.log', 'N55_D8_NC7_q3.m', 'N55_D8_NC7_q3.log',
        'replay_rr_kernels.py', 'replay_cartier_projection.py', 'replay_fi_jets.py')],
    '../shared/mdmagma/v2/mdmagma.spec', '../shared/fast_hecke.m',
    '../shared/fast_hecke_11.m', '../shared/fast_hecke_cyclic_commonzero.m',
]
packet = {
    'degree': 8, 'excluded_cyclic_order': 55, 'status': 'complete', 'complete': True,
    'Jacobian_rank_overQ': 0, 'quotient': 'X1(55)/diamond16_order5',
    'quotient_genus': 17, 'quotient_rank_overQ': 0,
    'good_reduction': {'noncuspidal_degrees_through8': nc,
                       'degree7_closed_places': 20, 'diamond_orbits': 1,
                       'A3_principal': [False], 'rational_cusp_complement_lifts_Q': True},
    'cuspidal_reduction': {'admissible_cusp_signature': {'1': 4, '4': 1},
                          'divisors': 201, 'single_cusp_kernels': 34,
                          'Cartier_image_dimension': 15, 'all_restriction_ranks': 8,
                          'exceptions': 0},
    'independent_exact_model_reproduction': True,
    'reproduced_kernel_export_identical': True,
    'files': [{'path': name, 'sha256': digest(P / name)} for name in sorted(files)],
}
(P / 'C55_manifest.json').write_text(json.dumps(packet, indent=2) + '\n')
print('OCTIC55_COMPLETE_CERTIFICATE', digest(P / 'C55_manifest.json'))
