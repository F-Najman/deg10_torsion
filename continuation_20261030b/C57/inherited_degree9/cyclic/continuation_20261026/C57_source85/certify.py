"""Read-only complete degree-nine C57 exclusion certificate.

--deep reruns the independent source census, finite quotient, cusp arcs,
formal-immersion enumeration, exact isogenies and genuine jets, canonical
multiplication identities, and the correctly normalized final RR tests.
No Magma process is launched or resumed.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
P=Path(__file__).resolve().parent;D=P.parents[2];Q=P.parent/'C57_quotient9'
def js(f):return json.loads((P/f).read_text())
m=js('manifest.json');assert m['degree']==9 and m['complete_exclusion'] and m['excluded_groups']==[[1,57]]
assert m['geometric_push_annihilator']=='T5-sigma-5' and m['necessary_divisor_source_coefficients']==[5,-4,-1]
for a in m['files']:
 p=D/a['path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==a['sha256'],p
checks_native=[
 (P,'completion.json','canonical_units.log','NONIC57_SOURCE85_CANONICAL_UNIT_STAGE1_COMPLETE'),
 (P,'operators_completion.json','canonical_operators.log','NONIC57_SOURCE85_CANONICAL_ANNIHILATOR_COMPLETE'),
 (P,'cusp_arc_v5_pivot_completion.json','cusp_quad_arc_v5_pivot.log','NONIC57_SOURCE85_QUADRATIC_CUSP_ARC_COMPLETE'),
 (P,'ring24_v2_completion.json','ring24_v2.log','NONIC57_SOURCE85_RING24_TENSOR_COMPLETE'),
 (P,'good6_v2_completion.json','good6_jets_v2.log','NONIC57_SOURCE85_GOOD_JETS_COMPLETE case=6'),
 (P,'good8_v3_completion.json','good8_jets_v3.log','NONIC57_SOURCE85_GOOD_JETS_COMPLETE case=8'),
 (P,'rr24_correct_completion.json','rr24_correct.log','NONIC57_SOURCE85_RR24_CORRECT_ALL_EXCLUSIONS_COMPLETE'),
 (Q,'completion.json','model57.log','NONIC57_MODEL_CERTIFIED'),
 (Q,'jmap_v2_completion.json','jmap57_v2.log','NONIC57_JMODEL_AND_BOUNDS_COMPLETE'),
 (Q,'action_v2_completion.json','action57_v2.log','NONIC57_DIAMOND_AND_PRYM_COMPLETE'),
 (Q,'finite_v8_export_completion.json','finite57_v8_export.log','NONIC57_QUOTIENT9_NC_EXPORT_COMPLETE'),
]
for base,cf,lf,marker in checks_native:
 a=json.loads((base/cf).read_text());s=(base/lf).read_text();assert a['complete']and a['actual_exit_code']==0 and marker in s and not any(w in s for w in['Runtime error','User error','syntax error','Assertion failed']),(base,cf)
for d in range(6,10):
 a=json.loads((D/f'parallel_expansion_20261026/N57_f{d}/status.json').read_text());assert a['actual_exit_code']==0 and a['completed_marker'] and a['status']=='COMPLETE'
a=json.loads((Q/'partial_sage_certificate.json').read_text());assert a['complete']and a['tested_conservative_non_cusp_cases']==23 and len(a['surviving_images'])==3 and a['non_global_entirely_cuspidal_reductions_excluded']
a=js('pole_degrees_sage_certificate.json');assert a['complete']and a['pole_degrees']==[37,39]
a=js('quadratic_arc_sage_certificate.json');assert a['complete']and a['verified_canonical_jets']==20 and a['independent_selected_minor_valuation']==95 and a['nonzero_tangent']
a=js('all_cusp_FI_sage_certificate.json');assert a['complete']and a['all_normalized_configurations']==796111 and a['all_FI_ranks_equal9'] and a['annihilating_differential_dimension']==82
a=js('good_jets_sage_certificate.json');assert a['complete']and a['independent_genuine_R1_R2_arcs']==21 and a['geometric_cyclic5_kernels']==18 and a['all_division_polynomial_norms_exact'] and a['all_canonical_tangents_nonzero']
a=js('rr_correct_divisors_sage_certificate.json');assert a['complete']and a['case_count']==3 and a['definition']=='(T5-delta5-5)(delta5-1)NC, correct geometric push'
a=js('ring_correct_RR_sage_certificate.json');assert a['complete']and a['geometric_push_annihilator']=='T5-sigma-5'and a['necessary_divisor_source_coefficients']==[5,-4,-1]
assert a['all_R1_products']==3655 and a['all_R2_products']==31878 and a['all_q_product_identities_exact']and a['physical_point5_is_canonical_diamond_inverse']
assert [r['positive_degree']for r in a['nonprincipal_divisors']]==[42,42,88] and all(r['final_criterion_rank']==252 and r['multiplication_by_F_rank']==252 for r in a['nonprincipal_divisors'])
w=(D/'parallel_expansion_20261026/character_probes/N57_winding_data.m').read_text();assert 'SourceGenus:=85;'in w and 'ZeroLFactors:=1;'in w and 'FiniteProjectiveDiamondOrder:=18;'in w and 'Complete:=true;'in w
replays=[
 (Q/'cusp57.py','NONIC57_SOURCE_AND_DUAL_CUSP_AUDIT_PASS'),
 (Q/'replay_partial_sage.py','C57_QUOTIENT9_PARTIAL_INDEPENDENT_SAGE_PASS'),
 (P/'replay_pole_degrees_sage.py','C57_SOURCE85_DIRECT_POLE_DEGREES_SAGE_PASS'),
 (P/'replay_quadratic_arc_sage.py','C57_QUADRATIC_CUSP_ARC_INDEPENDENT_SAGE_PASS'),
 (P/'replay_all_cusp_FI_sage.py','C57_ALL_GLOBAL_CUSP_DIVISORS_FI_INDEPENDENT_SAGE_PASS'),
 (P/'replay_good_jets_sage.py','C57_GOOD_HECKE_AND_GENUINE_JETS_INDEPENDENT_SAGE_PASS'),
 (P/'assemble_correct_RR_divisors_sage.py','C57_NECESSARY_RR_DIVISORS_ASSEMBLED'),
 (P/'replay_correct_ring_RR_sage.py','C57_CANONICAL_RING_AND_ALL_RR_OBSTRUCTIONS_INDEPENDENT_SAGE_PASS'),
]
if '--deep'in sys.argv:
 env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 for script,marker in replays:
  t=time.monotonic();r=subprocess.run([sys.executable,str(script),'--readonly'],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  assert r.returncode==0 and marker in r.stdout,(script,r.returncode,r.stdout)
  print('READ_ONLY_REPLAY_PASS',script.name,round(time.monotonic()-t,2),flush=True)
print('C57_D9_COMPLETE_READ_ONLY_PACKET_PASS excluded=C57 correct_geometric_push=true files='+str(len(m['files'])),flush=True)
