"""Read-only complete degree-nine C65 exclusion certificate.

--deep independently replays the scoped ordinary packet, supersingular
census, genus-nine cusp-class sieve, genuine quartic cusp arc, all 2250325
formal-immersion ranks, actual CM isogenies/jets, norm aggregation, and
all canonical multiplication identities and final RR matrices.
No Magma process is started or resumed.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
P=Path(__file__).resolve().parent;D=P.parents[2];C=P.parent/'C65_CM6';A=P.parent/'C65_allbad';Q=P.parent/'C65_quotient9';O=D/'prime_orders/continued_20261026';m=json.loads((P/'manifest.json').read_text())
assert m['degree']==9 and m['complete_exclusion']and m['excluded_groups']==[[1,65]]and m['necessary_point_coefficients']==[1,2,-3]
for r in m['files']:
 p=D/r['path'];assert p.is_file()and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],p
assert hashlib.sha256((O/'composite65/ordinary_manifest.json').read_bytes()).hexdigest()=='819088e6576329931d614fb1bb0ca43d9779c9ae98dd44c065eebd203d56172b'
ordinary_proof=(O/'point_convention_correction_20261027/case65_ordinary_corrected.tex').read_text();assert'A_3=T_3-1-3' in ordinary_proof and'336' in ordinary_proof
native=[
 (C,'good6_completion.json','good6_jets.log','NONIC65_CM6_GOOD_JETS_COMPLETE case=6'),
 (C,'rr12_completion.json','rr12.log','NONIC65_CM6_RR12_ALL_EXCLUSIONS_COMPLETE'),
 (A,'completion.json','canonical_operators.log','NONIC65_ALLBAD_CANONICAL_ANNIHILATOR_COMPLETE'),
 (A,'cusp_quartic_arc_completion.json','cusp_quartic_arc.log','NONIC65_ALLBAD_QUARTIC_CUSP_ARC_COMPLETE'),
 (Q,'preflight_completion.json','preflight.log','NONIC65_QUOTIENT_PREFLIGHT_COMPLETE'),
 (Q,'model_completion.json','model65.log','NONIC65_MODEL_CERTIFIED'),
 (Q,'jmap_completion.json','jmap65.log','NONIC65_JMODEL_AND_BOUNDS_COMPLETE'),
 (Q,'action_completion.json','action65.log','NONIC65_QUOTIENT_DIAMOND3_COMPLETE'),
 (Q,'finite_completion.json','finite65.log','NONIC65_FINITE_A3_COMPLETE'),
]
for B,cf,lf,marker in native:
 r=json.loads((B/cf).read_text());s=(B/lf).read_text();assert r['actual_exit_code']==0 and r['complete']and not r['wall_timeout']and marker in s and not any(w in s for w in['Runtime error','User error','syntax error','Assertion failed']),(B,cf,marker)
for d in range(6,10):
 s=json.loads((D/f'parallel_expansion_20261026/N65_f{d}/status.json').read_text());assert s['status']=='COMPLETE'and s['actual_exit_code']==0 and s['completed_marker']
w=(D/'parallel_expansion_20261026/character_probes/N65_winding_data.m').read_text();assert'SourceGenus:=121;'in w and'ZeroLFactors:=1;'in w and'FiniteProjectiveDiamondOrder:=24;'in w and'Complete:=true;'in w
s=json.loads((A/'quartic_arc_sage_certificate.json').read_text());assert s['complete']and s['verified_canonical_jets']==20 and s['independent_selected_minor_valuation']==98 and s['nonzero_tangent']
s=json.loads((A/'all_cusp_FI_sage_certificate.json').read_text());assert s['complete']and s['all_normalized_configurations']==2250325 and s['annihilating_differential_dimension']==106 and s['all_FI_ranks_equal9']
s=json.loads((Q/'cusp_class_sage_certificate.json').read_text());assert s['complete']and s['finite_J_invariants']==[4,16,336]and s['finite_alpha_cusp_image_order']==32 and s['mu_A3_kernel_order']==8 and s['mu_A3_boundary_and_base_corrections_independently_recomputed']and s['all_non_global_allbad_reductions_excluded']
s=json.loads((C/'supersingular_census_sage_certificate.json').read_text());assert s['complete']and s['unique_eligible_F64_isomorphism_class_matches_native']and s['eligible_F64_group_cyclic65']
s=json.loads((C/'cm6_sage_certificate.json').read_text());assert s['complete']and s['geometric_cyclic3_kernels']==4 and s['independent_genuine_R1_arcs']==7 and s['all_division_polynomial_norms_exact']and s['all_unit_recovery_ranks120']and s['all_canonical_tangents_nonzero']
s=json.loads((C/'ring_RR_sage_certificate.json').read_text());assert s['complete']and s['necessary_point_coefficients']==[1,2,-3]and s['all_ordered_R1_products']==14641 and s['all_q_product_identities_exact']and len(s['nonprincipal_divisors'])==2 and all(r['positive_degree']==36 and r['negative_degree']==36 and r['final_criterion_rank']==121 and r['multiplication_by_F_rank']==121 for r in s['nonprincipal_divisors'])
replays=[
 (O/'certify_65_ordinary_readonly.py',[],'NONIC65_ORDINARY_FROZEN_CERTIFICATE_PASS'),
 (C/'replay_supersingular_census_sage.py',['--readonly'],'C65_SUPERSINGULAR_COMPLETE_CENSUS_INDEPENDENT_PASS'),
 (Q/'replay_cusps_and_classes.py',['--readonly'],'C65_NONGLOBAL_CUSP_CLASS_INDEPENDENT_PASS'),
 (A/'replay_quartic_arc_sage.py',['--readonly'],'C65_QUARTIC_CUSP_ARC_INDEPENDENT_SAGE_PASS'),
 (A/'replay_all_cusp_FI_sage.py',['--readonly'],'C65_ALL_GLOBAL_CUSP_DIVISORS_FI_INDEPENDENT_SAGE_PASS'),
 (C/'replay_cm6_jets_sage.py',['--readonly'],'C65_CM6_GOOD_HECKE_AND_GENUINE_JETS_INDEPENDENT_SAGE_PASS'),
 (C/'assemble_RR_divisors_sage.py',['--readonly'],'C65_NECESSARY_RR_DIVISORS_ASSEMBLED'),
 (C/'replay_ring_RR_sage.py',['--readonly'],'C65_CM6_RING_RR_INDEPENDENT_SAGE_COMPLETE_PASS'),
]
if '--deep'in sys.argv:
 env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 for script,args,marker in replays:
  tt=time.monotonic();r=subprocess.run([sys.executable,str(script)]+args,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);assert r.returncode==0 and marker in r.stdout,(script,r.returncode,r.stdout);print('C65_READ_ONLY_REPLAY_PASS',script.name,round(time.monotonic()-tt,2),flush=True)
print('C65_D9_COMPLETE_READ_ONLY_PACKET_PASS excluded=C65 native_new_actual_zero='+str(len(native))+' ordinary_matrices=336 allbad_FI=2250325 CM_RR_rank=121 files='+str(len(m['files'])),flush=True)
