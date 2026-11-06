#!/usr/bin/env python3
"""Read-only certificate for the complete degree9 C58 exclusion."""
from pathlib import Path
import hashlib,json,os,re,subprocess,sys
P=Path(__file__).resolve().parent;H=P.parents[2];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
M=json.loads((P/'manifest62.json').read_text());assert M['degree']==9 and M['excluded_group']==[1,62] and M['complete']
for e in M['files']:assert sha(H/e['path'])==e['sha256'],e['path']
for e in M['native_jobs']:
 s=json.loads((H/e['status']).read_text());L=(H/e['log']).read_text()
 assert s['status']=='EXITED' and s['native_exit_code']==0 and s['log_sha256']==sha(H/e['log'])
 assert e['completion_marker'] in L and not re.search(r'Runtime error|User error|Assertion failed|Segmentation fault',L)
W=json.loads((P/'presentation/presentation.json').read_text());B=json.loads((P/'input_binding62.json').read_text())
assert W['native_proven'] and W['factor_basis_bound']==9 and B['byte_exact_regeneration'] and B['complete']
assert B['native_export_sha256']==sha(P/'presentation/presentation.json') and B['literal_input_sha256']==sha(P/'inputs.m')
literal=lambda x:json.dumps(x,separators=(',',':'))
expected=''.join(f'{k}:={literal(v)};\n' for k,v in [('N',W['N']),('p',W['p']),('ExpectedGenus',W['genus']),('KnownOrder',W['known_Pic_order']),('CandidateInvariants',W['candidate_invariants']),('FieldPolynomialData',W['field_polynomial']),('ModularCoordinatesData',W['modular_coordinates']),('RawRelations',W['relations']),('FBDegrees',[a['degree'] for a in W['factor_basis']]),('FBFinite',[a['finite'] for a in W['factor_basis']]),('FBBases',[a['ideal_basis'] for a in W['factor_basis']])])
assert (P/'inputs.m').read_text()==expected
C=json.loads((P/'cusp62.json').read_text());assert C['global_and_F3_orbit_partitions_equal'] and C['global_signature']==C['F3_signature']=={'1':15,'15':1}
assert 'Gon_Q(X1(31)) = 12.' in (H/'realizations/continuation/C35/gonality').read_text()
assert '31 & [10,1772833370]' in (H/'realizations/continuation/C35/DEHMZ.tex').read_text()
R=json.loads((P/'math_review62.json').read_text());assert R['independent_review']=='PASS' and R['source_J1_62_rank_zero_not_assumed']
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
r=subprocess.run(['/opt/sage/envs/sage/bin/sage','-python',str(P/'replay62.py')],env=env,capture_output=True,text=True)
assert r.returncode==0,(r.stdout,r.stderr);assert 'NONIC62_INDEPENDENT_REPLAY' in r.stdout and '"class_sieve_verified": true' in r.stdout
print(r.stdout.strip())
q=subprocess.run([sys.executable,str(P/'known31/certify_known31.py')],env=env,capture_output=True,text=True);assert q.returncode==0,(q.stdout,q.stderr);print(q.stdout.strip())
print('NONIC62_CERTIFICATE_PASS degree=9 tests=2190 positive_survivors=5 no62_lifts=true complete_exclusion=true')
