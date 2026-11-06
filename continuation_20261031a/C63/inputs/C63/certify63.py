#!/usr/bin/env python3
"""Read-only complete degree9 C63 packet audit."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
P=Path(__file__).resolve().parent
M=json.loads((P/'manifest63.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[63] and M['complete']
for n,h in M['files'].items():assert hashlib.sha256((P/n).read_bytes()).hexdigest()==h,n
for n,h in M['dependencies'].items():assert hashlib.sha256((P/n).read_bytes()).hexdigest()==h,n
log=(P/'cartier_rr_D9.log').read_text()
assert not re.search('Runtime error|User error|Assertion failed',log)
assert 'RANK0_CARTIER_DONE N=63 D=9 tested=105040 total=105040 exceptions=0 complete=true' in log
assert len(re.findall('C63_RR_KERNEL_DONE',log))==121
for e,k,pairs in [(7,3,189),(8,3,216),(9,4,324)]:
 log=(P/(f'good{e}_T5_fresh.log' if e<9 else 'good9_T5.log')).read_text()
 assert not re.search('Runtime error|User error|Assertion failed',log)
 assert f'Y63_NC_EXHAUSTIVE e={e} geometric_pairs={pairs} closed_places={9*k} diamond_orbits={k}' in log
 assert re.findall(r'Y63_T5_TEST index=(\d+) principal=(\w+)',log)==[(str(i),'false') for i in range(1,k+1)]
 assert f'Y63_NC_T5_DONE e={e} tests={k} complete=true' in log
log=(P/'winding63.log').read_text()
assert 'Y63_RANK0_DONE invariant_character_orbits=5 dimension_sum=49 all_winding_full=true' in log
assert 'closed=[ 0, 0, 0, 0, 0, 0, 54, 54, 72 ] complete=true' in (P/'source_counts.log').read_text()
C=json.loads((P/'cusp63_quotient.json').read_text())
assert C['involution_free_on_all_cusps'] and C['all_degree_le2_F2_cusps_lift_to_same_global_degree']
assert C['F2_signature']=={'1':9,'2':3,'3':9,'6':1}
assert C['global_signature']=={'1':9,'2':3,'3':6,'6':1,'9':1}
for rec in json.loads((P/'native_exit_records.json').read_text()):
 assert rec['actual_exit_code']==0
 assert hashlib.sha256((P/rec['source']).read_bytes()).hexdigest()==rec['source_sha256']
 assert hashlib.sha256((P/rec['log']).read_bytes()).hexdigest()==rec['log_sha256']
for rec in json.loads((P/'good7_8_exit_records.json').read_text()):assert rec['actual_exit_code']==0
for name in ['replay_rr_kernels.py','replay_cartier_projection.py','replay_fi_jets.py']:
 r=subprocess.run([sys.executable,str(P/name),str(P/'cartier_kernels_D9.m')],text=True,capture_output=True,check=True)
 if name=='replay_fi_jets.py':assert '"covered": 105040' in r.stdout and '"exceptions": 0' in r.stdout and '"complete": true' in r.stdout
print('NONIC_C63_READONLY_CERTIFICATE_PASS degree=9 quotient_genus=49 rankQ=0 U_dimension=45 exact_RR_blocks=121 configurations=105040 all_rank9=true good_Hecke_tests=10 all_nonprincipal=true complete=true')
