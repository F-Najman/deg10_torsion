#!/usr/bin/env python3
"""Bind the independently reconstructed C91 data to completed original producers."""
from pathlib import Path
import json,hashlib
from unit_arithmetic import lit
H=Path(__file__).resolve().parent;R=H.parents[2];U=R/'parallel_expansion_20261026/N91_source_units'
assert(H/'canonical_unit_stage1_data.m').read_bytes()==(U/'canonical_unit_data.m').read_bytes()
e=json.loads((U/'status_v2.json').read_text());assert e['actual_exit_code']==0 and e['status']=='COMPLETE'and e['completion_marker'];assert hashlib.sha256((U/'geometry_v2.m').read_bytes()).hexdigest()==e['source_sha256'];assert hashlib.sha256((U/'run_v2.log').read_bytes()).hexdigest()==e['log_sha256']
a=(H/'elliptic_hecke_jets_stage2_data.m').read_text();ss=[];tt=[];jj=[]
for i in range(1,5):
 P=U/f'hecke3_source{i}';e=json.loads((P/'run.exit.json').read_text());assert e['actual_exit_code']==0 and e['completed_marker'];assert'NONIC91_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE'in(P/'run.log').read_text();b=(P/'elliptic_hecke3_data.m').read_text();ss+=lit(b,'SourceRecords');tt+=lit(b,'SourceHeckeRecords');jj+=lit(b,'SourceJPolynomials')
assert ss==lit(a,'SourceRecords')and tt==lit(a,'SourceHeckeRecords')and jj==lit(a,'SourceJPolynomials')and len(ss)==4
for f in range(6,10):
 P=R/'parallel_expansion_20261026'/f'N91_f{f}';e=json.loads((P/'status.json').read_text());assert e['actual_exit_code']==0 and e['status']=='COMPLETE'and e['completed_marker'];assert hashlib.sha256((R/'parallel_expansion_20261026/odd_census.m').read_bytes()).hexdigest()==e['source_sha256'];assert hashlib.sha256((P/'run.log').read_bytes()).hexdigest()==e['log_sha256']
P=R/'cyclic/continuation_20261026';e=json.loads((P/'winding_N91_status.json').read_text());assert e['native_exit_code']==0 and e['status']=='EXITED';assert hashlib.sha256((P/'winding_N91.log').read_bytes()).hexdigest()==e['log_sha256']
s=(R/'parallel_expansion_20261026/character_probes/N91_winding_data.m').read_text();flags=lit(s,'ConstituentWindingFlags');HH={u for u in range(91)if __import__('math').gcd(u,91)==1}
for d,lev,w,o,k in flags:
 if not w:HH.intersection_update(k)
assert HH==set(lit(s,'FiniteDiamondUnits'))and 3 in HH and 27 in HH and len(HH)==12
print('NONIC91_INPUT_BINDING_PASS completed_native_exits10 exact_combined_records=true finite_diamond3_and27=true')
