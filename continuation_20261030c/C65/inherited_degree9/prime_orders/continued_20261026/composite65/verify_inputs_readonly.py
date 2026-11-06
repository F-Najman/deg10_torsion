#!/usr/bin/env python3
"""Bind completed C65 ordinary batches and exact winding/census provenance."""
from pathlib import Path
import json,hashlib,math
from unit_arithmetic import lit
H=Path(__file__).resolve().parent;R=H.parents[2];U=R/'parallel_expansion_20261026/N65_source_units';assert(H/'canonical_unit_stage1_data.m').read_bytes()==(U/'canonical_unit_data.m').read_bytes();e=json.loads((U/'status_v2.json').read_text());assert e['actual_exit_code']==0 and e['status']=='COMPLETE'and e['completion_marker'];assert hashlib.sha256((U/'geometry_v2.m').read_bytes()).hexdigest()==e['source_sha256'];assert hashlib.sha256((U/'run_v2.log').read_bytes()).hexdigest()==e['log_sha256'];ss=[];tt=[];jj=[];fds=[]
for f,k in[(7,2),(8,2),(9,3)]:
 for i in range(1,k+1):
  P=H/f'hecke3_f{f}_source{i}';e=json.loads((P/'run.exit.json').read_text());assert e['native_exit_code']==0 and e['termination_reason']is None and e['status']=='COMPLETED';assert'NONIC65_ELLIPTIC_HECKE_JETS_STAGE2_COMPLETE'in(P/'run.log').read_text();s=(P/'elliptic_hecke3_data.m').read_text();ss+=lit(s,'SourceRecords');tt+=lit(s,'SourceHeckeRecords');jj+=lit(s,'SourceJPolynomials');fds.append(f)
s=(H/'elliptic_hecke_jets_stage2_data.m').read_text();assert ss==lit(s,'SourceRecords')and tt==lit(s,'SourceHeckeRecords')and jj==lit(s,'SourceJPolynomials')and fds==lit(s,'SourceDegrees')
for f in range(6,10):
 P=R/'parallel_expansion_20261026'/f'N65_f{f}';e=json.loads((P/'status.json').read_text());assert e['actual_exit_code']==0 and e['status']=='COMPLETE'and e['completed_marker'];assert hashlib.sha256((R/'parallel_expansion_20261026/odd_census.m').read_bytes()).hexdigest()==e['source_sha256'];assert hashlib.sha256((P/'run.log').read_bytes()).hexdigest()==e['log_sha256']
P=R/'cyclic/continuation_20261026';e=json.loads((P/'winding_N65_status.json').read_text());assert e['native_exit_code']==0 and e['status']=='EXITED';assert hashlib.sha256((P/'winding_N65.log').read_bytes()).hexdigest()==e['log_sha256'];s=(R/'parallel_expansion_20261026/character_probes/N65_winding_data.m').read_text();HH={u for u in range(65)if math.gcd(u,65)==1}
for d,lev,w,o,k in lit(s,'ConstituentWindingFlags'):
 if not w:HH.intersection_update(k)
assert HH==set(lit(s,'FiniteDiamondUnits'))and len(HH)==48
print('NONIC65_ORDINARY_INPUT_BINDING_PASS completed_native_exits13 exact_combined_records=true finite_fullG24=true')
