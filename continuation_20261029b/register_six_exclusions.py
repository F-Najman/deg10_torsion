from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parents[1];P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def job(stem,marker=None):
 s=json.loads((H/'jobs'/(stem+'.status.json')).read_text())
 assert s['state']=='EXITED' and s['actual_exit_code']==0,(stem,s)
 assert sha(H/s['script'])==s['script_sha256']
 if marker:assert marker in (H/s['log']).read_text(),stem
 return [s['script'],s['log'],'jobs/'+stem+'.status.json']
def paths(ps):return [str(p.relative_to(H)) for p in ps]
def evidence(fs):return [{'path':s,'sha256':sha(H/s)} for s in sorted(set(fs))]
pre='continuation_20261029b/'
common=['exclusions/continuation_20261030.tex','exclusions/isogeny_transfer.tex',pre+'literature_inputs.json',pre+'six_exclusions_replay.json',pre+'four_exclusions_replay.json',pre+'independent_sources_all.json',pre+'independent_sources.json',pre+'cusp_tate.json',pre+'reuse_manifest.json']
for stem in ['late_winding_ranks','late_cusp_tate','late_replay_sources','late_replay_sources_all','late_replay_filters_v2','late_replay_six']:common+=job(stem)
common+=paths(P.glob('winding_N*.m'))+paths((P/'deps').glob('*.m'))
# The independent all-source replay also checks the C226 census, not its exclusion.
common+=paths((P/'C226').glob('source*.m'))+[pre+'C226/frobenius.m']
common+=['continuation_20261029/deps/fast_hecke.m','shared/deps/class_rational_reduction.m']
C38=common+paths((P/'C238').glob('source*.m'))+[pre+'C238/'+s for s in ['frobenius.m','filter19_v2.m','target19_v2_RR.log','target19_v2_RR.log.witnesses.py','hecke_template.m','native_control.m','hecke_completion.json','native_control_completion.json']]
C38+=job('late238_source')+job('late238_filter19_v2')
for k in range(4):
 C38+=job(f'late238_hecke_s{k}','C238_HECKE_SHARD_COMPLETE')+job(f'late238_control_s{k}','C238_NATIVE_CONTROL_SHARD_COMPLETE')
 C38+=[pre+f'C238/hecke_shard{k}.m',pre+f'C238/control_shard{k}.m']
C40=common+paths((P/'C240').glob('M2N40_D10*'))+[pre+'C240/'+s for s in ['cusp_degrees.txt','source_frobenius.txt','preflight.m','enumerate_direct.m','filter_quotient20_linear.m','filter_quotient20_isogeny.m','filter_quotient20_linear_arithmetic.txt','filter_quotient20_isogeny_arithmetic.txt','quotient20_linear_RR_identity.log','quotient20_linear_RR_identity.log.witnesses.py','quotient20_isogeny_RR_identity.log','quotient20_isogeny_RR_identity.log.witnesses.py','shard_completion.json','native_replay_v2_completion.json','replay_ideal_identities_v2.m']]
for stem in ['late240_preflight','late240_enum10','late240_linear','late240_isogeny']:C40+=job(stem)
C40+=paths((P/'C240/portable/imported').rglob('*'))
C40=[s for s in C40 if (H/s).is_file()]
C40+=[pre+'C240/portable/'+s for s in ['portable_inputs.m','portable_forward.m','import_provenance.json']]
for k in range(8):
 C40+=job(f'late240_picard_s{k}','complete=true')+job(f'late240_replay_v2_s{k}','PICARD_IDEAL_REPLAY_COMPLETE')
 C40+=[pre+f'C240/filter_cyclic40_shard{k}.m']
 C40+=[pre+f'C240/portable/{s}' for s in [f'shard{k}_arithmetic.txt',f'shard{k}_RR_identity.log',f'shard{k}_RR_identity.log.witnesses.py',f'replay{k}_data.m',f'replay{k}_result.txt']]
C50=common+paths((P/'C250').glob('source*.m'))+[pre+'C250/'+s for s in ['frobenius.m','filter25.m','target25_RR.log','target25_RR.log.witnesses.py','hecke_template.m','hecke_completion.json','kernel_certificates.m','independent_kernel_replay.json']]
C50+=job('late250_source')+job('late250_filter25')+job('late250_kernel_certificate')+job('late250_replay_kernels')
for k in range(3):C50+=job(f'late250_hecke_s{k}','C250_HECKE_SHARD_COMPLETE')+[pre+f'C250/hecke_shard{k}.m']
for f in ['six_exclusions_replay.json','four_exclusions_replay.json']:
 assert json.loads((P/f).read_text())['status']=='PASS'
assert json.loads((P/'C250/independent_kernel_replay.json').read_text())['status']=='PASS'
ledger=json.loads((H/'certified_results.json').read_text());present={tuple(r['group']) for r in ledger['exclusions']}
rows=[]
for n,fs,method in [(38,C38,'Four X1(19) Picard filters and complete T5-1-5<5> tests on X1(38); all native Hecke controls pass'),(40,C40,'Two X1(2,20) maps and complete proven Picard sieve on X1(40), with independent exact ideal-identity replay'),(50,C50,'Four X1(25) Picard filters and complete T7-1-7<7> tests on X1(50); all cyclic kernels independently checked in Sage')]:
 for group in [[2,n],[1,2*n]]:
  assert tuple(group) not in present
  desc=method if group[0]==2 else 'Isogeny-transfer consequence of the certified exclusion of C2 x C'+str(n)
  rows.append({'group':group,'verified':True,'method':desc,'evidence':evidence(fs)})
ledger['exclusions']+=rows
(H/'certified_results.json').write_text(json.dumps(ledger,indent=2)+'\n')
print('REGISTERED_SIX_EXCLUSIONS',[r['group'] for r in rows])
