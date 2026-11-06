from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent;H=P.parents[1];REPO=H.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
rows=[r for s in range(3)for r in rd(f'final_Hecke7_v2_s{s}.m')[2]if r[4]==0];assert len(rows)==43 and sum(len(r[3])for r in rows)==55
native=[rd(f'sieve_native_A7_s{s}.m')for s in[0,1]];ids=[r[0]for a in native for r in a[0]];assert sorted(ids)==sorted(r[0]for r in rows)and len(set(ids))==43;assert sum(a[1][1]for a in native)==292 and all(r[3]==0 for a in native for r in a[0]);sg=json.loads((P/'sieve_independent_genuine_A7.json').read_text());assert sg['status']=='PASS'and sg['excluded_supports']==43 and sg['physical_source_context_uses']==55 and sg['distinct_closed_marked_sources']==45 and sg['geometric7_kernels_over_distinct_sources']==360 and sg['native_prime_references_bound']==292;assert len(sg['observations'])==43
names=['d31b_m318_sieve_native_s0','d31b_m318_sieve_native_s1','d31b_m318_sieve_genuine_A7'];jobs=[]
for n in names:
 sp=H/'jobs'/(n+'.status.json');s=json.loads(sp.read_text());assert s['state']=='EXITED'and s['actual_exit_code']==0 and sha(H/s['script'])==s['script_sha256'];lp=H/s['log'];text=lp.read_text();marker='D31B_M318_SIEVE_NATIVE_COMPLETE'if'native'in n else'D31B_M318_SIEVE_GENUINE_ALL_COMPLETE';assert marker in text and not any(t in text for t in['Runtime error','User error:','Traceback (most recent call last):']);jobs.append({'purpose':n,'script':str((H/s['script']).relative_to(REPO)),'log':str(lp.relative_to(REPO)),'status':str(sp.relative_to(REPO)),'actual_exit_code':0})
active=[]
for x in Path('/proc').iterdir():
 if not x.name.isdigit():continue
 try:
  cmd=[v.decode(errors='replace')for v in(x/'cmdline').read_bytes().split(b'\0')if v]
  if cmd and(Path(cmd[0]).name.startswith('magma')or Path(cmd[0]).name.startswith('python'))and any(Path(v).stem.startswith('d31b_m318_sieve_')and Path(v).stem!='d31b_m318_sieve_build_A7_controls'for v in cmd[1:]):active.append(int(x.name))
 except OSError:pass
assert not active,active
summary={'status':'PASS','scope':'Independent certificate for43 zero A7 tests, not a complete exclusion of C3 x C18','group':'C3 x C18','prime':13,'hecke':7,'excluded_NC_supports':43,'physical_source_summands':55,'distinct_marked_sources':45,'complete_geometric7_kernels_over_distinct_sources':360,'explicit_Frobenius_kernel_orbit_payloads':sg['explicit_kernel_orbit_payloads_checked'],'actual_native_prime_references':292,'full_marked_source_and_diamond_Weil_component_checked':True,'both_independent_torsion_points_checked_under_every_isogeny':True,'all_stored_kernel_polynomials_source_embeddings_and_orbits_checked':True,'independent_manual_Velu_on_both_marked_points':True,'all_norms_and_physical_prime_keys_match':True,'all43_fresh_native_RR_dimensions_zero':True,'active_own_workers':active,'paused_own_workers':0,'own_automatic_launchers':0};(P/'sieve_A7_controls_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={P/'d31b_m318_sieve_build_A7_controls.py',P/'d31b_m318_sieve_native_common.m',P/'sieve_native_A7_s0.m',P/'sieve_native_A7_s1.m',P/'sieve_independent_genuine_A7.json',P/'sieve_A7_controls_summary.json',P/'paired_final_Hecke_inputs_v2.m',P/'final_Hecke_common_v2.m',P/'isogenies_witness.m',P/'inputs/fast_full_aut_hash_m3_18_root.m',H/'run_job.py'};paths|={P/f'final_Hecke7_v2_s{s}.m'for s in range(3)}
for j in jobs:
 for k in['script','log','status']:paths.add(REPO/j[k])
paths.update(z for z in(H/'shared/mdmagma/v2').rglob('*')if z.is_file());m={'status':'PASS','group':'C3 x C18','scope':'43 negative A7 support tests only','summary_sha256':sha(P/'sieve_A7_controls_summary.json'),'files':[{'path':str(z.relative_to(REPO)),'sha256':sha(z),'bytes':z.stat().st_size}for z in sorted(paths)],'jobs':jobs};(P/'sieve_A7_controls_manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('D31B_M318_SIEVE_A7_CONTROLS_FROZEN',len(paths),sha(P/'sieve_A7_controls_manifest.json'))
