"""Freeze C99's complete exclusion and distinguish non-excluding diagnostics."""
from pathlib import Path
import hashlib,json,shutil
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;C=H/'continuation_20261030c/C99'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old=json.loads((ROOT/'work/degree9/certified_results.json').read_text());r=next(r for r in old['exclusions']if r['group']==[1,99]);snaps=[]
for e in r['evidence']:
 src=ROOT/'work/degree9'/e['path'];assert sha(src)==e['sha256'],src
 dst=C/'inherited_degree9'/e['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);assert sha(dst)==e['sha256']
 snaps.append({'source':str(src.relative_to(ROOT)),'snapshot':str(dst.relative_to(ROOT)),'sha256':e['sha256']})
(C/'inherited_snapshots.json').write_text(json.dumps({'scope':'Exact degree-independent winding and finite-image evidence reused. Old degree9 inequalities, exclusion and CM tables are historical provenance only, not degree10 premises. All earlier-degree originals are unchanged.','files':snaps},indent=2)+'\n')
accepted=set();diagnostic=set()
def add(p,where=accepted):
 assert p.is_file();where.add(p)
for p in C.rglob('*'):
 if p.is_file():
  (diagnostic if any(s in p.name for s in ['quotient33_A7','quotient33_3isog'])else accepted).add(p)
for st in sorted((H/'jobs').glob('contc99*.status.json')):
 d=json.loads(st.read_text());assert d['state']=='EXITED';isdiag=d['actual_exit_code']!=0 or any(s in d['purpose']for s in['quotient33_A7','quotient33_3isog']);where=diagnostic if isdiag else accepted
 for p in[st,H/d['log'],H/d['script']]:add(p,where)
for st in[H/'jobs/contc_parent_source_census.status.json']:
 d=json.loads(st.read_text());assert d['state']=='EXITED'and d['actual_exit_code']==0
 for p in[st,H/d['log'],H/d['script'],H/'continuation_20261030c/parent_source_census.json']:add(p)
for p in(H/'shared/mdmagma/v2').rglob('*'):
 if p.is_file()and p.suffix in['.m','.sig','.spec']:add(p)
add(H/'shared/mdmagma_manifest.json');add(Path(__file__).resolve())
summary={'status':'PASS','degree':10,'group':[1,99],'complete_exclusion':True,'parent_mathematical_review':'PASS','source_closed_counts':[0]*9+[120],'source_diamond_orbits':4,'initial_target_genus':21,'initial_A5_h0':[0,0,1,0],'excluded_initial_orbits':[1,2,4],'final_target_genus':29,'final_operator':'(delta2-1)(T5-1-5delta5)','final_orbit':3,'final_h0':0,'signed_degrees':[110,110],'fresh_native_prime_reference_checks':23,'independent_marked_5_and_11_isogenies':'PASS; every source/diamond/paired Hecke key and norm bound','all_recorded_jobs_exited':True,'pure_cusp_argument':'finite diamond image, unramified2 torsion exponent2, gonality>40, nonintegral j excludes nontrivial projective stabilizers, free order15 cannot divide10','diagnostic_note':'A7 and extra3-isogeny X1(33) maps retained h0=1 and do not exclude residual3. Four initial Y29 model failures are not certificates. Only the final verified genus29 test closes residual3.'}
p=C/'certificate_summary.json';p.write_text(json.dumps(summary,indent=2)+'\n');add(p)
def rows(ps):return[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in sorted(ps)]
manifest={'status':'PASS','scope':'Complete C99 degree10 exclusion, independent native and genuine marked-isogeny checks; only exact winding input transferred from degree9','files':rows(accepted),'diagnostic_scope':'Preserved unsuccessful models and genuinely surviving smaller-quotient tests; these are not exclusions','diagnostic_files':rows(diagnostic)}
p=C/'certificate_manifest.json';p.write_text(json.dumps(manifest,indent=2)+'\n')
for e in manifest['files']+manifest['diagnostic_files']:assert sha(ROOT/e['path'])==e['sha256']
print(json.dumps({'status':'FROZEN','accepted_files':len(accepted),'diagnostic_files':len(diagnostic),'manifest_sha256':sha(p)},indent=2))
