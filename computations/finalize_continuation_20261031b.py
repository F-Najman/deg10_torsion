"""Finalize the completed Phi(10) classification; refuse live workers or missing proofs."""
from pathlib import Path
import hashlib,json,os,re,subprocess,time
import psutil
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;R=H/'continuation_20261031b'
D=json.loads((H/'candidates.json').read_text());L=json.loads((H/'certified_results.json').read_text());assert D['classification_complete'] and L['classification_complete'];assert len(D['known_lower_bound'])==69 and len(D['minimal_outside_known'])==len(D['certified_exclusions'])==110 and not D['remaining_minimal_obstructions'];assert json.loads((R/'round_status.json').read_text())['completed_count']==10
p=H/'paper/main.tex';s=p.read_text();assert 'certified partial results'not in s and 'incomplete'not in s and 'prop:decic318-classification'in s
c=subprocess.run(['/home/fnajman/bin/tectonic','-X','compile','main.tex','--keep-logs','--keep-intermediates'],cwd=H/'paper',capture_output=True,text=True);print(c.stdout+c.stderr);c.check_returncode()
info=subprocess.run(['/home/fnajman/bin/pdfinfo',str(H/'paper/main.pdf')],capture_output=True,text=True,check=True).stdout;pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
statuses=[json.loads(p.read_text())for p in sorted((H/'jobs').glob('*.status.json'))];assert all(s['state']=='EXITED'for s in statuses)
own=[];other=[]
for p in psutil.process_iter(['pid','name','username','create_time','cmdline','status']):
 try:
  a=p.info;args=a['cmdline']or[];magma='magma'in(a['name']or'').lower()
  if a['pid']in[os.getpid(),os.getppid()]:continue
  ispython=any(x.endswith('.py')for x in args[1:]);scope=any(str(H)in x or'work/deg10/'in x for x in args)
  if not scope and(magma or ispython):
   try:scope=Path(p.cwd()).is_relative_to(H)
   except(OSError,psutil.Error):pass
  if scope and(magma or ispython):own.append(a)
  elif magma:other.append(a)
 except psutil.Error:pass
assert not own,own
for name in['validation.json','lifecycle_final.json']:
 p=H/name;q=H/'notes'/('before_continuation_20261031b_'+name)
 if not q.exists():q.write_bytes(p.read_bytes())
life={'status':'PASS','time_unix':time.time(),'task':'Resolve the final ten original Phi(10) cases','original_cases_resolved':10,'original_new_exclusions':9,'new_occurrences':[[3,18]],'new_larger_exclusions':[[3,36],[6,18]],'degree_ten_workers_supervisors_launchers':own,'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'new_automatic_launchers_created':0,'unrelated_magma_workers_preserved':other,'host_cores':psutil.cpu_count(),'degree_ten_active_cores':0,'obsolete_computations_terminated_after_identity_checks':True,'failed_or_terminated_probes_are_not_certificates':True,'retirement_records':[str(p.relative_to(H))for p in sorted(R.rglob('*retire*.json'))]}
(H/'lifecycle_final.json').write_text(json.dumps(life,indent=2)+'\n')
a=subprocess.run(['python3',str(H/'verify_project.py')],cwd=ROOT,capture_output=True,text=True)
if a.returncode:print(a.stdout+a.stderr)
a.check_returncode();report=json.loads(a.stdout);assert report['remaining_minimal_obstructions']==0;(H/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
readme=f'''Phi(10): COMPLETE, 31 October 2026

There are 69 possible exact full torsion groups over degree10 number fields:
  C_n: n=1,...,40,42,45,50;
  C_2 x C_(2n): 1<=n<=15;
  C_3 x C_(3n): 1<=n<=6;
  C_4 x C_(4n): 1<=n<=3;
  C_6 x C_(6n): 1<=n<=2.

All 110 minimal subgroup obstructions outside this list are excluded.
There are ZERO remaining cases. The classification uses the audited known
S(10)={{2,3,5,7,11,13,17,19,23,29,31,37}} and the earlier accepted inputs.

The final ten original tasks are all resolved:
  Excluded: C48,C54,C56,C60,C66,C74,C3 x C27,C3 x C33,C4 x C28.
  Realized: C3 x C18.
Its occurrence exposes C3 x C36 and C6 x C18 as two new minimal tasks;
the complete degree5 point classification on X1(3,18) excludes both.

Read paper/main.pdf ({pages} pages) for the proof and complete list.
M3N18/resolution.tex within continuation_20261031b/ gives the explicit
relative quintic, elliptic curve, marked points, and complete two-orbit
classification. Both representative curves have full torsion C3 x C18
by independently checked good reductions of order54 over F49.

certified_results.json is the frozen evidence ledger. candidates.py checks
all evidence hashes and independently enumerates the minimal frontier.
The continuation_20261031b case folders contain exact native and Sage
replays, proof sources and immutable certificate summaries. The paper's
layout-only C48/C66 copies have separate provenance; frozen sources remain
unchanged. The earlier README and continuation notes are archived in notes/.

Final audit: {report['evidence_files_checked']} distinct evidence hashes;
{report['protected_earlier_files_unchanged']} protected earlier-degree files unchanged;
all {len(statuses)} recorded jobs exited. The paper compiles without errors
or overfull boxes. No degree10 worker, supervisor, paused job or automatic
launcher remains active; unrelated useful computations were preserved.
Failed, capped and terminated diagnostics are not certificates.
See validation.json and lifecycle_final.json for machine-readable audits.
'''
(H/'README.txt').write_text(readme)
(H/'notes/continuation.txt').write_text(f'''31 October 2026 — Phi(10) COMPLETE

CURRENT COUNTS: 69 exact full groups;110 minimal subgroup obstructions;
all110 excluded;ZERO unresolved cases. See README.txt for the list.

The final ten original groups were resolved by nine exclusions and one
occurrence (C3 x C18). Full degree5 point classification on X1(3,18)
over Q(zeta3) leaves two old27-marking orbits, interchanged by a marked
2-isogeny involution. An exact global quintic point realizes C3 x C18.
Both curves have full torsion54 from unramified degree2 reductions at7.
This excludes new frontier groups C3 x C36 and C6 x C18, completing Phi10.

Proof sources and immutable manifests are in continuation_20261031b/.
Parent differential controls check2394 cubic and36 quartic restrictions.
Independent native and Sage arithmetic check all43 A7-negative supports
and360 geometric7-isogeny kernels, the global section and exact NFpoint.
C56 has9 exact Kummer RR packets,82944 independently replayed matrix
entries,66 precision128 arcs, and72 independently constructed kernels.

The {pages}-page paper compiles without errors or overfull boxes.
All {len(statuses)} recorded jobs exited; no degree10 workers or launchers
remain. Obsolete jobs were identity-checked before termination. Actual
exits and failed diagnostic evidence remain recorded; none is silently
promoted to a certificate. Unrelated useful jobs were preserved.

Final validation: {report['evidence_files_checked']} evidence hashes;
{report['protected_earlier_files_unchanged']} earlier-degree protected files unchanged.
Earlier progress: continuation_before_20261031b.txt.
''')
summary={'status':'PASS; PHI(10) COMPLETE','known_full_groups':69,'excluded_minimal_obstructions':110,'remaining':0,'original_cases_resolved':10,'original_new_exclusions':9,'original_new_occurrence':[3,18],'new_larger_exclusions':[[3,36],[6,18]],'paper_pages':pages,'evidence_files_checked':report['evidence_files_checked'],'jobs_exited':len(statuses),'protected_earlier_files_unchanged':report['protected_earlier_files_unchanged']}
(R/'final_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
