#!/usr/bin/env python3
"""Finalize the six audited cyclic exclusions without rerunning computations."""
from pathlib import Path
import json, os, re, subprocess, time
import psutil
H=Path(__file__).resolve().parents[1]; ROOT=H.parent.parent
cases=[55,57,75,77,81,87]
L=json.loads((H/'certified_results.json').read_text())
assert all(any(r['group']==[1,n] and r['verified'] for r in L['exclusions']) for n in cases)
subprocess.run(['python3',str(H/'computations/sync_proof_counts.py')],cwd=ROOT,check=True)
D=json.loads((H/'candidates.json').read_text()); rem=D['remaining_minimal_obstructions']
cyc=[n for m,n in rem if m==1]; non=[(m,n) for m,n in rem if m!=1]
ne=len(D['certified_exclusions']); nr=len(rem)
assert (ne,nr,len(cyc),len(non))==(83,26,21,5)
checks={
 'jobs/contb77_parent_native.log':'C77_PARENT_NATIVE_COMPLETE 24',
 'jobs/contb77_parent_moduli_binding.log':'C77_PARENT_T3_MODULI_BINDING_COMPLETE 24',
}
for p,marker in checks.items():assert marker in (H/p).read_text(),(p,marker)
c=subprocess.run(['/home/fnajman/bin/tectonic','-X','compile','main.tex','--keep-logs','--keep-intermediates'],cwd=H/'paper',capture_output=True,text=True)
print(c.stdout+c.stderr);c.check_returncode()
info=subprocess.run(['/home/fnajman/bin/pdfinfo',str(H/'paper/main.pdf')],capture_output=True,text=True,check=True).stdout
pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
statuses=[json.loads(p.read_text()) for p in sorted((H/'jobs').glob('*.status.json'))]
assert all(s['state']=='EXITED' for s in statuses)
archive=H/'notes/lifecycle_before_continuation_20261030b.json'
if not archive.exists():archive.write_text((H/'lifecycle_final.json').read_text())
base=json.loads(archive.read_text())['time_unix']; own=[]; other=[]
for proc in psutil.process_iter(['pid','name','username','create_time','cmdline','status']):
 try:
  a=proc.info;args=a['cmdline'] or [];name=(a['name'] or '').lower();magma='magma' in name
  if a['pid'] in [os.getpid(),os.getppid()]:continue
  scope=any('work/deg10/' in x or str(H) in x for x in args)
  if not scope and (magma or any(x.endswith('.py') for x in args[1:])):
   try:scope=Path(proc.cwd()).is_relative_to(H)
   except (OSError,psutil.Error):pass
  if scope and (magma or any(x.endswith('.py') for x in args[1:])):own.append(a)
  elif magma:other.append(a)
 except psutil.Error:pass
assert not own,own
new=[s for s in statuses if s['start_time_unix']>base]
lifecycle={
 'status':'PASS','case':'C55, C57, C75, C77, C81 and C87 exclusions completed',
 'time_unix':time.time(),'degree_ten_workers_supervisors_launchers':own,
 'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'new_job_count':len(new),
 'new_jobs':[{'purpose':s['purpose'],'actual_exit_code':s.get('actual_exit_code'),'termination_reason':s.get('termination_reason')} for s in new],
 'new_automatic_launchers_created':0,'unrelated_magma_workers_preserved':other,
 'host_cores':psutil.cpu_count(),'new_case_active_cores':0,
 'native_alarm_or_failed_probes_are_not_completed_certificates':True,
 'retirement_records':['continuation_20261030b/C77/absolute_fallback_retired.json','continuation_20261030b/C77/absolute_model_trial_retirement.json'],
 'note':'All primary, replacement and independent computations have ended. Superseded absolute models and failed startup probes are retained as diagnostics and must not restart.'}
(H/'lifecycle_final.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
audit=subprocess.run(['python3',str(H/'verify_project.py')],cwd=ROOT,capture_output=True,text=True)
if audit.returncode: print(audit.stdout+audit.stderr)
audit.check_returncode(); report=json.loads(audit.stdout)
(H/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
p=H/'README.txt';s=p.read_text();a=H/'notes/README_before_continuation_20261030b.txt'
if not a.exists():a.write_text(s)
s=re.sub(r'Of the 109 minimal subgroup obstructions outside this lower list, \d+ are\nexcluded and \d+ remain unresolved \(\d+ cyclic, \d+ noncyclic\)\. These \d+ are',f'Of the 109 minimal subgroup obstructions outside this lower list, {ne} are\nexcluded and {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). These {nr} are',s)
s=re.sub(r'paper/main\.pdf \(\d+ pages\)',f'paper/main.pdf ({pages} pages)',s)
high='''  - C55, C57, C75, C77, C81 and C87 are excluded. Every final test has
    independent exact checks, including the physical isogeny markings and
    residue norm weights. The C57 and C77 proofs use finite Hecke images
    despite positive full-Jacobian rank. C81 uses unique Hecke divisors on
    a genus46 Kummer cover and an impossible finite 2-isogeny-closed j-set.
    The source cusp boundary is covered separately in every case.
    See continuation_20261030b/C{55,57,75,77,81,87}/resolution.tex.
'''
if high not in s:s=s.replace('Highlights:\n','Highlights:\n'+high)
if '\nFinal audit:' in s:s=s[:s.index('\nFinal audit:')]
s+=f'\nFinal audit: {report["evidence_files_checked"]} evidence hashes checked; all {len(statuses)} recorded jobs\nexited; 25 protected earlier-degree files unchanged. The {pages}-page paper\ncompiles without errors or overfull boxes. No degree-ten worker, supervisor\nor launcher remains active. See validation.json and lifecycle_final.json.\n';p.write_text(s)
p=H/'notes/continuation.txt';a=H/'notes/continuation_before_20261030b.txt'
if not a.exists():a.write_text(p.read_text())
s=f'''30 October 2026 — Phi(10), after excluding C55, C57, C75, C77, C81 and C87

CURRENT COUNTS: 68 exact full groups occur; {ne} minimal subgroup obstructions
are excluded; {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). The frontier
still has 109 groups. The classification is INCOMPLETE.
Cyclic orders remaining: {', '.join(map(str,cyc))}.
Noncyclic groups remaining: {', '.join(f'C{m} x C{n}' for m,n in non)}.
These are minimal proof tasks, not asserted occurrences or an exhaustive
upper list. New occurrences can expose new minimal obstructions.

All six new proofs, snapshots and independent checks are in
continuation_20261030b/C55, C57, C75, C77, C81 and C87, respectively.
Their accepted immutable evidence is pinned in certified_results.json.

C55: Source genus81, rank zero. The odd degree5 quotient of genus17
proves injection of the full rational 2-primary torsion at2. All source
cusps lift globally intact, so pure cusps give a forbidden pencil. The
old degree7/9 tests and five new degree10 A3 tests exhaust noncuspidal
reduction. Genuine Velu maps and a fresh native RR replay verify them.

C57: The full source Jacobian has positive rank. Products of Hecke
differences have finite image and kill the rational classes used in the
full-Jacobian dyadic formal-immersion proof. All2383311 normalized global
cusp divisors have rank10; the non-global9+1 cusp case has its own filter.
All required new degree10 Hecke divisors are nonprincipal by canonical
ring multiplication, with genuine local jets and independent Sage tests.
There are12 physical cyclic-subgroup orbits, not6 naive j-orbits.

C75: The genus73 quotient has rank zero. Six new A7 divisors and the
accepted degree7 test exclude noncuspidal reduction. All871630 admissible
cusp images are covered; fresh high-multiplicity ideal restrictions handle
the110 new boundary patterns. Source degree8 cusps split at2, so only the
proved intact quotient cusp lifts are used. Native replay checks all six
new RR tests, the point identifications and all residue norm weights.

C77: The source has positive rank; T2 has finite rational image.
Use the degree5 quotient of genus37 and the full-Jacobian dyadic criterion
for all4992 allowed cusp images. There are90 noncuspidal degree10 places,
three diamond orbits and24 possible quadratic Hecke-lift allocations.
A3T2 gives B=T+2delta(T)-7P-2delta(P)-3delta^-1(P), T=T3P.
Every h0(B) is zero. A separate native replay rebuilds all24 RR problems,
including exact two-generator prime identities and both diamond actions.
A separate moduli-binding replay identifies every positive Hecke prime
with independently verified marked isogeny coordinates and norm weight.
No rank-zero assertion for J1(77) is used.

C81: Rational2-torsion vanishes by exact integral Betti operators, and all
source cusps lift intact. The genus4 X0(81) filter leaves five cases.
On the odd degree3 quotient of genus46, each of the29 possible degree10
reductions E has h0(E+2delta(E))=1. A Kummer cover over genus13 makes this
exact computation feasible. Native completions rebuild87 base RR spaces
and all29 matrices; independent Sage jets agree entrywise. Uniqueness
and semicontinuity force an impossible finite 2-isogeny-closed j-set.
This supersedes the earlier unusable gonality>30 proposal; the argument
does not assert that inequality.

C87: All source cusps are intact. Doubling and gonality exclude pure cusp
reduction. A complete source census, including all cyclic3 subgroups on
curves with full3-torsion, gives12 relevant images on genus9 X0(87).
All twelve A5 tests are nonprincipal; independent Velu maps and native
RR tests check the exact norms and marked input points. No earlier
canonical-ring degree9 exclusion is transferred.

Inherited geometry and accepted old tests are copied into the new degree10
folders and hashed. No earlier-degree source was edited. TwoGenerators
exports are used on positive-genus relative fields; Basis(Ideal(P)) alone
would omit coefficient ideals. Magma exit0 is not sufficient evidence:
explicit success markers and independent proof checks are required.
Failed parsers, failed allocations and superseded models remain diagnostics.

Earlier accepted work is described in continuation_before_20261030b.txt.
Do not assume cuspidal saturation for C3 x C18. Relative quintic work
remains for C3 x C27, C3 x C33, C3 x C39 and C4 x C28.

LIFECYCLE: All {len(statuses)} recorded jobs have exited. No degree10 worker,
supervisor, automatic launcher or obsolete paused job remains. Superseded
C77 absolute models were explicitly retired after verifying owner, PID,
start time and command; their actual exits are retained. Failed startup
probes were checked for orphan workers and cleaned up. Unrelated Magma
processes were preserved. Continue leaving two host CPU cores idle.

FINAL AUDIT: {report['evidence_files_checked']} evidence hashes checked;
25 protected earlier-degree files unchanged. The {pages}-page paper compiles
without errors or overfull boxes. See validation.json and lifecycle_final.json.
''';p.write_text(s)
print(json.dumps({'known':68,'excluded':ne,'remaining':nr,'cyclic_remaining':len(cyc),'noncyclic_remaining':len(non),'pages':pages,'evidence_files':report['evidence_files_checked'],'jobs_exited':len(statuses)},indent=2))
