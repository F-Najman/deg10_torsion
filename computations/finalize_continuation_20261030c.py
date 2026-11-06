#!/usr/bin/env python3
"""Finalize six further independently verified cyclic decic exclusions."""
from pathlib import Path
import json,os,re,time,subprocess
import psutil
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;newcases=[65,85,91,93,99,105]
L=json.loads((H/'certified_results.json').read_text())
assert all(any(r['group']==[1,n] and r['verified'] for r in L['exclusions'])for n in newcases)
subprocess.run(['python3',str(H/'computations/sync_proof_counts.py')],cwd=ROOT,check=True)
D=json.loads((H/'candidates.json').read_text());rem=D['remaining_minimal_obstructions'];cyc=[n for m,n in rem if m==1];non=[(m,n)for m,n in rem if m!=1];ne=len(D['certified_exclusions']);nr=len(rem)
assert(ne,nr,len(cyc),len(non))==(89,20,15,5)
c=subprocess.run(['/home/fnajman/bin/tectonic','-X','compile','main.tex','--keep-logs','--keep-intermediates'],cwd=H/'paper',capture_output=True,text=True);print(c.stdout+c.stderr);c.check_returncode()
info=subprocess.run(['/home/fnajman/bin/pdfinfo',str(H/'paper/main.pdf')],capture_output=True,text=True,check=True).stdout;pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
statuses=[json.loads(p.read_text())for p in sorted((H/'jobs').glob('*.status.json'))];assert all(s['state']=='EXITED'for s in statuses)
archive=H/'notes/lifecycle_before_continuation_20261030c.json'
if not archive.exists():archive.write_text((H/'lifecycle_final.json').read_text())
base=json.loads(archive.read_text())['time_unix'];own=[];other=[]
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
new=[s for s in statuses if s['start_time_unix']>base]
lifecycle={'status':'PASS','case':'C65, C85, C91, C93, C99 and C105 exclusions completed','time_unix':time.time(),'degree_ten_workers_supervisors_launchers':own,'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'new_job_count':len(new),'new_jobs':[{'purpose':s['purpose'],'actual_exit_code':s.get('actual_exit_code'),'termination_reason':s.get('termination_reason')}for s in new],'new_automatic_launchers_created':0,'unrelated_magma_workers_preserved':other,'host_cores':psutil.cpu_count(),'new_case_active_cores':0,'native_alarm_or_failed_probes_are_not_completed_certificates':True,'retirement_note':'The superseded C91 order4 quotient fallback was identity-checked and terminated after the genus23 replacement passed. Its actual exit is recorded. Root source-level backups all exited. Failed model and parser trials are preserved as diagnostics.'}
(H/'lifecycle_final.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
a=subprocess.run(['python3',str(H/'verify_project.py')],cwd=ROOT,capture_output=True,text=True)
if a.returncode:print(a.stdout+a.stderr)
a.check_returncode();report=json.loads(a.stdout);(H/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
p=H/'README.txt';s=p.read_text();old=H/'notes/README_before_continuation_20261030c.txt'
if not old.exists():old.write_text(s)
s=re.sub(r'Of the 109 minimal subgroup obstructions outside this lower list, \d+ are\nexcluded and \d+ remain unresolved \(\d+ cyclic, \d+ noncyclic\)\. These \d+ are',f'Of the 109 minimal subgroup obstructions outside this lower list, {ne} are\nexcluded and {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). These {nr} are',s)
s=re.sub(r'paper/main\.pdf \(\d+ pages\)',f'paper/main.pdf ({pages} pages)',s)
high='''  - Six further exclusions: C65, C85, C91, C93, C99 and C105.
    C91 uses a relative formal-fibre argument on a genus23 quotient:
    the exact fibre consists of a fixed residual divisor plus pullbacks
    from X0(91), and the fixed part contains a cusp. All120 cusp and12
    mixed CM ranks were replayed independently. C99 needs a genus29
    quotient retaining the9-level to exclude its final Hecke class.
    C65 has7170048 independently checked cusp ranks. C85/C93 use the
    convention-correct full-Jacobian local operator T2mu-1, with genuine
    nonrational cusp arcs. C105 is resolved by six genus25 tests.
    See continuation_20261030c/C{65,85,91,93,99,105}/resolution.tex.
'''
if high not in s:s=s.replace('Highlights:\n','Highlights:\n'+high)
if'\nFinal audit:'in s:s=s[:s.index('\nFinal audit:')]
s+=f'\nFinal audit: {report["evidence_files_checked"]} evidence hashes checked; all {len(statuses)} recorded jobs\nexited; 25 protected earlier-degree files unchanged. The {pages}-page paper\ncompiles without errors or overfull boxes. No degree-ten worker, supervisor\nor launcher remains active. See validation.json and lifecycle_final.json.\n';p.write_text(s)
p=H/'notes/continuation.txt';old=H/'notes/continuation_before_20261030c.txt'
if not old.exists():old.write_text(p.read_text())
s=f'''30 October 2026 — Phi(10), after six further exclusions

CURRENT COUNTS: 68 exact full groups occur; {ne} minimal subgroup obstructions
are excluded; {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic).
The frontier still has109 groups. The classification is INCOMPLETE.
Cyclic orders remaining: {', '.join(map(str,cyc))}.
Noncyclic groups remaining: {', '.join(f'C{m} x C{n}'for m,n in non)}.
These are minimal proof tasks, not asserted occurrences or an exhaustive
upper list. New occurrences can expose new minimal obstructions.

NEW EXCLUSIONS: C65, C85, C91, C93, C99 and C105. Complete proofs and local
immutable snapshots are in continuation_20261030c/C65, C85, C91, C93, C99
and C105. All accepted evidence is pinned in certified_results.json.

C65: Five new ordinary decic diamond orbits, sizes24,24,24,12,24; one has
stabilizer2, handled explicitly. All240 canonical multiplication tests
have rank121, using genuine65-division deformations and all20 actual
3-isogenies. All7170048 normalized global cusp prefixes have rank10;
a separate C++ full-rank-at-leaf replay with another F16 functional basis
confirms every result. The split-sextic-cusp genus9 filter extends because
the degree4 complement, including quartic cusps, is globally intact.
All old degree6-9 tests and the entire inherited packet were replayed.

C85: Four new decic source orbits are excluded by independently verified
canonical multiplication. Exact integral k2(T2),k3(T3) have finite image;
the correct multiplicative-level local killer is T2mu-1. Rational prefix
relations have minimum cost32. All2330064 mixed quartic/octic cusp
patterns have rank10. Genuine Tate arcs, alternate-minor reconstruction,
Frobenius descent and independent full matrix checks agree.

C91: The genus23 normalization X1(7) x X0(13) retains an order3 diamond.
Its finite-image dyadic space U has dimension14. The120 admissible cusp
images have rank10-k, with k exactly the degree of their maximal
pullback from X0(91). A rigorous formal implicit-function argument proves
the local fibre is the embedded smooth family C0+pi*Sym^k. The fixed
residual C0 is a nonempty cusp divisor, contradicting noncuspidality.
The CM8 images are4P2 at one of two tame quadratic fixed places; all12
mixed tests have rank8, and the same lemma forces their degree2 cusp
complement. Parent native replay rebuilds all6 prime ideals,45 kernels,
Cartier/deck matrices and132 rank checks; Phi13 verifies the j-map.
The old288 canonical nonprincipality tests exclude ordinary degree9+1.
The discarded order4 quotient and full-source backups are not premises.

C93: No new decic noncuspidal reduction occurs; the old degree9 test
extends with a rational cusp complement. In the boundary, h32-stability
would descend31 to degree5, impossible by the accepted Phi5 list. This
excludes all quadratic/decic cusp supports after overlap cancellation.
Remaining rational/quintic supports use full-Jacobian T2mu-1 diamond
operators. Rational relation cost30; all278262 mixed patterns rank10.
Independent genuine arcs, Frobenius and full matrix controls pass.

C99: Pure cusp reduction forces delta2-stability by a doubled degree40
comparison below gonality. Integral j=0/1728 cannot specialize to a cusp,
so the projective action is free of order15, impossible in total degree10.
The four decic noncuspidal orbits first push to rank-zero X1(33), genus21.
Three A5 tests are nonprincipal; residual3 is genuinely principal there,
and neither A7 nor the additional3-isogeny map resolved it. These are
retained residual diagnostics. The genus29 quotient X1(9) x X0(11) retains
the necessary9-level. Its exact (delta2-1)(T5-1-5delta5) divisor has signed
degrees110 and h0=0. Fresh native ideal/RR replay and independent actual
5/11-isogeny keys verify every marked point and residue norm.

C105: Only six decic source orbits. All push to X1(35), genus25; every A3
divisor is nonprincipal, with fresh native RR and genuine Velu controls.
The finite-image order4 diamond22 excludes pure cusps: Abramovich gives
gonality>=40.32>40 and its action is free, whereas4 does not divide10.

CAUTIONS: The complete accepted Phi5 list includes BOTH13 and17. Early
memory-based suggested quintic obstructions for C85/C91 were rejected
and are NOT premises. Only the valid quintic exclusion of31 is used.
On the canonical multiplicative-level model the dyadic kernel killer is
T2mu-1, not T2mu-delta2; finite-image polynomials must be justified overZ
before reducing modulo2. For relative formal fibres, tangent rank alone
is insufficient: the embedded smooth pullback family and its closed
formal immersion are essential. Do not assume cuspidal saturation for
C3 x C18. Remaining noncyclic relative-quintic work is still unresolved.

LIFECYCLE: All {len(statuses)} recorded jobs have exited. No degree10 worker,
supervisor, automatic launcher or obsolete paused job remains. The C91
order4 fallback was identity-checked and terminated when the replacement
passed. Root backups completed; failed parsers/models are preserved and
are not certificates, even when Magma returned0. Unrelated Magma workers
were preserved. Continue leaving two host CPU cores idle.

Earlier results are described in continuation_before_20261030c.txt.
FINAL AUDIT: {report['evidence_files_checked']} evidence hashes checked;
25 protected earlier-degree files unchanged. The {pages}-page paper compiles
without errors or overfull boxes. See validation.json and lifecycle_final.json.
''';p.write_text(s)
print(json.dumps({'known':68,'excluded':ne,'remaining':nr,'cyclic_remaining':len(cyc),'noncyclic_remaining':len(non),'pages':pages,'evidence_files':report['evidence_files_checked'],'jobs_exited':len(statuses)},indent=2))
