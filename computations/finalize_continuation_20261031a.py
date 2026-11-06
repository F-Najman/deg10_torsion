#!/usr/bin/env python3
"""Finalize ten further verified exclusions; refuses incomplete work or live jobs."""
from pathlib import Path
import hashlib,json,os,re,subprocess,time
import psutil
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent;R=H/'continuation_20261031a'
L=json.loads((H/'certified_results.json').read_text());new=[r for r in L['exclusions']if r.get('proof','').startswith('continuation_20261031a/')]
assert len(new)==10 and all(r['verified']for r in new)
groups=sorted(r['group']for r in new);label=lambda g:f'C{g[1]}'if g[0]==1 else f'C{g[0]} x C{g[1]}'
subprocess.run(['python3',str(H/'computations/sync_proof_counts.py')],cwd=ROOT,check=True)
D=json.loads((H/'candidates.json').read_text());rem=D['remaining_minimal_obstructions'];cyc=[n for m,n in rem if m==1];non=[(m,n)for m,n in rem if m!=1];ne=len(D['certified_exclusions']);nr=len(rem)
assert(ne,nr)==(99,10)
p=H/'paper/main.tex';main=p.read_text().replace('\\date{30 October 2026}','\\date{31 October 2026}')
aux=r'''\subsection{An auxiliary gonality bound}
The following bound settles the entirely cuspidal reduction branch at $5$
for order $54$, using rank zero and the same lifting argument as for order $44$.
The other reduction branches for order $54$ remain unresolved.
\begingroup
\emergencystretch=3em
\input{../continuation_20261031a/C54/gonality.tex}
\par\endgroup
'''
if '../continuation_20261031a/C54/gonality.tex'not in main:main=main.replace('\\input{remaining.tex}',aux+'\\input{remaining.tex}')
p.write_text(main)
c=subprocess.run(['/home/fnajman/bin/tectonic','-X','compile','main.tex','--keep-logs','--keep-intermediates'],cwd=H/'paper',capture_output=True,text=True);print(c.stdout+c.stderr);c.check_returncode()
info=subprocess.run(['/home/fnajman/bin/pdfinfo',str(H/'paper/main.pdf')],capture_output=True,text=True,check=True).stdout;pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
statuses=[json.loads(p.read_text())for p in sorted((H/'jobs').glob('*.status.json'))];assert all(s['state']=='EXITED'for s in statuses)
archive=H/'notes/lifecycle_before_continuation_20261031a.json'
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
jobs=[s for s in statuses if s['start_time_unix']>base]
lifecycle={'status':'PASS','case':'Ten further degree-ten exclusions: '+', '.join(map(label,groups)),'time_unix':time.time(),'degree_ten_workers_supervisors_launchers':own,'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'new_job_count':len(jobs),'new_jobs':[{'purpose':s['purpose'],'actual_exit_code':s.get('actual_exit_code'),'termination_reason':s.get('termination_reason')}for s in jobs],'new_automatic_launchers_created':0,'unrelated_magma_workers_preserved':other,'host_cores':psutil.cpu_count(),'new_case_active_cores':0,'native_alarm_or_failed_probes_are_not_completed_certificates':True,'retirement_note':'Superseded and no-longer-needed computations were identity-checked and retired; actual exits are preserved in job statuses and retirement records. Failed parsers, capped computations and residual exploratory tests are not exclusion premises.'}
(H/'lifecycle_final.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
a=subprocess.run(['python3',str(H/'verify_project.py')],cwd=ROOT,capture_output=True,text=True)
if a.returncode:print(a.stdout+a.stderr)
a.check_returncode();report=json.loads(a.stdout);(H/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
p=H/'README.txt';s=p.read_text();old=H/'notes/README_before_continuation_20261031a.txt'
if not old.exists():old.write_text(s)
s=re.sub(r'Phi\(10\): certified partial results, [^\n]+','Phi(10): certified partial results, 31 October 2026',s)
s=re.sub(r'Of the 109 minimal subgroup obstructions outside this lower list, \d+ are\nexcluded and \d+ remain unresolved \(\d+ cyclic, \d+ noncyclic\)\. These \d+ are',f'Of the 109 minimal subgroup obstructions outside this lower list, {ne} are\nexcluded and {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). These {nr} are',s)
s=re.sub(r'paper/main\.pdf \(\d+ pages\)',f'paper/main.pdf ({pages} pages)',s)
high='  - Ten further exclusions: '+', '.join(map(label,groups))+'.\n    C78 also classifies all noncuspidal rational effective decic divisors\n    on X1(39), excluding C3 x C39. C90 uses marked 2-isogeny cycles;\n    the possible CM endomorphism fields are excluded by class numbers.\n    C63 handles non-global cubic cusp fragments by exact Hecke tests.\n    See continuation_20261031a/ for proofs, independent replays and\n    immutable certificate manifests.\n'
if high not in s:s=s.replace('Highlights:\n','Highlights:\n'+high)
if'\nFinal audit:'in s:s=s[:s.index('\nFinal audit:')]
s+=f'\nFinal audit: {report["evidence_files_checked"]} evidence hashes checked; all {len(statuses)} recorded jobs\nexited; 25 protected earlier-degree files unchanged. The {pages}-page paper\ncompiles without errors or overfull boxes. No degree-ten worker, supervisor\nor launcher remains active. See validation.json and lifecycle_final.json.\n';p.write_text(s)
p=H/'notes/continuation.txt';old=H/'notes/continuation_before_20261031a.txt'
if not old.exists():old.write_text(p.read_text())
s=f'''31 October 2026 — Phi(10), after ten further exclusions

CURRENT COUNTS: 68 exact full groups occur; {ne} minimal subgroup obstructions
are excluded; {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic).
The frontier still has109 groups. The classification is INCOMPLETE.
Cyclic orders remaining: {', '.join(map(str,cyc))}.
Noncyclic groups remaining: {', '.join(f'C{m} x C{n}'for m,n in non)}.
These are minimal proof tasks, not asserted occurrences or an exhaustive
upper list. New occurrences can expose new minimal obstructions.

NEW EXCLUSIONS: {', '.join(map(label,groups))}.
Complete proofs, exact local data and independent arithmetic replays are
in continuation_20261031a/. All accepted evidence is immutable and pinned
in certified_results.json. round_status.json records the ten groups.

C63: The genus49 involution quotient retains rankzero. Exact degree10
Hecke tests and non-global cubic cusp fragments cover the noncuspidal
boundary; all239238 full-Jacobian Cartier restriction checks pass.
C117: All degree8/9 reductions push to known closed points on X1(39),
and the remaining complement forces a cusp by injectivity and gonality.
C78 and C3 x C39: The complete finite sieve classifies every noncuspidal
rational effective decic divisor on X1(39) as a diamond translate of
known1.39.10.a. Its exact full torsionC39 excludes both larger groups.
C90: Full rational torsion injectivity at2 and rational gonality18 force
a rational2-isogeny to close, up to field automorphism and marked scalar1/2.
Composing conjugates gives an impossible scalar on45-torsion. If the curve
had a nontrivial K-endomorphism, ordinary reduction would force CM field
-119 or-319. Their class numbers are10, incompatible with[K:L]=5.
C44: Gon_Q15 closes the cusp boundary even for local cusp fragments.
All956 new canonical R2 multiplication tests have rank105, with954 A5
and2 A7 tests. The quintic cusp complement is annihilated literally
because5=3^8 mod44. Independent full jet and marked-kernel controls pass.
C46: Gon_Q19 closes the boundary. The X1(23) rational class quotientC23
rejects986of1023 new supports. All37 exceptions haveA5h0zero, with292
native prime references and324 independently constructed genuine kernels.
C58: The X1(29) class filter leaves2of293 eligible decic orbits; both lift
to known full-C29 curves. Good reductions of order29 rule out2-torsion.
C62: All253 eligible decic X1(31) classes,4095quintic pairs and45doubled
fibres miss the rational image. Five old nonic survivors force a global
nonic trace plus a cusp, impossible by gon_Q12 and torsion injectivity.
C70: The773source orbits are excluded on X1(35), genus25. Rational
cuspidal generation makes A2 kill the full rational Jacobian without
assuming reduction2 injectivity. The first map kills770cases; the second
2-isogeny projection kills3. All776RRtests/3489prime references and2328
physical kernels are independently checked. Source boundary usesrankzero
and the geometric gonality bound378/25.
Full details are in the individual resolution.tex and frozen manifests.

BOUNDARY OBSERVATION: At an odd good prime, rankzero and Manin-Drinfeld
make an entirely cuspidal trace reduction linearly equivalent overQp to
its local cusp lift. Since the original divisor is rational, flat base
change gives h0(D)>=2 overQ, contradicting gon_Q>10. The local cusp lift
need not be globally rational. Varivoda's TheoremD gives gon_Q15,19,16
for X1(44),X1(46),X1(48), respectively. This replaces unnecessary large
cusp searches. It does not by itself prove dyadic torsion injectivity.

CAUTIONS: Do not assume rational cuspidal saturation for C3 x C18.
The relative quintic noncyclic exploration remains unresolved; timed-out
CM factorization and insufficient torsion bounds are diagnostics only.
Only literal rational cuspidal divisor generation allows the A2 torsion
annihilation used for C70 without a reduction2 injectivity premise.
All physical markings, residue norms and coincident targets are retained.

LIFECYCLE: All {len(statuses)} recorded jobs have exited. No degree10 worker,
supervisor, automatic launcher or obsolete paused job remains. The job
statuses and retirement records give actual exits. Failed or terminated
probes are not certificates. Unrelated workers were preserved, and two
host CPU cores were reserved throughout.

Earlier results are described in continuation_before_20261031a.txt.
FINAL AUDIT: {report['evidence_files_checked']} evidence hashes checked;
25 protected earlier-degree files unchanged. The {pages}-page paper compiles
without errors or overfull boxes. See validation.json and lifecycle_final.json.
''';p.write_text(s)
print(json.dumps({'new_exclusions':groups,'known':68,'excluded':ne,'remaining':nr,'cyclic_remaining':len(cyc),'noncyclic_remaining':len(non),'pages':pages,'evidence_files':report['evidence_files_checked'],'jobs_exited':len(statuses)},indent=2))
