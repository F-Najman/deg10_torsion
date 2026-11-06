#!/usr/bin/env python3
"""Synchronize the degree-ten handoff after the four audited exclusions."""
from pathlib import Path
import hashlib,json,os,re,subprocess,time
import psutil
H=Path(__file__).resolve().parents[1];ROOT=H.parent.parent
L=json.loads((H/'certified_results.json').read_text())
assert all(any(r['group']==g and r['verified'] for r in L['exclusions']) for g in [[1,51],[1,69],[1,84],[2,42]])
subprocess.run(['python3',str(H/'candidates.py')],cwd=ROOT,check=True,stdout=subprocess.PIPE,text=True)
d=json.loads((H/'candidates.json').read_text());rem=d['remaining_minimal_obstructions'];cyc=[n for m,n in rem if m==1];non=[(m,n) for m,n in rem if m!=1];ne=len(d['certified_exclusions']);nr=len(rem)
assert (ne,nr,len(cyc),len(non))==(77,32,27,5)
p=H/'paper/main.tex';s=p.read_text();s=re.sub(r'\$\d+\$ minimal subgroup obstructions',lambda m:f'${ne}$ minimal subgroup obstructions',s);s=re.sub(r'\$\d+\$ minimal obstructions\nremain',lambda m:f'${nr}$ minimal obstructions\nremain',s);p.write_text(s)
p=H/'paper/remaining.tex';s=p.read_text();tail=s[s.index('These are minimal unresolved subgroup obstructions.'):]
prefix='\\section{The remaining proof tasks}\nThe certified lower list has $68$ groups. Its minimal subgroup obstructions,\nsubject to $S(10)$ and the Weil pairing, form a list of $109$ groups.\n'+f'The present note excludes ${ne}$ of them, leaving ${nr}$: ${len(cyc)}$ cyclic and\n${len(non)}$ noncyclic. The cyclic orders are\n'
prefix+='\\[\n\\begin{gathered}\n'+'\\\\\n'.join(','.join(map(str,cyc[i:i+10])) for i in range(0,len(cyc),10))+'\n\\end{gathered}\n\\]\nThe noncyclic groups are\n\\[\n\\begin{array}{c|l}\nm&n\\text{ in }C_m\\times C_n\\\\\\hline\n'
for m in sorted({m for m,n in non}):prefix+=str(m)+'&'+','.join(str(n) for mm,n in non if mm==m)+'\\\\\n'
prefix+='\\end{array}\n\\]\n';p.write_text(prefix+tail)
c=subprocess.run(['/home/fnajman/bin/tectonic','-X','compile','main.tex','--keep-logs','--keep-intermediates'],cwd=H/'paper',capture_output=True,text=True);print(c.stdout+c.stderr);c.check_returncode()
info=subprocess.run(['/home/fnajman/bin/pdfinfo',str(H/'paper/main.pdf')],capture_output=True,text=True,check=True).stdout;pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
statuses=[json.loads(p.read_text()) for p in sorted((H/'jobs').glob('*.status.json'))];assert all(s['state']=='EXITED' for s in statuses)
old=json.loads((H/'lifecycle_final.json').read_text());archive=H/'notes/lifecycle_before_C51_C242_C69.json'
if not archive.exists():archive.write_text(json.dumps(old,indent=2)+'\n')
base=json.loads(archive.read_text())['time_unix'];own=[];other=[]
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
lifecycle={'status':'PASS','case':'C51, C69, C2 x C42 and C84 exclusions completed','time_unix':time.time(),'degree_ten_workers_supervisors_launchers':own,'all_recorded_jobs_exited':True,'recorded_job_count':len(statuses),'new_job_count':len(new),'new_jobs':[{'purpose':s['purpose'],'actual_exit_code':s.get('actual_exit_code'),'termination_reason':s.get('termination_reason')} for s in new],'new_automatic_launchers_created':0,'unrelated_magma_workers_preserved':other,'host_cores':psutil.cpu_count(),'new_case_active_cores':0,'native_alarm_or_failed_probes_are_not_completed_certificates':True}
(H/'lifecycle_final.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
audit=subprocess.run(['python3',str(H/'verify_project.py')],cwd=ROOT,capture_output=True,text=True);audit.check_returncode();report=json.loads(audit.stdout);(H/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
p=H/'README.txt';s=p.read_text().replace('73 are\nexcluded and 36 remain unresolved (30 cyclic, 6 noncyclic). These 36 are',f'{ne} are\nexcluded and {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). These {nr} are').replace('paper/main.pdf (38 pages)',f'paper/main.pdf ({pages} pages)')
newhigh='''  - C51, C69, C2 x C42 and C84 are now excluded. The source data,
    quotient maps, Hecke norms and final Riemann-Roch tests were replayed
    independently. C51 uses full-Jacobian dyadic formal immersion on all
    9656144 cusp divisors. C242 ends with a unique-section argument on
    its diagonal genus17 quotient (24 patterns, 72 character spaces).
    C69 retains the 23-level generator on a genus67 quotient; its exact
    Hecke divisor is nonprincipal. See next_cases_20261030/{cyclic_further/C51,
    C242/geometry,C69}/resolution.tex.
'''
s=s.replace('Highlights:\n','Highlights:\n'+newhigh)
if '\nFinal audit:' in s:s=s[:s.index('\nFinal audit:')]
s+=f'\nFinal audit: {report["evidence_files_checked"]} evidence hashes checked; all {len(statuses)} recorded jobs\nexited; 25 protected earlier-degree files unchanged. The {pages}-page paper\ncompiles without errors or overfull boxes. No degree-ten worker, supervisor\nor launcher remains active. See validation.json and lifecycle_final.json.\n';p.write_text(s)
p=H/'notes/continuation.txt';archive=H/'notes/continuation_before_C51_C242_C69.txt'
if not archive.exists():archive.write_text(p.read_text())
s=f'''30 October 2026 — Phi(10), after excluding C51, C69, C2 x C42 and C84

CURRENT COUNTS: 68 exact full groups occur; {ne} minimal subgroup obstructions
are excluded; {nr} remain unresolved ({len(cyc)} cyclic, {len(non)} noncyclic). The frontier
still has 109 groups. The classification is INCOMPLETE.
Use candidates.json/txt and paper/remaining.tex for the current lists.

C51 IS EXCLUDED. The complete finite A5 sieve excludes all nine noncuspidal
orbits. The full-Jacobian Q2 formal-immersion criterion handles all 9656144
cuspidal degree-ten divisors. Fresh native cusp jets and a uniquely determined
Cartier intertwiner with an 8P kernel identify the two differential bases.
The exact proof and evidence are in next_cases_20261030/cyclic_further/C51/.

C2 x C42 AND C84 ARE EXCLUDED. A fresh source census covers all 60971
normalized support patterns. The filters leave 4647, 752, 111, 108, 10,
then one pattern. Nine final Hecke exclusions have independent kernel,
Velu, norm and native ideal/RR checks. Three order-three stabilizer cases
use the CM lemma that j=0 with C2 x C42 requires field degree at least12.
The last support has noncuspidal degree6 and a degree4 cusp complement.
An exact parity test forces cubic plus rational cusps. Both quadratic
candidates on the genus17 quotient and all12 cusp types give24 patterns;
all72 Kummer character spaces have dimensions [1,0,0]. Unique sections
force an invariant degree4 cusp divisor for a free order-three action,
which is impossible. Generator h has Frobenius phase2; the old model's
generator is h^2 with phase4. The complete proof/evidence packet is
next_cases_20261030/C242/geometry/. All556 manifest hashes were checked.

C69 IS EXCLUDED. Noncuspidal closed points first occur in degree7; the
counts through10 are 0,0,0,0,0,0,66,44,44,44. Every relevant small cusp
orbit is globally intact. Entirely cuspidal reduction gives a degree20
pencil, below gonality10725/512. Complete genus7 and genus13 Hecke tests
exclude every other branch except one degree10 diamond orbit. The
normalization of X1(23) x X0(3), a genus67 quotient, resolves this orbit.
All four actual input/diamond/Hecke lifts have unique same-degree linear
factors in a 1+3 quartic factorization. The two Hecke image degrees are
10 and50, with orbit sizes1 and5 and norm weights1 and1. The correct
T5-1-5<5> divisor has h0=0 in both primary and fresh native calculations.
See next_cases_20261030/C69/resolution.tex and parent_final_checks.json.
The small-cover file small_cover_verified_partial.json is historical
stage evidence; its residual is resolved by the new genus67 test.

Important serialization point: for a function field relative to a base
of positive genus, Basis(Ideal(P)) omits coefficient ideals. The accepted
C69 exports use TwoGenerators(P), exact ideal equality and serialization
roundtrips. Fresh native RR replays verify the corrected exports. A known
Sprint metadata tag is declared explicitly by the final reader; no numeric
proof data were altered. Failed diagnostic scripts and their actual exits
are preserved, and are not accepted as completed certificates.

All four exclusions are in certified_results.json and the paper. The
C49/C95/C111 and earlier noncyclic exclusions remain unchanged. Earlier
handoff details are in continuation_before_C51_C242_C69.txt.

Remaining cautions: C81's recorded Hecke method needs gonality greater
than30; the verified bound236925/8192 is only about28.92. Do not assume
cuspidal saturation for C3 x C18. Relative quintic work remains for
C3 x C27, C3 x C33, C3 x C39 and C4 x C28. The remaining obstructions
are proof tasks, not asserted occurrences or a full upper list.

LIFECYCLE: All {len(statuses)} recorded jobs have exited. No new automatic
launcher or paused worker remains. The full genus133 C69 attempt and the
first genus67 discovery stopped at their native alarms; two V9 quotient
profiles were bounded timeouts. These were superseded by complete verified
computations and must not be restarted automatically. Unrelated processes
were preserved. Continue leaving two host CPU cores idle.

FINAL AUDIT: {report['evidence_files_checked']} evidence hashes checked;
25 protected earlier-degree files unchanged. The {pages}-page paper compiles
without errors or overfull boxes. See validation.json and lifecycle_final.json.
'''
p.write_text(s)
print(json.dumps({'known':68,'excluded':ne,'remaining':nr,'cyclic_remaining':len(cyc),'noncyclic_remaining':len(non),'pages':pages,'evidence_files':report['evidence_files_checked'],'jobs_exited':len(statuses)},indent=2))
