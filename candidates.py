#!/usr/bin/env python3
"""Exact minimal-obstruction enumeration for Phi(10); no search is a proof."""
from pathlib import Path
from math import gcd
import json,hashlib
H=Path(__file__).resolve().parent

def divisors(n): return [d for d in range(1,n+1) if n%d==0]
def factors(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   out.append(p)
   while n%p==0:n//=p
  p+=1
 if n>1:out.append(n)
 return out

def subgroups(g):return {(a,b) for a in divisors(g[0]) for b in divisors(g[1]) if b%a==0}
def subgroup(a,b):return b[0]%a[0]==0 and b[1]%a[1]==0

def main():
 ledger=json.loads((H/'certified_results.json').read_text())
 primes=set(ledger['prime_orders']['set']);ms=ledger['first_invariant_factors']
 known=set(map(tuple,ledger['inherited_occurrences']['groups']))
 for r in ledger['realizations']:
  assert r['verified'] is True
  for e in r['evidence']:
   assert hashlib.sha256((H/e['path']).read_bytes()).hexdigest()==e['sha256'],e['path']
  known.add(tuple(r['group']))
 assert all(subgroups(g)<=known for g in known),'Known list must include independently established full-torsion subgroups.'
 for r in ledger['exclusions']:
  if r.get('verified') is True:
   for e in r['evidence']:
    assert hashlib.sha256((H/e['path']).read_bytes()).hexdigest()==e['sha256'],e['path']
 excluded={tuple(r['group']) for r in ledger['exclusions'] if r.get('verified') is True}
 assert not any(subgroup(e,g) for e in excluded for g in known)
 cyclic={n for m,n in known if m==1}
 trial=primes|{n*p for n in cyclic for p in primes}
 frontier={(1,n) for n in trial if (1,n) not in known and all((1,n//p) in known for p in factors(n))}
 for m in ms:
  if m==1:continue
  for n in cyclic:
   g=(m,n)
   if n%m or g in known:continue
   if subgroups(g)-{g}<=known:frontier.add(g)
 # Independent bounded enumeration; a minimal missing cyclic n has n/p in known.
 independent=[]
 for m in ms:
  for n in range(1,max(cyclic)*max(primes)+1) if m==1 else sorted(cyclic):
   g=(m,n)
   if n%m or g in known or not set(factors(n))<=primes:continue
   if subgroups(g)-{g}<=known:independent.append(g)
 assert sorted(frontier)==sorted(independent)
 covered=sorted(g for g in frontier if any(subgroup(e,g) for e in excluded))
 remaining=sorted(frontier-set(covered))
 result={'degree':10,'status':'COMPLETE' if not remaining else 'INCOMPLETE RESEARCH','classification_complete':not remaining,'prime_orders':sorted(primes),'first_invariant_factors':ms,'known_lower_bound':sorted(known),'minimal_outside_known':sorted(frontier),'certified_exclusions':sorted(excluded),'eliminated_minimal_obstructions':covered,'remaining_minimal_obstructions':remaining,'classification':sorted(known) if not remaining else None,'independent_frontier_enumeration':'PASS','warning':'Unresolved minimal obstructions are proof tasks, not asserted occurrences or a complete upper list. A new occurrence can expose further obstructions.'}
 (H/'candidates.json').write_text(json.dumps(result,indent=2)+'\n')
 label=lambda g:f'C{g[1]}' if g[0]==1 else f'C{g[0]} x C{g[1]}'
 text='Phi(10): '+result['status']+'\nKnown exact full torsion groups: '+str(len(known))+'\nMinimal obstructions: '+str(len(frontier))+'\nEliminated: '+str(len(covered))+'\nRemaining: '+str(len(remaining))+'\n\nKnown groups:\n'+', '.join(map(label,sorted(known)))+'\n\nRemaining proof tasks:\n'+', '.join(map(label,remaining))+'\n\n'+result['warning']+'\n'
 (H/'candidates.txt').write_text(text)
 print(text)
if __name__=='__main__':main()
