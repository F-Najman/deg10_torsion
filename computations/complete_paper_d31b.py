"""Rewrite the conclusion only after the certified ledger proves completeness."""
from pathlib import Path
import json
H=Path(__file__).resolve().parents[1];d=json.loads((H/'candidates.json').read_text())
assert d['classification_complete'] and len(d['known_lower_bound'])==69 and len(d['certified_exclusions'])==110 and not d['remaining_minimal_obstructions']
p=H/'paper/main.tex';s=p.read_text()
if '\\newtheorem{theorem}[proposition]{Theorem}' not in s:s=s.replace('\\newtheorem{proposition}{Proposition}[section]','\\newtheorem{proposition}{Proposition}[section]\n\\newtheorem{theorem}[proposition]{Theorem}')
s=s.replace('Torsion over decic number fields: certified partial results','Torsion of elliptic curves over decic number fields')
a=s.index('\\begin{abstract}');b=s.index('\\end{abstract}',a)+len('\\end{abstract}')
s=s[:a]+r'''\begin{abstract}
We determine the possible torsion subgroups of elliptic curves over
number fields of degree $10$, using the known set $S(10)$.
There are $69$ groups. In particular, $C_3\times C_{18}$ occurs
at isolated degree-$5$ points of $X_1(3,18)$ over $\mathbb Q(\zeta_3)$.
We give exact certificates for the occurrences and exclude all $110$
minimal subgroup obstructions outside the list.
\end{abstract}'''+s[b:]
s=s.replace('\\section{The certified lower list}','\\section{The classification}')
s=s.replace('\\begin{proposition}\\label{prop:lower}\nEvery group in the following list occurs as full torsion over a decic field:','\\begin{theorem}\\label{thm:classification}\\label{prop:lower}\nThe set $\\Phi(10)$ consists of the following groups:')
s=s.replace('&C_3\\times C_{3n} &&1\\leq n\\leq5;','&C_3\\times C_{3n} &&1\\leq n\\leq6;')
s=s.replace('There are $68$ groups in this list.\n\\end{proposition}','There are $69$ groups in this list.\n\\end{theorem}')
s=s.replace('These disjoint contributions have total size $38+15+10+2+1+1+1=68$.',r'''Proposition~\ref{prop:decic318-classification} supplies $C_3\times C_{18}$.
These disjoint contributions have total size $38+15+10+2+1+1+1+1=69$.
The exclusions below eliminate every minimal subgroup obstruction outside
this list. The independent enumeration in Section~\ref{sec:complete-enumeration}
then proves that no further group can occur.''')
completion_block=r'''\subsection{The degree-$5$ points on $X_1(3,18)$}
\begingroup
\emergencystretch=3em
\input{../continuation_20261031b/M3N18/resolution.tex}
\par\endgroup
\input{enumeration.tex}'''
if '../continuation_20261031b/M3N18/resolution.tex' in s:s=s.replace('\\input{remaining.tex}','\\input{enumeration.tex}')
else:s=s.replace('\\input{remaining.tex}',completion_block)
p.write_text(s)
(H/'paper/enumeration.tex').write_text(r'''\section{Completeness of the enumeration}\label{sec:complete-enumeration}
Let $\mathcal L$ be the list in Theorem~\ref{thm:classification}.
If a finite torsion group $G$ lies outside $\mathcal L$, choose a subgroup
$H\leq G$ minimal outside $\mathcal L$.
By the Weil pairing, its first invariant factor belongs to
$\{1,2,3,4,6,11,22\}$, and its prime divisors
belong to the known set $S(10)$.

There are exactly $110$ such minimal groups. All are excluded above.
Here is a finite bound for the independent enumeration.
If $C_n$ is minimal outside $\mathcal L$, then $C_{n/p}\in\mathcal L$
for every prime $p\mid n$. Since the largest known cyclic order is $50$,
we have $n\leq50\cdot37$. For a minimal noncyclic group $C_m\times C_n$,
its cyclic subgroup $C_n$ already belongs to $\mathcal L$, so $n\leq50$.
The program checks every group within these bounds and every proper
subgroup. Its output agrees with a second enumeration obtained by
multiplying the known cyclic orders by the allowed primes.

The occurrence of $C_3\times C_{18}$ adds two minimal groups,
$C_3\times C_{36}$ and $C_6\times C_{18}$, to the earlier frontier.
The complete degree-$5$ point calculation on $X_1(3,18)$ shows that all
resulting elliptic curves have exact full torsion $C_3\times C_{18}$.
It therefore excludes both larger groups. The final frontier has
$110$ groups, all with verified exclusion certificates.
''')
print('Completed paper statement:69 groups,110 exclusions')
