#!/usr/bin/env python3
"""Write the proof table from the exact accepted numerical packet."""
from pathlib import Path
import json
H=Path(__file__).resolve().parent;D=json.loads((H/'large_global_certificate.json').read_text());rr=D['rows']
names=','.join(str(r['N']) for r in rr)
text=r'''\begin{proposition}\label{prop:nonic-large-global}
No elliptic curve over a number field of degree $9$ has a point of order
\[
'''+names+r'''.
\]
\end{proposition}
\begin{proof}
Put $I_N=[\operatorname{PSL}_2(\mathbb Z):\overline{\Gamma_1(N)}]$
and $G_N=325I_N/2^{15}$. The Abramovich--Kim--Sarnak bound gives
$\operatorname{gon}_{\mathbb C}X_1(N)\geq G_N$.
For each row below we prove that
$(\langle a\rangle-1)J_1(N)(\mathbb Q)$ is finite.
There are two sources of rank information.

For the rows using archived newforms, every factor of level dividing
$N$ is included with its full oldform multiplicity. The dimension sum
is exactly $g(X_1(N))$, and every analytic-rank assertion is marked
proved. Every positive-rank character takes the value $1$ at $a$;
these values are recomputed exactly in
\texttt{realizations/large\_global/check\_nonic\_characters.m}.
Every remaining simple factor has nonzero central $L$-value, so its
rational points form a finite group by Kato's theorem.

For the other rows we compute the winding submodule for every
Galois orbit of even Dirichlet characters. The sum of the dimensions,
weighted by the character degrees, is the complete genus, including
oldforms. A character component with full winding submodule has
nonzero central $L$-value on every simple factor. For every component
whose winding submodule is not full, we check that its character is
trivial at $a$. This gives the same finite-image assertion without
requiring the rank of the full Jacobian to be zero.
The degree-independent winding records are copied without alteration;
all additional winding records are computed in the degree-nine folder.
At level $529$, a full-rank $T_2$ Krylov matrix modulo $2531$ gives
the characteristic-zero rank certificate for the character of order
$253$. Every denominator is a unit at $2531$, and the rank-$45$
matrix certificate is independently replayed.

We apply \cite[Proposition~4.4]{DN} with $q=3$ and the conservative
constant $k_3=7$ when $a=3$ or $3^{-1}$, and $k_3=8$ otherwise.
The table gives the exact data. The last column is the least degree
of a CM point, extracted from the accompanying tables of
Clark--Genao--Pollack--Saia.
\[
\begin{array}{r|r|r|r|r|r}
 N&a&\operatorname{ord}\langle a\rangle&k_3&I_N&d_{\rm CM}X_1(N)\\\hline
'''
for r in rr:text+='&'.join(str(r[k]) for k in ['N','a','diamond_order','kq','index','minimum_CM_degree'])+r'\\'+'\n'
text+=r'''\end{array}
\]
In every row $3\nmid N$, $9k_3<G_N$, and $d_{\rm CM}X_1(N)>9$.
These inequalities are recomputed for degree $9$; no exclusion in
another degree is transferred as such.

A hypothetical $K$-rational point of order $N$ gives an effective
rational trace divisor $D$ of degree $9$ on $X_1(N)$, supported away
from the cusps. The global criterion gives
$\langle a\rangle D=D$. No CM point can occur in its support by the
last column. The diamond action on non-CM noncuspidal points is free,
so every orbit has length $\operatorname{ord}\langle a\rangle$.
That order must divide $\deg D=9$, contrary to every row of the table.
\end{proof}
'''
# Split the displayed order list into short aligned rows for portable typesetting.
nums=[str(r['N']) for r in rr]
text=text.replace('\\[\n'+names+'.\n\\]','\\[\n\\begin{gathered}\n'+',\\ '.join(nums[:len(nums)//2])+',\\\\\n'+',\\ '.join(nums[len(nums)//2:])+'.\n\\end{gathered}\n\\]')
(H/'nonic_large_global.tex').write_text(text)
print('NONIC_LARGE_GLOBAL_PROOF_GENERATED cases='+str(len(rr)))
