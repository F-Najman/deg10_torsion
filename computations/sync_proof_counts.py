from pathlib import Path
import json,re,subprocess
H=Path(__file__).resolve().parents[1]
subprocess.run(['python3',str(H/'candidates.py')],check=True,stdout=subprocess.PIPE)
d=json.loads((H/'candidates.json').read_text());rem=d['remaining_minimal_obstructions'];cyc=[n for m,n in rem if m==1];non=[(m,n)for m,n in rem if m!=1];ne=len(d['certified_exclusions']);nr=len(rem)
p=H/'paper/main.tex';s=p.read_text();s=re.sub(r'\$\d+\$ minimal subgroup obstructions',lambda _:f'${ne}$ minimal subgroup obstructions',s);s=re.sub(r'\$\d+\$ minimal obstructions\nremain',lambda _:f'${nr}$ minimal obstructions\nremain',s);p.write_text(s)
p=H/'paper/remaining.tex';s=p.read_text();tail=s[s.index('These are minimal unresolved subgroup obstructions.'):]
prefix='\\section{The remaining proof tasks}\nThe certified lower list has $68$ groups. Its minimal subgroup obstructions,\nsubject to $S(10)$ and the Weil pairing, form a list of $109$ groups.\n'+f'The present note excludes ${ne}$ of them, leaving ${nr}$.\n'
if cyc:
 prefix+='The remaining cyclic orders are\n\\[\n\\begin{gathered}\n'+'\\\\\n'.join(','.join(map(str,cyc[i:i+10])) for i in range(0,len(cyc),10))+'\n\\end{gathered}\n\\]\n'
if non:
 prefix+='The remaining noncyclic groups are\n\\[\n\\begin{array}{c|l}\nm&n\\text{ in }C_m\\times C_n\\\\\\hline\n'
 for m in sorted({m for m,n in non}):prefix+=str(m)+'&'+','.join(str(n)for mm,n in non if mm==m)+'\\\\\n'
 prefix+='\\end{array}\n\\]\n'
p.write_text(prefix+tail)
print('Proof counts synchronized:',ne,'excluded,',nr,'remaining')
