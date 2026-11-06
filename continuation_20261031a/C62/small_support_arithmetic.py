"""Exact new degree-ten support types on the proven X1(31)/F3 quotient."""
from pathlib import Path
import json,itertools,hashlib
H=Path(__file__).resolve().parent;D=H/'inputs/C62_X131'
c=json.loads((D/'characters62.json').read_text());cl=json.loads((D/'classification62.json').read_text());mods=[16,31,61]
def add(a,b,m):return tuple((x+y)%m for x,y in zip(a,b))
Hs=[]
for m in mods:
 B=c['matrix'+str(m)];zero=(0,)*len(B);seen={zero};gens=[tuple(row[i-1]for row in B)for i in c['cusp_FB_indices']]
 todo=[zero]
 while todo:
  a=todo.pop()
  for b in gens:
   z=add(a,b,m)
   if z not in seen:seen.add(z);todo.append(z)
 Hs.append(seen)
assert list(map(len,Hs))==[4,1,1]
rows=[a for a in cl['tests']if a['degree']==5];assert len(rows)==90
point={a['FB_index']:a for a in map(json.loads,(D/'finite_points62.jsonseq').read_text().splitlines())if a['type']=='point'}
def val(a,m):return tuple(a['image'+str(m)])
def testpair(a,b):return all(add(val(a,m),val(b,m),m)in hs for m,hs in zip(mods,Hs))
pairs=[]
for a,b in itertools.combinations_with_replacement(rows,2):
 yes=testpair(a,b);pairs.append(dict(indices=[a['FB_index'],b['FB_index']],survives=yes))
doubles=[]
for a in rows:
 if point[a['FB_index']]['rational_two_points']==1:doubles.append(dict(index=a['FB_index'],survives=testpair(a,a)))
x=dict(status='PASS',pair_tests=pairs,quadratic_fibre_tests=doubles,pair_survivors=[a for a in pairs if a['survives']],quadratic_survivors=[a for a in doubles if a['survives']],
 method='Exact finite-primary character arithmetic; all repetitions allowed; rational cusps vanish in the quotient',
 input_sha256={q:hashlib.sha256((D/q).read_bytes()).hexdigest()for q in ['characters62.json','classification62.json','finite_points62.jsonseq']})
(H/'small_support_arithmetic.json').write_text(json.dumps(x,indent=2)+'\n')
print('D31A_C62_SMALL_SUPPORT_ARITHMETIC',len(pairs),len(doubles),len(x['pair_survivors']),len(x['quadratic_survivors']))
