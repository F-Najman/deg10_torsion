from pathlib import Path
import ast,json
from sage.all import QQ,NumberField,PolynomialRing,matrix
P=Path('work/deg10/C226_resolution_20261030')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
def q(v):return QQ(v[0])/QQ(v[1])
out=[]
for idx in [30,24]:
 fc,bc,cc,m=read(P/f'exact_points_{idx}.m')[0];R=PolynomialRing(QQ,'x');f=R([q(v) for v in fc]);K=NumberField(f,'a');b=K([q(v) for v in bc]);c=K([q(v) for v in cc]);elts=[K(1)]
 while True:
  cand=elts+[b*v for v in elts]+[c*v for v in elts];V=matrix(QQ,[v.list() for v in cand]).row_space();new=[K(v) for v in V.basis()]
  if len(new)==len(elts):break
  elts=new
 assert len(elts)==10;out.append({'class':idx,'coordinate_field_degree':10,'b_minpoly_degree':int(b.minpoly().degree()),'c_minpoly_degree':int(c.minpoly().degree())})
(P/'coordinate_fields.json').write_text(json.dumps({'status':'PASS','cases':out},indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
