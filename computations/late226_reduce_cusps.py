from pathlib import Path
import ast,json
from sage.all import ZZ,QQ,matrix,vector
P=Path('work/deg10/continuation_20261029b/C226')
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
labels,base,hgens,rels,ker,rows=read(P/'cusp_preimages.m');B=matrix(ZZ,ker).LLL();BS,_=B.gram_schmidt();out=[]
for idx,desc,cl,old in rows:
 v=vector(QQ,old)
 for b,s in reversed(list(zip(B.rows(),BS.rows()))):v-=round(v.dot_product(s)/s.dot_product(s))*b
 v=vector(ZZ,v);delta=vector(ZZ,old)-v;assert delta in matrix(ZZ,ker).row_module()
 actual=v*matrix(ZZ,hgens)-vector(ZZ,cl);assert actual in matrix(ZZ,rels).row_module()
 co=list(v);co[labels.index(base)]+=10-sum(co);assert sum(co)==10
 out.append([idx,desc,co]);print(idx,desc,'degree 10 signed cusp coefficients',co,'L1',sum(abs(c) for c in co),flush=True)
(P/'small_cusp_divisors.m').write_text(str([labels,out])+'\n')
(P/'cusp_lattice_replay.json').write_text(json.dumps({'status':'PASS','labels':labels,'base_label':base,'divisors':out,'kernel_LLL_rows':[list(map(int,b)) for b in B.rows()]},indent=2,default=int)+'\n')
