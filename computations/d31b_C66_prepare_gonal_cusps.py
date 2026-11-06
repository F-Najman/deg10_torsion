from pathlib import Path
import ast,json
C=Path(__file__).resolve().parents[1]/'continuation_20261031b/C66';r=lambda n:ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'));f=r('F2_proven_frame.m');g=json.loads((C/'F2_degree10_classification.json').read_text())['gonal_classes'];want={tuple(z['class']) for z in g};N=42373650;vs={}
cs=f[7]
def rec(i,left,a,b,v):
 if i==len(cs):
  z=(a%15,b%N,0)
  if left==0 and z in want:vs.setdefault(z,list(v))
  return
 d,_,cl,_=cs[i]
 for m in range(left//d+1):rec(i+1,left-m*d,a+m*cl[0],b+m*cl[1],v+[m])
rec(0,10,0,0,[]);assert set(vs)==want
rows=[vs[tuple(z['class'])]for z in g]
s='GonalLabels:='+repr([z['class']for z in g])+';GonalCuspVectors:='+repr(rows)+';';(C/'gonal_cusp_input.m').write_text(s+'\n');print('D31B_C66_GONAL_CUSP_INPUT_COMPLETE',len(rows))
