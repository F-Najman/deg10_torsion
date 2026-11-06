from sage.all import *
from pathlib import Path
import json,time
N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');exec(compile((N/'framed_prefix.snapshot.py').read_text(),str(N/'framed_prefix.snapshot.py'),'exec'));N=Path('work/deg10/continuation_20261031b/M3N18/endomorphism_gluing');START=time.monotonic()
a=json.loads((N/'actual_cusp_gluing_report.json').read_text());assert a['remaining_pairs']==2
inds=[i for i in range(24)if i!=base];CM=matrix(ZZ,[c7[i]for i in inds]);S=matrix(ZZ,list(CM.rows())+list(R7.rows()));kk=S.left_kernel().basis_matrix().matrix_from_columns(range(23)).hermite_form();assert kk.nrows()==kk.ncols()==23;assert abs(kk.det())==20253807
TT=zero_matrix(ZZ,23,24)
for j,i in enumerate(inds):TT[j,i]=1;TT[j,base]=-ds[i]
W=diagonal_matrix(ZZ,ds);B=(kk*TT*W).LLL();assert B.nrows()==23 and B.rank()==23
GS=[]
for row in B:
 x=vector(QQ,row)
 for y in GS:x-=row.dot_product(y)/y.dot_product(y)*y
 GS.append(x)
def babai(v):
 r=vector(QQ,v)
 for row,star in reversed(list(zip(B,GS))):
  fac=round(r.dot_product(star)/star.dot_product(star));r-=fac*row
 assert all(z in ZZ for z in r)
 return vector(ZZ,r)
out=[]
for row in a['observations']:
 co=row['cusp_coefficients'];v=vector(ZZ,[co[i]for i in inds])*TT;v[base]+=5;smallw=babai(v*W);small=smallw*W.inverse();assert all(z in ZZ for z in small);small=vector(ZZ,small);assert small.dot_product(vector(ZZ,ds))==5;cl=small*matrix(ZZ,c7);assert kr7(list(cl-vector(ZZ,row['source_class'])))==kr7([0]*n7);pred=small*matrix(ZZ,c13);assert key13(list(pred-vector(ZZ,row['target_class'])))==key13([0]*n13)
 out.append({'target_id':row['target_id'],'source_pattern':row['source_pattern'],'source_normalizing_Aut':row['source_normalizing_Aut'],'target_descriptor':row['target_descriptor'],'actual_degree5_cusp_divisor_coefficients':list(map(int,small)),'degree':5,'positive_degree':int(sum(ds[i]*max(0,small[i])for i in range(24))),'negative_degree':int(sum(ds[i]*max(0,-small[i])for i in range(24))),'source_class':row['source_class'],'target_class':row['target_class'],'native_frames_7_and_13_exact':True});print('M318_SHORT_CUSP_DIVISOR',row['target_id'],list(small),out[-1]['positive_degree'],out[-1]['negative_degree'],flush=True)
rep={'status':'PASS','meaning':'A=5c0+sum ci(Ci-degi*c0), simplified toactual24cuspcoefficients; degree5; sameunmultipliedclassasresidualcandidateatboth7and13','cusp_degrees':list(map(int,ds)),'cusp_unit_signatures':rd('cusp_automorphisms.m')[0],'base_cusp_index_1based':base+1,'data':out}
(N/'short_cusp_divisors.json').write_text(json.dumps(rep,indent=2)+'\n')
def mag(v):
 if isinstance(v,(int,Integer)):return str(v)
 flat=all(isinstance(z,(int,Integer))for z in v);return('['if flat else'[*')+','.join(mag(z)for z in v)+(']'if flat else'*]')
(N/'short_cusp_divisors.m').write_text(mag([rep['cusp_unit_signatures'],ds,[[r['target_id'],r['actual_degree5_cusp_divisor_coefficients']]for r in out]])+'\n');print('D31B_M318_SHORT_CUSP_DIVISORS_COMPLETE',len(out),time.monotonic()-START,flush=True)
