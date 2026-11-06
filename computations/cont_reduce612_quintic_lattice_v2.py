from pathlib import Path
import ast,json
from sage.all import ZZ,matrix,vector,MixedIntegerLinearProgram
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
a=ast.literal_eval((D/'M6N12_D5_all_source_lifts.m').read_text().replace('[*','[').replace('*]',']'));sig,ds,lifts=a
A=matrix(ZZ,sig).transpose();L=A.row_module();U=L.basis_matrix().LLL();assert U.nrows()==23 and U*vector(ZZ,ds)==0
out=[]
for k in [2791]:
 v=lifts[k-1][1];c=vector(ZZ,v);c[ds.index(1)]+=10-sum(v[i]*ds[i] for i in range(24));assert c*vector(ZZ,ds)==10
 p=MixedIntegerLinearProgram(maximization=False,solver='GLPK');t=p.new_variable(integer=True,nonnegative=False)
 for i in range(24):p.add_constraint(c[i]+sum(t[j]*U[j,i] for j in range(23)),min=0,max=10//ds[i])
 p.set_objective(0);p.solve();z=vector(ZZ,[ZZ(round(p.get_values(t)[j])) for j in range(23)]);eff=c+z*U
 assert all(x>=0 for x in eff) and eff*vector(ZZ,ds)==10 and eff-c in L
 out.append({'index':k,'signed':list(c),'effective':list(eff),'unit_lattice_basis_coefficients':list(z)})
 print('DECIC612_EFFECTIVE_L10',k,list(eff),flush=True)
(D/'M6N12_effective_quintic_lifts_v2.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
(D/'M6N12_effective_quintic_lifts_v2.m').write_text('[* '+', '.join('[* '+str(r['index'])+', '+str(r['effective'])+' *]' for r in out)+' *]\n')
print('DECIC612_EFFECTIVE_L10_LATTICE_COMPLETE',flush=True)
