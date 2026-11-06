from sage.all import *
from pathlib import Path
import time,json
T=time.monotonic();D=Path('work/deg10/continuation_20261031b/M3N18');Gamma=GammaH(54,[19]);M=ModularSymbols(Gamma,2,sign=0);C=M.cuspidal_submodule();E=M.eisenstein_submodule();B=matrix(QQ,list(C.basis_matrix().rows())+list(E.basis_matrix().rows()));pr=B.inverse()[:,:C.dimension()];L=C.integral_structure().basis_matrix();CL=L*pr;ML=M.integral_structure().basis_matrix();allperiods=ML*pr*CL.inverse();cs=M.cusps();assert len(cs)==36;perm=[]
for c in cs:
 a,b=ZZ(c.numerator()),ZZ(c.denominator());a*=7
 if b==0:b=54
 while gcd(a,b)!=1:a+=54
 n=Cusp(a,b);hits=[i for i,x in enumerate(cs)if n.is_gamma_h_equiv(x,Gamma)[0]];assert len(hits)==1;perm.append(hits[0])
assert all(perm[perm[perm[i]]]==i for i in range(36));base=next(i for i in range(36)if perm[i]==i)
cc=matrix(QQ,[M.modular_symbol([cs[base],c]).element()*pr*CL.inverse()for c in cs]);den=cc.denominator();H=matrix(ZZ,list((den*cc).rows())+list((den*identity_matrix(ZZ,20)).rows()));S,U=H.hermite_form(include_zero_rows=False,transformation=True);assert U*H==S and S.nrows()==20;A=den*S.inverse();assert A.denominator()==1;A=matrix(ZZ,A);order=abs(A.det());X=matrix(QQ,list((cc.matrix_from_rows(perm)-matrix(QQ,[cc.row(perm[base])]*36)).rows())+[[0]*20]*20);target=U*X*den*S.inverse();assert target.denominator()==1;target=matrix(ZZ,target)
# Every relation maps to a relation; the induced Galois action has order three.
def zero(v):return (vector(ZZ,v)*A.inverse()).denominator()==1
assert all(zero(v*target)for v in A.rows());assert all(zero(v)for v in (target**3-1).rows())
image_lattice=matrix(ZZ,list(A.rows())+list((target-1).rows())).row_module(ZZ).basis_matrix();fixedorder=abs(image_lattice.det());seen=set();orbs=[]
for i in range(36):
 if i in seen:continue
 orb=[];j=i
 while j not in orb:orb.append(j);seen.add(j);j=perm[j]
 orbs.append(orb)
assert sorted(map(len,orbs))==[1]*18+[3]*6
norm=matrix(QQ,[sum((cc.row(i)for i in o),vector(QQ,20))for o in orbs]);coords=norm*den*S.inverse();assert coords.denominator()==1;NL=matrix(ZZ,list(A.rows())+list(coords.rows())).row_module(ZZ).basis_matrix();normorder=order/abs(NL.det());assert normorder==20253807
periodden=allperiods.denominator();S2=matrix(ZZ,list((periodden*allperiods).rows())+list((periodden*identity_matrix(ZZ,20)).rows())).row_module(ZZ).basis_matrix();assert (S2/periodden*den*S.inverse()).denominator()==1 and (S/den*periodden*S2.inverse()).denominator()==1
rep={'status':'PASS','GammaH':[54,[19]],'relative_dimension':int(M.dimension()),'cuspidal_dimension':int(C.dimension()),'geometric_cuspidal_order':int(order),'geometric_invariants':[int(v)for v in A.smith_form()[0].diagonal()if abs(v)>1],'Frob7_permutation_1based':[int(i+1)for i in perm],'Frob_fixed_cusp_signature':{'1':18,'3':6},'Frob_fixed_geometric_class_order':int(fixedorder),'global_cusp_orbit_norm_subgroup_order':int(normorder),'fixed_classes_over_orbit_norms_index':int(fixedorder/normorder),'elapsed_seconds':time.monotonic()-T}
(D/'sage_geometric_cuspidal_group.json').write_text(json.dumps(rep,indent=2)+'\n')
def rat(v):return [[int(a.numerator()),int(a.denominator())]for a in v]
packet={'group_relations':[list(map(int,r))for r in A.rows()],'galois_matrix':[list(map(int,r))for r in target.rows()],'cusp_numerators_denominators':[[int(c.numerator()),int(c.denominator())]for c in cs],'base_cusp_index':base+1,'cusp_class_coordinates':[list(map(int,r))for r in (cc*den*S.inverse()).rows()],'integral_cuspidal_lattice':[rat(r)for r in L.rows()],'projection':[rat(r)for r in pr.rows()],'lattice_basis':[list(map(int,r))for r in S.rows()],'lattice_denominator':int(den)};(D/'sage_geometric_cuspidal_presentation.json').write_text(json.dumps(packet,separators=(',',':'))+'\n');print('D31B_M318_SAGE_GEOMETRIC_CUSPIDAL_COMPLETE',rep,flush=True)
