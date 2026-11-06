from sage.all import *
from pathlib import Path
from collections import defaultdict,Counter
import json,time,sys,hashlib
P=Path('work/deg10/continuation_20261031b/C56');D=P/'kummer_rr';START=time.monotonic();CID=int(sys.argv[1]);QELL=11 if CID==65 else 5;TAG=f'direct_case{CID}_A{QELL}';Fp=GF(3);RP=PolynomialRing(Fp,'t')
def rd(s):return json.loads((D/(TAG+'_'+s+'.m')).read_text().replace('[*','[').replace('*]',']'))
def freeze(v):return tuple(map(freeze,v))if isinstance(v,list)else v
md=rd('model_snapshot');op=rd('operator_snapshot');bs=rd('basis');jets=rd('jet_tables');result=rd('result');assert len(md[0])==16 and result[:2]==[CID,QELL]and result[2:5]==bs[:3];N0,N1=bs[1:3];NB=N0+N1;assert all(v[1]==[1]for v in md[0]);CURVE=[v[0]for v in md[0]]
def pointbc(E,Q):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=Q.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;s=b4/b3;n1=2*s+a1;n2=b2-s*s-a1*s;return[-n2**3/b3**2,1-n1*n2/b3]
def key(vals):
 K=vals[0].parent();nn=K.degree();dd=lcm([v.minpoly().degree()for v in vals]);bits=next(z for z in[[1,0],[0,1],[1,1],[1,2]]if(z[0]*vals[0]+z[1]*vals[1]).minpoly().degree()==dd);g=bits[0]*vals[0]+bits[1]*vals[1];M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(nn)]for i in range(dd)]);co=[list(map(int,M.transpose().solve_right(vector(Fp,[K(v).polynomial()[j]for j in range(nn)]))))for v in vals];return[int(dd),bits,list(map(int,g.minpoly().list())),co]
def physical(ph):
 K=GF(3**(len(ph[0])-1),'a',modulus=RP(ph[0]));aa=K.gen();decode=lambda v:RP(v)(aa);E=EllipticCurve(K,[decode(v)for v in ph[1]]);Q=E([decode(v)for v in ph[2]]);assert 56*Q==E(0)and 28*Q!=E(0)and 8*Q!=E(0);assert key(pointbc(E,Q))==ph[3];return K,E,Q
physcache={}
def physical_group(desc,ph):
 kk=freeze(ph)
 if kk not in physcache:
  K,E,Q=physical(ph);T=28*Q;Z=Q+T;a1,a2,a3,a4,a6=E.a_invariants();r=T[0];A2=a2+a1*a1/4;xx=Q[0]+Z[0]-T[0];yy=Q[1]+Z[1]-T[1]+(a1*xx+a3)/2;r2=-A2-2*r;delta=xx-r2;bb=delta**3/yy**2;cc=(3*r2+A2)*delta**2/yy**2;yk=key([bb,cc]);sk=ph[3];assert sk[0]in[yk[0],2*yk[0]];physcache[kk]=(freeze(yk),freeze(sk),int(yk[0]),int(sk[0]))
 yk,sk,e,dq=physcache[kk];assert desc[0]==dq and desc[1]==e and freeze(desc[4])==sk;return yk,sk,e,dq
groups={};ddeg=0
for ix in op[2]:
 row=op[0][ix-1];assert row[6]==QELL
 terms=[(row[3],-1,row[8]),(row[5],-QELL,row[9])]+[(h[2],h[1],h[3])for h in row[7]]
 for desc,n,ph in terms:
  yk,sk,e,dq=physical_group(desc,ph);gr=groups.setdefault(yk,{'degree':e,'coeffs':defaultdict(int)});assert gr['degree']==e;gr['coeffs'][sk]+=n;ddeg+=n*dq
assert ddeg==0
Adeg=sum(g['degree']*max([0]+list(g['coeffs'].values()))for g in groups.values());assert Adeg==bs[0]
class Eval:
 def __init__(self,S,U,V):self.S=S;self.U=U;self.V=V;self.vp=[S(1)];self.up=[S(1)];self.pc={};self.rc={};self.fc={}
 def pol(self,cs):
  kk=tuple(cs)
  if kk in self.pc:return self.pc[kk]
  while len(self.vp)<len(cs):self.vp.append(self.vp[-1]*self.V)
  v=self.S(0)
  for i,c in enumerate(cs):
   if c%3==1:v+=self.vp[i]
   elif c%3==2:v-=self.vp[i]
  self.pc[kk]=v;return v
 def rat(self,v):
  kk=freeze(v)
  if kk not in self.rc:self.rc[kk]=self.pol(v[0])/self.pol(v[1])
  return self.rc[kk]
 def ef(self,v):
  kk=freeze(v)
  if kk in self.fc:return self.fc[kk]
  while len(self.up)<len(v):self.up.append(self.up[-1]*self.U)
  ans=self.S(0)
  for i,z in enumerate(v):
   if z[0]:ans+=self.rat(z)*self.up[i]
  self.fc[kk]=ans;return ans
 def equation(self):
  ans=self.S(0)
  while len(self.up)<len(CURVE):self.up.append(self.up[-1]*self.U)
  for i,z in enumerate(CURVE):ans+=self.pol(z)*self.up[i]
  return ans
 def derivative(self,which):
  ans=self.S(0)
  while len(self.up)<len(CURVE):self.up.append(self.up[-1]*self.U)
  for i,z in enumerate(CURVE):
   if which==2:
    if i%3:ans+=i*self.pol(z)*self.up[i-1]
   else:ans+=self.pol([(j*z[j])%3 for j in range(1,len(z))])*self.up[i]
  return ans

def curve_scalar(u,v,der=None):
 if der=='u':return sum((i*RP(z)(v)*u**(i-1)for i,z in enumerate(CURVE)if i),u.parent()(0))
 if der=='v':return sum((RP(z).derivative()(v)*u**i for i,z in enumerate(CURVE)),u.parent()(0))
 return sum((RP(z)(v)*u**i for i,z in enumerate(CURVE)),u.parent()(0))
allmat=[[]for _ in range(NB)];seen=set();reports=[];denoms={tuple(c[1])for f in bs[3]+bs[4]+md[1]+md[2]for c in f}
for ji,jt in enumerate(jets,1):
 pd,a,e,modpol,bcvec,z0vec,Bnative,H0native,H1native,blocks,geometry=jt;K=GF(3**(len(modpol)-1),'a',modulus=RP(modpol));alpha=K.gen();decode=lambda v:RP(v)(alpha);bb,cc=map(decode,bcvec);z0=decode(z0vec);assert z0*z0==bb and bb!=0;yk=freeze(key([bb,cc]));assert yk in groups and yk not in seen;seen.add(yk);gr=groups[yk];assert e==gr['degree']==pd[0]and a==max([0]+list(gr['coeffs'].values()));which,prec,Udata,Vdata=geometry;assert prec==128;Ucoeff=list(map(decode,Udata));Vcoeff=list(map(decode,Vdata));u0,v0=Ucoeff[0],Vcoeff[0];assert curve_scalar(u0,v0)==0;du=curve_scalar(u0,v0,'u');dv=curve_scalar(u0,v0,'v');assert du!=0 or dv!=0;assert which==(2 if du!=0 else 1);assert lcm([u0.minpoly().degree(),v0.minpoly().degree()])==e
 # Recompute the actual implicit arc by precision-doubling Hensel lifting.
 uc=[u0];vc=[v0]
 for pr in[2,4,8,16,32,64,128]:
  S=LaurentSeriesRing(K,'s',default_prec=pr);s=S.gen();U=sum((S(c)*s**i for i,c in enumerate(uc)),S(0))+O(s**pr);V=sum((S(c)*s**i for i,c in enumerate(vc)),S(0))+O(s**pr)
  if which==2:V=S(v0)+s+O(s**pr)
  else:U=S(u0)+s+O(s**pr)
  ev=Eval(S,U,V);ff=ev.equation();dd=ev.derivative(which);assert dd.valuation()==0
  if which==2:U-=ff/dd
  else:V-=ff/dd
  assert Eval(S,U,V).equation().valuation()>=pr;uc=[U[i]for i in range(pr)];vc=[V[i]for i in range(pr)]
 assert uc==Ucoeff and vc==Vcoeff
 # The canonical source inverse determines the two cover branches physically.
 EE=EllipticCurve(K,[0,-2*cc,0,cc*cc-4*bb*(1-bb-cc),0]);sourcekeys=[]
 for sg in[1,-1]:
  zz=sg*z0;xx=2*bb+cc+2*zz;PP=EE(xx,2*zz*xx);assert 56*PP==EE(0)and 28*PP!=EE(0)and 8*PP!=EE(0);sourcekeys.append(key(pointbc(EE,PP)))
 dq=sourcekeys[0][0];assert dq in[e,2*e]and sourcekeys[1][0]==dq
 if dq==2*e:assert sourcekeys[0]==sourcekeys[1];sourcekeys=sourcekeys[:1];signs=[1]
 else:assert sourcekeys[0]!=sourcekeys[1];signs=[1,-1]
 assert set(gr['coeffs'])<={freeze(ky)for ky in sourcekeys};coeffs=[gr['coeffs'].get(freeze(ky),0)for ky in sourcekeys];needs=[a-n for n in coeffs];assert all(n>=0 for n in needs);maxneed=max(needs);assert maxneed==len(Bnative)and maxneed>0
 # Select enough working precision using exact denominator valuations.
 p0=v0.minpoly();vvord=next(i for i in range(1,128)if vc[i]!=0);maxpole=0
 for den in denoms|{tuple(c[1])for f in pd[2]for c in f}:
  poly=RP(den);mult=0
  while poly.degree()>=p0.degree()and poly%p0==0:poly//=p0;mult+=1
  maxpole=max(maxpole,mult*vvord)
 work=maxneed+2*maxpole+a+6;assert work<=128;S=LaurentSeriesRing(K,'s',default_prec=work);s=S.gen();U=sum((S(c)*s**i for i,c in enumerate(uc[:work])),S(0))+O(s**work);V=sum((S(c)*s**i for i,c in enumerate(vc[:work])),S(0))+O(s**work);ev=Eval(S,U,V)
 assert all(ev.ef(g).valuation()>0 for g in pd[2]);assert(ev.ef(md[2][0])-U).valuation()>=maxneed and(ev.ef(md[2][1])-V).valuation()>=maxneed;Bser=ev.ef(md[1][0]);Cser=ev.ef(md[1][1]);assert Bser.valuation()==0 and Bser[0]==bb and Cser[0]==cc;assert[Bser[i]for i in range(maxneed)]==list(map(decode,Bnative))
 Z=S(z0)+O(s**work)
 for it in range(8):Z=(Z+Bser/Z)/2
 assert(Z*Z-Bser).valuation()>=maxneed;H0=[];H1=[]
 for funcs,store,out in[(bs[3],H0native,H0),(bs[4],H1native,H1)]:
  for f,stored in zip(funcs,store):
   val=s**a*ev.ef(f);assert val.valuation()>=0 and val.precision_absolute()>=maxneed;assert[val[i]for i in range(maxneed)]==list(map(decode,stored));out.append(val)
 assert len(H0)==N0 and len(H1)==N1
 expectedblocks=[(ky,co,sg,need)for ky,co,sg,need in zip(sourcekeys,coeffs,signs,needs)if need];assert len(expectedblocks)==len(blocks);blockcols=0
 for block,(ky,co,sg,need)in zip(blocks,expectedblocks):
  kkey,n,sign,deg,nn,gvec,gpol,mat=block;assert kkey==ky and n==co and sign==sg and deg==dq and nn==need;valsbc=pointbc(EE,EE(2*bb+cc+2*sg*z0,2*sg*z0*(2*bb+cc+2*sg*z0)));g=ky[1][0]*valsbc[0]+ky[1][1]*valsbc[1];assert g==decode(gvec)and list(map(int,g.minpoly().list()))==gpol;M=matrix(Fp,[[K(g**i).polynomial()[j]for j in range(K.degree())]for i in range(dq)]);assert M.rank()==dq;values=[[h[i]for i in range(need)]for h in H0]+[[sg*(Z*h)[i]for i in range(need)]for h in H1];VV=matrix(Fp,[[K(z).polynomial()[j]for j in range(K.degree())]for row in values for z in row]);CC=M.transpose().solve_right(VV.transpose()).transpose();assert CC*M==VV;rows=[[int(CC[ii*need+ni,j])for ni in range(need)for j in range(dq)]for ii in range(NB)];nr,nc,flat=mat;assert nr==NB and nc==dq*need and sum(rows,[])==flat
  for ii,row in enumerate(rows):allmat[ii].extend(row)
  blockcols+=nc
 reports.append({'local_row':ji,'base_degree':e,'physical_ambient_degree':K.degree(),'pole_bound':a,'max_needed_jet_order':maxneed,'work_precision':work,'geometry_precision_independently_rebuilt':128,'uniformizer':'V'if which==2 else'U','matrix_columns':blockcols,'base_prime_generators_vanish':True,'all_basis_jets_and_residue_power_coordinates_entrywise_equal':True});print('D31B_C56_SIEVE_RR_LOCAL_PASS',CID,ji,e,K.degree(),work,blockcols,time.monotonic()-START,flush=True)
# Every nontrivial local constraint is represented; omitted groups impose none.
assert all(yk in seen or(a:=max([0]+list(gr['coeffs'].values())))==0 and all(n==0 for n in gr['coeffs'].values())or len(gr['coeffs'])==2 and len(set(gr['coeffs'].values()))==1 and min(gr['coeffs'].values())>0 for yk,gr in groups.items());M=matrix(Fp,allmat);assert M.nrows()==NB and M.ncols()==2*Adeg==result[5]and M.rank()==result[6]==NB and result[7]==0
out={'status':'PASS','case':CID,'hecke':QELL,'base_curve_genus':25,'source_level':56,'pole_bound_degree':Adeg,'basis_dimensions':[N0,N1],'matrix_rows':NB,'matrix_columns':M.ncols(),'rank':M.rank(),'h0':0,'all_geometry_series_coefficients_128_independently_recomputed':True,'all_operator_signed_closed_source_groups_bound_to_actual_physical_levels':True,'all_basis_rational_functions_reexpanded_in_genuine_local_arcs':True,'all_sqrt_branches_bound_to_actual_56_level':True,'all_matrix_entries_and_residue_power_coordinates_match':True,'local_rows':reports,'input_hashes':{s:hashlib.sha256((D/(TAG+'_'+s+'.m')).read_bytes()).hexdigest()for s in['model_snapshot','operator_snapshot','basis','jet_tables','result']},'elapsed_seconds':time.monotonic()-START};(P/f'sieve_independent_RR_case{CID}_A{QELL}.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C56_SIEVE_RR_COMPLETE',CID,QELL,M.nrows(),M.ncols(),M.rank(),0,time.monotonic()-START,flush=True)
