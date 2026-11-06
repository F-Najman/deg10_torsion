from pathlib import Path
from collections import Counter
from functools import lru_cache
import json,time,hashlib
def safe_dump(obj,*args,**kw):
 kw.setdefault("default",int);return json.dumps(obj,*args,**kw)
from sage.all import GF,PolynomialRing,EllipticCurve,matrix,identity_matrix,lcm,gcd
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C44';start=time.monotonic();k=GF(3);PF=PolynomialRing(k,'x');cache={};audits=[]
def load(n):return[json.loads(s)for s in(C/n).read_text().splitlines()if s.strip()]
def pointbc(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P.xy();b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;r=b4/b3;n2=b2-r*r-a1*r
 return[-n2**3/b3**2,-((2*r+a1)*n2-b3)/b3]
def unitxy(E,P):
 b,c=pointbc(E,P);r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def key2(vals):
 L=vals[0].parent();d=lcm([v.minpoly().degree()for v in vals]);found=False
 for bits in[[1,0],[0,1],[1,1],[1,2]]:
  g=bits[0]*vals[0]+bits[1]*vals[1];mp=g.minpoly()
  if mp.degree()==d:found=True;break
 assert found;n=L.degree();vv=lambda v:__import__('sage.all',fromlist=['vector']).vector(k,[L(v).polynomial()[i]for i in range(n)]);M=matrix(k,[vv(g**i)for i in range(d)]);assert M.rank()==d
 return[int(d),bits,list(map(int,mp.list())),[list(map(int,M.transpose().solve_right(vv(v))))for v in vals]]
def tag(vals):return safe_dump(key2(vals),separators=(',',':'))
@lru_cache(maxsize=30000)
def jetkey(pol,xy):
 d=len(pol)-1;L=GF(3**d,'a',modulus=PF(list(pol)));a=L.gen();v=[PF(list(w))(a)for w in xy];z=tag(v);assert json.loads(z)[0]==d;return z
@lru_cache(maxsize=5000)
def field(poly):return GF(3**(len(poly)-1),'a',modulus=PF(list(poly)))
def source(r):
 L=field(tuple(r['field_polynomial']));a=L.gen();aa=[PF(v)(a)for v in r['E_ai']];E=EllipticCurve(L,aa);P=E(*[PF(v)(a)for v in r['P_xy']]);d=r['source_degree'];assert L.degree()==d and 44*P==E(0)and 22*P!=E(0)and 4*P!=E(0);ky=tag(unitxy(E,P));assert json.loads(ky)[0]==d;return E,P,ky

def audit(E0,P0,ell,ky):
 ck=(ell,ky)
 if ck in cache:return cache[ck]
 d=E0.base_ring().degree();q=3**d;card=E0.cardinality();trace=q+1-card;M=matrix(GF(ell),[[0,-q],[1,trace]]);A=identity_matrix(GF(ell),2)
 for m in range(1,ell*(ell*ell-1)+1):
  A*=M
  if A.is_one():break
 else:raise AssertionError('Frobenius torsion extension bound')
 L=GF(3**(d*m),'w');RL=PolynomialRing(L,'z');z=RL.gen();pol=E0.base_ring().modulus();root=RL(list(pol)).roots(multiplicities=False)[0];embed=lambda v:PF(v.polynomial().list())(root);E=EllipticCurve(L,[embed(v)for v in E0.a_invariants()]);P=E(*[embed(v)for v in P0.xy()]);assert tag(unitxy(E,P))==ky
 assert 44*P==E(0)and 22*P!=E(0)and 4*P!=E(0)
 xs=E.division_polynomial(ell).roots(multiplicities=False);pts=[T for x in xs for T in E.lift_x(x,all=True)];assert len(xs)==(ell*ell-1)//2 and len(pts)==ell*ell-1
 kernels={tuple(__import__('sage.all',fromlist=['prod']).prod(z-(j*T)[0]for j in range(1,(ell+1)//2)).list()):T for T in pts};assert len(kernels)==ell+1;counts=Counter()
 for co,T in kernels.items():
  assert ell*T==E(0)and T!=E(0);phi=E.isogeny(RL(list(co)));EE=phi.codomain();Q=phi(P);assert phi.degree()==ell and 44*Q==EE(0)and 22*Q!=EE(0)and 4*Q!=EE(0)
  nz=[j*T for j in range(1,ell)];v=(P[0]+sum(((P+S)[0]-S[0]for S in nz),L(0)),P[1]+sum(((P+S)[1]-S[1]for S in nz),L(0)));assert v==Q.xy();counts[tag(unitxy(EE,Q))]+=1
 target=Counter()
 for zz,c in counts.items():
  deg=json.loads(zz)[0];assert(d*c)%deg==0;target[zz]=d*c//deg
 u=Counter({ky:1});u[tag(unitxy(E,ell*P))]+=ell
 result={'ell':ell,'source_degree':d,'source_key':json.loads(ky),'diamond_key':json.loads(tag(unitxy(E,ell*P))),'curve_cardinality':int(card),'torsion_extension_factor':m,'target_norms':[[json.loads(z),int(n)]for z,n in sorted(target.items())],'geometric_kernel_target_counts':[[json.loads(z),int(n)]for z,n in sorted(counts.items())]}
 cache[ck]=(target,u);audits.append(result);checkpoint.write(safe_dump(dict(type="unique_source_check",**result),separators=(",",":"))+"\n");checkpoint.flush()
 if len(audits)%5==0:print('D31A_C44_GENUINE_UNIQUE_SOURCE',len(audits),'ell',ell,'degree',d,'factor',m,'seconds',time.monotonic()-start,flush=True)
 return target,u


assert SHARD in range(6)
checkpoint=(C/('independent_genuine_shard%d_v2.checkpoint.jsonl'%SHARD)).open('x');checkpoint.write(safe_dump({'type':'header','shard':SHARD,'shards':6,'partition':'newA5 case(ci-1)mod6; q7 onlyshard0'})+'\n');checkpoint.flush()
def run_shard(name,ell,path,original=None):
 rows=load(path);assert rows[0]['type']=='header'and rows[-1]['type']=='complete';expected=rows[-1]['cases'];selected=[i for i in range(1,expected+1)if(name=='newA5'and(i-1)%6==SHARD)or(name!='newA5'and SHARD==0)]
 sources={i:[]for i in selected};jets={i:{q:Counter()for q in['T','U']}for i in selected};meta={};rep10=[];rep5=[];pair_keys=[];units=[u for u in range(1,23)if gcd(u,44)==1];coverage10=set();coverage5=set();source_case_keys={}
 for r in rows:
  ci=r.get('case')
  if ci not in sources:continue
  if r['type']=='source':sources[ci].append(r)
  elif r['type']=='case':meta[ci]=r
  elif r['type']=='jet':ky=jetkey(tuple(r['field_polynomial']),tuple(tuple(v)for v in r['xy']));jets[ci][r['side']][ky]+=r['multiplicity']
 for ci in selected:
  actualT=Counter();actualU=Counter();source_keys=[]
  for r in sources[ci]:
   E,P,ky=source(r);source_keys.append(ky);t,u=audit(E,P,ell,ky);actualT.update(t);actualU.update(u)
   if name=='newA5'and ci<=596:
    ks={tag(unitxy(E,a*P))for a in units}
    if r['source_degree']==10:coverage10.update(ks);rep10.append(ky)
    else:coverage5.update(ks);rep5.append(ky)
  assert actualT==jets[ci]['T']and actualU==jets[ci]['U'],(SHARD,name,ci,'physical target/norm binding')
  assert sum(json.loads(z)[0]*n for z,n in actualT.items())==meta[ci]['positive_degree']and sum(json.loads(z)[0]*n for z,n in actualU.items())==meta[ci]['negative_degree']
  source_case_keys[ci]=source_keys
  if name=='newA5'and ci>=597:assert len(source_keys)==2;pair_keys.append(tuple(source_keys))
  if original:assert source_keys==original[ci]
  checkpoint.write(safe_dump({'type':'case_check','stage':name,'case':ci,'source_keys':source_keys,'positive_degree':meta[ci]['positive_degree'],'all_target_norms_and_negative_terms_bound':True},separators=(',',':'))+'\n');checkpoint.flush()
  if len(source_case_keys)%5==0 or ci==selected[-1]:print('D31A_C44_GENUINE_SHARD_CASE',SHARD,name,ci,len(source_case_keys),len(selected),'seconds',time.monotonic()-start,flush=True)
 return{'stage':name,'case_ids':selected,'representative_source10_keys':rep10,'representative_source5_keys':rep5,'diamond_coverage10_keys':sorted(coverage10),'diamond_coverage5_keys':sorted(coverage5),'normalized_pair_keys':pair_keys,'all_source_and_diamond_Hecke_targets_and_norms_exact':True}
report=run_shard('newA5',5,'new_finite_jets.normalized.jsonl');reports=[report]
if SHARD==0:
 rows=load('new_finite_jets.normalized.jsonl');original={j:[source(next(r for r in rows if r['type']=='source'and r['case']==ci))[2]]for j,ci in[(1,510),(2,511)]};reports.append(run_shard('residualA7',7,'q7_finite_jets.normalized.jsonl',original))
res={'status':'PASS','shard':SHARD,'shards':6,'stages':reports,'unique_physical_sources_checked':len(audits),'genuine_geometric_kernels':sum(r['ell']+1 for r in audits),'manual_Velu_coordinate_checks':sum(r['ell']+1 for r in audits),'all_marked_point_orders44_verified':True,'all_frobenius_orbit_and_norm_multiplicities_verified':True,'all_bound_moduli_keys':audits,'elapsed_seconds':time.monotonic()-start}
(C/('independent_genuine_shard%d_v2.json'%SHARD)).write_text(safe_dump(res,indent=2)+'\n');checkpoint.write(safe_dump({'type':'complete','shard':SHARD,'unique_sources':len(audits),'cases':len(report['case_ids'])})+'\n');checkpoint.close();print('D31A_C44_INDEPENDENT_GENUINE_SHARD_COMPLETE',SHARD,len(report['case_ids']),len(audits),time.monotonic()-start,flush=True)
