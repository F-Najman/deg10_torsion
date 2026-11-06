from pathlib import Path
from collections import Counter
import ast,json,itertools,struct,subprocess,hashlib
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
def part(keys,labels,reps):
 assert len(keys)==len(labels)==len(set(map(tuple,keys)))
 assert set(labels)==set(range(1,len(reps)+1)) and all(labels[j-1]==i for i,j in enumerate(reps,1))
 return {i:frozenset(tuple(k) for k,l in zip(keys,labels) if l==i) for i in range(1,len(reps)+1)}
m,n,q,cs,nc,fp,counts=read(D/'cont_local_4_20_q9_census.m');assert (m,n,q)==(4,20,9)
tr=[0];cl=[0]
for k in range(1,11):tr.append(-k*fp[50-k]-sum(fp[50-j]*tr[k-j] for j in range(1,k)))
for k in range(1,6):
 a=9**k+1-tr[2*k]-sum(j*cl[j] for j in range(1,k) if k%j==0);assert a%k==0;cl.append(a//k)
assert cl[1:]==counts==[24,108,320,888,11904] and nc==[0,96,320,888,11904] and Counter(cs)=={1:24,2:12}
quintic={};reference=None
for sh in range(12):
 d,keys,labels,reps,tests,left=read(D/f'cont_sieve_4_20_quintic_shard{sh}_arithmetic.txt');assert d==5 and not left and len(keys)==11904 and len(reps)==124
 pp=part(keys,labels,reps);full=frozenset(pp.values())
 if reference is None:reference=full
 else:assert reference==full
 for i,q,pr in tests:assert q==13 and pr==0 and pp[i] not in quintic;quintic[pp[i]]=sh
assert set(quintic)==set(reference) and len(quintic)==124
mixed={};quadreference=None
for sh in range(8):
 keys,labels,reps,ykeys,tests,left=read(D/f'cont_sieve_4_20_mixed_shard{sh}_arithmetic.txt');assert len(keys)==96 and len(ykeys)==len(set(map(tuple,ykeys)))==320 and len(reps)==2
 pp=part(keys,labels,reps)
 if quadreference is None:quadreference=(keys,labels,reps,ykeys)
 else:assert (keys,labels,reps,ykeys)==quadreference
 truth=[]
 for typ,i,j,q,pr in tests:
  assert typ in (22,23) and 1<=i<=2 and 1<=j<=(96 if typ==22 else 320) and q==(29 if typ==22 else 13)
  ident=(typ,i,j);assert ident not in mixed;mixed[ident]=pr
  if pr:truth.append([typ,i,j])
 assert truth==left
assert set(mixed)=={(t,i,j) for t,lim in [(22,96),(23,320)] for i in (1,2) for j in range(1,lim+1)}
assert not any(pr for (t,i,j),pr in mixed.items() if t==23)
left0={(i,j) for (t,i,j),pr in mixed.items() if pr};assert len(left0)==8
keys,labels,reps,tests,left=read(D/'M4N20_pair_followup_arithmetic.txt');assert (keys,labels,reps)==quadreference[:3]
assert {tuple(t[:2]) for t in tests}==left0
pairleft={tuple(t) for t in left};assert pairleft=={(1,82),(2,84)}
assert {p for p in left0 if all(pr for i,j,q,pr in tests if (i,j)==p)}==pairleft
lower,left=read(D/'M4N20_lower_arithmetic.txt');assert len(lower)==3
for d,keys,labels,reps,tests in lower:
 pp=part(keys,labels,reps);assert len(keys)==nc[d-1] and len(reps)=={2:2,3:4,4:12}[d]
 if d<4:assert {i for i,q,pr in tests if q==29 and pr==0}==set(pp)
 else:q4part=pp;q4k=keys;q4l=labels;q4r=reps;q4bad={pp[i] for i,q,pr in tests if pr==0};q4needed=set(pp.values())-q4bad
_,keys,labels,reps,tests,left=read(D/'M4N20_quartic_followup_arithmetic.txt');pp=part(keys,labels,reps);assert set(pp.values())==set(q4part.values())
assert {pp[i] for i,q,pr in tests}==q4needed
q4left={pp[i] for i in left};assert len(q4left)==2
assert q4left=={p for p in q4needed if all(pr for i,q,pr in tests if pp[i]==p)}
q4keys,q4labels,q4reps=keys,labels,reps
pars,stabs=read(D/'M4N20_D5_stabilizers_arithmetic.txt');assert len(pars)==96 and len(stabs)==4
norm=lambda v:min(tuple(x%k for x,k in zip(v,[4,4,4,4,5])),tuple(-x%k for x,k in zip(v,[4,4,4,4,5])))
def five(p):a,b,c,d=p;return [a,b,c,d%4,d%5]
P=[five(p) for p in pars];assert len(set(map(norm,P)))==96
assert all((a*d-b*c)%4==1 and u in (1,2,3,4) for a,b,c,d,u in P)
def product(x,y):
 a,b,c,d,u=x;e,f,g,h,v=y;return norm([a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h,u*v])
id0=norm([1,0,0,1,1]);covered_pairs=set();covered_quartics=set()
for typ,descr,ks,images,stab,free in stabs:
 assert len(images)==96 and stab==[i+1 for i,im in enumerate(images) if im==ks] and free
 group={norm(P[i-1]) for i in stab};assert len(group)==len(stab) and id0 in group
 assert all(product(P[i-1],P[j-1]) in group for i,j in itertools.product(stab,repeat=2))
 for i in free:assert i in stab and P[i-1][-1] in (1,4) and norm(P[i-1])!=id0 and product(P[i-1],P[i-1])==id0
 if typ==22:
  i,j=descr;assert tuple(descr) in pairleft and ks==sorted([quadreference[0][quadreference[2][i-1]-1],quadreference[0][j-1]]) and len(stab)==4;covered_pairs.add(tuple(descr))
 else:
  assert typ==4 and len(stab)==8;i=descr[0];ky=q4keys[q4reps[i-1]-1];assert ks==[ky] and pp[i] in q4left;covered_quartics.add(pp[i])
assert covered_pairs==pairleft and covered_quartics==q4left
# Rebuild every projected canonical constraint space independently over F9 in C++.
poly,records,tested,bad,rjets,qjets,cjets=read(D/'M4N20_gonality6_arithmetic.txt')
assert len(poly)==3 and poly[-1]==1 and len(records)==5512 and tested==11156288 and not bad
expected={(1,i,j,k) for i in range(1,25) for j in range(i,25) for k in range(j,25)}|{(2,i,j) for i in range(1,109) for j in range(1,25)}|{(3,i) for i in range(1,321)}
assert {tuple(r[0]) for r in records}==expected
binary=D/'gonality420_replay_input.bin'
with binary.open('wb') as f:
 f.write(struct.pack('<III',poly[0],poly[1],len(records)))
 for place in rjets:
  for jet in place:f.write(bytes(row[0] for row in jet))
 for collection in (qjets,cjets):
  for mat in collection:f.write(bytes(x for row in mat for x in row))
 for desc,mat in records:
  dd=[desc[0]]+[x-1 for x in desc[1:]];dd+=[0]*(4-len(dd));f.write(struct.pack('<HHHH',*dd));assert len(mat)==22 and all(len(row)==24 for row in mat);f.write(bytes(x for row in mat for x in row))
exe=D/'replay_gonality420';subprocess.run(['g++','-O3','-std=c++17',str(D/'replay_gonality420.cpp'),'-o',str(exe)],check=True)
r=subprocess.run([str(exe),str(binary)],check=True,capture_output=True,text=True);print(r.stdout,flush=True)
tate=json.loads((D/'cusp_tate420_replay.json').read_text())[0];assert tate['all_residue_rational_cusps_globally_rational'] and tate['geometric_cusps']==48
report={'group':[4,20],'absolute_degree':10,'relative_degree':5,'genus':25,'native_NC_counts':nc,'quintic_orbits_excluded':124,'mixed_patterns_checked':832,'quartic_residual_orbits':2,'quadratic_pair_residual_orbits':2,'free_involutions_verified':4,'gonality_F9_lower_bound':7,'gonality_degree3_complements':5512,'canonical_rank_tests_independently_replayed':11156288,'all_5512_projected_kernel_spaces_recomputed':'PASS','complete_Hecke_orbit_and_mixed_support_coverage':'PASS','stabilizer_and_fine_5_level_freeness':'PASS','classification_case':'EXCLUDED'}
(D/'M4N20_independent_replay.json').write_text(json.dumps(report,indent=2)+'\n');print('DECIC420_INDEPENDENT_REPLAY_COMPLETE',report,flush=True)
