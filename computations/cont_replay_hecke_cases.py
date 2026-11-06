from pathlib import Path
import ast,json,hashlib,re,math
from collections import Counter
H=Path(__file__).resolve().parents[1];D=H/'continuation_20261029'
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'))
def census(stem,g):
 m,n,q,cusps,nc,f,reported=read(D/(stem+'_census.m'));e=2 if q in (4,9) else 1;p=2 if q==4 else 3 if q==9 else q
 assert len(f)==2*g+1 and f[-1]==1
 tr=[0]
 for k in range(1,5*e+1):tr.append(-k*f[2*g-k]-sum(f[2*g-j]*tr[k-j] for j in range(1,k)))
 cl=[0]
 for k in range(1,6):
  a=q**k+1-tr[e*k]-sum(d*cl[d] for d in range(1,k) if k%d==0);assert a%k==0;cl.append(a//k)
 assert cl[1:]==reported==[nc[k-1]+cusps.count(k) for k in range(1,6)]
 return nc
reports=[]
for m,n,genus,stem,patterns,expected in [
 (4,24,33,'cont_local_4_24_q5',['cont_sieve_4_24_D5_d*_arithmetic.txt'],{4:4,5:11}),
 (3,24,25,'cont_local_3_24_q7',['cont_sieve_3_24_D5_shard*_arithmetic.txt','cont_sieve_3_24_D5_lower_supports_arithmetic.txt'],{3:1,4:13,5:72})]:
 nc=census(stem,genus);reference={};covered={d:set() for d in expected};checked=[]
 files=sorted({p for pattern in patterns for p in D.glob(pattern)})
 for p in files:
  record,remaining=read(p);assert not remaining
  jobstem=p.name.removesuffix('_arithmetic.txt');status=json.loads((H/'jobs'/(jobstem+'.status.json')).read_text());assert status['state']=='EXITED' and status['actual_exit_code']==0
  for dd,keys,labels,reps,tests in record:
   assert dd in expected and len(keys)==len(set(map(tuple,keys)))==len(labels)==nc[dd-1]
   assert len(reps)==expected[dd]
   assert {labels[i-1] for i in reps}==set(range(1,len(reps)+1))
   orbits={i:frozenset(tuple(k) for k,lab in zip(keys,labels) if lab==i) for i in range(1,len(reps)+1)}
   assert all(orbits.values()) and sum(map(len,orbits.values()))==len(keys)
   part=frozenset(orbits.values())
   if dd in reference:assert part==reference[dd]
   else:reference[dd]=part
   for i,q,pr in tests:
    assert q==13 and pr==0;covered[dd].add(orbits[i])
  checked.append(p.name)
 assert all(covered[d]==set(reference[d]) for d in expected)
 if m==3:
  rec,left=read(D/'cont_sieve_3_24_cubic31_arithmetic.txt');assert not left and rec[0][-1]==[[1,31,0]]
  assert 31%24==7 and 31%3==1
 tate=json.loads((D/'cusp_tate_replay.json').read_text());tate=next(r for r in tate if (r['m'],r['n'])==(m,n));assert tate['all_residue_rational_cusps_globally_rational']
 reports.append({'group':[m,n],'absolute_degree':10,'relative_degree':5,'genus':genus,'noncuspidal_counts':nc,'complete_Hecke_orbit_counts':expected,'all_images_nonprincipal':True,'independent_Newton_and_exact_orbit_coverage':'PASS','files':checked})
 print('DECIC_HECKE_CASE_REPLAY_PASS',m,n,expected,flush=True)
# The descent in the C155 argument is exact, and every possible good residue-field trace is checked.
assert pow(16,5,155)==1 and all(pow(16,i,155) not in (1,154) for i in range(1,5)) and 16%5==1
assert {pow(16,i,31) for i in range(5)}=={1,2,4,8,16}
for q in [2,4]:
 traces=[t for t in range(-math.isqrt(4*q),math.isqrt(4*q)+1) if (q+1-t)%5==0]
 assert traces==([-2] if q==2 else [0])
 for t in traces:assert pow((t*t-4*q)%31,15,31)==30
for N,gg,bad in [(31,26,[]),(155,841,[1,12,3])]:
 log=(H/f'jobs/cont_winding{N}.log').read_text();rows=re.findall(fr'WINDING N={N} order=(\d+) conductor=(\d+) relative_dimension=(\d+) winding_dimension=(\d+) coefficient_degree=(\d+)',log)
 phi=lambda n:sum(math.gcd(i,n)==1 for i in range(1,n+1))
 assert sum(phi(int(o))*int(d) for o,c,d,w,k in rows)==gg
 assert [int(o) for o,c,d,w,k in rows if d!=w]==bad
 if N==155:assert all(math.gcd(o,5)==1 for o in bad)
 for p in D.glob(f'krylov_N{N}_*.json'):
  d=json.loads(p.read_text());mod=d['p'];n=d['d'];A=[d['A'][i*n:(i+1)*n] for i in range(n)];v=d['v'];mat=[]
  for _ in range(n):mat.append(v);v=[sum(v[i]*A[i][j] for i in range(n))%mod for j in range(n)]
  rank=0
  for j in range(n):
   r=next((r for r in range(rank,n) if mat[r][j]%mod),None)
   if r is None:continue
   mat[rank],mat[r]=mat[r],mat[rank];a=pow(mat[rank][j],-1,mod);mat[rank]=[x*a%mod for x in mat[rank]]
   for r in range(rank+1,n):
    a=mat[r][j];mat[r]=[(x-a*y)%mod for x,y in zip(mat[r],mat[rank])]
   rank+=1
  if rank<n:
   marker=fr'WINDING N={N} order={d["order"]} conductor=\d+ relative_dimension={n} winding_dimension={n}'
   assert re.search(marker,log),'Deficient Krylov rank without a completed exact replacement'
 primary=json.loads((D/'Orlic_source.json').read_text());assert hashlib.sha256((D/'Orlic_2407.14512.pdf').read_bytes()).hexdigest()==primary['sha256']
assert 'DECIC155_LOCAL_TWO_COMPLETE' in (H/'jobs/cont_C155_local2.log').read_text()
reports.append({'group':[1,155],'absolute_degree':10,'descent_degree':2,'diamond':16,'diamond_order':5,'quotient':'X1(31)/<2>','quotient_genus':6,'quotient_Q_gonality_lower_bound':5,'gonality_source':'Orlic, Tetragonal intermediate modular curves, Proposition 3.6, arXiv:2407.14512','local_q2_trace':-2,'local_q4_trace':0,'no_residue_31_isogeny_with_5_torsion':True,'independent_arithmetic_and_input_replay':'PASS'})
(D/'new_cases_independent_replay.json').write_text(json.dumps(reports,indent=2)+'\n')
print('DECIC_NEW_HECkE_AND_C155_REPLAY_COMPLETE',flush=True)
