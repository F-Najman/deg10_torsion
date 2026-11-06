from pathlib import Path
from collections import Counter
from math import gcd
import ast,json,hashlib,os
ROOT=Path(__file__).resolve().parent;HOME=ROOT.parents[1];REPO=HOME.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
src=rd(ROOT/'source70_p3.m');assert len(src[2])==773 and src[1]==[0,0,0,0,0,48,408,735,2144,5880]
p=rd(ROOT/'frobenius_p3.m');assert len(p)==195;tr=[];closed=[]
for k in range(1,11):
 z=-k*p[194-k]-sum(p[194-j]*tr[k-j-1]for j in range(1,k));tr.append(z);v=3**k+1-z-sum(j*closed[j-1]for j in range(1,k)if k%j==0);assert v%k==0;closed.append(v//k)
# Reconstruct every primitive Tate cusp, every global orbit and Frobenius3.
n=70;labels={(a,b)for a in range(n)for b in range(n)if gcd(gcd(a,b),n)==1};unseen=set(labels);cusps=[]
while unseen:
 a,b=min(unseen);o={(s*(a+k*b)%n,s*b%n)for k in range(n)for s in[-1,1]};assert o<=labels;unseen-=o;cusps.append(sorted(o))
ix={x:i for i,C in enumerate(cusps)for x in C}
def perm(u):return[ix[(u*C[0][0]%n,C[0][1])]for C in cusps]
def orbits(perms):
 unseen=set(range(len(cusps)));out=[]
 while unseen:
  x=min(unseen);todo=[x];o={x}
  while todo:
   a=todo.pop()
   for pp in perms:
    b=pp[a]
    if b not in o:o.add(b);todo.append(b)
  unseen-=o;out.append(tuple(sorted(o)))
 return set(out)
glob=orbits([perm(u)for u in range(1,n)if gcd(u,n)==1]);fin=orbits([perm(3)]);assert glob==fin
sig=Counter(map(len,fin));assert dict(sig)=={1:24,4:6,6:4,12:2}
assert[closed[k-1]-sig[k]for k in range(1,11)]==src[1]
primary=json.loads((ROOT/'A2_all.json').read_text())[2];assert len(primary)==773 and[v[0]for v in primary]==list(range(1,774))
native={}
for sh in range(4):
 for row in rd(ROOT/('parent_native_s%d.m'%sh)):assert row[0]not in native;native[row[0]]=row
assert sorted(native)==list(range(1,774))
actual=json.loads((ROOT/'independent_genuine_isogenies.json').read_text());assert actual['status']=='PASS'and actual['actual_cyclic2_kernels']==2328
nrefs=0;nmaps=0;firstfails=[]
for producer,sage in zip(primary,actual['cases']):
 ci,d,maps=producer;row=native[ci];assert row[1]==d==sage['degree']and len(row[2])==len(maps)==len(sage['maps'])
 if maps[0][5]:firstfails.append(ci)
 assert any(m[5]==0 for m in maps)
 for mi,(pr,na,ss)in enumerate(zip(maps,row[2],sage['maps']),1):
  assert na[0]==mi==ss['map']and na[1]==ss['source_key']and na[3]==ss['diamond2_key']and na[2]==pr[1]and na[4]==pr[3]and na[6]==pr[5]==ss['h0'];assert na[2]*na[1][0]==d and na[4]*na[3][0]==2*d
  expected={}
  for e,w,key in na[5]:
   assert e*d==w*key[0];s=json.dumps(key);v=expected.setdefault(s,[0,0,key[0],key]);v[0]+=e;v[1]+=w
  assert sorted(expected.values(),key=str)==sorted(ss['targets'],key=str)
  nrefs+=2+len(na[5]);nmaps+=1
assert firstfails==[53,78,84]and nmaps==776 and nrefs==3489
for v in json.loads((ROOT/'reuse_manifest.json').read_text()):assert sha(REPO/v['source'])==sha(REPO/v['copy'])==v['sha256']
markers={'d31a_C70_A2_pilot_v2':'D31A_C70_A2_COMPLETE','d31a_C70_independent_genuine_v2':'D31A_C70_ALL_GENUINE_ISOGENY_CONTROLS_COMPLETE','d31a_even_p3_low':'D31A_P3_ALL_SOURCE_CENSUS_COMPLETE','d31a_even_p3_d8':'D31A_P3_ALL_SOURCE_CENSUS_COMPLETE','d31a_even_p3_d9':'D31A_P3_ALL_SOURCE_CENSUS_COMPLETE','d31a_even_p3_d10':'D31A_P3_ALL_SOURCE_CENSUS_COMPLETE','d31a_even_p3_counts':'D31A_SOURCE_FROBENIUS_COMPLETE'}
for sh in range(4):markers['d31a_C70_A2_s'+str(sh)]='D31A_C70_A2_COMPLETE';markers['d31a_70_parent_native_s'+str(sh)]='C70_PARENT_NATIVE_COMPLETE'
diagnostics={'d31a_C70_A2_pilot':'Initial heterogeneous source serialization rejected; replaced by pilot_v2','d31a_C70_independent_genuine':'Initial comparison did not aggregate coincident closed targets; replaced by genuine_v2'}
jobs=[]
for name in sorted(set(markers)|set(diagnostics)):
 sp=HOME/'jobs'/(name+'.status.json');a=json.loads(sp.read_text());assert a['state']=='EXITED'and sha(HOME/a['script'])==a['script_sha256'];text=(HOME/a['log']).read_text()
 if name in markers:
  assert a['actual_exit_code']==0 and markers[name]in text,name
  assert not any(v in text for v in['User error:','Runtime error','Traceback (most recent call last):']),name
 jobs.append({'purpose':name,'state':a['state'],'actual_exit_code':a['actual_exit_code'],'successful_certificate':name in markers,'diagnostic':diagnostics.get(name),'script':a['script'],'log':a['log'],'status':str(sp.relative_to(REPO))})
active=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  aa=[v.decode(errors='replace')for v in(p/'cmdline').read_bytes().split(b'\0')if v]
  if aa and(Path(aa[0]).name.startswith('magma')or Path(aa[0]).name.startswith('python'))and any(Path(v).name.startswith(('d31a_C70_','d31a_70_parent_'))for v in aa[1:]):active.append(int(p.name))
 except OSError:pass
assert not active,active
summary={'status':'COMPLETE EXCLUSION','degree':10,'group':'C70','source_genus':97,'source_rank':0,'gonality_bound':[378,25],'finite_NC_closed_counts':src[1],'source_diamond_orbits':773,'source_orbits_by_degree':[4,34,65,180,490],'first_projection_excludes':770,'second_projection_case_ids':firstfails,'RR_target':'X1(35)','RR_target_genus':25,'operator':'T2-1-2diamond2','global_annihilator_justification':'DEHMZ rational cuspidal generation and literal point-level cusp boundary; no dyadic torsion-injectivity claim','tested_projection_maps':776,'independent_native_prime_references':3489,'independent_actual_cyclic2_kernels':2328,'all_source_diamond_and_target_keys_match':True,'remaining_cases':0,'lifecycle':{'all_recorded_C70_jobs_exited':True,'active_workers':active,'paused_workers':0,'automatic_launchers':0}}
(ROOT/'certificate_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths={p for p in ROOT.rglob('*')if p.is_file()and p.name!='certificate_manifest.json'and'__pycache__'not in p.parts}
paths.update((HOME/'computations').glob('d31a_C70_*'));paths.update((HOME/'computations').glob('d31a_even_p3_*'));paths.add(HOME/'run_job.py')
for job in jobs:paths.add(REPO/job['status']);paths.add(HOME/job['log']);paths.add(HOME/job['script'])
paths.update(p for p in(HOME/'shared/mdmagma/v2').rglob('*')if p.is_file())
files=[{'path':str(p.relative_to(REPO)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
out={'status':'PASS','group':'C70','degree':10,'summary_sha256':sha(ROOT/'certificate_summary.json'),'files':files,'jobs':jobs}
(ROOT/'certificate_manifest.json').write_text(json.dumps(out,indent=2)+'\n')
print('C70_FINAL_PACKET_ALL_HASHES_PASS',len(files),len(markers),len(diagnostics));print('manifest_sha256',sha(ROOT/'certificate_manifest.json'))
