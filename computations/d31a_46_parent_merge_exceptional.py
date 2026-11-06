from pathlib import Path
import json,ast,hashlib
H=Path(__file__).resolve().parents[1];C=H/'continuation_20261031a/C46'
def read(name):return ast.literal_eval((C/name).read_text().replace('[*','[').replace('*]',']'))
s=(C/'source_exceptional_cases.m').read_text();ex=ast.literal_eval(s.split(':=',1)[1].strip().rstrip(';').replace('[*','[').replace('*]',']'));want={(tuple(t),int(i))for t,i,parts in ex};assert len(want)==37
seen=set();nc=0;nk=0;npr=0;rows=[];jobs=[]
for sh in range(3):
 raw=read(f'parent_exceptional_hecke_s{sh}.m');nat=read(f'parent_exceptional_native_s{sh}.m');gen=json.loads((C/f'parent_exceptional_genuine_s{sh}.json').read_text());assert gen['status']=='PASS'and len(raw)==len(nat)==gen['case_count']
 for kind in ['hecke','native','genuine']:
  p=H/f'jobs/d31a_46_parent_exceptional_{kind}_s{sh}.status.json';d=json.loads(p.read_text());assert d['state']=='EXITED'and d['actual_exit_code']==0
  script=H/d['script'];assert hashlib.sha256(script.read_bytes()).hexdigest()==d['script_sha256'];jobs.append(str(p.relative_to(H)))
 for a,b,g in zip(raw,nat,gen['cases']):
  typ,idx,attempts=a;key=(tuple(typ),int(idx));assert key in want and key not in seen;seen.add(key);assert b[:2]==a[:2]and g['type']==typ and g['index']==idx
  assert len(attempts)==1 and attempts[0][0]==5 and attempts[0][2]==0 and b[2][0][2]==0
  for part in attempts[0][1]:nc+=1;npr+=2+len(part[7])
  rows.append({'type':typ,'index':idx,'q':5,'h0':0})
 nk+=gen['geometric_kernels']
assert seen==want and nc==54 and nk==324
out={'status':'PASS','scope':'Independent Hecke exclusion for the37 source supports not rejected by the separate quotient-class sieve; this is not by itself the full C46 proof','source_supports':37,'support_type_counts':{str(list(t)):sum(tuple(r['type'])==t for r in rows)for t in [(10,),(6,4),(5,5)]},'source_components':nc,'native_two_generator_prime_references':npr,'genuine_geometric5_kernels':nk,'manual_Velu_checks':nk,'all_h0_zero':True,'all_source_diamond_and_target_keys_and_norms_verified':True,'all9jobs_exited0':True,'jobs':jobs,'cases':sorted(rows,key=lambda r:(r['type'],r['index']))}
(C/'parent_exceptional_summary.json').write_text(json.dumps(out,indent=2)+'\n');print('C46_PARENT_MERGED_EXCEPTIONAL_COMPLETE',37,nc,npr,nk)
