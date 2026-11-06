from pathlib import Path
H=Path(__file__).resolve().parents[1]
jobs=[]
for d,count in [(4,4),(5,11)]:
 for i in range(1,count+1):
  stem=f'cont_sieve_4_24_D5_d{d}_r{i}'
  s=(H/'computations/cont_sieve_4_24_D5.m').read_text()
  s=s.replace('for dd in [4, 5] do',f'for dd in [{d}] do').replace('for i in [1..#rs] do',f'assert #rs eq {count};\n for i in [{i}] do')
  s=s.replace('cont_sieve_4_24_D5_arithmetic.txt',stem+'_arithmetic.txt')
  (H/'computations'/(stem+'.m')).write_text(s);jobs.append((stem,600))
for j in range(6):
 stem=f'cont_sieve_3_24_D5_shard{j}'
 s=(H/'computations/cont_sieve_3_24_D5.m').read_text().replace('for dd in [3, 4, 5] do','for dd in [5] do')
 s=s.replace('for i in [1..#rs] do',f'assert #rs eq 72;\n for i in [k:k in [1..#rs]|(k-1) mod 6 eq {j}] do')
 s=s.replace('cont_sieve_3_24_D5_arithmetic.txt',stem+'_arithmetic.txt')
 (H/'computations'/(stem+'.m')).write_text(s);jobs.append((stem,600))
stem='cont_sieve_3_24_D5_lower_supports'
s=(H/'computations/cont_sieve_3_24_D5.m').read_text().replace('for dd in [3, 4, 5] do','for dd in [3,4] do').replace('cont_sieve_3_24_D5_arithmetic.txt',stem+'_arithmetic.txt')
(H/'computations'/(stem+'.m')).write_text(s);jobs.append((stem,600))
import json
(H/'continuation_20261029/shard_jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
print('Prepared',len(jobs),'bounded independent Hecke shards.')
