from pathlib import Path
import ast,json,re,subprocess,time
T=time.monotonic();R=Path('work/deg10/continuation_20261031b/M3N33');p=R/'analytic_cusp_jets_v4.m';data=ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'));rows,wits,perm,idx=data;rr=[i for i,r in enumerate(rows)if r[2]==1];assert len(rr)==40;positions={j:i for i,j in enumerate(rr)};text=['40 36']
for i in rr:
 r=rows[i];assert len(r[8])==36 and all(len(v)==5 for v in r[8]);text.append(str(positions[perm[i]-1]));text.extend(' '.join(map(str,v))for v in r[8])
(R/'prefix_rank_input.txt').write_text('\n'.join(text)+'\n');subprocess.run(['g++','-O3','-std=c++17',str(R/'prefix_rank_gf4.cpp'),'-o',str(R/'prefix_rank_gf4')],check=True);subprocess.run([str(R/'prefix_rank_gf4'),str(R/'prefix_rank_input.txt'),str(R/'prefix_ranks_primary.json')],check=True);ans=json.loads((R/'prefix_ranks_primary.json').read_text());assert ans['status']=='PASS';print('D31B_M333_ALL_PREFIX_RANKS_COMPLETE',time.monotonic()-T,flush=True)
