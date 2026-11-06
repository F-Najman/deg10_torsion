from pathlib import Path
import ast,json,subprocess,time,hashlib
D=Path(__file__).resolve().parents[1];C=D/'continuation_20261031b/C74';read=lambda n:ast.literal_eval((C/n).read_text().replace('[*','[').replace('*]',']'))
def cols(m):
 r,c,v=m;assert len(v)==r*c;return[sum(v[i*c+j]<<i for i in range(r))for j in range(c)]
cs=read('rational_cusp_jets.m');j=[cols(m)for m in cs[2]];assert len(j)==18 and all(len(v)==10 for v in j)
ss=cols(read('sextic_restriction.m')[2]);assert len(ss)==6
exe=C/'prefix_ranks';src=D/'computations/d31b_C74_prefix_ranks.cpp';subprocess.run(['g++','-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
reports=[]
for name,total,base,want in [('rational',10,[],8436285),('sextic',4,ss,5985)]:
 p=C/(name+'_prefix_rank_input.txt');p.write_text(f'18 10 {len(base)} {total}\n'+'\n'.join(' '.join(map(str,v))for v in j)+'\n'+' '.join(map(str,base))+'\n');t=time.monotonic();r=subprocess.run([str(exe),str(p)],text=True,capture_output=True,check=True);print(r.stdout,flush=True);vals=r.stdout.strip().split();assert int(vals[-3])==want and vals[-2]=='0';reports.append({'case':name,'patterns':want,'failures':0,'rank':int(vals[-1]),'seconds':time.monotonic()-t,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(C/'prefix_ranks_verified.json').write_text(json.dumps({'status':'PASS','algorithm':'incremental binary pivot basis along complete multiplicity recursion','reports':reports},indent=2)+'\n');print('D31B_C74_PREFIX_RANK_DRIVER_COMPLETE',flush=True)
